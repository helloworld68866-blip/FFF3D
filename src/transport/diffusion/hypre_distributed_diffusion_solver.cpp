#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include "mesh/boundary/spherical_scalar_remap.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "physics/units/physical_constants.hpp"

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dec3d::transport {

namespace {

constexpr const char* kHypreBackend = "hypre_parcsr_gmres_boomeramg";
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kBoundaryTolerance = 1.0e-14;
constexpr int kBoomerAmgPreconditionerMaxIterations = 1;
constexpr double kBoomerAmgPreconditionerTolerance = 0.0;
constexpr int kBoomerAmgPreconditionerPrintLevel = 0;
constexpr double kBoomerAmgStrongThreshold = 0.5;
constexpr std::size_t kExpectedDistributedStencilEntriesPerRow = 8u;

using RowAccumulator = std::vector<std::vector<std::pair<std::size_t, double>>>;

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] const char* ToString(DiffusionPoleMetricMode mode) noexcept {
  switch (mode) {
    case DiffusionPoleMetricMode::cell_centered_spherical:
      return "cell_centered_spherical";
    case DiffusionPoleMetricMode::axis_regular_polar_phi:
      return "axis_regular_polar_phi";
  }
  return "unknown";
}

[[nodiscard]] bool OptionsAreValid(const GenericDiffusionHypreSolveOptions& options) noexcept {
  return std::isfinite(options.relative_tolerance) &&
         options.relative_tolerance > 0.0 &&
         options.max_iterations > 0 &&
         options.krylov_dimension > 0 &&
         options.print_level >= 0 &&
         options.require_boomeramg_preconditioner;
}

[[nodiscard]] bool FitsHypreInt(std::size_t value) noexcept {
  return value <= static_cast<std::size_t>(std::numeric_limits<HYPRE_Int>::max());
}

[[nodiscard]] int EnsureHypreInitialized() noexcept {
  const HYPRE_Int initialized = HYPRE_Initialized();
  if (initialized != 0) {
    return 0;
  }
  return HYPRE_Initialize();
}

[[nodiscard]] HYPRE_Int NoOpParSolverSetup(
    HYPRE_Solver,
    HYPRE_ParCSRMatrix,
    HYPRE_ParVector,
    HYPRE_ParVector) {
  return 0;
}

struct HyprePreconditionerSetupStatus {
  HYPRE_Int status{0};
  const char* reason{"none"};
  const char* phase{"none"};
};

[[nodiscard]] HyprePreconditionerSetupStatus ConfigureBoomerAmg(
    HYPRE_Solver boomeramg) {
  auto status = HYPRE_BoomerAMGSetPrintLevel(
      boomeramg,
      kBoomerAmgPreconditionerPrintLevel);
  if (status != 0) {
    return {status,
            "HYPRE distributed BoomerAMG print-level setup failed",
            "BoomerAMGSetPrintLevel"};
  }
  status = HYPRE_BoomerAMGSetMaxIter(
      boomeramg,
      kBoomerAmgPreconditionerMaxIterations);
  if (status != 0) {
    return {status,
            "HYPRE distributed BoomerAMG max-iteration setup failed",
            "BoomerAMGSetMaxIter"};
  }
  status = HYPRE_BoomerAMGSetTol(
      boomeramg,
      kBoomerAmgPreconditionerTolerance);
  if (status != 0) {
    return {status,
            "HYPRE distributed BoomerAMG tolerance setup failed",
            "BoomerAMGSetTol"};
  }
  status = HYPRE_BoomerAMGSetStrongThreshold(
      boomeramg,
      kBoomerAmgStrongThreshold);
  if (status != 0) {
    return {status,
            "HYPRE distributed BoomerAMG strong-threshold setup failed",
            "BoomerAMGSetStrongThreshold"};
  }
  return {};
}

struct HypreDistributedObjects {
  HYPRE_IJMatrix ij_matrix{nullptr};
  HYPRE_IJVector ij_rhs{nullptr};
  HYPRE_IJVector ij_solution{nullptr};
  HYPRE_Solver gmres{nullptr};
  HYPRE_Solver boomeramg{nullptr};

  ~HypreDistributedObjects() {
    if (gmres != nullptr) {
      HYPRE_ParCSRGMRESDestroy(gmres);
    }
    if (boomeramg != nullptr) {
      HYPRE_BoomerAMGDestroy(boomeramg);
    }
    if (ij_solution != nullptr) {
      HYPRE_IJVectorDestroy(ij_solution);
    }
    if (ij_rhs != nullptr) {
      HYPRE_IJVectorDestroy(ij_rhs);
    }
    if (ij_matrix != nullptr) {
      HYPRE_IJMatrixDestroy(ij_matrix);
    }
  }
};

[[nodiscard]] std::string DistributedFailureLine(
    const char* reason,
    int rank,
    int rank_count,
    const char* mpi_phase = "none",
    const char* hypre_phase = "none",
    int hypre_status = 0,
    int hypre_error_code = 0,
    int gmres_iterations = 0,
    double gmres_final_relative_residual = 0.0) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.distributed_failure"
      << "; failure_reason=" << reason
      << "; mpi_rank=" << rank
      << "; mpi_rank_count=" << rank_count
      << "; mpi_phase=" << mpi_phase
      << "; hypre_phase=" << hypre_phase
      << "; hypre_status=" << hypre_status
      << "; hypre_error_code=" << hypre_error_code
      << "; gmres_iterations=" << gmres_iterations
      << "; gmres_final_relative_residual=" << gmres_final_relative_residual
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  return out.str();
}

[[nodiscard]] DistributedGenericDiffusionAssemblyResult AssemblyFailure(
    const char* reason,
    const DistributedDiffusionRowOwnership& ownership) {
  DistributedGenericDiffusionAssemblyResult result;
  result.ownership = ownership;
  result.failure_reason = reason;
  result.failure_diagnostics =
      DistributedFailureLine(reason, ownership.rank, ownership.rank_count);
  return result;
}

[[nodiscard]] DistributedGenericDiffusionHypreSolveResult SolveFailure(
    const char* reason,
    const DistributedDiffusionRowOwnership& ownership,
    const char* hypre_phase = "none",
    int hypre_status = 0,
    int gmres_iterations = 0,
    double gmres_final_relative_residual = 0.0) {
  DistributedGenericDiffusionHypreSolveResult result;
  result.backend = kHypreBackend;
  result.hypre_enabled = true;
  result.ownership = ownership;
  result.gmres_iterations = gmres_iterations;
  result.gmres_final_relative_residual = gmres_final_relative_residual;
  result.failure_reason = reason;
  result.failure_diagnostics =
      DistributedFailureLine(reason,
                             ownership.rank,
                             ownership.rank_count,
                             "none",
                             hypre_phase,
                             hypre_status,
                             HYPRE_GetError(),
                             gmres_iterations,
                             gmres_final_relative_residual);
  return result;
}

[[nodiscard]] std::string BuildOwnershipReport(const DistributedDiffusionRowOwnership& ownership) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.distributed_row_ownership"
      << "; row_ownership=distributed_radial_slab"
      << "; collective_row_ownership_valid=true"
      << "; global_layout_consistent=true"
      << "; ownership_gap_free=true"
      << "; ownership_overlap_free=true"
      << "; mpi_rank=" << ownership.rank
      << "; mpi_rank_count=" << ownership.rank_count
      << "; global_radial_cells=" << ownership.global_radial_cells
      << "; global_theta_cells=" << ownership.global_theta_cells
      << "; global_phi_cells=" << ownership.global_phi_cells
      << "; global_row_count=" << ownership.global_row_count
      << "; local_global_radial_begin=" << ownership.global_radial_begin
      << "; local_global_radial_end=" << ownership.global_radial_end
      << "; local_row_begin=" << ownership.local_row_begin
      << "; local_row_end=" << ownership.local_row_end
      << "; local_row_count=" << ownership.local_row_count
      << "; canonical_state_mutated=false";
  return out.str();
}

[[nodiscard]] std::string BuildOwnershipFailureReport(
    const DistributedDiffusionRowOwnership& ownership,
    const char* reason,
    const char* mpi_phase = "MPI_Allgather") {
  std::ostringstream out;
  out << "diagnostic_id=p2.diffusion.distributed_row_ownership_failure"
      << "; failure_reason=" << reason
      << "; mpi_rank=" << ownership.rank
      << "; mpi_rank_count=" << ownership.rank_count
      << "; mpi_phase=" << mpi_phase
      << "; collective_row_ownership_valid=false"
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  return out.str();
}

struct CollectiveOwnershipValidation {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] CollectiveOwnershipValidation ValidateCollectiveGlobalLayout(
    const DistributedDiffusionRowOwnership& ownership) {
  constexpr int kFieldCount = 3;
  unsigned long long local[kFieldCount] = {
      static_cast<unsigned long long>(ownership.global_radial_cells),
      static_cast<unsigned long long>(ownership.global_theta_cells),
      static_cast<unsigned long long>(ownership.global_phi_cells)};
  std::vector<unsigned long long> gathered(
      static_cast<std::size_t>(ownership.rank_count) * static_cast<std::size_t>(kFieldCount),
      0u);
  const int status = MPI_Allgather(
      local,
      kFieldCount,
      MPI_UNSIGNED_LONG_LONG,
      gathered.data(),
      kFieldCount,
      MPI_UNSIGNED_LONG_LONG,
      ownership.communicator);

  CollectiveOwnershipValidation validation;
  if (status != MPI_SUCCESS) {
    validation.failure_reason = "collective global layout validation MPI_Allgather failed";
    validation.failure_diagnostics =
        BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
    return validation;
  }

  for (int peer = 0; peer < ownership.rank_count; ++peer) {
    const std::size_t offset =
        static_cast<std::size_t>(peer) * static_cast<std::size_t>(kFieldCount);
    if (gathered[offset] != local[0] ||
        gathered[offset + 1u] != local[1] ||
        gathered[offset + 2u] != local[2]) {
      validation.failure_reason = "inconsistent global layout across ranks";
      validation.failure_diagnostics =
          BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
      return validation;
    }
  }

  if (ownership.global_radial_cells == 0u ||
      ownership.global_theta_cells == 0u ||
      ownership.global_phi_cells == 0u) {
    validation.failure_reason = "global radial/theta/phi cells must be positive";
    validation.failure_diagnostics =
        BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
    return validation;
  }

  validation.success = true;
  return validation;
}

[[nodiscard]] CollectiveOwnershipValidation ValidateCollectiveRowOwnership(
    const DistributedDiffusionRowOwnership& ownership) {
  constexpr int kFieldCount = 9;
  std::vector<unsigned long long> local(kFieldCount, 0u);
  local[0] = static_cast<unsigned long long>(ownership.global_radial_cells);
  local[1] = static_cast<unsigned long long>(ownership.global_theta_cells);
  local[2] = static_cast<unsigned long long>(ownership.global_phi_cells);
  local[3] = static_cast<unsigned long long>(ownership.global_row_count);
  local[4] = static_cast<unsigned long long>(ownership.global_radial_begin);
  local[5] = static_cast<unsigned long long>(ownership.global_radial_end);
  local[6] = static_cast<unsigned long long>(ownership.local_row_begin);
  local[7] = static_cast<unsigned long long>(ownership.local_row_end);
  local[8] = static_cast<unsigned long long>(ownership.local_row_count);

  std::vector<unsigned long long> gathered(
      static_cast<std::size_t>(ownership.rank_count) * static_cast<std::size_t>(kFieldCount),
      0u);
  const int status = MPI_Allgather(
      local.data(),
      kFieldCount,
      MPI_UNSIGNED_LONG_LONG,
      gathered.data(),
      kFieldCount,
      MPI_UNSIGNED_LONG_LONG,
      ownership.communicator);
  CollectiveOwnershipValidation validation;
  if (status != MPI_SUCCESS) {
    validation.failure_reason = "collective row ownership validation MPI_Allgather failed";
    validation.failure_diagnostics =
        BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
    return validation;
  }

  const auto field = [&](int rank, int slot) -> unsigned long long {
    return gathered[static_cast<std::size_t>(rank) * static_cast<std::size_t>(kFieldCount) +
                    static_cast<std::size_t>(slot)];
  };

  for (int peer = 0; peer < ownership.rank_count; ++peer) {
    if (field(peer, 0) != local[0] ||
        field(peer, 1) != local[1] ||
        field(peer, 2) != local[2] ||
        field(peer, 3) != local[3]) {
      validation.failure_reason = "inconsistent global layout across ranks";
      validation.failure_diagnostics =
          BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
      return validation;
    }
  }

  unsigned long long expected_radial_begin = 0u;
  unsigned long long expected_row_begin = 0u;
  for (int peer = 0; peer < ownership.rank_count; ++peer) {
    const unsigned long long radial_begin = field(peer, 4);
    const unsigned long long radial_end = field(peer, 5);
    const unsigned long long row_begin = field(peer, 6);
    const unsigned long long row_end = field(peer, 7);
    const unsigned long long row_count = field(peer, 8);
    if (radial_begin != expected_radial_begin ||
        row_begin != expected_row_begin ||
        radial_end <= radial_begin ||
        row_end <= row_begin ||
        row_count != row_end - row_begin) {
      validation.failure_reason = "distributed row ownership has gap or overlap";
      validation.failure_diagnostics =
          BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
      return validation;
    }
    expected_radial_begin = radial_end;
    expected_row_begin = row_end;
  }

  if (expected_radial_begin != local[0] || expected_row_begin != local[3]) {
    validation.failure_reason = "distributed row ownership does not cover global layout";
    validation.failure_diagnostics =
        BuildOwnershipFailureReport(ownership, validation.failure_reason.c_str());
    return validation;
  }

  validation.success = true;
  return validation;
}

[[nodiscard]] std::size_t LocalLinearIndex(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t local_radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((local_radial * ownership.global_theta_cells) + theta) * ownership.global_phi_cells + phi;
}

[[nodiscard]] std::size_t GlobalRow(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t global_radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((global_radial * ownership.global_theta_cells) + theta) * ownership.global_phi_cells + phi;
}

[[nodiscard]] bool LocalArrayShapeMatches(
    const dec3d::core::Array3D<double>& values,
    const DistributedDiffusionRowOwnership& ownership) noexcept {
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  return values.extent_r() == local_radial_count &&
         values.extent_theta() == ownership.global_theta_cells &&
         values.extent_phi() == ownership.global_phi_cells;
}

[[nodiscard]] bool LocalFaceEffectiveShapeMatches(
    const GenericDiffusionFaceEffectiveCoefficients& face,
    const DistributedDiffusionRowOwnership& ownership) noexcept {
  if (!face.enabled) {
    return true;
  }
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  return face.radial_face_D.extent_r() == local_radial_count + 1u &&
         face.radial_face_D.extent_theta() == ownership.global_theta_cells &&
         face.radial_face_D.extent_phi() == ownership.global_phi_cells &&
         face.theta_face_D.extent_r() == local_radial_count &&
         face.theta_face_D.extent_theta() == ownership.global_theta_cells + 1u &&
         face.theta_face_D.extent_phi() == ownership.global_phi_cells &&
         face.phi_face_D.extent_r() == local_radial_count &&
         face.phi_face_D.extent_theta() == ownership.global_theta_cells &&
         face.phi_face_D.extent_phi() == ownership.global_phi_cells;
}

[[nodiscard]] double RadialCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

[[nodiscard]] double ThetaCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] double RadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  const double polar_factor =
      std::cos(geometry.theta_faces[theta]) - std::cos(geometry.theta_faces[theta + 1u]);
  const double dphi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return radius * radius * polar_factor * dphi;
}

[[nodiscard]] double MarshakGFace(
    double diffusion_coefficient,
    double center_to_boundary_distance_cm) noexcept {
  if (diffusion_coefficient <= 0.0) {
    return 0.0;
  }
  const double h = 0.5 * dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double d_over_distance = diffusion_coefficient / center_to_boundary_distance_cm;
  return (h * d_over_distance) / (h + d_over_distance);
}

[[nodiscard]] double ThetaFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  const double radial_factor =
      0.5 * (std::pow(geometry.radial_faces[radial + 1u], 2) -
             std::pow(geometry.radial_faces[radial], 2));
  const double dphi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return radial_factor * std::sin(geometry.theta_faces[theta_face]) * dphi;
}

[[nodiscard]] double PhiFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta) noexcept {
  const double radial_factor =
      0.5 * (std::pow(geometry.radial_faces[radial + 1u], 2) -
             std::pow(geometry.radial_faces[radial], 2));
  const double dtheta = geometry.theta_faces[theta + 1u] - geometry.theta_faces[theta];
  return radial_factor * dtheta;
}

[[nodiscard]] double SolidAngleWeightedMeanSinTheta(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  const double theta0 = geometry.theta_faces[theta];
  const double theta1 = geometry.theta_faces[theta + 1u];
  const double solid_angle_factor = std::cos(theta0) - std::cos(theta1);
  if (solid_angle_factor == 0.0) {
    return std::sin(ThetaCenter(geometry, theta));
  }
  const double sin_squared_integral =
      0.5 * (theta1 - theta0) -
      0.25 * (std::sin(2.0 * theta1) - std::sin(2.0 * theta0));
  return sin_squared_integral / solid_angle_factor;
}

[[nodiscard]] double PhiMetricSinTheta(
    const DistributedGenericDiffusionProblem& problem,
    std::size_t theta) noexcept {
  if (problem.pole_metric_mode == DiffusionPoleMetricMode::axis_regular_polar_phi &&
      (theta == 0u || theta + 1u == problem.ownership.global_theta_cells)) {
    return SolidAngleWeightedMeanSinTheta(problem.global_geometry, theta);
  }
  return std::sin(ThetaCenter(problem.global_geometry, theta));
}

struct DistributedGeometryMetrics {
  std::vector<double> radial_center;
  std::vector<double> theta_center;
  std::vector<double> radial_center_distance;
  std::vector<double> radial_face_area;
  std::vector<double> theta_face_area;
  std::vector<double> phi_face_area;
  std::vector<double> phi_metric_sin_theta;
  double dphi{0.0};
};

[[nodiscard]] std::size_t RadialFaceMetricIndex(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((radial_face * ownership.global_theta_cells) + theta) *
             ownership.global_phi_cells +
         phi;
}

[[nodiscard]] std::size_t ThetaFaceMetricIndex(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  return ((radial * (ownership.global_theta_cells + 1u)) + theta_face) *
             ownership.global_phi_cells +
         phi;
}

[[nodiscard]] std::size_t PhiFaceMetricIndex(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t radial,
    std::size_t theta) noexcept {
  return radial * ownership.global_theta_cells + theta;
}

[[nodiscard]] DistributedGeometryMetrics BuildDistributedGeometryMetrics(
    const DistributedGenericDiffusionProblem& problem) {
  const auto& ownership = problem.ownership;
  DistributedGeometryMetrics metrics;
  metrics.radial_center.assign(ownership.global_radial_cells, 0.0);
  metrics.theta_center.assign(ownership.global_theta_cells, 0.0);
  metrics.radial_center_distance.assign(
      ownership.global_radial_cells > 0u ? ownership.global_radial_cells - 1u : 0u, 0.0);
  metrics.radial_face_area.assign(
      (ownership.global_radial_cells + 1u) * ownership.global_theta_cells *
          ownership.global_phi_cells,
      0.0);
  metrics.theta_face_area.assign(
      ownership.global_radial_cells * (ownership.global_theta_cells + 1u) *
          ownership.global_phi_cells,
      0.0);
  metrics.phi_face_area.assign(
      ownership.global_radial_cells * ownership.global_theta_cells, 0.0);
  metrics.phi_metric_sin_theta.assign(ownership.global_theta_cells, 0.0);

  for (std::size_t gr = 0u; gr < ownership.global_radial_cells; ++gr) {
    metrics.radial_center[gr] = RadialCenter(problem.global_geometry, gr);
  }
  for (std::size_t gr = 0u; gr + 1u < ownership.global_radial_cells; ++gr) {
    metrics.radial_center_distance[gr] = metrics.radial_center[gr + 1u] -
                                         metrics.radial_center[gr];
  }
  for (std::size_t theta = 0u; theta < ownership.global_theta_cells; ++theta) {
    metrics.theta_center[theta] = ThetaCenter(problem.global_geometry, theta);
    metrics.phi_metric_sin_theta[theta] = PhiMetricSinTheta(problem, theta);
  }
  metrics.dphi = problem.global_geometry.phi_faces[1u] -
                 problem.global_geometry.phi_faces[0u];

  std::vector<double> radial_face_radius_squared(ownership.global_radial_cells + 1u, 0.0);
  for (std::size_t rf = 0u; rf <= ownership.global_radial_cells; ++rf) {
    const double radius = problem.global_geometry.radial_faces[rf];
    radial_face_radius_squared[rf] = radius * radius;
  }
  std::vector<double> radial_cell_area_factor(ownership.global_radial_cells, 0.0);
  for (std::size_t gr = 0u; gr < ownership.global_radial_cells; ++gr) {
    radial_cell_area_factor[gr] =
        0.5 * (std::pow(problem.global_geometry.radial_faces[gr + 1u], 2) -
               std::pow(problem.global_geometry.radial_faces[gr], 2));
  }
  std::vector<double> polar_face_factor(ownership.global_theta_cells, 0.0);
  std::vector<double> theta_face_sin(ownership.global_theta_cells + 1u, 0.0);
  std::vector<double> dtheta(ownership.global_theta_cells, 0.0);
  for (std::size_t t = 0u; t < ownership.global_theta_cells; ++t) {
    polar_face_factor[t] =
        std::cos(problem.global_geometry.theta_faces[t]) -
        std::cos(problem.global_geometry.theta_faces[t + 1u]);
    dtheta[t] =
        problem.global_geometry.theta_faces[t + 1u] -
        problem.global_geometry.theta_faces[t];
  }
  for (std::size_t tf = 0u; tf <= ownership.global_theta_cells; ++tf) {
    theta_face_sin[tf] = std::sin(problem.global_geometry.theta_faces[tf]);
  }
  std::vector<double> dphi_by_cell(ownership.global_phi_cells, 0.0);
  for (std::size_t p = 0u; p < ownership.global_phi_cells; ++p) {
    dphi_by_cell[p] =
        problem.global_geometry.phi_faces[p + 1u] -
        problem.global_geometry.phi_faces[p];
  }

  const std::size_t radial_face_fill_begin = ownership.global_radial_begin;
  const std::size_t radial_face_fill_end = ownership.global_radial_end;
  for (std::size_t rf = radial_face_fill_begin; rf <= radial_face_fill_end; ++rf) {
    for (std::size_t t = 0u; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0u; p < ownership.global_phi_cells; ++p) {
        metrics.radial_face_area[RadialFaceMetricIndex(ownership, rf, t, p)] =
            radial_face_radius_squared[rf] * polar_face_factor[t] * dphi_by_cell[p];
      }
    }
  }
  for (std::size_t gr = ownership.global_radial_begin; gr < ownership.global_radial_end; ++gr) {
    for (std::size_t tf = 0u; tf <= ownership.global_theta_cells; ++tf) {
      for (std::size_t p = 0u; p < ownership.global_phi_cells; ++p) {
        metrics.theta_face_area[ThetaFaceMetricIndex(ownership, gr, tf, p)] =
            radial_cell_area_factor[gr] * theta_face_sin[tf] * dphi_by_cell[p];
      }
    }
    for (std::size_t t = 0u; t < ownership.global_theta_cells; ++t) {
      metrics.phi_face_area[PhiFaceMetricIndex(ownership, gr, t)] =
          radial_cell_area_factor[gr] * dtheta[t];
    }
  }
  return metrics;
}

[[nodiscard]] double ArithmeticMean(double lhs, double rhs) noexcept {
  return 0.5 * (lhs + rhs);
}

[[nodiscard]] double RadialInterfaceD(
    const DistributedGenericDiffusionProblem& problem,
    std::size_t global_radial_face,
    std::size_t theta,
    std::size_t phi,
    double fallback_left_d,
    double fallback_right_d) noexcept {
  if (problem.face_effective_coefficients.enabled &&
      global_radial_face >= problem.ownership.global_radial_begin &&
      global_radial_face <= problem.ownership.global_radial_end) {
    return problem.face_effective_coefficients.radial_face_D(
        global_radial_face - problem.ownership.global_radial_begin,
        theta,
        phi);
  }
  return ArithmeticMean(fallback_left_d, fallback_right_d);
}

[[nodiscard]] double ThetaInterfaceD(
    const DistributedGenericDiffusionProblem& problem,
    std::size_t local_radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.theta_face_D(local_radial, theta_face, phi);
  }
  return ArithmeticMean(
      problem.local_coefficient_D(local_radial, theta_face - 1u, phi),
      problem.local_coefficient_D(local_radial, theta_face, phi));
}

[[nodiscard]] double PhiInterfaceD(
    const DistributedGenericDiffusionProblem& problem,
    std::size_t local_radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.phi_face_D(local_radial, theta, phi);
  }
  const std::size_t next_phi = (phi + 1u) % problem.ownership.global_phi_cells;
  return ArithmeticMean(
      problem.local_coefficient_D(local_radial, theta, phi),
      problem.local_coefficient_D(local_radial, theta, next_phi));
}

struct CartesianCenter {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

[[nodiscard]] CartesianCenter CellCenterCartesian(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double radius = RadialCenter(geometry, radial);
  const double polar = ThetaCenter(geometry, theta);
  const double azimuth = 0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
  const double sin_theta = std::sin(polar);
  return CartesianCenter{
      radius * sin_theta * std::cos(azimuth),
      radius * sin_theta * std::sin(azimuth),
      radius * std::cos(polar)};
}

[[nodiscard]] double SquaredDistance(CartesianCenter lhs, CartesianCenter rhs) noexcept {
  const double dx = lhs.x - rhs.x;
  const double dy = lhs.y - rhs.y;
  const double dz = lhs.z - rhs.z;
  return (dx * dx) + (dy * dy) + (dz * dz);
}

[[nodiscard]] double MappedCenterConductance(
    const DistributedGenericDiffusionProblem& problem,
    std::size_t lhs_global_radial,
    std::size_t lhs_local_radial,
    std::size_t lhs_theta,
    std::size_t lhs_phi,
    std::size_t rhs_global_radial,
    std::size_t rhs_local_radial,
    std::size_t rhs_theta,
    std::size_t rhs_phi) noexcept {
  const std::size_t lhs_global_row =
      GlobalRow(problem.ownership, lhs_global_radial, lhs_theta, lhs_phi);
  const std::size_t rhs_global_row =
      GlobalRow(problem.ownership, rhs_global_radial, rhs_theta, rhs_phi);
  const double lhs_volume = problem.global_geometry.cell_volumes[lhs_global_row];
  const double rhs_volume = problem.global_geometry.cell_volumes[rhs_global_row];
  const double d_face = ArithmeticMean(
      problem.local_coefficient_D(lhs_local_radial, lhs_theta, lhs_phi),
      problem.local_coefficient_D(rhs_local_radial, rhs_theta, rhs_phi));
  const double distance_squared = SquaredDistance(
      CellCenterCartesian(problem.global_geometry, lhs_global_radial, lhs_theta, lhs_phi),
      CellCenterCartesian(problem.global_geometry, rhs_global_radial, rhs_theta, rhs_phi));
  if (!std::isfinite(lhs_volume) || !std::isfinite(rhs_volume) ||
      !std::isfinite(d_face) || !std::isfinite(distance_squared) ||
      lhs_volume <= 0.0 || rhs_volume <= 0.0 || distance_squared <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return d_face * (0.5 * (lhs_volume + rhs_volume)) / distance_squared;
}

void AddMatrixEntry(
    RowAccumulator& rows,
    std::size_t local_row,
    std::size_t global_column,
    double value) {
  auto& entries = rows[local_row];
  for (auto& entry : entries) {
    if (entry.first == global_column) {
      entry.second += value;
      return;
    }
  }
  entries.emplace_back(global_column, value);
}

[[nodiscard]] bool OwnsGlobalRow(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t global_row) noexcept {
  return global_row >= ownership.local_row_begin &&
         global_row < ownership.local_row_end;
}

[[nodiscard]] std::size_t LocalRowFromGlobalRow(
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t global_row) noexcept {
  return global_row - ownership.local_row_begin;
}

void AddOwnedPairConductance(
    RowAccumulator& rows,
    const DistributedDiffusionRowOwnership& ownership,
    std::size_t lhs_global_row,
    std::size_t rhs_global_row,
    double conductance,
    bool count_radial_seam,
    DistributedGenericDiffusionAssemblyResult& result) {
  if (conductance == 0.0 || lhs_global_row == rhs_global_row) {
    return;
  }
  const bool lhs_owned = OwnsGlobalRow(ownership, lhs_global_row);
  const bool rhs_owned = OwnsGlobalRow(ownership, rhs_global_row);

  if (lhs_owned) {
    const std::size_t local_row = LocalRowFromGlobalRow(ownership, lhs_global_row);
    AddMatrixEntry(rows, local_row, lhs_global_row, conductance);
    AddMatrixEntry(rows, local_row, rhs_global_row, -conductance);
    if (count_radial_seam && !rhs_owned) {
      ++result.local_radial_seam_coupling_count;
      ++result.local_off_rank_column_count;
    }
  }

  if (rhs_owned) {
    const std::size_t local_row = LocalRowFromGlobalRow(ownership, rhs_global_row);
    AddMatrixEntry(rows, local_row, rhs_global_row, conductance);
    AddMatrixEntry(rows, local_row, lhs_global_row, -conductance);
    if (count_radial_seam && !lhs_owned) {
      ++result.local_radial_seam_coupling_count;
      ++result.local_off_rank_column_count;
    }
  }
}

struct DHalo {
  bool lower_received{false};
  bool upper_received{false};
  std::vector<double> lower;
  std::vector<double> upper;
};

[[nodiscard]] DHalo ExchangeCoefficientDHalo(const DistributedGenericDiffusionProblem& problem) {
  const auto& ownership = problem.ownership;
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  DHalo halo;
  halo.lower.assign(plane_size, 0.0);
  halo.upper.assign(plane_size, 0.0);

  std::vector<double> send_lower(plane_size, 0.0);
  std::vector<double> send_upper(plane_size, 0.0);
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
    for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
      const std::size_t slot = t * ownership.global_phi_cells + p;
      send_lower[slot] = problem.local_coefficient_D(0, t, p);
      send_upper[slot] = problem.local_coefficient_D(local_radial_count - 1u, t, p);
    }
  }

  MPI_Status status{};
  const int lower_rank = ownership.rank - 1;
  const int upper_rank = ownership.rank + 1;
  if (ownership.global_radial_begin > 0u) {
    MPI_Sendrecv(
        send_lower.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        lower_rank,
        1101,
        halo.lower.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        lower_rank,
        1102,
        ownership.communicator,
        &status);
    halo.lower_received = true;
  }
  if (ownership.global_radial_end < ownership.global_radial_cells) {
    MPI_Sendrecv(
        send_upper.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        upper_rank,
        1102,
        halo.upper.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        upper_rank,
        1101,
        ownership.communicator,
        &status);
    halo.upper_received = true;
  }
  return halo;
}

[[nodiscard]] double CoefficientDAt(
    const DistributedGenericDiffusionProblem& problem,
    const DHalo& halo,
    std::size_t global_radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const auto& ownership = problem.ownership;
  if (global_radial >= ownership.global_radial_begin &&
      global_radial < ownership.global_radial_end) {
    return problem.local_coefficient_D(global_radial - ownership.global_radial_begin, theta, phi);
  }

  const std::size_t slot = theta * ownership.global_phi_cells + phi;
  if (global_radial + 1u == ownership.global_radial_begin && halo.lower_received) {
    return halo.lower[slot];
  }
  if (global_radial == ownership.global_radial_end && halo.upper_received) {
    return halo.upper[slot];
  }
  return std::numeric_limits<double>::quiet_NaN();
}

[[nodiscard]] DistributedLocalCsrMatrix BuildLocalCsr(
    RowAccumulator& rows,
    const DistributedDiffusionRowOwnership& ownership) {
  DistributedLocalCsrMatrix matrix;
  matrix.global_row_count = ownership.global_row_count;
  matrix.local_row_begin = ownership.local_row_begin;
  matrix.local_row_end = ownership.local_row_end;
  matrix.row_offsets.reserve(rows.size() + 1u);
  matrix.row_offsets.push_back(0u);
  for (auto& row : rows) {
    std::sort(row.begin(), row.end(), [](const auto& lhs, const auto& rhs) {
      return lhs.first < rhs.first;
    });
    for (const auto& entry : row) {
      if (entry.second != 0.0) {
        matrix.column_indices.push_back(entry.first);
        matrix.values.push_back(entry.second);
      }
    }
    matrix.row_offsets.push_back(matrix.values.size());
  }
  return matrix;
}

[[nodiscard]] std::string BuildDistributedAssemblyReport(
    const DistributedGenericDiffusionAssemblyResult& result) {
  const auto& ownership = result.ownership;
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.distributed_assembly"
      << "; implementation_id=p2.diffusion.distributed_radial_slab_assembly_v1"
      << "; equation_form=A_dT_dt_div_D_grad_T_plus_C_T_plus_B"
      << "; unit_system=cgs"
      << "; unknown=generic_scalar"
      << "; direction_splitting=false"
      << "; interface_diffusion_mean=arithmetic"
      << "; pole_metric_mode=" << ToString(result.pole_metric_mode)
      << "; row_ownership=distributed_radial_slab"
      << "; mpi_rank=" << ownership.rank
      << "; mpi_rank_count=" << ownership.rank_count
      << "; global_radial_cells=" << ownership.global_radial_cells
      << "; global_theta_cells=" << ownership.global_theta_cells
      << "; global_phi_cells=" << ownership.global_phi_cells
      << "; global_row_count=" << ownership.global_row_count
      << "; local_global_radial_begin=" << ownership.global_radial_begin
      << "; local_global_radial_end=" << ownership.global_radial_end
      << "; local_row_begin=" << ownership.local_row_begin
      << "; local_row_end=" << ownership.local_row_end
      << "; local_row_count=" << ownership.local_row_count
      << "; local_nonzero_count=" << result.local_nonzero_count
      << "; global_nonzero_count=" << result.global_nonzero_count
      << "; local_off_rank_column_count=" << result.local_off_rank_column_count
      << "; global_off_rank_column_count=" << result.global_off_rank_column_count
      << "; local_radial_seam_coupling_count=" << result.local_radial_seam_coupling_count
      << "; global_radial_seam_coupling_count=" << result.global_radial_seam_coupling_count
      << "; coefficient_D_halo_lower_received="
      << (result.coefficient_D_halo_lower_received ? "true" : "false")
      << "; coefficient_D_halo_upper_received="
      << (result.coefficient_D_halo_upper_received ? "true" : "false")
      << "; origin_remap_used=" << (result.global_origin_remap_used ? "true" : "false")
      << "; pole_remap_used=" << (result.global_pole_remap_used ? "true" : "false")
      << "; origin_remap_coupling_mode="
      << (result.global_origin_remap_used ? "metric_regular_origin" : "not_used")
      << "; pole_remap_coupling_mode="
      << (result.global_pole_remap_used ? "metric_regular_axis" : "not_used")
      << "; singular_boundary_face_conductance_used=false"
      << "; canonical_state_mutated=false";
  return out.str();
}

struct SolutionHalo {
  bool lower_received{false};
  bool upper_received{false};
  std::vector<double> lower;
  std::vector<double> upper;
};

[[nodiscard]] SolutionHalo ExchangeSolutionHalo(
    const DistributedDiffusionRowOwnership& ownership,
    const std::vector<double>& local_solution) {
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  SolutionHalo halo;
  halo.lower.assign(plane_size, 0.0);
  halo.upper.assign(plane_size, 0.0);

  std::vector<double> send_lower(plane_size, 0.0);
  std::vector<double> send_upper(plane_size, 0.0);
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
    for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
      const std::size_t slot = t * ownership.global_phi_cells + p;
      send_lower[slot] = local_solution[LocalLinearIndex(ownership, 0u, t, p)];
      send_upper[slot] =
          local_solution[LocalLinearIndex(ownership, local_radial_count - 1u, t, p)];
    }
  }

  MPI_Status status{};
  const int lower_rank = ownership.rank - 1;
  const int upper_rank = ownership.rank + 1;
  if (ownership.global_radial_begin > 0u) {
    MPI_Sendrecv(
        send_lower.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        lower_rank,
        1201,
        halo.lower.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        lower_rank,
        1202,
        ownership.communicator,
        &status);
    halo.lower_received = true;
  }
  if (ownership.global_radial_end < ownership.global_radial_cells) {
    MPI_Sendrecv(
        send_upper.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        upper_rank,
        1202,
        halo.upper.data(),
        static_cast<int>(plane_size),
        MPI_DOUBLE,
        upper_rank,
        1201,
        ownership.communicator,
        &status);
    halo.upper_received = true;
  }
  return halo;
}

[[nodiscard]] double SolutionAtGlobalColumn(
    const DistributedDiffusionRowOwnership& ownership,
    const std::vector<double>& local_solution,
    const SolutionHalo& halo,
    std::size_t global_column) noexcept {
  if (OwnsGlobalRow(ownership, global_column)) {
    return local_solution[global_column - ownership.local_row_begin];
  }
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  const std::size_t global_radial = global_column / plane_size;
  const std::size_t plane_slot = global_column % plane_size;
  if (global_radial + 1u == ownership.global_radial_begin && halo.lower_received) {
    return halo.lower[plane_slot];
  }
  if (global_radial == ownership.global_radial_end && halo.upper_received) {
    return halo.upper[plane_slot];
  }
  return std::numeric_limits<double>::quiet_NaN();
}

void ComputeDistributedResiduals(
    DistributedGenericDiffusionHypreSolveResult& result,
    const DistributedGenericDiffusionAssemblyResult& assembly) {
  const auto halo = ExchangeSolutionHalo(assembly.ownership, result.local_scalar_new);
  double local_l2_sq = 0.0;
  double local_linf = 0.0;
  double local_rhs_linf = 0.0;
  double local_delta = 0.0;

  for (std::size_t local_row = 0; local_row < assembly.ownership.local_row_count; ++local_row) {
    double ax = 0.0;
    for (std::size_t slot = assembly.local_matrix.row_offsets[local_row];
         slot < assembly.local_matrix.row_offsets[local_row + 1u];
         ++slot) {
      const double value = SolutionAtGlobalColumn(
          assembly.ownership,
          result.local_scalar_new,
          halo,
          assembly.local_matrix.column_indices[slot]);
      ax += assembly.local_matrix.values[slot] * value;
    }
    const double residual = ax - assembly.local_rhs[local_row];
    local_l2_sq += residual * residual;
    local_linf = std::max(local_linf, std::abs(residual));
    local_rhs_linf = std::max(local_rhs_linf, std::abs(assembly.local_rhs[local_row]));
    local_delta = std::max(
        local_delta,
        std::abs(result.local_scalar_new[local_row] - assembly.local_scalar_old_flat[local_row]));
  }

  double global_l2_sq = 0.0;
  MPI_Allreduce(&local_l2_sq, &global_l2_sq, 1, MPI_DOUBLE, MPI_SUM, assembly.ownership.communicator);
  MPI_Allreduce(&local_linf, &result.global_residual_linf, 1, MPI_DOUBLE, MPI_MAX, assembly.ownership.communicator);
  MPI_Allreduce(&local_rhs_linf, &result.global_rhs_linf, 1, MPI_DOUBLE, MPI_MAX, assembly.ownership.communicator);
  MPI_Allreduce(&local_delta, &result.max_abs_delta_global, 1, MPI_DOUBLE, MPI_MAX, assembly.ownership.communicator);
  result.local_residual_l2 = std::sqrt(local_l2_sq);
  result.global_residual_l2 = std::sqrt(global_l2_sq);
  result.local_residual_linf = local_linf;
  result.local_rhs_linf = local_rhs_linf;
  result.max_abs_delta_local = local_delta;
  result.global_residual_linf_relative =
      result.global_residual_linf / std::max(1.0, result.global_rhs_linf);
}

[[nodiscard]] double LocalMatrixRelativeChange(
    const DistributedGenericDiffusionAssemblyResult& assembly,
    const DistributedLaggedBoomerAmgCacheEntry& cache) noexcept {
  if (!cache.valid ||
      cache.row_offsets != assembly.local_matrix.row_offsets ||
      cache.column_indices != assembly.local_matrix.column_indices ||
      cache.values.size() != assembly.local_matrix.values.size()) {
    return std::numeric_limits<double>::infinity();
  }
  double max_relative_change = 0.0;
  for (std::size_t i = 0; i < assembly.local_matrix.values.size(); ++i) {
    const double old_abs = std::abs(cache.values[i]);
    const double new_abs = std::abs(assembly.local_matrix.values[i]);
    const double scale = std::max({old_abs, new_abs, 1.0e-300});
    const double relative_change = std::abs(assembly.local_matrix.values[i] - cache.values[i]) /
                                   scale;
    max_relative_change = std::max(max_relative_change, relative_change);
  }
  return max_relative_change;
}

void StoreLaggedAmgSetupSnapshot(
    DistributedLaggedBoomerAmgCacheEntry& cache,
    const DistributedGenericDiffusionAssemblyResult& assembly,
    int iterations) {
  cache.valid = true;
  cache.last_iterations = iterations;
  cache.reuse_count_since_setup = 0;
  cache.row_offsets = assembly.local_matrix.row_offsets;
  cache.column_indices = assembly.local_matrix.column_indices;
  cache.values = assembly.local_matrix.values;
}

void RecordLaggedAmgReuse(
    DistributedLaggedBoomerAmgCacheEntry& cache,
    int iterations) noexcept {
  cache.last_iterations = iterations;
  ++cache.reuse_count_since_setup;
}

[[nodiscard]] std::string BuildDistributedSolveReport(
    const DistributedGenericDiffusionHypreSolveResult& result,
    const GenericDiffusionHypreSolveOptions& options) {
  const auto& ownership = result.ownership;
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.hypre_distributed_solve"
      << "; implementation_id=p2.diffusion.hypre_distributed_parcsr_gmres_boomeramg_v1"
      << "; backend=" << kHypreBackend
      << "; hypre_enabled=true"
      << "; row_ownership=distributed_radial_slab"
      << "; matrix_source=p2.diffusion.distributed_radial_slab_assembly_v1"
      << "; matrix_format=hypre_ij_parcsr"
      << "; direction_splitting=false"
      << "; mpi_rank=" << ownership.rank
      << "; mpi_rank_count=" << ownership.rank_count
      << "; global_row_count=" << ownership.global_row_count
      << "; local_row_begin=" << ownership.local_row_begin
      << "; local_row_end=" << ownership.local_row_end
      << "; local_row_count=" << ownership.local_row_count
      << "; local_nonzero_count=" << result.local_nonzero_count
      << "; global_nonzero_count=" << result.global_nonzero_count
      << "; local_off_rank_column_count=" << result.local_off_rank_column_count
      << "; global_off_rank_column_count=" << result.global_off_rank_column_count
      << "; solver=ParCSRGMRES"
      << "; preconditioner=BoomerAMG"
      << "; boomeramg_parameters=preconditioner_max_iter_1_tol_0"
      << "; boomeramg_max_iterations=" << kBoomerAmgPreconditionerMaxIterations
      << "; boomeramg_tolerance=" << kBoomerAmgPreconditionerTolerance
      << "; boomeramg_print_level=" << kBoomerAmgPreconditionerPrintLevel
      << "; boomeramg_strong_threshold=" << kBoomerAmgStrongThreshold
      << "; lagged_amg_enabled=" << (result.lagged_amg_enabled ? "true" : "false")
      << "; lagged_amg_candidate=" << (result.lagged_amg_candidate ? "true" : "false")
      << "; lagged_amg_reuse_attempted="
      << (result.lagged_amg_reuse_attempted ? "true" : "false")
      << "; lagged_amg_reuse_accepted="
      << (result.lagged_amg_reuse_accepted ? "true" : "false")
      << "; lagged_amg_rebuild_used="
      << (result.lagged_amg_rebuild_used ? "true" : "false")
      << "; lagged_amg_fallback_rebuild_used="
      << (result.lagged_amg_fallback_rebuild_used ? "true" : "false")
      << "; lagged_amg_local_matrix_rel_change="
      << result.lagged_amg_local_matrix_rel_change
      << "; lagged_amg_global_matrix_rel_change="
      << result.lagged_amg_global_matrix_rel_change
      << "; lagged_amg_rebuild_every=" << result.lagged_amg_rebuild_every
      << "; lagged_amg_cached_reuse_count="
      << result.lagged_amg_cached_reuse_count
      << "; lagged_amg_cached_iterations=" << result.lagged_amg_cached_iterations
      << "; lagged_amg_reuse_iterations=" << result.lagged_amg_reuse_iterations
      << "; lagged_amg_reuse_count_after="
      << result.lagged_amg_reuse_count_after
      << "; lagged_amg_reuse_final_relative_residual="
      << result.lagged_amg_reuse_final_relative_residual
      << "; lagged_amg_reuse_setup_wall_s=" << result.lagged_amg_reuse_setup_wall_s
      << "; lagged_amg_reuse_solve_wall_s=" << result.lagged_amg_reuse_solve_wall_s
      << "; gmres_relative_tolerance=" << options.relative_tolerance
      << "; gmres_max_iterations=" << options.max_iterations
      << "; gmres_krylov_dimension=" << options.krylov_dimension
      << "; gmres_iterations=" << result.gmres_iterations
      << "; gmres_final_relative_residual=" << result.gmres_final_relative_residual
      << "; hypre_setup_wall_s=" << result.hypre_setup_wall_s
      << "; hypre_solve_wall_s=" << result.hypre_solve_wall_s
      << "; local_residual_l2=" << result.local_residual_l2
      << "; global_residual_l2=" << result.global_residual_l2
      << "; local_residual_linf=" << result.local_residual_linf
      << "; global_residual_linf=" << result.global_residual_linf
      << "; local_rhs_linf=" << result.local_rhs_linf
      << "; global_rhs_linf=" << result.global_rhs_linf
      << "; global_residual_linf_relative=" << result.global_residual_linf_relative
      << "; max_abs_delta_local=" << result.max_abs_delta_local
      << "; max_abs_delta_global=" << result.max_abs_delta_global
      << "; canonical_state_mutated=false"
      << "; fallback_used=false"
      << "; gathered_solve_used=false";
  return out.str();
}

}  // namespace

DistributedLaggedBoomerAmgCache::~DistributedLaggedBoomerAmgCache() {
  for (auto& entry : entries) {
    if (entry.boomeramg != nullptr) {
      HYPRE_BoomerAMGDestroy(static_cast<HYPRE_Solver>(entry.boomeramg));
      entry.boomeramg = nullptr;
    }
  }
}

bool DistributedLocalCsrMatrix::is_shape_complete() const noexcept {
  return local_row_end >= local_row_begin &&
         row_offsets.size() == (local_row_end - local_row_begin) + 1u &&
         column_indices.size() == values.size() &&
         (row_offsets.empty() || row_offsets.back() == values.size());
}

DistributedDiffusionRowOwnership BuildDistributedDiffusionRowOwnership(
    MPI_Comm communicator,
    std::size_t global_radial_cells,
    std::size_t global_theta_cells,
    std::size_t global_phi_cells) noexcept {
  DistributedDiffusionRowOwnership ownership;
  ownership.communicator = communicator;
  ownership.global_radial_cells = global_radial_cells;
  ownership.global_theta_cells = global_theta_cells;
  ownership.global_phi_cells = global_phi_cells;

  int rank = 0;
  int rank_count = 1;
  if (MPI_Comm_rank(communicator, &rank) != MPI_SUCCESS ||
      MPI_Comm_size(communicator, &rank_count) != MPI_SUCCESS) {
    ownership.failure_reason = "MPI communicator query failed";
    ownership.failure_diagnostics = DistributedFailureLine(
        ownership.failure_reason.c_str(),
        rank,
        rank_count,
        "MPI_Comm_rank_or_size");
    return ownership;
  }
  ownership.rank = rank;
  ownership.rank_count = rank_count;

  const auto global_layout_validation = ValidateCollectiveGlobalLayout(ownership);
  if (!global_layout_validation.success) {
    ownership.failure_reason = global_layout_validation.failure_reason;
    ownership.failure_diagnostics = global_layout_validation.failure_diagnostics;
    return ownership;
  }

  const auto radial = dec3d::mesh::BuildRadialOwnership(
      global_radial_cells,
      static_cast<std::size_t>(rank_count));
  if (!radial.is_valid()) {
    ownership.failure_reason = radial.failure_reason;
    ownership.failure_diagnostics =
        DistributedFailureLine(ownership.failure_reason.c_str(), rank, rank_count);
    return ownership;
  }

  const auto& slice = radial.slices[static_cast<std::size_t>(rank)];
  ownership.global_row_count = global_radial_cells * global_theta_cells * global_phi_cells;
  ownership.global_radial_begin = slice.begin_index;
  ownership.global_radial_end = slice.end_index;
  ownership.local_row_begin = slice.begin_index * global_theta_cells * global_phi_cells;
  ownership.local_row_end = slice.end_index * global_theta_cells * global_phi_cells;
  ownership.local_row_count = ownership.local_row_end - ownership.local_row_begin;

  const auto collective = ValidateCollectiveRowOwnership(ownership);
  if (!collective.success) {
    ownership.failure_reason = collective.failure_reason;
    ownership.failure_diagnostics = collective.failure_diagnostics;
    return ownership;
  }

  ownership.success = true;
  ownership.report_line = BuildOwnershipReport(ownership);
  return ownership;
}

DistributedGenericDiffusionAssemblyResult AssembleDistributedGenericDiffusionSystem(
    const DistributedGenericDiffusionProblem& problem) noexcept {
  if (!problem.ownership.success) {
    return AssemblyFailure("distributed row ownership is invalid", problem.ownership);
  }
  const auto& ownership = problem.ownership;
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  if (!LocalArrayShapeMatches(problem.local_coefficient_A, ownership) ||
      !LocalArrayShapeMatches(problem.local_coefficient_D, ownership) ||
      !LocalArrayShapeMatches(problem.local_coefficient_C, ownership) ||
      !LocalArrayShapeMatches(problem.local_coefficient_B, ownership) ||
      !LocalArrayShapeMatches(problem.local_scalar_old, ownership)) {
    return AssemblyFailure("local arrays do not match distributed row ownership", ownership);
  }
  if (!LocalFaceEffectiveShapeMatches(problem.face_effective_coefficients, ownership)) {
    return AssemblyFailure("local face-effective arrays do not match distributed row ownership", ownership);
  }
  if (!std::isfinite(problem.dt_s) || problem.dt_s <= 0.0) {
    return AssemblyFailure("dt_s must be positive for distributed diffusion assembly", ownership);
  }
  if (!problem.global_geometry.valid ||
      problem.global_geometry.radial_faces.size() != ownership.global_radial_cells + 1u ||
      problem.global_geometry.theta_faces.size() != ownership.global_theta_cells + 1u ||
      problem.global_geometry.phi_faces.size() != ownership.global_phi_cells + 1u ||
      problem.global_geometry.cell_volumes.size() != ownership.global_row_count) {
    return AssemblyFailure("global geometry does not match distributed row ownership", ownership);
  }
  const auto geometry_metrics = BuildDistributedGeometryMetrics(problem);

  DistributedGenericDiffusionAssemblyResult result;
  result.success = true;
  result.ownership = ownership;
  result.singular_boundary_face_conductance_used = false;
  RowAccumulator rows(ownership.local_row_count);
  for (auto& row : rows) {
    row.reserve(kExpectedDistributedStencilEntriesPerRow);
  }
  result.local_rhs.assign(ownership.local_row_count, 0.0);
  result.local_scalar_old_flat.assign(ownership.local_row_count, 0.0);

  const auto halo = ExchangeCoefficientDHalo(problem);
  result.coefficient_D_halo_lower_received = halo.lower_received;
  result.coefficient_D_halo_upper_received = halo.upper_received;

  for (std::size_t lr = 0; lr < local_radial_count; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const std::size_t local_row = LocalLinearIndex(ownership, lr, t, p);
        const std::size_t global_row = GlobalRow(ownership, gr, t, p);
        const double volume = problem.global_geometry.cell_volumes[global_row];
        if (!std::isfinite(volume) || volume <= 0.0) {
          return AssemblyFailure("distributed geometry volume is not positive", ownership);
        }
        const double time_coefficient = volume * problem.local_coefficient_A(lr, t, p) / problem.dt_s;
        const double c_coefficient = volume * problem.local_coefficient_C(lr, t, p);
        const double b_source = volume * problem.local_coefficient_B(lr, t, p);
        if (!std::isfinite(time_coefficient) || !std::isfinite(c_coefficient) ||
            !std::isfinite(b_source) || !std::isfinite(problem.local_scalar_old(lr, t, p))) {
          return AssemblyFailure("distributed coefficient is non-finite", ownership);
        }
        AddMatrixEntry(rows, local_row, global_row, time_coefficient - c_coefficient);
        result.local_rhs[local_row] = time_coefficient * problem.local_scalar_old(lr, t, p) + b_source;
        result.local_scalar_old_flat[local_row] = problem.local_scalar_old(lr, t, p);
      }
    }
  }

  if (problem.face_effective_coefficients.enabled) {
    for (double value : problem.face_effective_coefficients.radial_face_D.storage()) {
      if (!std::isfinite(value) || value < 0.0) {
        return AssemblyFailure("distributed radial face-effective D must be non-negative", ownership);
      }
    }
    for (double value : problem.face_effective_coefficients.theta_face_D.storage()) {
      if (!std::isfinite(value) || value < 0.0) {
        return AssemblyFailure("distributed theta face-effective D must be non-negative", ownership);
      }
    }
    for (double value : problem.face_effective_coefficients.phi_face_D.storage()) {
      if (!std::isfinite(value) || value < 0.0) {
        return AssemblyFailure("distributed phi face-effective D must be non-negative", ownership);
      }
    }
  }

  const std::size_t radial_pair_begin =
      ownership.global_radial_begin == 0u ? 0u : ownership.global_radial_begin - 1u;
  const std::size_t radial_pair_end = std::min(
      ownership.global_radial_end,
      ownership.global_radial_cells == 0u ? 0u : ownership.global_radial_cells - 1u);
  for (std::size_t gr = radial_pair_begin; gr < radial_pair_end; ++gr) {
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const double left_d = CoefficientDAt(problem, halo, gr, t, p);
        const double right_d = CoefficientDAt(problem, halo, gr + 1u, t, p);
        if (!std::isfinite(left_d) || !std::isfinite(right_d)) {
          return AssemblyFailure("distributed radial D halo is incomplete", ownership);
        }
        const double area =
            geometry_metrics.radial_face_area[
                RadialFaceMetricIndex(ownership, gr + 1u, t, p)];
        const double distance = geometry_metrics.radial_center_distance[gr];
        const double conductance = area *
            RadialInterfaceD(problem, gr + 1u, t, p, left_d, right_d) /
            distance;
        AddOwnedPairConductance(
            rows,
            ownership,
            GlobalRow(ownership, gr, t, p),
            GlobalRow(ownership, gr + 1u, t, p),
            conductance,
            true,
            result);
      }
    }
  }

  for (std::size_t lr = 0; lr < local_radial_count; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;

    for (std::size_t theta = 0; theta + 1u < ownership.global_theta_cells; ++theta) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const double area =
            geometry_metrics.theta_face_area[
                ThetaFaceMetricIndex(ownership, gr, theta + 1u, p)];
        const double radius = geometry_metrics.radial_center[gr];
        const double distance =
            radius * (geometry_metrics.theta_center[theta + 1u] -
                      geometry_metrics.theta_center[theta]);
        const double conductance = area * ThetaInterfaceD(problem, lr, theta + 1u, p) /
            distance;
        AddOwnedPairConductance(
            rows,
            ownership,
            GlobalRow(ownership, gr, theta, p),
            GlobalRow(ownership, gr, theta + 1u, p),
            conductance,
            false,
            result);
      }
    }

    if (ownership.global_phi_cells > 1u) {
      for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
        for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
          if (ownership.global_phi_cells == 2u && p == 1u) {
            continue;
          }
          const std::size_t next_phi = (p + 1u) % ownership.global_phi_cells;
          const double radius = geometry_metrics.radial_center[gr];
          const double sin_theta = geometry_metrics.phi_metric_sin_theta[t];
          const double distance = radius * sin_theta * geometry_metrics.dphi;
          const double area =
              geometry_metrics.phi_face_area[PhiFaceMetricIndex(ownership, gr, t)];
          const double conductance = area * PhiInterfaceD(problem, lr, t, p) / distance;
          AddOwnedPairConductance(
              rows,
              ownership,
              GlobalRow(ownership, gr, t, p),
              GlobalRow(ownership, gr, t, next_phi),
              conductance,
              false,
              result);
        }
      }
    }
  }

  const auto remap_layout = dec3d::mesh::ScalarRemapLayout{
      ownership.global_radial_cells,
      ownership.global_theta_cells,
      ownership.global_phi_cells};

  if (problem.boundary_policy.inner_radial == DiffusionBoundaryKind::scalar_origin_remap_required &&
      ownership.global_radial_begin == 0u &&
      std::abs(problem.global_geometry.radial_faces.front()) <= kBoundaryTolerance) {
    const auto validation = dec3d::mesh::ValidateScalarHalfTurnTopology(remap_layout);
    if (!validation.success) {
      return AssemblyFailure(validation.failure_reason.c_str(), ownership);
    }
    result.local_origin_remap_used = true;
  }

  const bool lower_pole_remap_requested =
      problem.boundary_policy.theta_lower == DiffusionBoundaryKind::scalar_pole_remap_required &&
      std::abs(problem.global_geometry.theta_faces.front()) <= kBoundaryTolerance;
  const bool upper_pole_remap_requested =
      problem.boundary_policy.theta_upper == DiffusionBoundaryKind::scalar_pole_remap_required &&
      std::abs(problem.global_geometry.theta_faces.back() - kPi) <= kBoundaryTolerance;
  if (lower_pole_remap_requested || upper_pole_remap_requested) {
    const auto validation = dec3d::mesh::ValidateScalarHalfTurnTopology(remap_layout);
    if (!validation.success) {
      return AssemblyFailure(validation.failure_reason.c_str(), ownership);
    }
    for (std::size_t lr = 0; lr < local_radial_count; ++lr) {
      const std::size_t gr = ownership.global_radial_begin + lr;
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const auto lower = dec3d::mesh::MapScalarLowerPoleNeighbor(1u, gr, p, remap_layout);
        if (lower_pole_remap_requested && lower.success) {
          result.local_pole_remap_used = true;
        }
        const auto upper = dec3d::mesh::MapScalarUpperPoleNeighbor(1u, gr, p, remap_layout);
        if (upper_pole_remap_requested && upper.success) {
          result.local_pole_remap_used = true;
        }
      }
    }
  }

  if (problem.boundary_policy.outer_radial == DiffusionBoundaryKind::radiation_marshak_vacuum &&
      ownership.global_radial_end == ownership.global_radial_cells) {
    const std::size_t global_radial = ownership.global_radial_cells - 1u;
    const std::size_t local_radial = local_radial_count - 1u;
    const double r_center = geometry_metrics.radial_center[global_radial];
    const double distance = problem.global_geometry.radial_faces.back() - r_center;
    if (!std::isfinite(distance) || distance <= 0.0) {
      return AssemblyFailure("distributed Marshak center-to-boundary distance must be positive", ownership);
    }
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const double d_face =
            problem.face_effective_coefficients.enabled
                ? problem.face_effective_coefficients.radial_face_D(
                      local_radial_count, t, p)
                : problem.local_coefficient_D(local_radial, t, p);
        const double g_face = MarshakGFace(d_face, distance);
        const double area =
            geometry_metrics.radial_face_area[
                RadialFaceMetricIndex(ownership, ownership.global_radial_cells, t, p)];
        const double loss = area * g_face;
        if (!std::isfinite(d_face) || !std::isfinite(g_face) ||
            !std::isfinite(area) || !std::isfinite(loss) ||
            d_face < 0.0 || g_face < 0.0 || area < 0.0 || loss < 0.0) {
          return AssemblyFailure("distributed Marshak boundary conductance is not physical", ownership);
        }
        if (loss > 0.0) {
          const std::size_t local_row = LocalLinearIndex(ownership, local_radial, t, p);
          const std::size_t global_row = GlobalRow(ownership, global_radial, t, p);
          AddMatrixEntry(rows, local_row, global_row, loss);
        }
      }
    }
  }

  result.local_matrix = BuildLocalCsr(rows, ownership);
  result.pole_metric_mode = problem.pole_metric_mode;
  result.local_nonzero_count = result.local_matrix.values.size();
  unsigned long long local_counts[3] = {
      static_cast<unsigned long long>(result.local_nonzero_count),
      static_cast<unsigned long long>(result.local_off_rank_column_count),
      static_cast<unsigned long long>(result.local_radial_seam_coupling_count)};
  unsigned long long global_counts[3] = {0ull, 0ull, 0ull};
  MPI_Allreduce(
      local_counts,
      global_counts,
      3,
      MPI_UNSIGNED_LONG_LONG,
      MPI_SUM,
      ownership.communicator);
  result.global_nonzero_count = static_cast<std::size_t>(global_counts[0]);
  result.global_off_rank_column_count = static_cast<std::size_t>(global_counts[1]);
  result.global_radial_seam_coupling_count = static_cast<std::size_t>(global_counts[2]);
  int local_origin_used = result.local_origin_remap_used ? 1 : 0;
  int local_pole_used = result.local_pole_remap_used ? 1 : 0;
  int global_origin_used = 0;
  int global_pole_used = 0;
  MPI_Allreduce(&local_origin_used, &global_origin_used, 1, MPI_INT, MPI_MAX, ownership.communicator);
  MPI_Allreduce(&local_pole_used, &global_pole_used, 1, MPI_INT, MPI_MAX, ownership.communicator);
  result.global_origin_remap_used = global_origin_used != 0;
  result.global_pole_remap_used = global_pole_used != 0;
  result.report_line = BuildDistributedAssemblyReport(result);
  return result;
}

DistributedGenericDiffusionHypreSolveResult SolveDistributedGenericDiffusionHypre(
    const DistributedGenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionHypreSolveOptions& options) noexcept {
  return SolveDistributedGenericDiffusionHypre(
      assembly,
      options,
      nullptr,
      DistributedLaggedBoomerAmgSolveOptions{});
}

DistributedGenericDiffusionHypreSolveResult SolveDistributedGenericDiffusionHypre(
    const DistributedGenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionHypreSolveOptions& options,
    DistributedLaggedBoomerAmgCache* lagged_cache,
    const DistributedLaggedBoomerAmgSolveOptions& lagged_options) noexcept {
  if (!assembly.success || !assembly.is_complete()) {
    return SolveFailure("distributed assembly is incomplete or failed", assembly.ownership);
  }
  if (!assembly.local_matrix.is_shape_complete()) {
    return SolveFailure("distributed local matrix shape is incomplete", assembly.ownership);
  }
  if (assembly.local_rhs.size() != assembly.ownership.local_row_count ||
      assembly.local_scalar_old_flat.size() != assembly.ownership.local_row_count) {
    return SolveFailure("distributed RHS or scalar_old size is invalid", assembly.ownership);
  }
  if (!OptionsAreValid(options)) {
    return SolveFailure("HYPRE solver options are invalid", assembly.ownership);
  }
  if (!FitsHypreInt(assembly.ownership.local_row_count) ||
      !FitsHypreInt(assembly.ownership.global_row_count)) {
    return SolveFailure("row count exceeds HYPRE_Int range", assembly.ownership);
  }

  int status = EnsureHypreInitialized();
  if (status != 0) {
    return SolveFailure("HYPRE library initialization failed", assembly.ownership, "HYPRE_Initialize");
  }

  HypreDistributedObjects objects;
  HYPRE_ParCSRMatrix parcsr_matrix = nullptr;
  HYPRE_ParVector par_rhs = nullptr;
  HYPRE_ParVector par_solution = nullptr;
  const auto row_begin = static_cast<HYPRE_BigInt>(assembly.ownership.local_row_begin);
  const auto row_end = static_cast<HYPRE_BigInt>(assembly.ownership.local_row_end - 1u);

  status = HYPRE_IJMatrixCreate(
      assembly.ownership.communicator,
      row_begin,
      row_end,
      row_begin,
      row_end,
      &objects.ij_matrix);
  if (status != 0) {
    return SolveFailure("HYPRE distributed matrix creation failed", assembly.ownership, "IJMatrixCreate");
  }
  status = HYPRE_IJMatrixSetObjectType(objects.ij_matrix, HYPRE_PARCSR);
  if (status != 0) {
    return SolveFailure("HYPRE matrix object type failed", assembly.ownership, "IJMatrixSetObjectType");
  }
  status = HYPRE_IJMatrixInitialize(objects.ij_matrix);
  if (status != 0) {
    return SolveFailure("HYPRE matrix initialize failed", assembly.ownership, "IJMatrixInitialize");
  }

  const HYPRE_Int local_row_count_int =
      static_cast<HYPRE_Int>(assembly.ownership.local_row_count);
  std::vector<HYPRE_BigInt> rows(assembly.ownership.local_row_count);
  std::vector<HYPRE_Int> row_nonzero_counts(assembly.ownership.local_row_count);
  std::vector<HYPRE_BigInt> matrix_columns;
  std::vector<HYPRE_Complex> matrix_values;
  matrix_columns.reserve(assembly.local_matrix.column_indices.size());
  matrix_values.reserve(assembly.local_matrix.values.size());
  for (std::size_t local_row = 0; local_row < assembly.ownership.local_row_count; ++local_row) {
    const auto begin = assembly.local_matrix.row_offsets[local_row];
    const auto end = assembly.local_matrix.row_offsets[local_row + 1u];
    if (!FitsHypreInt(end - begin)) {
      return SolveFailure("row nonzero count exceeds HYPRE_Int range", assembly.ownership);
    }
    rows[local_row] =
        static_cast<HYPRE_BigInt>(assembly.ownership.local_row_begin + local_row);
    row_nonzero_counts[local_row] = static_cast<HYPRE_Int>(end - begin);
    for (std::size_t slot = begin; slot < end; ++slot) {
      matrix_columns.push_back(
          static_cast<HYPRE_BigInt>(assembly.local_matrix.column_indices[slot]));
      matrix_values.push_back(static_cast<HYPRE_Complex>(assembly.local_matrix.values[slot]));
    }
  }
  status = HYPRE_IJMatrixSetValues(
      objects.ij_matrix,
      local_row_count_int,
      row_nonzero_counts.data(),
      rows.data(),
      matrix_columns.data(),
      matrix_values.data());
  if (status != 0) {
    return SolveFailure("HYPRE distributed matrix row insertion failed",
                        assembly.ownership,
                        "IJMatrixSetValues");
  }
  status = HYPRE_IJMatrixAssemble(objects.ij_matrix);
  if (status != 0) {
    return SolveFailure("HYPRE distributed matrix assemble failed", assembly.ownership, "IJMatrixAssemble");
  }
  status = HYPRE_IJMatrixGetObject(objects.ij_matrix, reinterpret_cast<void**>(&parcsr_matrix));
  if (status != 0) {
    return SolveFailure("HYPRE distributed matrix object extraction failed", assembly.ownership, "IJMatrixGetObject");
  }

  std::vector<HYPRE_Complex> rhs_values(assembly.ownership.local_row_count);
  std::vector<HYPRE_Complex> initial_values(assembly.ownership.local_row_count);
  for (std::size_t i = 0; i < assembly.ownership.local_row_count; ++i) {
    rhs_values[i] = static_cast<HYPRE_Complex>(assembly.local_rhs[i]);
    initial_values[i] = static_cast<HYPRE_Complex>(assembly.local_scalar_old_flat[i]);
  }

  status = HYPRE_IJVectorCreate(assembly.ownership.communicator, row_begin, row_end, &objects.ij_rhs);
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS vector creation failed", assembly.ownership, "IJVectorCreate");
  }
  status = HYPRE_IJVectorSetObjectType(objects.ij_rhs, HYPRE_PARCSR);
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS object type failed", assembly.ownership, "IJVectorSetObjectType");
  }
  status = HYPRE_IJVectorInitialize(objects.ij_rhs);
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS initialize failed", assembly.ownership, "IJVectorInitialize");
  }
  status = HYPRE_IJVectorSetValues(objects.ij_rhs, local_row_count_int, rows.data(), rhs_values.data());
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS insertion failed", assembly.ownership, "IJVectorSetValues");
  }
  status = HYPRE_IJVectorAssemble(objects.ij_rhs);
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS assemble failed", assembly.ownership, "IJVectorAssemble");
  }
  status = HYPRE_IJVectorGetObject(objects.ij_rhs, reinterpret_cast<void**>(&par_rhs));
  if (status != 0) {
    return SolveFailure("HYPRE distributed RHS object extraction failed", assembly.ownership, "IJVectorGetObject");
  }

  status = HYPRE_IJVectorCreate(assembly.ownership.communicator, row_begin, row_end, &objects.ij_solution);
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution vector creation failed", assembly.ownership, "IJVectorCreate");
  }
  status = HYPRE_IJVectorSetObjectType(objects.ij_solution, HYPRE_PARCSR);
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution object type failed", assembly.ownership, "IJVectorSetObjectType");
  }
  status = HYPRE_IJVectorInitialize(objects.ij_solution);
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution initialize failed", assembly.ownership, "IJVectorInitialize");
  }
  status = HYPRE_IJVectorSetValues(
      objects.ij_solution,
      local_row_count_int,
      rows.data(),
      initial_values.data());
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution insertion failed", assembly.ownership, "IJVectorSetValues");
  }
  status = HYPRE_IJVectorAssemble(objects.ij_solution);
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution assemble failed", assembly.ownership, "IJVectorAssemble");
  }
  status = HYPRE_IJVectorGetObject(objects.ij_solution, reinterpret_cast<void**>(&par_solution));
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution object extraction failed", assembly.ownership, "IJVectorGetObject");
  }

  DistributedGenericDiffusionHypreSolveResult result;
  result.backend = kHypreBackend;
  result.hypre_enabled = true;
  result.ownership = assembly.ownership;
  result.local_nonzero_count = assembly.local_nonzero_count;
  result.global_nonzero_count = assembly.global_nonzero_count;
  result.local_off_rank_column_count = assembly.local_off_rank_column_count;
  result.global_off_rank_column_count = assembly.global_off_rank_column_count;

  DistributedLaggedBoomerAmgCacheEntry* lagged_entry = nullptr;
  const bool lagged_enabled =
      lagged_cache != nullptr &&
      lagged_options.enabled &&
      lagged_options.rebuild_every > 1 &&
      std::isfinite(lagged_options.max_matrix_relative_change) &&
      lagged_options.max_matrix_relative_change >= 0.0 &&
      std::isfinite(lagged_options.max_iteration_growth) &&
      lagged_options.max_iteration_growth >= 1.0;
  if (lagged_enabled) {
    result.lagged_amg_enabled = true;
    result.lagged_amg_rebuild_every = lagged_options.rebuild_every;
    if (lagged_cache->entries.size() <= lagged_options.group_index) {
      lagged_cache->entries.resize(lagged_options.group_index + 1u);
    }
    lagged_entry = &lagged_cache->entries[lagged_options.group_index];
    result.lagged_amg_local_matrix_rel_change =
        LocalMatrixRelativeChange(assembly, *lagged_entry);
    MPI_Allreduce(&result.lagged_amg_local_matrix_rel_change,
                  &result.lagged_amg_global_matrix_rel_change,
                  1,
                  MPI_DOUBLE,
                  MPI_MAX,
                  assembly.ownership.communicator);
    result.lagged_amg_cached_reuse_count =
        lagged_entry->reuse_count_since_setup;
    result.lagged_amg_cached_iterations = lagged_entry->last_iterations;
    result.lagged_amg_candidate =
        lagged_entry->valid &&
        lagged_entry->boomeramg != nullptr &&
        lagged_entry->reuse_count_since_setup < lagged_options.rebuild_every - 1 &&
        result.lagged_amg_global_matrix_rel_change <=
            lagged_options.max_matrix_relative_change;
  }

  if (result.lagged_amg_candidate && lagged_entry != nullptr) {
    result.lagged_amg_reuse_attempted = true;
    HYPRE_IJVector lagged_ij_solution = nullptr;
    HYPRE_ParVector lagged_par_solution = nullptr;
    HYPRE_Solver lagged_gmres = nullptr;
    bool lagged_ok = true;
    status = HYPRE_IJVectorCreate(
        assembly.ownership.communicator,
        row_begin,
        row_end,
        &lagged_ij_solution);
    lagged_ok = status == 0;
    if (lagged_ok) {
      status = HYPRE_IJVectorSetObjectType(lagged_ij_solution, HYPRE_PARCSR);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_IJVectorInitialize(lagged_ij_solution);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_IJVectorSetValues(
          lagged_ij_solution,
          local_row_count_int,
          rows.data(),
          initial_values.data());
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_IJVectorAssemble(lagged_ij_solution);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_IJVectorGetObject(
          lagged_ij_solution,
          reinterpret_cast<void**>(&lagged_par_solution));
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESCreate(assembly.ownership.communicator, &lagged_gmres);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetTol(lagged_gmres, options.relative_tolerance);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetMaxIter(lagged_gmres, options.max_iterations);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetKDim(lagged_gmres, options.krylov_dimension);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetLogging(
          lagged_gmres,
          options.logging_enabled ? 1 : 0);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetPrintLevel(lagged_gmres, options.print_level);
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      status = HYPRE_ParCSRGMRESSetPrecond(
          lagged_gmres,
          reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSolve),
          NoOpParSolverSetup,
          static_cast<HYPRE_Solver>(lagged_entry->boomeramg));
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      const auto lagged_setup_timer = std::chrono::steady_clock::now();
      status = HYPRE_ParCSRGMRESSetup(
          lagged_gmres,
          parcsr_matrix,
          par_rhs,
          lagged_par_solution);
      result.lagged_amg_reuse_setup_wall_s = ElapsedSecondsSince(lagged_setup_timer);
      result.hypre_setup_wall_s += result.lagged_amg_reuse_setup_wall_s;
      lagged_ok = status == 0;
    }
    if (lagged_ok) {
      const auto lagged_solve_timer = std::chrono::steady_clock::now();
      status = HYPRE_ParCSRGMRESSolve(
          lagged_gmres,
          parcsr_matrix,
          par_rhs,
          lagged_par_solution);
      result.lagged_amg_reuse_solve_wall_s = ElapsedSecondsSince(lagged_solve_timer);
      result.hypre_solve_wall_s += result.lagged_amg_reuse_solve_wall_s;
      (void)HYPRE_ParCSRGMRESGetNumIterations(
          lagged_gmres,
          &result.lagged_amg_reuse_iterations);
      (void)HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(
          lagged_gmres,
          &result.lagged_amg_reuse_final_relative_residual);
      lagged_ok = status == 0;
    }
    if (lagged_gmres != nullptr) {
      HYPRE_ParCSRGMRESDestroy(lagged_gmres);
    }
    if (lagged_ok) {
      std::vector<HYPRE_Complex> lagged_solution_values(assembly.ownership.local_row_count);
      status = HYPRE_IJVectorGetValues(
          lagged_ij_solution,
          local_row_count_int,
          rows.data(),
          lagged_solution_values.data());
      lagged_ok = status == 0;
      result.local_scalar_new.resize(assembly.ownership.local_row_count);
      for (std::size_t i = 0; lagged_ok && i < assembly.ownership.local_row_count; ++i) {
        result.local_scalar_new[i] = static_cast<double>(lagged_solution_values[i]);
        lagged_ok = std::isfinite(result.local_scalar_new[i]);
      }
    }
    if (lagged_ij_solution != nullptr) {
      HYPRE_IJVectorDestroy(lagged_ij_solution);
    }
    int local_lagged_ok = lagged_ok ? 1 : 0;
    int global_lagged_ok = 0;
    MPI_Allreduce(&local_lagged_ok,
                  &global_lagged_ok,
                  1,
                  MPI_INT,
                  MPI_MIN,
                  assembly.ownership.communicator);
    lagged_ok = global_lagged_ok != 0;
    if (lagged_ok) {
      ComputeDistributedResiduals(result, assembly);
      lagged_ok =
          result.global_residual_linf_relative <= options.relative_tolerance * 100.0 &&
          result.lagged_amg_reuse_iterations <=
              static_cast<int>(std::ceil(lagged_options.max_iteration_growth *
                                         std::max(1, lagged_entry->last_iterations)));
    }
    local_lagged_ok = lagged_ok ? 1 : 0;
    MPI_Allreduce(&local_lagged_ok,
                  &global_lagged_ok,
                  1,
                  MPI_INT,
                  MPI_MIN,
                  assembly.ownership.communicator);
    if (global_lagged_ok != 0) {
      result.lagged_amg_reuse_accepted = true;
      result.gmres_iterations = result.lagged_amg_reuse_iterations;
      result.gmres_final_relative_residual =
          result.lagged_amg_reuse_final_relative_residual;
      RecordLaggedAmgReuse(*lagged_entry, result.gmres_iterations);
      result.lagged_amg_reuse_count_after =
          lagged_entry->reuse_count_since_setup;
      result.success = true;
      result.report_line = BuildDistributedSolveReport(result, options);
      return result;
    }
    result.local_scalar_new.clear();
    result.lagged_amg_fallback_rebuild_used = true;
    (void)HYPRE_ClearAllErrors();
  }

  status = HYPRE_ParCSRGMRESCreate(assembly.ownership.communicator, &objects.gmres);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES creation failed", assembly.ownership, "ParCSRGMRESCreate");
  }
  HYPRE_Solver active_boomeramg = nullptr;
  if (lagged_enabled && lagged_entry != nullptr) {
    if (lagged_entry->boomeramg != nullptr) {
      HYPRE_BoomerAMGDestroy(static_cast<HYPRE_Solver>(lagged_entry->boomeramg));
      lagged_entry->boomeramg = nullptr;
    }
    status = HYPRE_BoomerAMGCreate(
        reinterpret_cast<HYPRE_Solver*>(&lagged_entry->boomeramg));
    if (status != 0) {
      return SolveFailure("HYPRE distributed BoomerAMG creation failed",
                          assembly.ownership,
                          "BoomerAMGCreate");
    }
    active_boomeramg = static_cast<HYPRE_Solver>(lagged_entry->boomeramg);
    result.lagged_amg_rebuild_used = true;
  } else {
    status = HYPRE_BoomerAMGCreate(&objects.boomeramg);
    if (status != 0) {
      return SolveFailure("HYPRE distributed BoomerAMG creation failed",
                          assembly.ownership,
                          "BoomerAMGCreate");
    }
    active_boomeramg = objects.boomeramg;
  }
  const auto amg_setup = ConfigureBoomerAmg(active_boomeramg);
  if (amg_setup.status != 0) {
    return SolveFailure(amg_setup.reason,
                        assembly.ownership,
                        amg_setup.phase,
                        amg_setup.status);
  }
  status = HYPRE_ParCSRGMRESSetTol(objects.gmres, options.relative_tolerance);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES tolerance setup failed", assembly.ownership, "ParCSRGMRESSetTol");
  }
  status = HYPRE_ParCSRGMRESSetMaxIter(objects.gmres, options.max_iterations);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES max-iteration setup failed", assembly.ownership, "ParCSRGMRESSetMaxIter");
  }
  status = HYPRE_ParCSRGMRESSetKDim(objects.gmres, options.krylov_dimension);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES Krylov setup failed", assembly.ownership, "ParCSRGMRESSetKDim");
  }
  status = HYPRE_ParCSRGMRESSetLogging(objects.gmres, options.logging_enabled ? 1 : 0);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES logging setup failed", assembly.ownership, "ParCSRGMRESSetLogging");
  }
  status = HYPRE_ParCSRGMRESSetPrintLevel(objects.gmres, options.print_level);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES print-level setup failed", assembly.ownership, "ParCSRGMRESSetPrintLevel");
  }
  status = HYPRE_ParCSRGMRESSetPrecond(
      objects.gmres,
      reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSolve),
      reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSetup),
      active_boomeramg);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES preconditioner setup failed", assembly.ownership, "ParCSRGMRESSetPrecond");
  }
  const auto setup_timer = std::chrono::steady_clock::now();
  status = HYPRE_ParCSRGMRESSetup(objects.gmres, parcsr_matrix, par_rhs, par_solution);
  result.hypre_setup_wall_s += ElapsedSecondsSince(setup_timer);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES setup failed", assembly.ownership, "ParCSRGMRESSetup");
  }
  const auto solve_timer = std::chrono::steady_clock::now();
  status = HYPRE_ParCSRGMRESSolve(objects.gmres, parcsr_matrix, par_rhs, par_solution);
  result.hypre_solve_wall_s += ElapsedSecondsSince(solve_timer);
  if (status != 0) {
    int failed_iterations = 0;
    double failed_residual = 0.0;
    (void)HYPRE_ParCSRGMRESGetNumIterations(objects.gmres, &failed_iterations);
    (void)HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(
        objects.gmres,
        &failed_residual);
    return SolveFailure("HYPRE distributed GMRES solve failed",
                        assembly.ownership,
                        "ParCSRGMRESSolve",
                        status,
                        failed_iterations,
                        failed_residual);
  }
  status = HYPRE_ParCSRGMRESGetNumIterations(objects.gmres, &result.gmres_iterations);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES iteration query failed", assembly.ownership, "ParCSRGMRESGetNumIterations");
  }
  status = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(
      objects.gmres,
      &result.gmres_final_relative_residual);
  if (status != 0) {
    return SolveFailure("HYPRE distributed GMRES residual query failed", assembly.ownership, "ParCSRGMRESGetFinalRelativeResidualNorm");
  }

  std::vector<HYPRE_Complex> solution_values(assembly.ownership.local_row_count);
  status = HYPRE_IJVectorGetValues(
      objects.ij_solution,
      local_row_count_int,
      rows.data(),
      solution_values.data());
  if (status != 0) {
    return SolveFailure("HYPRE distributed solution extraction failed", assembly.ownership, "IJVectorGetValues");
  }
  result.local_scalar_new.resize(assembly.ownership.local_row_count);
  for (std::size_t i = 0; i < assembly.ownership.local_row_count; ++i) {
    result.local_scalar_new[i] = static_cast<double>(solution_values[i]);
    if (!std::isfinite(result.local_scalar_new[i])) {
      return SolveFailure("HYPRE distributed solution contains non-finite values", assembly.ownership);
    }
  }

  ComputeDistributedResiduals(result, assembly);
  if (result.global_residual_linf_relative > options.relative_tolerance * 100.0) {
    return SolveFailure("HYPRE distributed recomputed residual exceeds tolerance", assembly.ownership);
  }
  if (lagged_enabled && lagged_entry != nullptr) {
    StoreLaggedAmgSetupSnapshot(*lagged_entry, assembly, result.gmres_iterations);
    result.lagged_amg_reuse_count_after =
        lagged_entry->reuse_count_since_setup;
  }
  result.success = true;
  result.report_line = BuildDistributedSolveReport(result, options);
  return result;
}

std::vector<double> GatherDistributedScalarForTest(
    const DistributedGenericDiffusionHypreSolveResult& result,
    MPI_Comm communicator) {
  int rank_count = 1;
  MPI_Comm_size(communicator, &rank_count);
  std::vector<int> counts(static_cast<std::size_t>(rank_count), 0);
  std::vector<int> displacements(static_cast<std::size_t>(rank_count), 0);
  const int local_count = static_cast<int>(result.local_scalar_new.size());
  MPI_Allgather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator);
  int total = 0;
  for (int i = 0; i < rank_count; ++i) {
    displacements[static_cast<std::size_t>(i)] = total;
    total += counts[static_cast<std::size_t>(i)];
  }
  std::vector<double> gathered(static_cast<std::size_t>(total), 0.0);
  MPI_Allgatherv(
      result.local_scalar_new.data(),
      local_count,
      MPI_DOUBLE,
      gathered.data(),
      counts.data(),
      displacements.data(),
      MPI_DOUBLE,
      communicator);
  return gathered;
}

bool DistributedGenericDiffusionAssemblyResult::is_complete() const noexcept {
  return success && ownership.success && local_matrix.is_shape_complete() &&
         local_rhs.size() == ownership.local_row_count &&
         local_scalar_old_flat.size() == ownership.local_row_count &&
         !report_line.empty() &&
         ValidateDistributedGenericDiffusionAssemblyDiagnostics(*this);
}

bool DistributedGenericDiffusionHypreSolveResult::is_complete() const noexcept {
  return success && backend == kHypreBackend && hypre_enabled && ownership.success &&
         local_scalar_new.size() == ownership.local_row_count &&
         !report_line.empty() &&
         ValidateDistributedGenericDiffusionHypreSolveDiagnostics(*this);
}

bool ValidateDistributedGenericDiffusionAssemblyDiagnostics(
    const DistributedGenericDiffusionAssemblyResult& result) noexcept {
  const auto& line = result.report_line;
  return result.success &&
         Contains(line, "diagnostic_id=p2.diffusion.distributed_assembly") &&
         Contains(line, "implementation_id=p2.diffusion.distributed_radial_slab_assembly_v1") &&
         Contains(line, "equation_form=A_dT_dt_div_D_grad_T_plus_C_T_plus_B") &&
         Contains(line, "unit_system=cgs") &&
         Contains(line, "unknown=generic_scalar") &&
         Contains(line, "direction_splitting=false") &&
         Contains(line, "interface_diffusion_mean=arithmetic") &&
         Contains(line, "pole_metric_mode=") &&
         Contains(line, "row_ownership=distributed_radial_slab") &&
         Contains(line, "mpi_rank=") &&
         Contains(line, "mpi_rank_count=") &&
         Contains(line, "global_radial_cells=") &&
         Contains(line, "global_theta_cells=") &&
         Contains(line, "global_phi_cells=") &&
         Contains(line, "global_row_count=") &&
         Contains(line, "local_global_radial_begin=") &&
         Contains(line, "local_global_radial_end=") &&
         Contains(line, "local_row_begin=") &&
         Contains(line, "local_row_end=") &&
         Contains(line, "local_row_count=") &&
         Contains(line, "local_nonzero_count=") &&
         Contains(line, "global_nonzero_count=") &&
         Contains(line, "local_off_rank_column_count=") &&
         Contains(line, "global_off_rank_column_count=") &&
         Contains(line, "local_radial_seam_coupling_count=") &&
         Contains(line, "global_radial_seam_coupling_count=") &&
         Contains(line, "coefficient_D_halo_lower_received=") &&
         Contains(line, "coefficient_D_halo_upper_received=") &&
         Contains(line, "origin_remap_used=") &&
         Contains(line, "pole_remap_used=") &&
         Contains(line, "origin_remap_coupling_mode=") &&
         Contains(line, "pole_remap_coupling_mode=") &&
         Contains(line, "singular_boundary_face_conductance_used=false") &&
         Contains(line, "canonical_state_mutated=false");
}

bool ValidateDistributedGenericDiffusionHypreSolveDiagnostics(
    const DistributedGenericDiffusionHypreSolveResult& result) noexcept {
  const auto& line = result.report_line;
  return result.success &&
         Contains(line, "diagnostic_id=p2.diffusion.hypre_distributed_solve") &&
         Contains(line, "implementation_id=p2.diffusion.hypre_distributed_parcsr_gmres_boomeramg_v1") &&
         Contains(line, "backend=hypre_parcsr_gmres_boomeramg") &&
         Contains(line, "hypre_enabled=true") &&
         Contains(line, "row_ownership=distributed_radial_slab") &&
         Contains(line, "matrix_source=p2.diffusion.distributed_radial_slab_assembly_v1") &&
         Contains(line, "matrix_format=hypre_ij_parcsr") &&
         Contains(line, "direction_splitting=false") &&
         Contains(line, "mpi_rank=") &&
         Contains(line, "mpi_rank_count=") &&
         Contains(line, "global_row_count=") &&
         Contains(line, "local_row_begin=") &&
         Contains(line, "local_row_end=") &&
         Contains(line, "local_row_count=") &&
         Contains(line, "local_nonzero_count=") &&
         Contains(line, "global_nonzero_count=") &&
         Contains(line, "local_off_rank_column_count=") &&
         Contains(line, "global_off_rank_column_count=") &&
         Contains(line, "solver=ParCSRGMRES") &&
         Contains(line, "preconditioner=BoomerAMG") &&
         Contains(line, "boomeramg_parameters=preconditioner_max_iter_1_tol_0") &&
         Contains(line, "boomeramg_max_iterations=1") &&
         Contains(line, "boomeramg_tolerance=0") &&
         Contains(line, "boomeramg_print_level=0") &&
         Contains(line, "boomeramg_strong_threshold=0.5") &&
         Contains(line, "lagged_amg_enabled=") &&
         Contains(line, "lagged_amg_candidate=") &&
         Contains(line, "lagged_amg_reuse_attempted=") &&
         Contains(line, "lagged_amg_reuse_accepted=") &&
         Contains(line, "lagged_amg_rebuild_used=") &&
         Contains(line, "lagged_amg_fallback_rebuild_used=") &&
         Contains(line, "lagged_amg_local_matrix_rel_change=") &&
         Contains(line, "lagged_amg_global_matrix_rel_change=") &&
         Contains(line, "lagged_amg_rebuild_every=") &&
         Contains(line, "lagged_amg_cached_reuse_count=") &&
         Contains(line, "lagged_amg_cached_iterations=") &&
         Contains(line, "lagged_amg_reuse_iterations=") &&
         Contains(line, "lagged_amg_reuse_count_after=") &&
         Contains(line, "lagged_amg_reuse_final_relative_residual=") &&
         Contains(line, "lagged_amg_reuse_setup_wall_s=") &&
         Contains(line, "lagged_amg_reuse_solve_wall_s=") &&
         Contains(line, "gmres_relative_tolerance=") &&
         Contains(line, "gmres_max_iterations=") &&
        Contains(line, "gmres_krylov_dimension=") &&
        Contains(line, "gmres_iterations=") &&
        Contains(line, "gmres_final_relative_residual=") &&
        Contains(line, "hypre_setup_wall_s=") &&
        Contains(line, "hypre_solve_wall_s=") &&
        Contains(line, "local_residual_l2=") &&
         Contains(line, "global_residual_l2=") &&
         Contains(line, "local_residual_linf=") &&
         Contains(line, "global_residual_linf=") &&
         Contains(line, "local_rhs_linf=") &&
         Contains(line, "global_rhs_linf=") &&
         Contains(line, "global_residual_linf_relative=") &&
         Contains(line, "max_abs_delta_local=") &&
         Contains(line, "max_abs_delta_global=") &&
         Contains(line, "canonical_state_mutated=false") &&
         Contains(line, "gathered_solve_used=false") &&
         Contains(line, "fallback_used=false");
}

}  // namespace dec3d::transport
