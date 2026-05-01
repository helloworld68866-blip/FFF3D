#include "transport/diffusion/generic_diffusion.hpp"

#include "mesh/boundary/spherical_scalar_remap.hpp"
#include "physics/units/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dec3d::transport {

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kBoundaryTolerance = 1.0e-14;

using RowAccumulator = std::vector<std::vector<std::pair<std::size_t, double>>>;

struct RemapAssemblyStats {
  std::size_t origin_unique_pair_count{0};
  std::size_t origin_row_coupling_count{0};
  std::size_t pole_unique_pair_count{0};
  std::size_t pole_row_coupling_count{0};
  std::size_t origin_nonzero_offdiag_count{0};
  std::size_t pole_nonzero_offdiag_count{0};
  double origin_min_abs_offdiag{0.0};
  double pole_min_abs_offdiag{0.0};
};

struct MarshakBoundaryStats {
  bool used{false};
  std::size_t outer_face_count{0};
  std::size_t nonzero_diagonal_loss_count{0};
  double min_g_face{0.0};
  double max_g_face{0.0};
  double min_distance_cm{0.0};
  double max_distance_cm{0.0};
  double min_abs_diagonal_loss{0.0};
  double max_abs_diagonal_loss{0.0};
};

enum class ScalarRemapCouplingKind {
  origin,
  pole,
};

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::size_t LinearIndex(
    const DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((radial * layout.theta_cells) + theta) * layout.phi_cells + phi;
}

[[nodiscard]] dec3d::mesh::ScalarRemapLayout ToScalarRemapLayout(
    const DiffusionGridLayout& layout) noexcept {
  return dec3d::mesh::ScalarRemapLayout{
      layout.radial_cells,
      layout.theta_cells,
      layout.phi_cells};
}

[[nodiscard]] bool ShapeMatches(
    const dec3d::core::Array3D<double>& values,
    const DiffusionGridLayout& layout) noexcept {
  return values.extent_r() == layout.radial_cells &&
         values.extent_theta() == layout.theta_cells &&
         values.extent_phi() == layout.phi_cells;
}

[[nodiscard]] bool FaceEffectiveShapeMatches(
    const GenericDiffusionFaceEffectiveCoefficients& face,
    const DiffusionGridLayout& layout) noexcept {
  if (!face.enabled) {
    return true;
  }
  return face.radial_face_D.extent_r() == layout.radial_cells + 1u &&
         face.radial_face_D.extent_theta() == layout.theta_cells &&
         face.radial_face_D.extent_phi() == layout.phi_cells &&
         face.theta_face_D.extent_r() == layout.radial_cells &&
         face.theta_face_D.extent_theta() == layout.theta_cells + 1u &&
         face.theta_face_D.extent_phi() == layout.phi_cells &&
         face.phi_face_D.extent_r() == layout.radial_cells &&
         face.phi_face_D.extent_theta() == layout.theta_cells &&
         face.phi_face_D.extent_phi() == layout.phi_cells;
}

[[nodiscard]] bool GeometryMatches(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const DiffusionGridLayout& layout) noexcept {
  return geometry.valid &&
         geometry.radial_faces.size() == layout.radial_cells + 1u &&
         geometry.theta_faces.size() == layout.theta_cells + 1u &&
         geometry.phi_faces.size() == layout.phi_cells + 1u &&
         geometry.cell_volumes.size() == layout.cell_count();
}

[[nodiscard]] bool BoundaryPolicyComplete(
    const GenericDiffusionBoundaryPolicy& policy) noexcept {
  return policy.inner_radial != DiffusionBoundaryKind::missing &&
         policy.outer_radial != DiffusionBoundaryKind::missing &&
         policy.theta_lower != DiffusionBoundaryKind::missing &&
         policy.theta_upper != DiffusionBoundaryKind::missing &&
         policy.phi != DiffusionBoundaryKind::missing;
}

[[nodiscard]] bool IsValidInnerRadialPolicy(DiffusionBoundaryKind kind) noexcept {
  return kind == DiffusionBoundaryKind::interior_patch_no_origin ||
         kind == DiffusionBoundaryKind::scalar_origin_remap_required;
}

[[nodiscard]] bool IsValidThetaPolicy(DiffusionBoundaryKind kind) noexcept {
  return kind == DiffusionBoundaryKind::interior_patch_no_pole ||
         kind == DiffusionBoundaryKind::scalar_pole_remap_required;
}

[[nodiscard]] const char* ToString(DiffusionBoundaryKind kind) noexcept {
  switch (kind) {
    case DiffusionBoundaryKind::missing:
      return "missing";
    case DiffusionBoundaryKind::neumann_zero_flux:
      return "neumann_zero_flux";
    case DiffusionBoundaryKind::periodic:
      return "periodic";
    case DiffusionBoundaryKind::scalar_origin_remap_required:
      return "scalar_origin_remap_required";
    case DiffusionBoundaryKind::scalar_pole_remap_required:
      return "scalar_pole_remap_required";
    case DiffusionBoundaryKind::interior_patch_no_origin:
      return "interior_patch_no_origin";
    case DiffusionBoundaryKind::interior_patch_no_pole:
      return "interior_patch_no_pole";
    case DiffusionBoundaryKind::radiation_marshak_vacuum:
      return "radiation_marshak_vacuum";
  }
  return "unknown";
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

[[nodiscard]] GenericDiffusionAssemblyResult AssemblyFailure(
    const std::string& reason,
    std::size_t radial = 0,
    std::size_t theta = 0,
    std::size_t phi = 0) {
  GenericDiffusionAssemblyResult result;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p2.diffusion.failure"
      << "; failure_reason=" << reason
      << "; first_bad_cell=" << radial << "," << theta << "," << phi
      << "; canonical_state_mutated=false";
  result.failure_diagnostics = out.str();
  return result;
}

[[nodiscard]] GenericDiffusionSolveResult SolveFailure(
    const std::string& reason) {
  GenericDiffusionSolveResult result;
  result.backend = "serial_dense_reference";
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p2.diffusion.failure"
      << "; failure_reason=" << reason
      << "; backend=serial_dense_reference"
      << "; canonical_state_mutated=false";
  result.failure_diagnostics = out.str();
  return result;
}

[[nodiscard]] double CellVolume(
    const GenericDiffusionProblem& problem,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return problem.geometry.cell_volumes[LinearIndex(problem.layout, radial, theta, phi)];
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

[[nodiscard]] double PhiCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t phi) noexcept {
  return 0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
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
    const GenericDiffusionProblem& problem,
    std::size_t theta) noexcept {
  if (problem.pole_metric_mode == DiffusionPoleMetricMode::axis_regular_polar_phi &&
      (theta == 0u || theta + 1u == problem.layout.theta_cells)) {
    return SolidAngleWeightedMeanSinTheta(problem.geometry, theta);
  }
  return std::sin(ThetaCenter(problem.geometry, theta));
}

[[nodiscard]] double ArithmeticMean(double lhs, double rhs) noexcept {
  return 0.5 * (lhs + rhs);
}

[[nodiscard]] double RadialInterfaceD(
    const GenericDiffusionProblem& problem,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.radial_face_D(radial_face, theta, phi);
  }
  return ArithmeticMean(
      problem.coefficient_D(radial_face - 1u, theta, phi),
      problem.coefficient_D(radial_face, theta, phi));
}

[[nodiscard]] double ThetaInterfaceD(
    const GenericDiffusionProblem& problem,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.theta_face_D(radial, theta_face, phi);
  }
  return ArithmeticMean(
      problem.coefficient_D(radial, theta_face - 1u, phi),
      problem.coefficient_D(radial, theta_face, phi));
}

[[nodiscard]] double PhiInterfaceD(
    const GenericDiffusionProblem& problem,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.phi_face_D(radial, theta, phi);
  }
  const std::size_t next_phi = (phi + 1u) % problem.layout.phi_cells;
  return ArithmeticMean(
      problem.coefficient_D(radial, theta, phi),
      problem.coefficient_D(radial, theta, next_phi));
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
  const double azimuth = PhiCenter(geometry, phi);
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
    const GenericDiffusionProblem& problem,
    std::size_t lhs_radial,
    std::size_t lhs_theta,
    std::size_t lhs_phi,
    std::size_t rhs_radial,
    std::size_t rhs_theta,
    std::size_t rhs_phi) noexcept {
  const double lhs_volume = CellVolume(problem, lhs_radial, lhs_theta, lhs_phi);
  const double rhs_volume = CellVolume(problem, rhs_radial, rhs_theta, rhs_phi);
  const double d_face = ArithmeticMean(
      problem.coefficient_D(lhs_radial, lhs_theta, lhs_phi),
      problem.coefficient_D(rhs_radial, rhs_theta, rhs_phi));
  const double distance_squared = SquaredDistance(
      CellCenterCartesian(problem.geometry, lhs_radial, lhs_theta, lhs_phi),
      CellCenterCartesian(problem.geometry, rhs_radial, rhs_theta, rhs_phi));
  if (!std::isfinite(lhs_volume) || !std::isfinite(rhs_volume) ||
      !std::isfinite(d_face) || !std::isfinite(distance_squared) ||
      lhs_volume <= 0.0 || rhs_volume <= 0.0 || distance_squared <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double effective_volume = 0.5 * (lhs_volume + rhs_volume);
  return d_face * effective_volume / distance_squared;
}

void AddMatrixEntry(
    RowAccumulator& rows,
    std::size_t row,
    std::size_t column,
    double value) {
  auto& entries = rows[row];
  for (auto& entry : entries) {
    if (entry.first == column) {
      entry.second += value;
      return;
    }
  }
  entries.emplace_back(column, value);
}

void AddConductance(
    RowAccumulator& rows,
    std::size_t lhs,
    std::size_t rhs,
    double conductance) {
  if (conductance == 0.0 || lhs == rhs) {
    return;
  }

  AddMatrixEntry(rows, lhs, lhs, conductance);
  AddMatrixEntry(rows, lhs, rhs, -conductance);
  AddMatrixEntry(rows, rhs, rhs, conductance);
  AddMatrixEntry(rows, rhs, lhs, -conductance);
}

void AddMappedConductance(
    const GenericDiffusionProblem& problem,
    RowAccumulator& rows,
    RemapAssemblyStats& stats,
    ScalarRemapCouplingKind kind,
    std::size_t lhs_radial,
    std::size_t lhs_theta,
    std::size_t lhs_phi,
    const dec3d::mesh::ScalarRemapIndex& rhs_index) {
  const auto& layout = problem.layout;
  const std::size_t lhs = LinearIndex(layout, lhs_radial, lhs_theta, lhs_phi);
  const std::size_t rhs =
      LinearIndex(layout, rhs_index.radial, rhs_index.theta, rhs_index.phi);
  if (lhs == rhs) {
    return;
  }

  const std::size_t first = std::min(lhs, rhs);
  const std::size_t second = std::max(lhs, rhs);
  if (lhs != first) {
    return;
  }

  const double conductance = MappedCenterConductance(
      problem,
      lhs_radial,
      lhs_theta,
      lhs_phi,
      rhs_index.radial,
      rhs_index.theta,
      rhs_index.phi);
  if (!std::isfinite(conductance) || conductance <= 0.0) {
    return;
  }

  AddConductance(rows, first, second, conductance);

  switch (kind) {
    case ScalarRemapCouplingKind::origin:
      ++stats.origin_unique_pair_count;
      stats.origin_row_coupling_count += 2u;
      stats.origin_nonzero_offdiag_count += 2u;
      stats.origin_min_abs_offdiag =
          stats.origin_min_abs_offdiag == 0.0
              ? conductance
              : std::min(stats.origin_min_abs_offdiag, conductance);
      break;
    case ScalarRemapCouplingKind::pole:
      ++stats.pole_unique_pair_count;
      stats.pole_row_coupling_count += 2u;
      stats.pole_nonzero_offdiag_count += 2u;
      stats.pole_min_abs_offdiag =
          stats.pole_min_abs_offdiag == 0.0
              ? conductance
              : std::min(stats.pole_min_abs_offdiag, conductance);
      break;
  }
}

[[nodiscard]] SparseMatrixCsr BuildCsrFromRows(RowAccumulator& rows) {
  SparseMatrixCsr matrix;
  matrix.row_count = rows.size();
  matrix.column_count = rows.size();
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

[[nodiscard]] std::size_t CountBoundaryRows(const DiffusionGridLayout& layout) noexcept {
  std::size_t count = 0;
  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const bool touches_boundary =
            radial == 0 || radial + 1u == layout.radial_cells ||
            theta == 0 || theta + 1u == layout.theta_cells ||
            layout.phi_cells > 1u;
        if (touches_boundary) {
          ++count;
        }
      }
    }
  }
  return count;
}

void FillScalarOldFlat(
    const GenericDiffusionProblem& problem,
    GenericDiffusionAssemblyResult& result) {
  result.scalar_old_flat.assign(problem.layout.cell_count(), 0.0);
  for (std::size_t radial = 0; radial < problem.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < problem.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
        result.scalar_old_flat[LinearIndex(problem.layout, radial, theta, phi)] =
            problem.scalar_old(radial, theta, phi);
      }
    }
  }
}

void MarkScalarRemapAvailability(
    GenericDiffusionAssemblyResult& result,
    bool origin_policy_available,
    bool pole_policy_available) noexcept {
  result.scalar_origin_remap_available = origin_policy_available;
  result.scalar_pole_remap_available = pole_policy_available;
  result.singular_boundary_face_conductance_used = false;
}

void PublishRemapStats(
    GenericDiffusionAssemblyResult& result,
    bool origin_touched,
    bool pole_touched,
    const RemapAssemblyStats& stats) noexcept {
  result.origin_remap_used = origin_touched;
  result.pole_remap_used = pole_touched;
  result.origin_remap_unique_pair_count = stats.origin_unique_pair_count;
  result.origin_remap_row_coupling_count = stats.origin_row_coupling_count;
  result.pole_remap_unique_pair_count = stats.pole_unique_pair_count;
  result.pole_remap_row_coupling_count = stats.pole_row_coupling_count;
  result.origin_remap_nonzero_offdiag_count = stats.origin_nonzero_offdiag_count;
  result.pole_remap_nonzero_offdiag_count = stats.pole_nonzero_offdiag_count;
  result.origin_remap_min_abs_offdiag = stats.origin_min_abs_offdiag;
  result.pole_remap_min_abs_offdiag = stats.pole_min_abs_offdiag;
  result.singular_boundary_face_conductance_used = false;
}

void PublishMarshakStats(
    GenericDiffusionAssemblyResult& result,
    const MarshakBoundaryStats& stats) noexcept {
  result.marshak_boundary_used = stats.used;
  result.marshak_outer_face_count = stats.outer_face_count;
  result.marshak_nonzero_diagonal_loss_count = stats.nonzero_diagonal_loss_count;
  result.marshak_min_G_face = stats.min_g_face;
  result.marshak_max_G_face = stats.max_g_face;
  result.marshak_min_center_to_boundary_distance_cm = stats.min_distance_cm;
  result.marshak_max_center_to_boundary_distance_cm = stats.max_distance_cm;
  result.marshak_min_abs_diagonal_loss = stats.min_abs_diagonal_loss;
  result.marshak_max_abs_diagonal_loss = stats.max_abs_diagonal_loss;
}

void FinalizeMatrixDiagnostics(GenericDiffusionAssemblyResult& result) {
  result.row_count = result.matrix.row_count;
  result.nonzero_count = result.matrix.values.size();
  result.min_diagonal = std::numeric_limits<double>::infinity();
  result.max_abs_offdiagonal_row_sum = 0.0;
  result.max_rhs_abs = 0.0;

  for (std::size_t row = 0; row < result.matrix.row_count; ++row) {
    double offdiagonal_sum = 0.0;
    double diagonal = 0.0;
    for (std::size_t slot = result.matrix.row_offsets[row];
         slot < result.matrix.row_offsets[row + 1u];
         ++slot) {
      if (result.matrix.column_indices[slot] == row) {
        diagonal += result.matrix.values[slot];
      } else {
        offdiagonal_sum += std::abs(result.matrix.values[slot]);
      }
    }
    result.min_diagonal = std::min(result.min_diagonal, diagonal);
    result.max_abs_offdiagonal_row_sum =
        std::max(result.max_abs_offdiagonal_row_sum, offdiagonal_sum);
  }

  if (!std::isfinite(result.min_diagonal)) {
    result.min_diagonal = 0.0;
  }

  for (double rhs_value : result.rhs) {
    result.max_rhs_abs = std::max(result.max_rhs_abs, std::abs(rhs_value));
  }
}

[[nodiscard]] std::string BuildAssemblyReport(
    const GenericDiffusionProblem& problem,
    const GenericDiffusionAssemblyResult& result,
    bool origin_touched,
    bool pole_touched,
    double min_a,
    double min_d,
    double max_d) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.assembly"
      << "; implementation_id=p2.diffusion.generic_implicit_assembly_v1"
      << "; equation_form=A_dT_dt_div_D_grad_T_plus_C_T_plus_B"
      << "; unit_system=cgs"
      << "; dt_s=" << problem.dt_s
      << "; unknown=generic_scalar"
      << "; direction_splitting=false"
      << "; matrix_format=csr"
      << "; interface_diffusion_mean=arithmetic"
      << "; pole_metric_mode=" << ToString(problem.pole_metric_mode)
      << "; row_count=" << result.row_count
      << "; nonzero_count=" << result.nonzero_count
      << "; boundary_row_count=" << result.boundary_row_count
      << "; origin_boundary_touched=" << (origin_touched ? "true" : "false")
      << "; pole_boundary_touched=" << (pole_touched ? "true" : "false")
      << "; outer_boundary_policy=" << ToString(problem.boundary_policy.outer_radial)
      << "; phi_boundary_policy=periodic"
      << "; scalar_origin_remap_available="
      << (result.scalar_origin_remap_available ? "true" : "false")
      << "; scalar_pole_remap_available="
      << (result.scalar_pole_remap_available ? "true" : "false")
      << "; scalar_remap_value_transform=identity"
      << "; origin_remap_used=" << (result.origin_remap_used ? "true" : "false")
      << "; pole_remap_used=" << (result.pole_remap_used ? "true" : "false")
      << "; origin_remap_coupling_mode="
      << (result.origin_remap_used ? "metric_regular_origin" : "not_used")
      << "; pole_remap_coupling_mode="
      << (result.pole_remap_used ? "metric_regular_axis" : "not_used")
      << "; origin_remap_unique_pair_count=" << result.origin_remap_unique_pair_count
      << "; origin_remap_row_coupling_count=" << result.origin_remap_row_coupling_count
      << "; pole_remap_unique_pair_count=" << result.pole_remap_unique_pair_count
      << "; pole_remap_row_coupling_count=" << result.pole_remap_row_coupling_count
      << "; origin_remap_nonzero_offdiag_count="
      << result.origin_remap_nonzero_offdiag_count
      << "; pole_remap_nonzero_offdiag_count=" << result.pole_remap_nonzero_offdiag_count
      << "; origin_remap_min_abs_offdiag=" << result.origin_remap_min_abs_offdiag
      << "; pole_remap_min_abs_offdiag=" << result.pole_remap_min_abs_offdiag
      << "; singular_boundary_face_conductance_used="
      << (result.singular_boundary_face_conductance_used ? "true" : "false")
      << "; marshak_boundary_used=" << (result.marshak_boundary_used ? "true" : "false")
      << "; marshak_outer_face_count=" << result.marshak_outer_face_count
      << "; marshak_nonzero_diagonal_loss_count="
      << result.marshak_nonzero_diagonal_loss_count
      << "; marshak_min_G_face=" << result.marshak_min_G_face
      << "; marshak_max_G_face=" << result.marshak_max_G_face
      << "; marshak_min_center_to_boundary_distance_cm="
      << result.marshak_min_center_to_boundary_distance_cm
      << "; marshak_max_center_to_boundary_distance_cm="
      << result.marshak_max_center_to_boundary_distance_cm
      << "; marshak_min_abs_diagonal_loss=" << result.marshak_min_abs_diagonal_loss
      << "; marshak_max_abs_diagonal_loss=" << result.marshak_max_abs_diagonal_loss
      << "; min_A=" << min_a
      << "; min_D=" << min_d
      << "; max_D=" << max_d
      << "; min_diagonal=" << result.min_diagonal
      << "; max_abs_offdiagonal_row_sum=" << result.max_abs_offdiagonal_row_sum
      << "; max_rhs_abs=" << result.max_rhs_abs;
  return out.str();
}

[[nodiscard]] bool ValidateInputs(
    const GenericDiffusionProblem& problem,
    double* min_a,
    double* min_d,
    double* max_d,
    bool* origin_touched,
    bool* pole_touched,
    GenericDiffusionAssemblyResult* failure) {
  if (!problem.layout.is_valid()) {
    *failure = AssemblyFailure("layout is invalid");
    return false;
  }
  if (!GeometryMatches(problem.geometry, problem.layout)) {
    *failure = AssemblyFailure("geometry is invalid or shape-mismatched");
    return false;
  }
  if (!std::isfinite(problem.dt_s) || problem.dt_s < 0.0) {
    *failure = AssemblyFailure("dt_s must be finite and non-negative");
    return false;
  }
  if (!ShapeMatches(problem.coefficient_A, problem.layout) ||
      !ShapeMatches(problem.coefficient_D, problem.layout) ||
      !ShapeMatches(problem.coefficient_C, problem.layout) ||
      !ShapeMatches(problem.coefficient_B, problem.layout) ||
      !ShapeMatches(problem.scalar_old, problem.layout)) {
    *failure = AssemblyFailure("coefficient or scalar array shape mismatch");
    return false;
  }
  if (!FaceEffectiveShapeMatches(problem.face_effective_coefficients, problem.layout)) {
    *failure = AssemblyFailure("face-effective diffusion coefficient shape mismatch");
    return false;
  }
  if (!BoundaryPolicyComplete(problem.boundary_policy)) {
    *failure = AssemblyFailure("boundary policy is missing");
    return false;
  }
  if (!IsValidInnerRadialPolicy(problem.boundary_policy.inner_radial)) {
    *failure = AssemblyFailure("inner radial boundary policy is unsupported");
    return false;
  }
  if (!IsValidThetaPolicy(problem.boundary_policy.theta_lower) ||
      !IsValidThetaPolicy(problem.boundary_policy.theta_upper)) {
    *failure = AssemblyFailure("theta boundary policy is unsupported");
    return false;
  }
  if (problem.boundary_policy.phi != DiffusionBoundaryKind::periodic) {
    *failure = AssemblyFailure("phi boundary policy must be periodic");
    return false;
  }
  if (problem.boundary_policy.outer_radial != DiffusionBoundaryKind::neumann_zero_flux &&
      problem.boundary_policy.outer_radial != DiffusionBoundaryKind::radiation_marshak_vacuum) {
    *failure = AssemblyFailure("outer radial boundary policy is unsupported");
    return false;
  }

  *origin_touched =
      !problem.geometry.radial_faces.empty() &&
      std::abs(problem.geometry.radial_faces.front()) <= kBoundaryTolerance;
  *pole_touched =
      !problem.geometry.theta_faces.empty() &&
      (std::abs(problem.geometry.theta_faces.front()) <= kBoundaryTolerance ||
       std::abs(problem.geometry.theta_faces.back() - kPi) <= kBoundaryTolerance);

  if (*origin_touched &&
      problem.boundary_policy.inner_radial == DiffusionBoundaryKind::interior_patch_no_origin) {
    *failure = AssemblyFailure("patch boundary policy touches origin");
    return false;
  }
  if (*pole_touched &&
      (problem.boundary_policy.theta_lower == DiffusionBoundaryKind::interior_patch_no_pole ||
       problem.boundary_policy.theta_upper == DiffusionBoundaryKind::interior_patch_no_pole)) {
    *failure = AssemblyFailure("patch boundary policy touches pole");
    return false;
  }
  if (*origin_touched &&
      problem.boundary_policy.inner_radial != DiffusionBoundaryKind::scalar_origin_remap_required) {
    *failure = AssemblyFailure("scalar origin remap boundary policy is required");
    return false;
  }
  if (*pole_touched &&
      (problem.boundary_policy.theta_lower != DiffusionBoundaryKind::scalar_pole_remap_required ||
       problem.boundary_policy.theta_upper != DiffusionBoundaryKind::scalar_pole_remap_required)) {
    *failure = AssemblyFailure("scalar pole remap boundary policy is required");
    return false;
  }
  if (*origin_touched || *pole_touched) {
    const auto validation =
        dec3d::mesh::ValidateScalarHalfTurnTopology(ToScalarRemapLayout(problem.layout));
    if (!validation.success) {
      *failure = AssemblyFailure(validation.failure_reason);
      return false;
    }
  }

  *min_a = std::numeric_limits<double>::infinity();
  *min_d = std::numeric_limits<double>::infinity();
  *max_d = 0.0;

  for (std::size_t radial = 0; radial < problem.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < problem.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
        const double a = problem.coefficient_A(radial, theta, phi);
        const double d = problem.coefficient_D(radial, theta, phi);
        const double c = problem.coefficient_C(radial, theta, phi);
        const double b = problem.coefficient_B(radial, theta, phi);
        const double old_scalar = problem.scalar_old(radial, theta, phi);
        if (!std::isfinite(a) || a <= 0.0) {
          *failure = AssemblyFailure("coefficient_A must be positive", radial, theta, phi);
          return false;
        }
        if (!std::isfinite(d) || d < 0.0) {
          *failure = AssemblyFailure("coefficient_D must be non-negative", radial, theta, phi);
          return false;
        }
        if (!std::isfinite(c)) {
          *failure = AssemblyFailure("coefficient_C must be finite", radial, theta, phi);
          return false;
        }
        if (!std::isfinite(b)) {
          *failure = AssemblyFailure("coefficient_B must be finite", radial, theta, phi);
          return false;
        }
        if (!std::isfinite(old_scalar)) {
          *failure = AssemblyFailure("scalar_old must be finite", radial, theta, phi);
          return false;
        }
        *min_a = std::min(*min_a, a);
        *min_d = std::min(*min_d, d);
        *max_d = std::max(*max_d, d);
      }
    }
  }

  if (problem.face_effective_coefficients.enabled) {
    for (std::size_t radial_face = 0; radial_face <= problem.layout.radial_cells; ++radial_face) {
      for (std::size_t theta = 0; theta < problem.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
          const double value = problem.face_effective_coefficients.radial_face_D(radial_face, theta, phi);
          if (!std::isfinite(value) || value < 0.0) {
            *failure = AssemblyFailure("radial face-effective D must be non-negative");
            return false;
          }
        }
      }
    }
    for (std::size_t radial = 0; radial < problem.layout.radial_cells; ++radial) {
      for (std::size_t theta_face = 0; theta_face <= problem.layout.theta_cells; ++theta_face) {
        for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
          const double value = problem.face_effective_coefficients.theta_face_D(radial, theta_face, phi);
          if (!std::isfinite(value) || value < 0.0) {
            *failure = AssemblyFailure("theta face-effective D must be non-negative");
            return false;
          }
        }
      }
    }
    for (std::size_t radial = 0; radial < problem.layout.radial_cells; ++radial) {
      for (std::size_t theta = 0; theta < problem.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
          const double value = problem.face_effective_coefficients.phi_face_D(radial, theta, phi);
          if (!std::isfinite(value) || value < 0.0) {
            *failure = AssemblyFailure("phi face-effective D must be non-negative");
            return false;
          }
        }
      }
    }
  }

  return true;
}

[[nodiscard]] GenericDiffusionAssemblyResult BuildZeroDtAssembly(
    const GenericDiffusionProblem& problem,
    bool origin_touched,
    bool pole_touched,
    double min_a,
    double min_d,
    double max_d) {
  GenericDiffusionAssemblyResult result;
  result.success = true;
  result.row_count = problem.layout.cell_count();
  result.matrix.row_count = result.row_count;
  result.matrix.column_count = result.row_count;
  result.matrix.row_offsets.reserve(result.row_count + 1u);
  result.matrix.row_offsets.push_back(0u);
  result.rhs.resize(result.row_count, 0.0);
  FillScalarOldFlat(problem, result);
  MarkScalarRemapAvailability(
      result,
      problem.boundary_policy.inner_radial == DiffusionBoundaryKind::scalar_origin_remap_required,
      problem.boundary_policy.theta_lower == DiffusionBoundaryKind::scalar_pole_remap_required &&
          problem.boundary_policy.theta_upper == DiffusionBoundaryKind::scalar_pole_remap_required);

  for (std::size_t row = 0; row < result.row_count; ++row) {
    result.matrix.column_indices.push_back(row);
    result.matrix.values.push_back(1.0);
    result.matrix.row_offsets.push_back(result.matrix.values.size());
    result.rhs[row] = result.scalar_old_flat[row];
  }

  result.boundary_row_count = result.row_count;
  FinalizeMatrixDiagnostics(result);
  result.report_line = BuildAssemblyReport(
      problem,
      result,
      origin_touched,
      pole_touched,
      min_a,
      min_d,
      max_d);
  return result;
}

void AssembleInteriorConductances(
    const GenericDiffusionProblem& problem,
    RowAccumulator& rows) {
  const auto& layout = problem.layout;

  for (std::size_t radial = 0; radial + 1u < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const std::size_t left = LinearIndex(layout, radial, theta, phi);
        const std::size_t right = LinearIndex(layout, radial + 1u, theta, phi);
        const double area = RadialFaceArea(problem.geometry, radial + 1u, theta, phi);
        const double distance =
            RadialCenter(problem.geometry, radial + 1u) -
            RadialCenter(problem.geometry, radial);
        const double d_face = RadialInterfaceD(problem, radial + 1u, theta, phi);
        AddConductance(rows, left, right, area * d_face / distance);
      }
    }
  }

  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta + 1u < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const std::size_t lower = LinearIndex(layout, radial, theta, phi);
        const std::size_t upper = LinearIndex(layout, radial, theta + 1u, phi);
        const double area = ThetaFaceArea(problem.geometry, radial, theta + 1u, phi);
        const double radius = RadialCenter(problem.geometry, radial);
        const double distance =
            radius * (ThetaCenter(problem.geometry, theta + 1u) -
                      ThetaCenter(problem.geometry, theta));
        const double d_face = ThetaInterfaceD(problem, radial, theta + 1u, phi);
        AddConductance(rows, lower, upper, area * d_face / distance);
      }
    }
  }

  if (layout.phi_cells <= 1u) {
    return;
  }

  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const std::size_t next_phi = (phi + 1u) % layout.phi_cells;
        if (layout.phi_cells == 2u && phi == 1u) {
          continue;
        }
        const std::size_t left = LinearIndex(layout, radial, theta, phi);
        const std::size_t right = LinearIndex(layout, radial, theta, next_phi);
        const double radius = RadialCenter(problem.geometry, radial);
        const double sin_theta = PhiMetricSinTheta(problem, theta);
        const double dphi = problem.geometry.phi_faces[1u] - problem.geometry.phi_faces[0u];
        const double distance = radius * sin_theta * dphi;
        const double area = PhiFaceArea(problem.geometry, radial, theta);
        const double d_face = PhiInterfaceD(problem, radial, theta, phi);
        AddConductance(rows, left, right, area * d_face / distance);
      }
    }
  }
}

void AssembleScalarRemapConductances(
    const GenericDiffusionProblem& problem,
    bool origin_touched,
    bool pole_touched,
    RowAccumulator& rows,
    RemapAssemblyStats& stats) {
  (void)problem;
  (void)origin_touched;
  (void)pole_touched;
  (void)rows;
  (void)stats;
}

[[nodiscard]] double OuterRadialBoundaryD(
    const GenericDiffusionProblem& problem,
    std::size_t theta,
    std::size_t phi) noexcept {
  if (problem.face_effective_coefficients.enabled) {
    return problem.face_effective_coefficients.radial_face_D(
        problem.layout.radial_cells, theta, phi);
  }
  return problem.coefficient_D(problem.layout.radial_cells - 1u, theta, phi);
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

void AccumulateMarshakRange(
    double value,
    double* min_value,
    double* max_value) noexcept {
  if (*min_value == 0.0) {
    *min_value = value;
  } else {
    *min_value = std::min(*min_value, value);
  }
  *max_value = std::max(*max_value, value);
}

[[nodiscard]] bool AssembleOuterRadialMarshakBoundary(
    const GenericDiffusionProblem& problem,
    RowAccumulator& rows,
    MarshakBoundaryStats* stats,
    GenericDiffusionAssemblyResult* failure) {
  if (problem.boundary_policy.outer_radial != DiffusionBoundaryKind::radiation_marshak_vacuum) {
    return true;
  }

  stats->used = true;
  const auto& layout = problem.layout;
  const std::size_t radial = layout.radial_cells - 1u;
  const double r_center = RadialCenter(problem.geometry, radial);
  const double distance = problem.geometry.radial_faces.back() - r_center;
  if (!std::isfinite(distance) || distance <= 0.0) {
    *failure = AssemblyFailure("Marshak center-to-boundary distance must be positive");
    return false;
  }

  for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
    for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
      ++stats->outer_face_count;
      const double d_face = OuterRadialBoundaryD(problem, theta, phi);
      const double g_face = MarshakGFace(d_face, distance);
      const double area = RadialFaceArea(problem.geometry, layout.radial_cells, theta, phi);
      const double loss = area * g_face;
      if (!std::isfinite(d_face) || !std::isfinite(g_face) ||
          !std::isfinite(area) || !std::isfinite(loss) ||
          d_face < 0.0 || g_face < 0.0 || area < 0.0 || loss < 0.0) {
        *failure = AssemblyFailure("Marshak boundary conductance is not physical");
        return false;
      }

      AccumulateMarshakRange(g_face, &stats->min_g_face, &stats->max_g_face);
      AccumulateMarshakRange(distance, &stats->min_distance_cm, &stats->max_distance_cm);
      if (loss > 0.0) {
        const std::size_t row = LinearIndex(layout, radial, theta, phi);
        AddMatrixEntry(rows, row, row, loss);
        ++stats->nonzero_diagonal_loss_count;
        AccumulateMarshakRange(
            std::abs(loss),
            &stats->min_abs_diagonal_loss,
            &stats->max_abs_diagonal_loss);
      }
    }
  }

  if (stats->max_g_face > 0.0 && stats->nonzero_diagonal_loss_count == 0u) {
    *failure = AssemblyFailure("Marshak boundary requested but no diagonal loss was assembled");
    return false;
  }
  return true;
}

void AssembleTimeSourceAndRhs(
    const GenericDiffusionProblem& problem,
    RowAccumulator& rows,
    GenericDiffusionAssemblyResult& result) {
  const auto& layout = problem.layout;
  result.rhs.assign(layout.cell_count(), 0.0);
  FillScalarOldFlat(problem, result);

  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const std::size_t row = LinearIndex(layout, radial, theta, phi);
        const double volume = CellVolume(problem, radial, theta, phi);
        const double time_coefficient =
            volume * problem.coefficient_A(radial, theta, phi) / problem.dt_s;
        const double c_coefficient = volume * problem.coefficient_C(radial, theta, phi);
        const double b_source = volume * problem.coefficient_B(radial, theta, phi);
        AddMatrixEntry(rows, row, row, time_coefficient - c_coefficient);
        result.rhs[row] =
            time_coefficient * problem.scalar_old(radial, theta, phi) + b_source;
      }
    }
  }
}

[[nodiscard]] bool SolveDense(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs,
    std::vector<double>* solution_out) {
  const std::size_t n = rhs.size();
  for (std::size_t pivot = 0; pivot < n; ++pivot) {
    std::size_t best = pivot;
    double best_abs = std::abs(matrix[pivot][pivot]);
    for (std::size_t row = pivot + 1u; row < n; ++row) {
      const double candidate = std::abs(matrix[row][pivot]);
      if (candidate > best_abs) {
        best = row;
        best_abs = candidate;
      }
    }
    if (best_abs <= 0.0 || !std::isfinite(best_abs)) {
      return false;
    }
    if (best != pivot) {
      std::swap(matrix[pivot], matrix[best]);
      std::swap(rhs[pivot], rhs[best]);
    }
    const double diagonal = matrix[pivot][pivot];
    for (std::size_t column = pivot; column < n; ++column) {
      matrix[pivot][column] /= diagonal;
    }
    rhs[pivot] /= diagonal;
    for (std::size_t row = 0; row < n; ++row) {
      if (row == pivot) {
        continue;
      }
      const double factor = matrix[row][pivot];
      if (factor == 0.0) {
        continue;
      }
      for (std::size_t column = pivot; column < n; ++column) {
        matrix[row][column] -= factor * matrix[pivot][column];
      }
      rhs[row] -= factor * rhs[pivot];
    }
  }
  *solution_out = std::move(rhs);
  return true;
}

[[nodiscard]] std::vector<std::vector<double>> CsrToDense(const SparseMatrixCsr& matrix) {
  std::vector<std::vector<double>> dense(
      matrix.row_count,
      std::vector<double>(matrix.column_count, 0.0));
  for (std::size_t row = 0; row < matrix.row_count; ++row) {
    for (std::size_t slot = matrix.row_offsets[row];
         slot < matrix.row_offsets[row + 1u];
         ++slot) {
      dense[row][matrix.column_indices[slot]] += matrix.values[slot];
    }
  }
  return dense;
}

void ComputeResiduals(
    const SparseMatrixCsr& matrix,
    const std::vector<double>& solution,
    const std::vector<double>& rhs,
    double* residual_l2,
    double* residual_linf,
    double* rhs_linf) {
  double l2_sum = 0.0;
  double linf = 0.0;
  double rhs_max = 0.0;

  for (std::size_t row = 0; row < matrix.row_count; ++row) {
    double ax = 0.0;
    for (std::size_t slot = matrix.row_offsets[row];
         slot < matrix.row_offsets[row + 1u];
         ++slot) {
      ax += matrix.values[slot] * solution[matrix.column_indices[slot]];
    }
    const double residual = ax - rhs[row];
    l2_sum += residual * residual;
    linf = std::max(linf, std::abs(residual));
    rhs_max = std::max(rhs_max, std::abs(rhs[row]));
  }

  *residual_l2 = std::sqrt(l2_sum);
  *residual_linf = linf;
  *rhs_linf = rhs_max;
}

[[nodiscard]] std::string BuildSolveReport(
    const GenericDiffusionSolveResult& result,
    const GenericDiffusionReferenceSolveOptions& options) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.reference_solve"
      << "; implementation_id=p2.diffusion.serial_dense_reference_v1"
      << "; backend=serial_dense_reference"
      << "; row_count=" << result.row_count
      << "; row_limit=" << options.row_limit
      << "; residual_l2=" << result.residual_l2
      << "; residual_linf=" << result.residual_linf
      << "; rhs_linf=" << result.rhs_linf
      << "; residual_linf_relative=" << result.residual_linf_relative
      << "; max_abs_delta=" << result.max_abs_delta
      << "; canonical_state_mutated=false";
  return out.str();
}

[[nodiscard]] GenericDiffusionSolveResult SolveResidualFailure(
    GenericDiffusionSolveResult result,
    const GenericDiffusionReferenceSolveOptions& options) {
  result.success = false;
  result.failure_reason = "reference solve residual exceeds locked tolerance";
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.failure"
      << "; failure_reason=" << result.failure_reason
      << "; backend=serial_dense_reference"
      << "; row_count=" << result.row_count
      << "; residual_l2=" << result.residual_l2
      << "; residual_linf=" << result.residual_linf
      << "; rhs_linf=" << result.rhs_linf
      << "; residual_linf_relative=" << result.residual_linf_relative
      << "; residual_tolerance=" << options.residual_tolerance
      << "; canonical_state_mutated=false";
  result.failure_diagnostics = out.str();
  return result;
}

}  // namespace

bool SparseMatrixCsr::is_shape_complete() const noexcept {
  return row_count == column_count &&
         row_offsets.size() == row_count + 1u &&
         column_indices.size() == values.size() &&
         (row_offsets.empty() || row_offsets.back() == values.size());
}

GenericDiffusionAssemblyResult AssembleGenericImplicitDiffusionSystem(
    const GenericDiffusionProblem& problem) noexcept {
  double min_a = 0.0;
  double min_d = 0.0;
  double max_d = 0.0;
  bool origin_touched = false;
  bool pole_touched = false;
  GenericDiffusionAssemblyResult failure;
  if (!ValidateInputs(
          problem,
          &min_a,
          &min_d,
          &max_d,
          &origin_touched,
          &pole_touched,
          &failure)) {
    return failure;
  }

  if (problem.dt_s == 0.0) {
    return BuildZeroDtAssembly(problem, origin_touched, pole_touched, min_a, min_d, max_d);
  }

  GenericDiffusionAssemblyResult result;
  result.success = true;
  MarkScalarRemapAvailability(
      result,
      problem.boundary_policy.inner_radial == DiffusionBoundaryKind::scalar_origin_remap_required,
      problem.boundary_policy.theta_lower == DiffusionBoundaryKind::scalar_pole_remap_required &&
          problem.boundary_policy.theta_upper == DiffusionBoundaryKind::scalar_pole_remap_required);
  RowAccumulator rows(problem.layout.cell_count());
  AssembleInteriorConductances(problem, rows);
  RemapAssemblyStats remap_stats;
  AssembleScalarRemapConductances(problem, origin_touched, pole_touched, rows, remap_stats);
  PublishRemapStats(result, origin_touched, pole_touched, remap_stats);
  MarshakBoundaryStats marshak_stats;
  GenericDiffusionAssemblyResult marshak_failure;
  if (!AssembleOuterRadialMarshakBoundary(problem, rows, &marshak_stats, &marshak_failure)) {
    return marshak_failure;
  }
  PublishMarshakStats(result, marshak_stats);
  AssembleTimeSourceAndRhs(problem, rows, result);
  result.matrix = BuildCsrFromRows(rows);
  result.boundary_row_count = CountBoundaryRows(problem.layout);

  if (!result.matrix.is_shape_complete()) {
    return AssemblyFailure("assembled matrix has missing rows");
  }
  for (double value : result.matrix.values) {
    if (!std::isfinite(value)) {
      return AssemblyFailure("assembled matrix has non-finite entries");
    }
  }
  for (double value : result.rhs) {
    if (!std::isfinite(value)) {
      return AssemblyFailure("rhs has non-finite entries");
    }
  }

  FinalizeMatrixDiagnostics(result);
  result.report_line = BuildAssemblyReport(
      problem,
      result,
      origin_touched,
      pole_touched,
      min_a,
      min_d,
      max_d);
  return result;
}

GenericDiffusionSolveResult SolveGenericDiffusionReference(
    const GenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionReferenceSolveOptions& options) noexcept {
  if (!assembly.success || !assembly.matrix.is_shape_complete() ||
      assembly.rhs.size() != assembly.matrix.row_count ||
      assembly.scalar_old_flat.size() != assembly.matrix.row_count) {
    return SolveFailure("assembly result is incomplete");
  }
  if (assembly.matrix.row_count > options.row_limit) {
    return SolveFailure("reference backend row_count exceeds limit");
  }
  if (!std::isfinite(options.residual_tolerance) || options.residual_tolerance < 0.0) {
    return SolveFailure("reference residual tolerance must be finite and non-negative");
  }

  auto dense = CsrToDense(assembly.matrix);
  std::vector<double> solution;
  if (!SolveDense(std::move(dense), assembly.rhs, &solution)) {
    return SolveFailure("serial dense reference solve failed");
  }

  for (double value : solution) {
    if (!std::isfinite(value)) {
      return SolveFailure("serial dense reference solution is non-finite");
    }
  }

  GenericDiffusionSolveResult result;
  result.success = true;
  result.backend = "serial_dense_reference";
  result.scalar_new = std::move(solution);
  result.row_count = assembly.matrix.row_count;

  ComputeResiduals(
      assembly.matrix,
      result.scalar_new,
      assembly.rhs,
      &result.residual_l2,
      &result.residual_linf,
      &result.rhs_linf);
  result.residual_linf_relative =
      result.residual_linf / std::max(1.0, result.rhs_linf);
  if (!std::isfinite(result.residual_l2) ||
      !std::isfinite(result.residual_linf) ||
      !std::isfinite(result.rhs_linf) ||
      !std::isfinite(result.residual_linf_relative) ||
      result.residual_linf_relative > options.residual_tolerance) {
    return SolveResidualFailure(std::move(result), options);
  }

  result.max_abs_delta = 0.0;
  for (std::size_t row = 0; row < result.scalar_new.size(); ++row) {
    result.max_abs_delta = std::max(
        result.max_abs_delta,
        std::abs(result.scalar_new[row] - assembly.scalar_old_flat[row]));
  }

  result.report_line = BuildSolveReport(result, options);
  return result;
}

double GetCsrValue(
    const SparseMatrixCsr& matrix,
    std::size_t row,
    std::size_t column) noexcept {
  if (!matrix.is_shape_complete() || row >= matrix.row_count || column >= matrix.column_count) {
    return 0.0;
  }

  for (std::size_t slot = matrix.row_offsets[row]; slot < matrix.row_offsets[row + 1u]; ++slot) {
    if (matrix.column_indices[slot] == column) {
      return matrix.values[slot];
    }
  }
  return 0.0;
}

bool GenericDiffusionAssemblyResult::is_complete() const noexcept {
  return success && matrix.is_shape_complete() && rhs.size() == row_count &&
         scalar_old_flat.size() == row_count && !report_line.empty() &&
         ValidateGenericDiffusionAssemblyDiagnostics(*this);
}

bool GenericDiffusionSolveResult::is_complete() const noexcept {
  return success && backend == "serial_dense_reference" && !scalar_new.empty() &&
         !report_line.empty() && ValidateGenericDiffusionSolveDiagnostics(*this);
}

bool ValidateGenericDiffusionAssemblyDiagnostics(
    const GenericDiffusionAssemblyResult& result) noexcept {
  if (!result.success) {
    return false;
  }

  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p2.diffusion.assembly") &&
         Contains(line, "implementation_id=p2.diffusion.generic_implicit_assembly_v1") &&
         Contains(line, "equation_form=A_dT_dt_div_D_grad_T_plus_C_T_plus_B") &&
         Contains(line, "unit_system=cgs") &&
         Contains(line, "dt_s=") &&
         Contains(line, "unknown=generic_scalar") &&
         Contains(line, "direction_splitting=false") &&
         Contains(line, "matrix_format=csr") &&
         Contains(line, "interface_diffusion_mean=arithmetic") &&
         Contains(line, "pole_metric_mode=") &&
         Contains(line, "row_count=") &&
         Contains(line, "nonzero_count=") &&
         Contains(line, "boundary_row_count=") &&
         Contains(line, "origin_boundary_touched=") &&
         Contains(line, "pole_boundary_touched=") &&
         Contains(line, "outer_boundary_policy=") &&
         Contains(line, "phi_boundary_policy=periodic") &&
         Contains(line, "scalar_origin_remap_available=") &&
         Contains(line, "scalar_pole_remap_available=") &&
         Contains(line, "scalar_remap_value_transform=identity") &&
         Contains(line, "origin_remap_used=") &&
         Contains(line, "pole_remap_used=") &&
         Contains(line, "origin_remap_coupling_mode=") &&
         Contains(line, "pole_remap_coupling_mode=") &&
         Contains(line, "origin_remap_unique_pair_count=") &&
         Contains(line, "origin_remap_row_coupling_count=") &&
         Contains(line, "pole_remap_unique_pair_count=") &&
         Contains(line, "pole_remap_row_coupling_count=") &&
         Contains(line, "origin_remap_nonzero_offdiag_count=") &&
         Contains(line, "pole_remap_nonzero_offdiag_count=") &&
         Contains(line, "origin_remap_min_abs_offdiag=") &&
         Contains(line, "pole_remap_min_abs_offdiag=") &&
         Contains(line, "singular_boundary_face_conductance_used=false") &&
         Contains(line, "marshak_boundary_used=") &&
         Contains(line, "marshak_outer_face_count=") &&
         Contains(line, "marshak_nonzero_diagonal_loss_count=") &&
         Contains(line, "marshak_min_G_face=") &&
         Contains(line, "marshak_max_G_face=") &&
         Contains(line, "marshak_min_center_to_boundary_distance_cm=") &&
         Contains(line, "marshak_max_center_to_boundary_distance_cm=") &&
         Contains(line, "marshak_min_abs_diagonal_loss=") &&
         Contains(line, "marshak_max_abs_diagonal_loss=") &&
         Contains(line, "min_A=") &&
         Contains(line, "min_D=") &&
         Contains(line, "max_D=") &&
         Contains(line, "min_diagonal=") &&
         Contains(line, "max_abs_offdiagonal_row_sum=") &&
         Contains(line, "max_rhs_abs=");
}

bool ValidateGenericDiffusionSolveDiagnostics(
    const GenericDiffusionSolveResult& result) noexcept {
  if (!result.success) {
    return false;
  }

  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p2.diffusion.reference_solve") &&
         Contains(line, "implementation_id=p2.diffusion.serial_dense_reference_v1") &&
         Contains(line, "backend=serial_dense_reference") &&
         Contains(line, "row_count=") &&
         Contains(line, "row_limit=") &&
         Contains(line, "residual_l2=") &&
         Contains(line, "residual_linf=") &&
         Contains(line, "rhs_linf=") &&
         Contains(line, "residual_linf_relative=") &&
         Contains(line, "max_abs_delta=") &&
         Contains(line, "canonical_state_mutated=false");
}

}  // namespace dec3d::transport
