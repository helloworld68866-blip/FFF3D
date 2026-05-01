#include "core/array/array3d.hpp"
#include "hydro_mpi_real_case_support.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/ale/radial_ale_mpi.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kRadialCells = 10u;
constexpr std::size_t kThetaCells = 3u;
constexpr std::size_t kPhiCells = 4u;
constexpr double kDt = 1.0e-4;
constexpr double kBeta = 0.9;
constexpr double kTolerance = 1.0e-14;

[[nodiscard]] bool Near(double lhs, double rhs) noexcept {
  return std::abs(lhs - rhs) <= kTolerance;
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildGeometry() {
  dec3d::mesh::SphericalMeshDescriptor descriptor;
  descriptor.radial_cells = kRadialCells;
  descriptor.theta_cells = kThetaCells;
  descriptor.phi_cells = kPhiCells;
  descriptor.inner_radius = 0.0;
  descriptor.outer_radius = 1.0;
  return dec3d::mesh::BuildSphericalGeometry(descriptor);
}

[[nodiscard]] dec3d::state::CanonicalState BuildGlobalState() {
  auto state = dec3d::state::CanonicalState::Create(
      {kRadialCells, kThetaCells, kPhiCells, 1u});
  for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
    const double velocity = -0.04 - 0.003 * static_cast<double>(radial);
    for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
      for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
        const double rho =
            1.0 + 0.02 * static_cast<double>(radial) +
            0.001 * static_cast<double>(theta) +
            0.0001 * static_cast<double>(phi);
        state.rho(radial, theta, phi) = rho;
        state.mom_r(radial, theta, phi) = rho * velocity;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = 1.0;
        state.e_electron(radial, theta, phi) = 0.5;
      }
    }
  }
  return state;
}

[[nodiscard]] double FaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t local_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double r = geometry.radial_faces[local_face];
  const double theta_factor =
      std::cos(geometry.theta_faces[theta]) -
      std::cos(geometry.theta_faces[theta + 1u]);
  const double phi_factor = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return r * r * theta_factor * phi_factor;
}

[[nodiscard]] bool CheckDenseSharedFaceEquality(
    const std::vector<double>& local_values,
    const std::vector<int>& local_has_values,
    MPI_Comm communicator,
    const char* label) {
  const std::size_t global_face_count = local_values.size();
  const double positive_infinity = std::numeric_limits<double>::infinity();
  std::vector<double> local_min(global_face_count, positive_infinity);
  std::vector<double> local_max(global_face_count, -positive_infinity);
  std::vector<int> local_count(global_face_count, 0);
  for (std::size_t face = 0; face < global_face_count; ++face) {
    if (local_has_values[face] != 0) {
      local_min[face] = local_values[face];
      local_max[face] = local_values[face];
      local_count[face] = 1;
    }
  }

  std::vector<double> global_min(global_face_count, positive_infinity);
  std::vector<double> global_max(global_face_count, -positive_infinity);
  std::vector<int> global_count(global_face_count, 0);
  MPI_Allreduce(
      local_min.data(),
      global_min.data(),
      static_cast<int>(global_face_count),
      MPI_DOUBLE,
      MPI_MIN,
      communicator);
  MPI_Allreduce(
      local_max.data(),
      global_max.data(),
      static_cast<int>(global_face_count),
      MPI_DOUBLE,
      MPI_MAX,
      communicator);
  MPI_Allreduce(
      local_count.data(),
      global_count.data(),
      static_cast<int>(global_face_count),
      MPI_INT,
      MPI_SUM,
      communicator);

  bool ok = true;
  for (std::size_t face = 1u; face + 1u < global_face_count; ++face) {
    if (global_count[face] == 2 && !Near(global_min[face], global_max[face])) {
      int rank = 0;
      MPI_Comm_rank(communicator, &rank);
      if (rank == 0) {
        std::cerr << label << " mismatch at face " << face
                  << ": min=" << global_min[face]
                  << " max=" << global_max[face] << '\n';
      }
      ok = false;
    }
  }
  return ok;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

  bool ok = true;
  const auto geometry = BuildGeometry();
  const auto decomposition =
      dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
  if (!decomposition.is_valid()) {
    if (rank == 0) {
      std::cerr << decomposition.failure_reason << '\n';
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  const auto global_state = BuildGlobalState();
  const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
  auto local_geometry =
      dec3d::testsupport::BuildLocalGeometrySlice(
          geometry,
          slice,
          kThetaCells,
          kPhiCells);
  auto local_state =
      dec3d::testsupport::BuildLocalStateSlice(global_state, slice);

  const auto mpi_result = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
      local_state.rho,
      local_state.mom_r,
      decomposition,
      geometry.radial_faces,
      kDt,
      kBeta,
      MPI_COMM_WORLD);
  ok &= mpi_result.is_complete();

  const auto single_shell =
      dec3d::mesh::BuildRadialAleShellVelocitySlice(global_state.rho, global_state.mom_r);
  const auto single_proposal =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          single_shell,
          geometry.radial_faces,
          kDt,
          kBeta);

  ok &= single_proposal.is_complete_global(geometry.radial_faces.size());
  ok &= mpi_result.global_shell_velocity.size() == single_shell.size();
  for (std::size_t i = 0; i < single_shell.size(); ++i) {
    if (!Near(mpi_result.global_shell_velocity[i], single_shell[i])) {
      if (rank == 0) {
        std::cerr << "global_proposal_parity: shell velocity mismatch at " << i << '\n';
      }
      ok = false;
    }
  }
  for (std::size_t face = 0; face < geometry.radial_faces.size(); ++face) {
    if (!Near(
            mpi_result.proposal.radial_face_velocities[face],
            single_proposal.radial_face_velocities[face])) {
      if (rank == 0) {
        std::cerr << "global_proposal_parity: face velocity mismatch at " << face << '\n';
      }
      ok = false;
    }
    if (!Near(
            mpi_result.proposal.proposed_radial_faces[face],
            single_proposal.proposed_radial_faces[face])) {
      if (rank == 0) {
        std::cerr << "global_proposal_parity: proposed face mismatch at " << face << '\n';
      }
      ok = false;
    }
  }

  const auto local_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
      mpi_result.proposal,
      slice.begin_index,
      local_geometry.radial_faces.size());
  ok &= local_window.is_complete();

  const auto commit = dec3d::mesh::ApplyRadialAleMpiGlobalProposalToLocalGeometry(
      mpi_result.proposal,
      local_geometry,
      decomposition,
      MPI_COMM_WORLD);
  ok &= commit.success;

  const std::size_t global_face_count = geometry.radial_faces.size();
  std::vector<double> local_velocity_values(global_face_count, 0.0);
  std::vector<double> local_proposed_values(global_face_count, 0.0);
  std::vector<double> local_committed_values(global_face_count, 0.0);
  std::vector<double> local_area_values(global_face_count, 0.0);
  std::vector<int> local_has_values(global_face_count, 0);
  for (std::size_t local_face = 0; local_face < local_geometry.radial_faces.size(); ++local_face) {
    const std::size_t global_face = slice.begin_index + local_face;
    const bool has_face = global_face > 0u && global_face + 1u < global_face_count;
    if (has_face) {
      local_has_values[global_face] = 1;
      local_velocity_values[global_face] = local_window.radial_face_velocities[local_face];
      local_proposed_values[global_face] = local_window.proposed_radial_faces[local_face];
      local_committed_values[global_face] = local_geometry.radial_faces[local_face];
      local_area_values[global_face] = FaceArea(local_geometry, local_face, 1u, 2u);
    }
  }
  ok &= CheckDenseSharedFaceEquality(
      local_velocity_values,
      local_has_values,
      MPI_COMM_WORLD,
      "seam_face_velocity_equality");
  ok &= CheckDenseSharedFaceEquality(
      local_proposed_values,
      local_has_values,
      MPI_COMM_WORLD,
      "seam_proposed_face_equality");
  ok &= CheckDenseSharedFaceEquality(
      local_committed_values,
      local_has_values,
      MPI_COMM_WORLD,
      "seam_committed_face_equality");
  ok &= CheckDenseSharedFaceEquality(
      local_area_values,
      local_has_values,
      MPI_COMM_WORLD,
      "seam_geometry_continuity");

  for (std::size_t face = 1u; face < mpi_result.proposal.proposed_radial_faces.size(); ++face) {
    ok &= mpi_result.proposal.proposed_radial_faces[face] >
          mpi_result.proposal.proposed_radial_faces[face - 1u];
  }
  for (std::size_t face = 1u; face < local_geometry.radial_faces.size(); ++face) {
    ok &= local_geometry.radial_faces[face] > local_geometry.radial_faces[face - 1u];
  }
  for (const double volume : local_geometry.cell_volumes) {
    ok &= volume > 0.0 && std::isfinite(volume);
  }

  dec3d::core::MeshUpdateProposal bad_proposal = mpi_result.proposal;
  bad_proposal.proposed_radial_faces.pop_back();
  auto bad_local_geometry = local_geometry;
  const auto before_bad_commit = bad_local_geometry.radial_faces;
  const auto bad_commit = dec3d::mesh::ApplyRadialAleMpiGlobalProposalToLocalGeometry(
      bad_proposal,
      bad_local_geometry,
      decomposition,
      MPI_COMM_WORLD);
  const int local_failed_cleanly =
      (!bad_commit.success &&
       !bad_commit.local_stage_attempted &&
       !bad_commit.published &&
       bad_local_geometry.radial_faces == before_bad_commit) ? 1 : 0;
  int global_failed_cleanly = 0;
  MPI_Allreduce(
      &local_failed_cleanly,
      &global_failed_cleanly,
      1,
      MPI_INT,
      MPI_SUM,
      MPI_COMM_WORLD);
  ok &= global_failed_cleanly == rank_count;
  if (rank == 0 && global_failed_cleanly != rank_count) {
    std::cerr << "transaction_failure_no_partial_commit: invalid proposal mutated geometry\n";
  }

  int local_ok = ok ? 1 : 0;
  int global_ok = 0;
  MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  if (rank == 0 && global_ok != 1) {
    std::cerr << "A1 radial ALE MPI global proposal test failed\n";
  }

  MPI_Finalize();
  return global_ok == 1 ? EXIT_SUCCESS : EXIT_FAILURE;
}
