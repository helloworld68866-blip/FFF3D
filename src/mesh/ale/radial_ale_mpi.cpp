#include "mesh/ale/radial_ale_mpi.hpp"

#include "mesh/ale/radial_ale.hpp"

#include <sstream>
#include <utility>

namespace dec3d::mesh {

bool RadialAleMpiGlobalProposalResult::is_complete() const noexcept {
  return success &&
         rank_count > 0 &&
         !global_shell_velocity.empty() &&
         proposal.is_complete_global(global_shell_velocity.size() + 1u) &&
         failure_reason.empty() &&
         !report_line.empty();
}

bool RadialAleMpiCommitResult::is_complete() const noexcept {
  return !report_line.empty();
}

RadialAleMpiGlobalProposalResult BuildRadialAleMpiGlobalProposal(
    const dec3d::core::Array3D<double>& local_rho,
    const dec3d::core::Array3D<double>& local_mom_r,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const std::vector<double>& global_radial_faces,
    double dt_s,
    double beta,
    MPI_Comm communicator) noexcept {
  RadialAleMpiGlobalProposalResult result;
  MPI_Comm_rank(communicator, &result.rank);
  MPI_Comm_size(communicator, &result.rank_count);

  const auto finish = [&]() {
    std::ostringstream report;
    report << "rank=" << result.rank
           << "; rank_count=" << result.rank_count
           << "; local_shell_velocity_count=" << result.local_shell_velocity.size()
           << "; global_shell_velocity_count=" << result.global_shell_velocity.size()
           << "; success=" << (result.success ? "true" : "false");
    if (!result.failure_reason.empty()) {
      report << "; failure_reason=" << result.failure_reason;
    }
    result.report_line = report.str();
    return result;
  };

  if (!decomposition.is_valid() ||
      decomposition.slices.size() != static_cast<std::size_t>(result.rank_count) ||
      result.rank < 0 ||
      static_cast<std::size_t>(result.rank) >= decomposition.slices.size()) {
    result.failure_reason =
        "radial ALE MPI proposal requires a valid rank-matched radial decomposition";
    return finish();
  }

  const auto& local_slice = decomposition.slices[static_cast<std::size_t>(result.rank)];
  result.local_shell_velocity = BuildRadialAleShellVelocitySlice(local_rho, local_mom_r);
  if (result.local_shell_velocity.size() != local_slice.local_cell_count()) {
    result.failure_reason = "local shell velocity count does not match radial ownership slice";
    return finish();
  }

  std::vector<int> recvcounts(decomposition.slices.size(), 0);
  std::vector<int> displacements(decomposition.slices.size(), 0);
  for (std::size_t slice_index = 0; slice_index < decomposition.slices.size(); ++slice_index) {
    const auto& slice = decomposition.slices[slice_index];
    recvcounts[slice_index] = static_cast<int>(slice.local_cell_count());
    displacements[slice_index] = static_cast<int>(slice.begin_index);
  }

  result.global_shell_velocity.assign(decomposition.global_radial_cells, 0.0);
  MPI_Allgatherv(
      result.local_shell_velocity.data(),
      static_cast<int>(result.local_shell_velocity.size()),
      MPI_DOUBLE,
      result.global_shell_velocity.data(),
      recvcounts.data(),
      displacements.data(),
      MPI_DOUBLE,
      communicator);

  result.proposal = BuildRadialAleProposalFromGlobalShellVelocity(
      result.global_shell_velocity,
      global_radial_faces,
      dt_s,
      beta);
  if (!result.proposal.is_complete_global(global_radial_faces.size())) {
    result.failure_reason = result.proposal.summary.empty()
                                ? "global radial ALE proposal construction failed"
                                : result.proposal.summary;
    return finish();
  }

  result.success = true;
  return finish();
}

RadialAleMpiCommitResult ApplyRadialAleMpiGlobalProposalToLocalGeometry(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    MPI_Comm communicator) noexcept {
  RadialAleMpiCommitResult result;
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &rank_count);

  const auto finish = [&]() {
    std::ostringstream report;
    report << "rank=" << rank
           << "; rank_count=" << rank_count
           << "; local_preflight_passed=" << (result.local_preflight_passed ? "true" : "false")
           << "; global_preflight_passed=" << (result.global_preflight_passed ? "true" : "false")
           << "; local_stage_attempted=" << (result.local_stage_attempted ? "true" : "false")
           << "; local_stage_success=" << (result.local_stage_success ? "true" : "false")
           << "; global_stage_success=" << (result.global_stage_success ? "true" : "false")
           << "; published=" << (result.published ? "true" : "false")
           << "; success=" << (result.success ? "true" : "false");
    if (!result.failure_reason.empty()) {
      report << "; failure_reason=" << result.failure_reason;
    }
    result.report_line = report.str();
    return result;
  };

  if (!decomposition.is_valid() ||
      decomposition.slices.size() != static_cast<std::size_t>(rank_count) ||
      rank < 0 ||
      static_cast<std::size_t>(rank) >= decomposition.slices.size()) {
    result.failure_reason = "radial ALE MPI commit requires a valid rank-matched decomposition";
  } else if (!proposal.is_complete_global(decomposition.global_radial_cells + 1u)) {
    result.failure_reason = "radial ALE MPI commit requires a complete global proposal";
  } else {
    bool monotonic = true;
    for (std::size_t face = 1u; face < proposal.proposed_radial_faces.size(); ++face) {
      if (!(proposal.proposed_radial_faces[face] > proposal.proposed_radial_faces[face - 1u])) {
        monotonic = false;
        break;
      }
    }

    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
    const bool has_window = proposal.has_radial_face_window(
        slice.begin_index,
        local_geometry.radial_faces.size());
    result.local_preflight_passed =
        monotonic && has_window && local_geometry.is_valid();
    if (!result.local_preflight_passed) {
      result.failure_reason = "radial ALE MPI commit preflight failed";
    }
  }

  const int local_ready = result.local_preflight_passed ? 1 : 0;
  int global_ready = 0;
  MPI_Allreduce(&local_ready, &global_ready, 1, MPI_INT, MPI_MIN, communicator);
  result.global_preflight_passed = global_ready == 1;
  if (!result.global_preflight_passed) {
    if (result.failure_reason.empty()) {
      result.failure_reason = "radial ALE MPI commit global preflight failed";
    }
    return finish();
  }

  const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
  auto staged_geometry = local_geometry;
  result.local_stage_attempted = true;
  const auto staged_commit = ApplyRadialAleMeshUpdateProposal(
      proposal,
      staged_geometry,
      slice.begin_index);
  result.local_stage_success = staged_commit.success;
  const int local_stage_ok = staged_commit.success ? 1 : 0;
  int global_stage_ok = 0;
  MPI_Allreduce(&local_stage_ok, &global_stage_ok, 1, MPI_INT, MPI_MIN, communicator);
  result.global_stage_success = global_stage_ok == 1;
  if (!result.global_stage_success) {
    result.failure_reason = staged_commit.failure_reason.empty()
                                ? "radial ALE MPI staged commit failed"
                                : staged_commit.failure_reason;
    return finish();
  }

  local_geometry = std::move(staged_geometry);
  result.published = true;
  result.success = true;
  return finish();
}

}  // namespace dec3d::mesh
