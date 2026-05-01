#include "hydro/driver/macro_ale_mpi_compatibility.hpp"

#include "mesh/ale/radial_ale.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace dec3d::hydro {
namespace {

constexpr double kA3SeamTolerance = 1.0e-14;
constexpr double kA3ConservationTolerance = 1.0e-10;
constexpr double kA3CoverageTolerance = 1.0e-10;

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] bool SameRadialFaces(
    const dec3d::mesh::SphericalGeometryMetadata& lhs,
    const dec3d::mesh::SphericalGeometryMetadata& rhs) noexcept {
  return lhs.radial_faces == rhs.radial_faces;
}

[[nodiscard]] std::string BuildReportLine(
    const MacroAleMpiCompatibilityDiagnostics& diagnostics) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "rank=" << diagnostics.rank
         << "; rank_count=" << diagnostics.rank_count
         << "; global_radial_cells=" << diagnostics.global_radial_cells
         << "; global_face_begin=" << diagnostics.global_face_begin
         << "; local_face_count=" << diagnostics.local_face_count
         << "; local_radial_cells=" << diagnostics.local_radial_cells
         << "; radial_face_indexing=global"
         << "; remap_order=first_order_proposal_mapped_overlap"
         << "; local_mass_residual=" << diagnostics.local_mass_residual
         << "; local_mom_r_residual=" << diagnostics.local_mom_r_residual
         << "; local_mom_theta_residual=" << diagnostics.local_mom_theta_residual
         << "; local_mom_phi_residual=" << diagnostics.local_mom_phi_residual
         << "; local_e_fluid_total_residual="
         << diagnostics.local_e_fluid_total_residual
         << "; local_chi_e_residual=" << diagnostics.local_chi_e_residual
         << "; global_mass_residual=" << diagnostics.global_mass_residual
         << "; global_mom_r_residual=" << diagnostics.global_mom_r_residual
         << "; global_mom_theta_residual=" << diagnostics.global_mom_theta_residual
         << "; global_mom_phi_residual=" << diagnostics.global_mom_phi_residual
         << "; global_e_fluid_total_residual="
         << diagnostics.global_e_fluid_total_residual
         << "; global_chi_e_residual=" << diagnostics.global_chi_e_residual
         << "; max_seam_face_velocity_delta="
         << diagnostics.max_seam_face_velocity_delta
         << "; max_seam_proposed_face_delta="
         << diagnostics.max_seam_proposed_face_delta
         << "; max_seam_preview_radius_delta="
         << diagnostics.max_seam_preview_radius_delta
         << "; global_stage_ok="
         << (diagnostics.global_stage_ok ? "true" : "false")
         << "; global_publish_ok="
         << (diagnostics.global_publish_ok ? "true" : "false")
         << "; canonical_writeback_published="
         << (diagnostics.canonical_writeback_published ? "true" : "false")
         << "; mesh_commit_published="
         << (diagnostics.mesh_commit_published ? "true" : "false")
         << "; published_from_staged_geometry="
         << (diagnostics.published_from_staged_geometry ? "true" : "false");
  return report.str();
}

void AppendResultDiagnostics(
    MacroAleMpiCompatibilityResult& result) {
  const auto report = result.diagnostics.report_line;
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_mpi.global_transaction.preflight",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_mpi.seam_face_velocity_equality",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_mpi.seam_proposed_face_equality",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_mpi.seam_preview_geometry_equality",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_mpi.global_conservation",
      report);
  if (result.diagnostics.canonical_writeback_published &&
      result.diagnostics.mesh_commit_published) {
    AppendDiagnostic(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_mpi.transaction.clean",
        report);
  }
}

[[nodiscard]] double MaxSharedFaceDelta(
    const std::vector<double>& local_values,
    const std::vector<int>& local_has_values,
    MPI_Comm communicator) {
  const std::size_t face_count = local_values.size();
  const double positive_infinity = std::numeric_limits<double>::infinity();
  std::vector<double> local_min(face_count, positive_infinity);
  std::vector<double> local_max(face_count, -positive_infinity);
  std::vector<int> local_count(face_count, 0);
  for (std::size_t face = 0; face < face_count; ++face) {
    if (local_has_values[face] != 0) {
      local_min[face] = local_values[face];
      local_max[face] = local_values[face];
      local_count[face] = 1;
    }
  }

  std::vector<double> global_min(face_count, positive_infinity);
  std::vector<double> global_max(face_count, -positive_infinity);
  std::vector<int> global_count(face_count, 0);
  MPI_Allreduce(
      local_min.data(),
      global_min.data(),
      static_cast<int>(face_count),
      MPI_DOUBLE,
      MPI_MIN,
      communicator);
  MPI_Allreduce(
      local_max.data(),
      global_max.data(),
      static_cast<int>(face_count),
      MPI_DOUBLE,
      MPI_MAX,
      communicator);
  MPI_Allreduce(
      local_count.data(),
      global_count.data(),
      static_cast<int>(face_count),
      MPI_INT,
      MPI_SUM,
      communicator);

  double max_delta = 0.0;
  for (std::size_t face = 1u; face + 1u < face_count; ++face) {
    if (global_count[face] == 2) {
      max_delta = std::max(
          max_delta,
          std::abs(global_max[face] - global_min[face]));
    }
  }
  return max_delta;
}

}  // namespace

bool MacroAleMpiCompatibilityDiagnostics::is_complete() const noexcept {
  return executed &&
         rank_count > 0 &&
         local_face_count == local_radial_cells + 1u &&
         Finite(local_mass_residual) &&
         Finite(local_mom_r_residual) &&
         Finite(local_mom_theta_residual) &&
         Finite(local_mom_phi_residual) &&
         Finite(local_e_fluid_total_residual) &&
         Finite(local_chi_e_residual) &&
         Finite(global_mass_residual) &&
         Finite(global_mom_r_residual) &&
         Finite(global_mom_theta_residual) &&
         Finite(global_mom_phi_residual) &&
         Finite(global_e_fluid_total_residual) &&
         Finite(global_chi_e_residual) &&
         Finite(max_seam_face_velocity_delta) &&
         Finite(max_seam_proposed_face_delta) &&
         Finite(max_seam_preview_radius_delta) &&
         !report_line.empty();
}

bool MacroAleMpiCompatibilityResult::is_complete() const noexcept {
  return success &&
         local_stage_ok &&
         global_stage_ok &&
         local_publishable &&
         global_publish_ok &&
         canonical_writeback_published &&
         mesh_commit_published &&
         staged_hydro_result.success &&
         diagnostics.is_complete() &&
         diagnostics_payload.has_entries() &&
         failure_reason.empty();
}

MacroAleMpiCompatibilityResult AdvanceMacroAleMpiCompatibilityStep(
    dec3d::state::HydroStateView& local_hydro_view,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const dec3d::core::MeshUpdateProposal& global_radial_ale_proposal,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    StaticGridHydroOptions options,
    MPI_Comm communicator) noexcept {
  MacroAleMpiCompatibilityResult result;
  MPI_Comm_rank(communicator, &result.diagnostics.rank);
  MPI_Comm_size(communicator, &result.diagnostics.rank_count);
  const int rank = result.diagnostics.rank;
  const int rank_count = result.diagnostics.rank_count;
  result.diagnostics.executed = true;
  result.diagnostics.global_radial_cells =
      global_radial_ale_proposal.global_radial_face_count > 0u
          ? global_radial_ale_proposal.global_radial_face_count - 1u
          : 0u;

  bool local_preflight_ok = false;
  std::size_t global_face_begin = 0u;
  std::size_t local_face_count = local_geometry.radial_faces.size();
  std::size_t local_radial_cells = 0u;
  std::string local_failure;

  if (!decomposition.is_valid()) {
    local_failure = decomposition.failure_reason.empty()
                        ? "macro ALE MPI requires a valid radial decomposition"
                        : decomposition.failure_reason;
  } else if (rank < 0 ||
             static_cast<std::size_t>(rank) >= decomposition.slices.size() ||
             rank_count != static_cast<int>(decomposition.slices.size())) {
    local_failure = "macro ALE MPI rank count does not match decomposition";
  } else {
    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
    global_face_begin = slice.begin_index;
    local_radial_cells = slice.local_cell_count();
    result.diagnostics.global_face_begin = global_face_begin;
    result.diagnostics.local_face_count = local_face_count;
    result.diagnostics.local_radial_cells = local_radial_cells;

    if (!local_hydro_view.is_complete() || !local_geometry.is_valid()) {
      local_failure = "macro ALE MPI requires complete local hydro and geometry";
    } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      local_failure = "macro ALE MPI requires a finite positive timestep";
    } else if (!options.use_macro_zoning ||
               !options.use_ppm_reconstruction ||
               !options.apply_radial_ale_flux_correction) {
      local_failure = "macro ALE MPI requires macro zoning, PPM, and ALE enabled";
    } else if (!global_radial_ale_proposal.radial_face_indexing_is_global ||
               !global_radial_ale_proposal.has_radial_face_window(
                   global_face_begin,
                   local_face_count)) {
      local_failure = "macro ALE MPI proposal/window failure";
    } else if (local_face_count != local_radial_cells + 1u) {
      local_failure = "macro ALE MPI local geometry is not the declared face window";
    } else {
      local_preflight_ok = true;
    }
  }

  dec3d::mesh::SphericalGeometryMetadata staged_geometry = local_geometry;
  bool preview_ok = false;
  if (local_preflight_ok) {
    const auto preview = dec3d::mesh::ApplyRadialAleMeshUpdateProposal(
        global_radial_ale_proposal,
        staged_geometry,
        global_face_begin);
    preview_ok = preview.success;
    if (!preview_ok) {
      local_failure = preview.failure_reason.empty()
                          ? "macro ALE MPI geometry preview failed"
                          : preview.failure_reason;
    }
  }

  auto local_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
      global_radial_ale_proposal,
      global_face_begin,
      local_face_count);
  if (local_window.is_complete()) {
    AppendDiagnostic(
        result.diagnostics_payload,
        "p1.hydro.macro_ale.local_window",
        local_window.report_line);
  }

  options.radial_ale_global_face_begin_index = global_face_begin;
  options.defer_macro_ale_compatibility_writeback = true;
  if (local_preflight_ok && preview_ok) {
    result.staged_hydro_result = AdvanceStaticGridHydro(
        local_hydro_view,
        local_geometry,
        dt_s,
        radial_ghost_override,
        options,
        &global_radial_ale_proposal);
    result.diagnostics_payload.entries.insert(
        result.diagnostics_payload.entries.end(),
        result.staged_hydro_result.diagnostics.entries.begin(),
        result.staged_hydro_result.diagnostics.entries.end());
    if (!result.staged_hydro_result.success) {
      local_failure =
          result.staged_hydro_result.failure_reason.empty()
              ? "macro ALE MPI staged hydro failed"
              : result.staged_hydro_result.failure_reason;
    }
  }

  const bool staged_geometry_matches =
      result.staged_hydro_result.macro_ale_staged_geometry_available &&
      SameRadialFaces(result.staged_hydro_result.macro_ale_staged_geometry, staged_geometry);

  if (result.staged_hydro_result.macro_ale_remap_diagnostics.is_complete()) {
    const auto& remap_diag =
        result.staged_hydro_result.macro_ale_remap_diagnostics;
    result.diagnostics.local_mass_residual = remap_diag.mass_residual;
    result.diagnostics.local_mom_r_residual = remap_diag.mom_r_residual;
    result.diagnostics.local_mom_theta_residual = remap_diag.mom_theta_residual;
    result.diagnostics.local_mom_phi_residual = remap_diag.mom_phi_residual;
    result.diagnostics.local_e_fluid_total_residual =
        remap_diag.e_fluid_total_residual;
    result.diagnostics.local_chi_e_residual = remap_diag.chi_e_residual;
  }

  const std::size_t global_face_count =
      std::max<std::size_t>(global_radial_ale_proposal.global_radial_face_count, 1u);
  std::vector<double> local_velocity_values(global_face_count, 0.0);
  std::vector<double> local_proposed_values(global_face_count, 0.0);
  std::vector<double> local_preview_values(global_face_count, 0.0);
  std::vector<int> local_has_values(global_face_count, 0);
  if (local_preflight_ok && preview_ok && local_window.is_complete()) {
    for (std::size_t local_face = 0u; local_face < local_face_count; ++local_face) {
      const std::size_t global_face = global_face_begin + local_face;
      if (global_face > 0u && global_face + 1u < global_face_count) {
        local_has_values[global_face] = 1;
        local_velocity_values[global_face] =
            global_radial_ale_proposal.radial_face_velocities[global_face];
        local_proposed_values[global_face] =
            global_radial_ale_proposal.proposed_radial_faces[global_face];
        local_preview_values[global_face] = staged_geometry.radial_faces[local_face];
      }
    }
  }
  result.diagnostics.max_seam_face_velocity_delta =
      MaxSharedFaceDelta(local_velocity_values, local_has_values, communicator);
  result.diagnostics.max_seam_proposed_face_delta =
      MaxSharedFaceDelta(local_proposed_values, local_has_values, communicator);
  result.diagnostics.max_seam_preview_radius_delta =
      MaxSharedFaceDelta(local_preview_values, local_has_values, communicator);

  const double local_abs_residuals[] = {
      std::abs(result.diagnostics.local_mass_residual),
      std::abs(result.diagnostics.local_mom_r_residual),
      std::abs(result.diagnostics.local_mom_theta_residual),
      std::abs(result.diagnostics.local_mom_phi_residual),
      std::abs(result.diagnostics.local_e_fluid_total_residual),
      std::abs(result.diagnostics.local_chi_e_residual)};
  double global_abs_residuals[6] = {};
  MPI_Allreduce(
      local_abs_residuals,
      global_abs_residuals,
      6,
      MPI_DOUBLE,
      MPI_SUM,
      communicator);
  result.diagnostics.global_mass_residual = global_abs_residuals[0];
  result.diagnostics.global_mom_r_residual = global_abs_residuals[1];
  result.diagnostics.global_mom_theta_residual = global_abs_residuals[2];
  result.diagnostics.global_mom_phi_residual = global_abs_residuals[3];
  result.diagnostics.global_e_fluid_total_residual = global_abs_residuals[4];
  result.diagnostics.global_chi_e_residual = global_abs_residuals[5];

  const bool seam_diagnostics_finite =
      Finite(result.diagnostics.max_seam_face_velocity_delta) &&
      Finite(result.diagnostics.max_seam_proposed_face_delta) &&
      Finite(result.diagnostics.max_seam_preview_radius_delta);
  const bool local_residuals_finite =
      Finite(result.diagnostics.local_mass_residual) &&
      Finite(result.diagnostics.local_mom_r_residual) &&
      Finite(result.diagnostics.local_mom_theta_residual) &&
      Finite(result.diagnostics.local_mom_phi_residual) &&
      Finite(result.diagnostics.local_e_fluid_total_residual) &&
      Finite(result.diagnostics.local_chi_e_residual);
  const bool seam_within_tolerance =
      result.diagnostics.max_seam_face_velocity_delta <= kA3SeamTolerance &&
      result.diagnostics.max_seam_proposed_face_delta <= kA3SeamTolerance &&
      result.diagnostics.max_seam_preview_radius_delta <= kA3SeamTolerance;
  const bool local_residuals_within_tolerance =
      local_abs_residuals[0] <= kA3ConservationTolerance &&
      local_abs_residuals[1] <= kA3ConservationTolerance &&
      local_abs_residuals[2] <= kA3ConservationTolerance &&
      local_abs_residuals[3] <= kA3ConservationTolerance &&
      local_abs_residuals[4] <= kA3ConservationTolerance &&
      local_abs_residuals[5] <= kA3ConservationTolerance;
  const bool coverage_within_tolerance =
      result.staged_hydro_result.macro_ale_remap_diagnostics.max_volume_coverage_error <=
      kA3CoverageTolerance;
  const bool global_residuals_within_tolerance =
      result.diagnostics.global_mass_residual <= kA3ConservationTolerance &&
      result.diagnostics.global_mom_r_residual <= kA3ConservationTolerance &&
      result.diagnostics.global_mom_theta_residual <= kA3ConservationTolerance &&
      result.diagnostics.global_mom_phi_residual <= kA3ConservationTolerance &&
      result.diagnostics.global_e_fluid_total_residual <= kA3ConservationTolerance &&
      result.diagnostics.global_chi_e_residual <= kA3ConservationTolerance;

  if (local_failure.empty()) {
    if (!seam_within_tolerance) {
      local_failure = "macro ALE MPI seam equality tolerance exceeded";
    } else if (!local_residuals_within_tolerance) {
      local_failure = "macro ALE MPI local conservation tolerance exceeded";
    } else if (!coverage_within_tolerance) {
      local_failure = "macro ALE MPI remap coverage tolerance exceeded";
    } else if (!global_residuals_within_tolerance) {
      local_failure = "macro ALE MPI global conservation tolerance exceeded";
    }
  }

  const int local_stage_ok =
      local_preflight_ok &&
              preview_ok &&
              result.staged_hydro_result.success &&
              result.staged_hydro_result.macro_ale_staged_writeback_available &&
              result.staged_hydro_result.macro_ale_staged_geometry_available &&
              staged_geometry_matches &&
              seam_diagnostics_finite &&
              local_residuals_finite &&
              seam_within_tolerance &&
              local_residuals_within_tolerance &&
              coverage_within_tolerance &&
              global_residuals_within_tolerance
          ? 1
          : 0;
  int global_stage_ok = 0;
  MPI_Allreduce(
      &local_stage_ok,
      &global_stage_ok,
      1,
      MPI_INT,
      MPI_MIN,
      communicator);
  result.local_stage_ok = local_stage_ok == 1;
  result.global_stage_ok = global_stage_ok == 1;
  result.diagnostics.global_stage_ok = result.global_stage_ok;

  const bool local_publishable =
      result.global_stage_ok &&
      local_hydro_view.is_complete() &&
      result.staged_hydro_result.success &&
      result.staged_hydro_result.macro_ale_staged_writeback_available &&
      result.staged_hydro_result.macro_ale_staged_geometry_available &&
      !result.staged_hydro_result.macro_ale_staged_cells.empty() &&
      result.staged_hydro_result.macro_ale_staged_cells.extent_r() ==
          local_hydro_view.rho->extent_r() &&
      result.staged_hydro_result.macro_ale_staged_cells.extent_theta() ==
          local_hydro_view.rho->extent_theta() &&
      result.staged_hydro_result.macro_ale_staged_cells.extent_phi() ==
          local_hydro_view.rho->extent_phi() &&
      staged_geometry_matches;
  const int local_publish_ok = local_publishable ? 1 : 0;
  int global_publish_ok = 0;
  MPI_Allreduce(
      &local_publish_ok,
      &global_publish_ok,
      1,
      MPI_INT,
      MPI_MIN,
      communicator);
  result.local_publishable = local_publishable;
  result.global_publish_ok = global_publish_ok == 1;
  result.diagnostics.global_publish_ok = result.global_publish_ok;

  bool hydro_published = false;
  if (result.global_stage_ok && result.global_publish_ok) {
    std::string publish_failure;
    hydro_published = PublishMacroAleStagedHydroState(
        local_hydro_view,
        result.staged_hydro_result,
        &publish_failure);
    if (hydro_published) {
      local_geometry = result.staged_hydro_result.macro_ale_staged_geometry;
    } else {
      local_failure = publish_failure.empty()
                          ? "macro ALE MPI publish failed after global publish preflight"
                          : publish_failure;
    }
  } else if (result.global_stage_ok && !result.global_publish_ok) {
    local_failure = "macro ALE MPI global publish preflight failed";
  }

  const int local_published_ok = hydro_published ? 1 : 0;
  int global_published_ok = 0;
  MPI_Allreduce(
      &local_published_ok,
      &global_published_ok,
      1,
      MPI_INT,
      MPI_MIN,
      communicator);
  result.canonical_writeback_published = hydro_published && global_published_ok == 1;
  result.mesh_commit_published = result.canonical_writeback_published;
  result.diagnostics.canonical_writeback_published =
      result.canonical_writeback_published;
  result.diagnostics.mesh_commit_published = result.mesh_commit_published;
  result.diagnostics.published_from_staged_geometry =
      result.canonical_writeback_published && result.mesh_commit_published;

  if (result.global_stage_ok && result.global_publish_ok && global_published_ok != 1) {
    local_failure = "macro ALE MPI publish failed after global publish preflight";
  }
  if (!result.global_stage_ok && local_failure.empty()) {
    local_failure = "macro ALE MPI global stage failed";
  }

  result.success =
      result.global_stage_ok &&
      result.global_publish_ok &&
      result.canonical_writeback_published &&
      result.mesh_commit_published &&
      local_failure.empty();
  result.failure_reason = result.success ? std::string{} : local_failure;
  result.diagnostics.report_line = BuildReportLine(result.diagnostics);
  AppendResultDiagnostics(result);
  return result;
}

}  // namespace dec3d::hydro
