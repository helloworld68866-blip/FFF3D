#include "hydro/driver/macro_ale_hllc_mpi.hpp"

#include "mesh/ale/radial_ale.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace dec3d::hydro {
namespace {

constexpr double kA4SeamGeometryTolerance = 1.0e-14;
constexpr double kA4FluxTolerance = 1.0e-10;
constexpr double kA4EnergyFluxTolerance = 1.0e-9;
constexpr double kA4FluxRelativeTolerance = 1.0e-12;
constexpr double kA4ConservationTolerance = 1.0e-10;
constexpr double kA4EnergyConservationTolerance = 1.0e-8;
constexpr double kA4EnergyConservationRelativeTolerance = 1.0e-12;

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

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  const std::size_t index = ((radial * theta_cells) + theta) * phi_cells + phi;
  if (index >= geometry.cell_volumes.size()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return geometry.cell_volumes[index];
}

[[nodiscard]] double SumAbsEnergyExtensive(
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  if (hydro_view.e_fluid_total == nullptr) {
    return 0.0;
  }
  double sum = 0.0;
  const std::size_t theta_cells = hydro_view.e_fluid_total->extent_theta();
  const std::size_t phi_cells = hydro_view.e_fluid_total->extent_phi();
  for (std::size_t radial = 0; radial < hydro_view.e_fluid_total->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const double volume = CellVolume(geometry, radial, theta, phi, theta_cells, phi_cells);
        const double density = (*hydro_view.e_fluid_total)(radial, theta, phi);
        if (Finite(volume) && Finite(density)) {
          sum += std::abs(density * volume);
        }
      }
    }
  }
  return sum;
}

[[nodiscard]] double SumAbsEnergyExtensive(
    const dec3d::core::Array3D<HydroConservativeState>& cells,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  if (cells.empty()) {
    return 0.0;
  }
  double sum = 0.0;
  const std::size_t theta_cells = cells.extent_theta();
  const std::size_t phi_cells = cells.extent_phi();
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const double volume = CellVolume(geometry, radial, theta, phi, theta_cells, phi_cells);
        const double density = cells(radial, theta, phi).e_fluid_total;
        if (Finite(volume) && Finite(density)) {
          sum += std::abs(density * volume);
        }
      }
    }
  }
  return sum;
}

[[nodiscard]] double RadiationEnergyDensityFromPressureScalar(double chi_rad) noexcept {
  return (chi_rad >= 0.0 && std::isfinite(chi_rad))
             ? 3.0 * std::pow(chi_rad, 4.0 / 3.0)
             : std::numeric_limits<double>::quiet_NaN();
}

[[nodiscard]] double AlphaEnergyDensityFromPressureScalar(double chi_alpha) noexcept {
  return (chi_alpha >= 0.0 && std::isfinite(chi_alpha))
             ? 1.5 * std::pow(chi_alpha, 5.0 / 3.0)
             : std::numeric_limits<double>::quiet_NaN();
}

[[nodiscard]] dec3d::state::CanonicalState CloneHydroViewToCanonicalState(
    const dec3d::state::HydroStateView& hydro_view) {
  auto clone = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{
          hydro_view.rho->extent_r(),
          hydro_view.rho->extent_theta(),
          hydro_view.rho->extent_phi(),
          hydro_view.radiation_chi.size()});
  for (std::size_t radial = 0; radial < hydro_view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < hydro_view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < hydro_view.rho->extent_phi(); ++phi) {
        clone.rho(radial, theta, phi) = (*hydro_view.rho)(radial, theta, phi);
        clone.mom_r(radial, theta, phi) = (*hydro_view.mom_r)(radial, theta, phi);
        clone.mom_theta(radial, theta, phi) =
            (*hydro_view.mom_theta)(radial, theta, phi);
        clone.mom_phi(radial, theta, phi) =
            (*hydro_view.mom_phi)(radial, theta, phi);
        clone.e_fluid_total(radial, theta, phi) =
            (*hydro_view.e_fluid_total)(radial, theta, phi);
        clone.e_electron(radial, theta, phi) =
            (*hydro_view.e_electron)(radial, theta, phi);
        if (hydro_view.operator_local_alpha_chi) {
          clone.alpha_state.storage(radial, theta, phi) =
              AlphaEnergyDensityFromPressureScalar(
                  hydro_view.alpha_chi(radial, theta, phi));
        }
        for (std::size_t group = 0; group < hydro_view.radiation_chi.size(); ++group) {
          clone.radiation_groups[group](radial, theta, phi) =
              RadiationEnergyDensityFromPressureScalar(
                  hydro_view.radiation_chi[group](radial, theta, phi));
        }
      }
    }
  }
  return clone;
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
  for (std::size_t face = 0u; face < face_count; ++face) {
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

[[nodiscard]] double MaxSharedFaceScale(
    const std::vector<double>& local_values,
    const std::vector<int>& local_has_values,
    MPI_Comm communicator) {
  const std::size_t face_count = local_values.size();
  std::vector<double> local_abs_max(face_count, 0.0);
  std::vector<int> local_count(face_count, 0);
  for (std::size_t face = 0u; face < face_count; ++face) {
    if (local_has_values[face] != 0) {
      local_abs_max[face] = std::abs(local_values[face]);
      local_count[face] = 1;
    }
  }

  std::vector<double> global_abs_max(face_count, 0.0);
  std::vector<int> global_count(face_count, 0);
  MPI_Allreduce(
      local_abs_max.data(),
      global_abs_max.data(),
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

  double max_scale = 1.0;
  for (std::size_t face = 1u; face + 1u < face_count; ++face) {
    if (global_count[face] == 2) {
      max_scale = std::max(max_scale, global_abs_max[face]);
    }
  }
  return max_scale;
}

[[nodiscard]] bool WithinScaledTolerance(
    double delta,
    double scale,
    double absolute_tolerance) noexcept {
  return delta <= absolute_tolerance ||
         delta / std::max(1.0, scale) <= kA4FluxRelativeTolerance;
}

[[nodiscard]] bool PublishableA4StagedResult(
    const dec3d::state::HydroStateView& hydro_view,
    const StaticGridHydroResult& staged_result,
    const dec3d::mesh::SphericalGeometryMetadata& staged_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& preview_geometry) noexcept {
  return hydro_view.is_complete() &&
         staged_result.success &&
         staged_result.macro_ale_hllc_executed &&
         staged_result.macro_ale_hllc_staged_writeback_available &&
         staged_result.macro_ale_hllc_staged_geometry_available &&
         !staged_result.macro_ale_hllc_staged_cells.empty() &&
         staged_result.macro_ale_hllc_staged_cells.extent_r() ==
             hydro_view.rho->extent_r() &&
         staged_result.macro_ale_hllc_staged_cells.extent_theta() ==
             hydro_view.rho->extent_theta() &&
         staged_result.macro_ale_hllc_staged_cells.extent_phi() ==
             hydro_view.rho->extent_phi() &&
         SameRadialFaces(staged_geometry, preview_geometry);
}

void AppendMpiDiagnostics(MacroAleHllcMpiResult& result) {
  const std::string report =
      result.diagnostics.report_line.empty()
          ? BuildMacroAleHllcReportLine(result.diagnostics)
          : result.diagnostics.report_line;
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_hllc.global_transaction.preflight",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_hllc.seam_flux_equality",
      report);
  AppendDiagnostic(
      result.diagnostics_payload,
      "p1.hydro.macro_ale_hllc.global_conservation",
      report);
  if (result.diagnostics.no_partial_canonical_writeback) {
    AppendDiagnostic(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_hllc.transaction.clean",
        report);
  }
}

}  // namespace

bool MacroAleHllcMpiResult::is_complete() const noexcept {
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

MacroAleHllcMpiResult AdvanceMacroAleHllcMpiStep(
    dec3d::state::HydroStateView& local_hydro_view,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const dec3d::core::MeshUpdateProposal& global_radial_ale_proposal,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    StaticGridHydroOptions options,
    MPI_Comm communicator) noexcept {
  MacroAleHllcMpiResult result;
  int rank = 0;
  int rank_count = 0;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &rank_count);
  result.diagnostics.executed = true;

  bool local_preflight_ok = false;
  bool preview_ok = false;
  std::size_t global_face_begin = 0u;
  std::size_t local_face_count = local_geometry.radial_faces.size();
  std::size_t local_radial_cells = 0u;
  std::string local_failure;

  if (!decomposition.is_valid()) {
    local_failure = decomposition.failure_reason.empty()
                        ? "direct moving-face HLLC MPI requires a valid radial decomposition"
                        : decomposition.failure_reason;
  } else if (rank < 0 ||
             static_cast<std::size_t>(rank) >= decomposition.slices.size() ||
             rank_count != static_cast<int>(decomposition.slices.size())) {
    local_failure =
        "direct moving-face HLLC MPI rank count does not match decomposition";
  } else {
    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
    global_face_begin = slice.begin_index;
    local_radial_cells = slice.local_cell_count();
    result.diagnostics.global_face_begin = global_face_begin;
    result.diagnostics.local_face_count = local_face_count;

    if (!local_hydro_view.is_complete() || !local_geometry.is_valid()) {
      local_failure =
          "direct moving-face HLLC MPI requires complete local hydro and geometry";
    } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      local_failure =
          "direct moving-face HLLC MPI requires a finite positive timestep";
    } else if (!global_radial_ale_proposal.radial_face_indexing_is_global ||
               !global_radial_ale_proposal.has_radial_face_window(
                   global_face_begin,
                   local_face_count)) {
      local_failure =
          "direct moving-face HLLC MPI proposal/window failure";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::proposal_window_missing;
    } else if (local_face_count != local_radial_cells + 1u) {
      local_failure =
          "direct moving-face HLLC MPI local geometry is not the declared face window";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::proposal_window_missing;
    } else {
      local_preflight_ok = true;
    }
  }

  dec3d::mesh::SphericalGeometryMetadata staged_geometry = local_geometry;
  if (local_preflight_ok) {
    const auto preview = dec3d::mesh::ApplyRadialAleMeshUpdateProposal(
        global_radial_ale_proposal,
        staged_geometry,
        global_face_begin);
    preview_ok = preview.success;
    if (!preview_ok) {
      local_failure = preview.failure_reason.empty()
                          ? "direct moving-face HLLC MPI geometry preview failed"
                          : preview.failure_reason;
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::geometry_preview_failed;
    }
  }

  options.use_macro_zoning = true;
  options.use_ppm_reconstruction = true;
  options.apply_radial_ale_flux_correction = true;
  options.use_macro_ale_direct_moving_face_hllc = true;
  options.radial_ale_global_face_begin_index = global_face_begin;

  if (local_preflight_ok && preview_ok) {
    auto staged_state = CloneHydroViewToCanonicalState(local_hydro_view);
    dec3d::state::HydroAlphaViewOptions alpha_options;
    alpha_options.enabled = local_hydro_view.operator_local_alpha_chi;
    auto staged_view = dec3d::state::BuildHydroWorkView(staged_state, alpha_options);
    result.staged_hydro_result = AdvanceStaticGridHydro(
        staged_view,
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
              ? "direct moving-face HLLC MPI staged hydro failed"
              : result.staged_hydro_result.failure_reason;
      result.diagnostics.failure_class =
          result.staged_hydro_result.macro_ale_hllc_diagnostics.failure_class;
    }
  }

  if (result.staged_hydro_result.macro_ale_hllc_diagnostics.executed) {
    const auto& local_diag = result.staged_hydro_result.macro_ale_hllc_diagnostics;
    result.diagnostics.branch_left_flux_count = local_diag.branch_left_flux_count;
    result.diagnostics.branch_left_star_count = local_diag.branch_left_star_count;
    result.diagnostics.branch_right_star_count = local_diag.branch_right_star_count;
    result.diagnostics.branch_right_flux_count = local_diag.branch_right_flux_count;
    result.diagnostics.branch_static_zero_mismatch_count =
        local_diag.branch_static_zero_mismatch_count;
    result.diagnostics.branch_invalid_count = local_diag.branch_invalid_count;
    result.diagnostics.max_abs_w_face = local_diag.max_abs_w_face;
  }

  const int local_radiation_group_count =
      static_cast<int>(
          result.staged_hydro_result.macro_ale_hllc_diagnostics.radiation_group_count);
  int global_radiation_group_count = 0;
  MPI_Allreduce(
      &local_radiation_group_count,
      &global_radiation_group_count,
      1,
      MPI_INT,
      MPI_MAX,
      communicator);
  result.diagnostics.radiation_group_count =
      static_cast<std::size_t>(std::max(global_radiation_group_count, 0));

  const bool staged_geometry_matches =
      result.staged_hydro_result.macro_ale_hllc_staged_geometry_available &&
      SameRadialFaces(
          result.staged_hydro_result.macro_ale_hllc_staged_geometry,
          staged_geometry);

  const double local_signed_residuals[8] = {
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_mass_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_mom_r_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_mom_theta_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_mom_phi_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_e_fluid_total_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_chi_e_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_alpha_chi_residual,
      result.staged_hydro_result.macro_ale_hllc_diagnostics.global_radiation_chi_residual};
  double global_signed_residuals[8] = {};
  MPI_Allreduce(
      local_signed_residuals,
      global_signed_residuals,
      8,
      MPI_DOUBLE,
      MPI_SUM,
      communicator);
  const double local_abs_residuals[8] = {
      std::abs(local_signed_residuals[0]),
      std::abs(local_signed_residuals[1]),
      std::abs(local_signed_residuals[2]),
      std::abs(local_signed_residuals[3]),
      std::abs(local_signed_residuals[4]),
      std::abs(local_signed_residuals[5]),
      std::abs(local_signed_residuals[6]),
      std::abs(local_signed_residuals[7])};
  double global_max_residuals[8] = {};
  MPI_Allreduce(
      local_abs_residuals,
      global_max_residuals,
      8,
      MPI_DOUBLE,
      MPI_MAX,
      communicator);
  result.diagnostics.global_mass_residual = global_signed_residuals[0];
  result.diagnostics.global_mom_r_residual = global_signed_residuals[1];
  result.diagnostics.global_mom_theta_residual = global_signed_residuals[2];
  result.diagnostics.global_mom_phi_residual = global_signed_residuals[3];
  result.diagnostics.global_e_fluid_total_residual = global_signed_residuals[4];
  result.diagnostics.global_chi_e_residual = global_signed_residuals[5];
  result.diagnostics.global_alpha_chi_residual = global_signed_residuals[6];
  result.diagnostics.global_radiation_chi_residual = global_signed_residuals[7];
  result.diagnostics.local_max_mass_residual = global_max_residuals[0];
  result.diagnostics.local_max_mom_r_residual = global_max_residuals[1];
  result.diagnostics.local_max_mom_theta_residual = global_max_residuals[2];
  result.diagnostics.local_max_mom_phi_residual = global_max_residuals[3];
  result.diagnostics.local_max_e_fluid_total_residual = global_max_residuals[4];
  result.diagnostics.local_max_chi_e_residual = global_max_residuals[5];
  result.diagnostics.local_max_alpha_chi_residual = global_max_residuals[6];
  result.diagnostics.local_max_radiation_chi_residual = global_max_residuals[7];

  const double local_energy_scale = std::max(
      SumAbsEnergyExtensive(local_hydro_view, local_geometry),
      SumAbsEnergyExtensive(
          result.staged_hydro_result.macro_ale_hllc_staged_cells,
          staged_geometry));
  double global_energy_scale = 0.0;
  MPI_Allreduce(
      &local_energy_scale,
      &global_energy_scale,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      communicator);
  result.diagnostics.global_e_fluid_total_budget_scale =
      (global_energy_scale > 0.0 && Finite(global_energy_scale))
          ? global_energy_scale
          : 1.0;
  result.diagnostics.global_e_fluid_total_relative_residual =
      std::abs(result.diagnostics.global_e_fluid_total_residual) /
      result.diagnostics.global_e_fluid_total_budget_scale;

  const std::size_t global_face_count =
      std::max<std::size_t>(
          global_radial_ale_proposal.global_radial_face_count,
          1u);
  std::vector<int> local_has_values(global_face_count, 0);
  std::vector<double> local_w_values(global_face_count, 0.0);
  std::vector<double> local_proposed_values(global_face_count, 0.0);
  std::vector<double> local_preview_values(global_face_count, 0.0);
  std::vector<double> local_flux_mass(global_face_count, 0.0);
  std::vector<double> local_flux_mom_r(global_face_count, 0.0);
  std::vector<double> local_flux_mom_theta(global_face_count, 0.0);
  std::vector<double> local_flux_mom_phi(global_face_count, 0.0);
  std::vector<double> local_flux_e(global_face_count, 0.0);
  std::vector<double> local_flux_chi_e(global_face_count, 0.0);
  std::vector<double> local_flux_alpha_chi(global_face_count, 0.0);
  std::vector<double> local_flux_radiation_chi(global_face_count, 0.0);
  const auto& face_window =
      result.staged_hydro_result.macro_ale_hllc_local_face_window;
  if (local_preflight_ok &&
      preview_ok &&
      face_window.is_complete() &&
      result.staged_hydro_result.macro_ale_hllc_staged_geometry_available) {
    for (std::size_t local_face = 0u; local_face < face_window.local_face_count; ++local_face) {
      const std::size_t global_face = face_window.global_face_begin + local_face;
      if (global_face > 0u && global_face + 1u < global_face_count) {
        local_has_values[global_face] = 1;
        local_w_values[global_face] = face_window.face_speeds[local_face];
        local_proposed_values[global_face] =
            face_window.proposed_face_radii[local_face];
        local_preview_values[global_face] =
            result.staged_hydro_result
                .macro_ale_hllc_staged_geometry
                .radial_faces[local_face];
        local_flux_mass[global_face] =
            face_window.face_fluxes[local_face].rho;
        local_flux_mom_r[global_face] =
            face_window.face_fluxes[local_face].mom_r;
        local_flux_mom_theta[global_face] =
            face_window.face_fluxes[local_face].mom_theta;
        local_flux_mom_phi[global_face] =
            face_window.face_fluxes[local_face].mom_phi;
        local_flux_e[global_face] =
            face_window.face_fluxes[local_face].e_fluid_total;
        local_flux_chi_e[global_face] =
            face_window.face_fluxes[local_face].chi_e;
        local_flux_alpha_chi[global_face] =
            face_window.face_fluxes[local_face].alpha_chi;
        for (const double flux : face_window.face_fluxes[local_face].radiation_chi) {
          local_flux_radiation_chi[global_face] =
              std::max(local_flux_radiation_chi[global_face], std::abs(flux));
        }
      }
    }
  }

  result.diagnostics.max_seam_w_face_delta =
      MaxSharedFaceDelta(local_w_values, local_has_values, communicator);
  result.diagnostics.max_seam_proposed_radius_delta =
      MaxSharedFaceDelta(local_proposed_values, local_has_values, communicator);
  result.diagnostics.max_seam_preview_radius_delta =
      MaxSharedFaceDelta(local_preview_values, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mass_delta =
      MaxSharedFaceDelta(local_flux_mass, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mass_scale =
      MaxSharedFaceScale(local_flux_mass, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mass_relative_delta =
      result.diagnostics.max_seam_flux_mass_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_mass_scale);
  result.diagnostics.max_seam_flux_mom_r_delta =
      MaxSharedFaceDelta(local_flux_mom_r, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_r_scale =
      MaxSharedFaceScale(local_flux_mom_r, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_r_relative_delta =
      result.diagnostics.max_seam_flux_mom_r_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_mom_r_scale);
  result.diagnostics.max_seam_flux_mom_theta_delta =
      MaxSharedFaceDelta(local_flux_mom_theta, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_theta_scale =
      MaxSharedFaceScale(local_flux_mom_theta, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_theta_relative_delta =
      result.diagnostics.max_seam_flux_mom_theta_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_mom_theta_scale);
  result.diagnostics.max_seam_flux_mom_phi_delta =
      MaxSharedFaceDelta(local_flux_mom_phi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_phi_scale =
      MaxSharedFaceScale(local_flux_mom_phi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_mom_phi_relative_delta =
      result.diagnostics.max_seam_flux_mom_phi_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_mom_phi_scale);
  result.diagnostics.max_seam_flux_e_fluid_total_delta =
      MaxSharedFaceDelta(local_flux_e, local_has_values, communicator);
  result.diagnostics.max_seam_flux_e_fluid_total_scale =
      MaxSharedFaceScale(local_flux_e, local_has_values, communicator);
  result.diagnostics.max_seam_flux_e_fluid_total_relative_delta =
      result.diagnostics.max_seam_flux_e_fluid_total_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_e_fluid_total_scale);
  result.diagnostics.max_seam_flux_chi_e_delta =
      MaxSharedFaceDelta(local_flux_chi_e, local_has_values, communicator);
  result.diagnostics.max_seam_flux_chi_e_scale =
      MaxSharedFaceScale(local_flux_chi_e, local_has_values, communicator);
  result.diagnostics.max_seam_flux_chi_e_relative_delta =
      result.diagnostics.max_seam_flux_chi_e_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_chi_e_scale);
  result.diagnostics.max_seam_flux_alpha_chi_delta =
      MaxSharedFaceDelta(local_flux_alpha_chi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_alpha_chi_scale =
      MaxSharedFaceScale(local_flux_alpha_chi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_alpha_chi_relative_delta =
      result.diagnostics.max_seam_flux_alpha_chi_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_alpha_chi_scale);
  result.diagnostics.max_seam_flux_radiation_chi_delta =
      MaxSharedFaceDelta(local_flux_radiation_chi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_radiation_chi_scale =
      MaxSharedFaceScale(local_flux_radiation_chi, local_has_values, communicator);
  result.diagnostics.max_seam_flux_radiation_chi_relative_delta =
      result.diagnostics.max_seam_flux_radiation_chi_delta /
      std::max(1.0, result.diagnostics.max_seam_flux_radiation_chi_scale);

  const bool seam_finite =
      Finite(result.diagnostics.max_seam_w_face_delta) &&
      Finite(result.diagnostics.max_seam_proposed_radius_delta) &&
      Finite(result.diagnostics.max_seam_preview_radius_delta) &&
      Finite(result.diagnostics.max_seam_flux_mass_delta) &&
      Finite(result.diagnostics.max_seam_flux_mass_scale) &&
      Finite(result.diagnostics.max_seam_flux_mass_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_r_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_r_scale) &&
      Finite(result.diagnostics.max_seam_flux_mom_r_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_theta_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_theta_scale) &&
      Finite(result.diagnostics.max_seam_flux_mom_theta_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_phi_delta) &&
      Finite(result.diagnostics.max_seam_flux_mom_phi_scale) &&
      Finite(result.diagnostics.max_seam_flux_mom_phi_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_e_fluid_total_delta) &&
      Finite(result.diagnostics.max_seam_flux_e_fluid_total_scale) &&
      Finite(result.diagnostics.max_seam_flux_e_fluid_total_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_chi_e_delta) &&
      Finite(result.diagnostics.max_seam_flux_chi_e_scale) &&
      Finite(result.diagnostics.max_seam_flux_chi_e_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_alpha_chi_delta) &&
      Finite(result.diagnostics.max_seam_flux_alpha_chi_scale) &&
      Finite(result.diagnostics.max_seam_flux_alpha_chi_relative_delta) &&
      Finite(result.diagnostics.max_seam_flux_radiation_chi_delta) &&
      Finite(result.diagnostics.max_seam_flux_radiation_chi_scale) &&
      Finite(result.diagnostics.max_seam_flux_radiation_chi_relative_delta);
  const bool residuals_finite =
      Finite(result.diagnostics.global_mass_residual) &&
      Finite(result.diagnostics.global_mom_r_residual) &&
      Finite(result.diagnostics.global_mom_theta_residual) &&
      Finite(result.diagnostics.global_mom_phi_residual) &&
      Finite(result.diagnostics.global_e_fluid_total_residual) &&
      Finite(result.diagnostics.global_e_fluid_total_budget_scale) &&
      Finite(result.diagnostics.global_e_fluid_total_relative_residual) &&
      Finite(result.diagnostics.global_chi_e_residual) &&
      Finite(result.diagnostics.global_alpha_chi_residual) &&
      Finite(result.diagnostics.global_radiation_chi_residual) &&
      Finite(result.diagnostics.local_max_mass_residual) &&
      Finite(result.diagnostics.local_max_mom_r_residual) &&
      Finite(result.diagnostics.local_max_mom_theta_residual) &&
      Finite(result.diagnostics.local_max_mom_phi_residual) &&
      Finite(result.diagnostics.local_max_e_fluid_total_residual) &&
      Finite(result.diagnostics.local_max_chi_e_residual) &&
      Finite(result.diagnostics.local_max_alpha_chi_residual) &&
      Finite(result.diagnostics.local_max_radiation_chi_residual);
  const bool seam_within_tolerance =
      result.diagnostics.max_seam_w_face_delta <= kA4SeamGeometryTolerance &&
      result.diagnostics.max_seam_proposed_radius_delta <= kA4SeamGeometryTolerance &&
      result.diagnostics.max_seam_preview_radius_delta <= kA4SeamGeometryTolerance &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_mass_delta,
          result.diagnostics.max_seam_flux_mass_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_mom_r_delta,
          result.diagnostics.max_seam_flux_mom_r_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_mom_theta_delta,
          result.diagnostics.max_seam_flux_mom_theta_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_mom_phi_delta,
          result.diagnostics.max_seam_flux_mom_phi_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_e_fluid_total_delta,
          result.diagnostics.max_seam_flux_e_fluid_total_scale,
          kA4EnergyFluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_chi_e_delta,
          result.diagnostics.max_seam_flux_chi_e_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_alpha_chi_delta,
          result.diagnostics.max_seam_flux_alpha_chi_scale,
          kA4FluxTolerance) &&
      WithinScaledTolerance(
          result.diagnostics.max_seam_flux_radiation_chi_delta,
          result.diagnostics.max_seam_flux_radiation_chi_scale,
          kA4FluxTolerance);
  const double energy_conservation_tolerance = std::max(
      kA4EnergyConservationTolerance,
      kA4EnergyConservationRelativeTolerance *
          result.diagnostics.global_e_fluid_total_budget_scale);
  const bool residuals_within_tolerance =
      std::abs(result.diagnostics.global_mass_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_mom_r_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_mom_theta_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_mom_phi_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_e_fluid_total_residual) <=
          energy_conservation_tolerance &&
      std::abs(result.diagnostics.global_chi_e_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_alpha_chi_residual) <= kA4ConservationTolerance &&
      std::abs(result.diagnostics.global_radiation_chi_residual) <= kA4ConservationTolerance &&
      result.diagnostics.local_max_mass_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_mom_r_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_mom_theta_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_mom_phi_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_e_fluid_total_residual <=
          energy_conservation_tolerance &&
      result.diagnostics.local_max_chi_e_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_alpha_chi_residual <= kA4ConservationTolerance &&
      result.diagnostics.local_max_radiation_chi_residual <= kA4ConservationTolerance;

  if (local_failure.empty()) {
    if (!seam_within_tolerance) {
      local_failure =
          "direct moving-face HLLC MPI seam flux equality tolerance exceeded";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::seam_flux_delta_exceeded;
    } else if (!residuals_within_tolerance) {
      local_failure =
          "direct moving-face HLLC MPI conservation tolerance exceeded";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
    }
  }

  const bool local_stage_bool =
      local_preflight_ok &&
      preview_ok &&
      result.staged_hydro_result.success &&
      result.staged_hydro_result.macro_ale_hllc_executed &&
      result.staged_hydro_result.macro_ale_hllc_staged_writeback_available &&
      result.staged_hydro_result.macro_ale_hllc_staged_geometry_available &&
      result.staged_hydro_result.macro_ale_hllc_local_face_window.is_complete() &&
      staged_geometry_matches &&
      seam_finite &&
      residuals_finite &&
      seam_within_tolerance &&
      residuals_within_tolerance &&
      result.diagnostics.branch_invalid_count == 0u;
  int mpi_rank = -1;
  MPI_Comm_rank(communicator, &mpi_rank);
  result.diagnostics.mpi_rank = mpi_rank;
  result.diagnostics.mpi_local_preflight_ok = local_preflight_ok;
  result.diagnostics.mpi_preview_ok = preview_ok;
  result.diagnostics.mpi_staged_hydro_success =
      result.staged_hydro_result.success;
  result.diagnostics.mpi_staged_hydro_executed =
      result.staged_hydro_result.macro_ale_hllc_executed;
  result.diagnostics.mpi_staged_writeback_available =
      result.staged_hydro_result.macro_ale_hllc_staged_writeback_available;
  result.diagnostics.mpi_staged_geometry_available =
      result.staged_hydro_result.macro_ale_hllc_staged_geometry_available;
  result.diagnostics.mpi_local_face_window_complete =
      result.staged_hydro_result.macro_ale_hllc_local_face_window.is_complete();
  result.diagnostics.mpi_staged_geometry_matches = staged_geometry_matches;
  result.diagnostics.mpi_seam_finite = seam_finite;
  result.diagnostics.mpi_residuals_finite = residuals_finite;
  result.diagnostics.mpi_seam_within_tolerance = seam_within_tolerance;
  result.diagnostics.mpi_residuals_within_tolerance = residuals_within_tolerance;
  result.diagnostics.mpi_local_stage_ok = local_stage_bool;
  if (!local_stage_bool && local_failure.empty()) {
    if (!local_preflight_ok) {
      local_failure = "direct moving-face HLLC MPI local preflight failed";
    } else if (!preview_ok) {
      local_failure = "direct moving-face HLLC MPI local geometry preview failed";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::geometry_preview_failed;
    } else if (!result.staged_hydro_result.success) {
      local_failure =
          result.staged_hydro_result.failure_reason.empty()
              ? "direct moving-face HLLC MPI staged hydro failed"
              : result.staged_hydro_result.failure_reason;
      result.diagnostics.failure_class =
          result.staged_hydro_result.macro_ale_hllc_diagnostics.failure_class;
    } else if (!result.staged_hydro_result.macro_ale_hllc_executed) {
      local_failure = "direct moving-face HLLC MPI staged hydro did not execute direct HLLC";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::transaction_aborted;
    } else if (!result.staged_hydro_result.macro_ale_hllc_staged_writeback_available) {
      local_failure =
          "direct moving-face HLLC MPI staged hydro writeback cells are missing";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::publish_preflight_failed;
    } else if (!result.staged_hydro_result.macro_ale_hllc_staged_geometry_available) {
      local_failure =
          "direct moving-face HLLC MPI staged hydro geometry is missing";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::geometry_preview_failed;
    } else if (!result.staged_hydro_result.macro_ale_hllc_local_face_window.is_complete()) {
      local_failure =
          "direct moving-face HLLC MPI local face window is incomplete";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::proposal_window_missing;
    } else if (!staged_geometry_matches) {
      local_failure =
          "direct moving-face HLLC MPI staged geometry does not match preview";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::geometry_preview_failed;
    } else if (!seam_finite || !residuals_finite) {
      local_failure =
          "direct moving-face HLLC MPI diagnostics contain non-finite values";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
    } else if (!seam_within_tolerance) {
      local_failure =
          "direct moving-face HLLC MPI seam flux equality tolerance exceeded";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::seam_flux_delta_exceeded;
    } else if (!residuals_within_tolerance) {
      local_failure =
          "direct moving-face HLLC MPI conservation tolerance exceeded";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
    } else if (result.diagnostics.branch_invalid_count != 0u) {
      local_failure =
          "direct moving-face HLLC MPI moving-interface branch classification failed";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::moving_hllc_branch_unclassified;
    } else {
      local_failure =
          "direct moving-face HLLC MPI local stage failed an unknown preflight gate";
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::transaction_aborted;
    }
  }
  const int local_stage_ok = local_stage_bool ? 1 : 0;
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
  if (!result.global_stage_ok && local_failure.empty()) {
    local_failure = "direct moving-face HLLC MPI global stage failed";
    result.diagnostics.failure_class =
        MacroAleHllcFailureClass::transaction_aborted;
  }

  const bool local_publishable = PublishableA4StagedResult(
      local_hydro_view,
      result.staged_hydro_result,
      result.staged_hydro_result.macro_ale_hllc_staged_geometry,
      staged_geometry);
  result.diagnostics.mpi_local_publishable = local_publishable;
  const int local_publish_ok =
      result.global_stage_ok && local_publishable ? 1 : 0;
  int global_publish_ok = 0;
  MPI_Allreduce(
      &local_publish_ok,
      &global_publish_ok,
      1,
      MPI_INT,
      MPI_MIN,
      communicator);
  result.local_publishable = result.global_stage_ok && local_publishable;
  result.global_publish_ok = global_publish_ok == 1;
  result.diagnostics.global_publish_ok = result.global_publish_ok;
  if (result.global_stage_ok && !result.global_publish_ok && local_failure.empty()) {
    local_failure = "direct moving-face HLLC MPI global publish preflight failed";
    result.diagnostics.failure_class =
        MacroAleHllcFailureClass::publish_preflight_failed;
  }

  bool hydro_published = false;
  if (result.global_stage_ok && result.global_publish_ok) {
    std::string publish_failure;
    hydro_published = PublishMacroAleStagedHydroState(
        local_hydro_view,
        result.staged_hydro_result,
        &publish_failure);
    if (hydro_published) {
      local_geometry = result.staged_hydro_result.macro_ale_hllc_staged_geometry;
    } else {
      local_failure = publish_failure.empty()
                          ? "direct moving-face HLLC MPI publish failed after global publish preflight"
                          : publish_failure;
      result.diagnostics.failure_class =
          MacroAleHllcFailureClass::publish_preflight_failed;
    }
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
  result.canonical_writeback_published =
      hydro_published && global_published_ok == 1;
  result.mesh_commit_published = result.canonical_writeback_published;
  result.diagnostics.no_partial_canonical_writeback =
      result.canonical_writeback_published && result.mesh_commit_published;
  if (result.global_stage_ok &&
      result.global_publish_ok &&
      global_published_ok != 1) {
    local_failure =
        "direct moving-face HLLC MPI publish failed after global publish preflight";
    result.diagnostics.failure_class =
        MacroAleHllcFailureClass::publish_preflight_failed;
  }

  result.success =
      result.global_stage_ok &&
      result.global_publish_ok &&
      result.canonical_writeback_published &&
      result.mesh_commit_published &&
      local_failure.empty();
  result.failure_reason = result.success ? std::string{} : local_failure;
  if (result.success) {
    result.diagnostics.failure_class = MacroAleHllcFailureClass::none;
  }
  result.diagnostics.report_line =
      BuildMacroAleHllcReportLine(result.diagnostics);
  AppendMpiDiagnostics(result);
  return result;
}

}  // namespace dec3d::hydro
