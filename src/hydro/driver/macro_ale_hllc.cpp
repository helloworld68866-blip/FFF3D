#include "hydro/driver/macro_ale_hllc.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <utility>

namespace dec3d::hydro {
namespace {

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& payload,
    std::string code,
    std::string message) {
  payload.entries.push_back({std::move(code), std::move(message)});
}

}  // namespace

bool MacroAleHllcLocalFaceWindow::is_complete() const noexcept {
  return available &&
         local_face_count > 0u &&
         face_fluxes.size() == local_face_count &&
         face_speeds.size() == local_face_count &&
         proposed_face_radii.size() == local_face_count;
}

bool MacroAleHllcDiagnostics::is_complete() const noexcept {
  return executed &&
         failure_class == MacroAleHllcFailureClass::none &&
         branch_invalid_count == 0u &&
         Finite(max_abs_w_face) &&
         Finite(global_mass_residual) &&
         Finite(global_mom_r_residual) &&
         Finite(global_mom_theta_residual) &&
         Finite(global_mom_phi_residual) &&
         Finite(global_e_fluid_total_residual) &&
         Finite(global_e_fluid_total_budget_scale) &&
         Finite(global_e_fluid_total_relative_residual) &&
         Finite(global_chi_e_residual) &&
         Finite(global_alpha_chi_residual) &&
         Finite(global_radiation_chi_residual) &&
         Finite(local_max_mass_residual) &&
         Finite(local_max_mom_r_residual) &&
         Finite(local_max_mom_theta_residual) &&
         Finite(local_max_mom_phi_residual) &&
         Finite(local_max_e_fluid_total_residual) &&
         Finite(local_max_chi_e_residual) &&
         Finite(local_max_alpha_chi_residual) &&
         Finite(local_max_radiation_chi_residual) &&
         Finite(max_seam_w_face_delta) &&
         Finite(max_seam_proposed_radius_delta) &&
         Finite(max_seam_preview_radius_delta) &&
         Finite(max_seam_flux_mass_delta) &&
         Finite(max_seam_flux_mass_scale) &&
         Finite(max_seam_flux_mass_relative_delta) &&
         Finite(max_seam_flux_mom_r_delta) &&
         Finite(max_seam_flux_mom_r_scale) &&
         Finite(max_seam_flux_mom_r_relative_delta) &&
         Finite(max_seam_flux_mom_theta_delta) &&
         Finite(max_seam_flux_mom_theta_scale) &&
         Finite(max_seam_flux_mom_theta_relative_delta) &&
         Finite(max_seam_flux_mom_phi_delta) &&
         Finite(max_seam_flux_mom_phi_scale) &&
         Finite(max_seam_flux_mom_phi_relative_delta) &&
         Finite(max_seam_flux_e_fluid_total_delta) &&
         Finite(max_seam_flux_e_fluid_total_scale) &&
         Finite(max_seam_flux_e_fluid_total_relative_delta) &&
         Finite(max_seam_flux_chi_e_delta) &&
         Finite(max_seam_flux_chi_e_scale) &&
         Finite(max_seam_flux_chi_e_relative_delta) &&
         Finite(max_seam_flux_alpha_chi_delta) &&
         Finite(max_seam_flux_alpha_chi_scale) &&
         Finite(max_seam_flux_alpha_chi_relative_delta) &&
         Finite(max_seam_flux_radiation_chi_delta) &&
         Finite(max_seam_flux_radiation_chi_scale) &&
         Finite(max_seam_flux_radiation_chi_relative_delta) &&
         global_stage_ok &&
         global_publish_ok &&
         no_partial_canonical_writeback &&
         !report_line.empty();
}

const char* MacroAleHllcModeName(MacroAleHllcMode mode) noexcept {
  switch (mode) {
    case MacroAleHllcMode::direct_moving_face_hllc:
      return "direct_moving_face_hllc";
  }
  return "direct_moving_face_hllc";
}

const char* MacroAleHllcFailureClassName(
    MacroAleHllcFailureClass failure_class) noexcept {
  switch (failure_class) {
    case MacroAleHllcFailureClass::none:
      return "none";
    case MacroAleHllcFailureClass::proposal_not_global_indexed:
      return "proposal_not_global_indexed";
    case MacroAleHllcFailureClass::proposal_window_missing:
      return "proposal_window_missing";
    case MacroAleHllcFailureClass::proposal_non_monotone:
      return "proposal_non_monotone";
    case MacroAleHllcFailureClass::geometry_preview_failed:
      return "geometry_preview_failed";
    case MacroAleHllcFailureClass::coarse_ghost_missing:
      return "coarse_ghost_missing";
    case MacroAleHllcFailureClass::moving_hllc_nonphysical:
      return "moving_hllc_nonphysical";
    case MacroAleHllcFailureClass::moving_hllc_branch_unclassified:
      return "moving_hllc_branch_unclassified";
    case MacroAleHllcFailureClass::extensive_budget_residual_exceeded:
      return "extensive_budget_residual_exceeded";
    case MacroAleHllcFailureClass::seam_flux_delta_exceeded:
      return "seam_flux_delta_exceeded";
    case MacroAleHllcFailureClass::physical_recovery_failed:
      return "physical_recovery_failed";
    case MacroAleHllcFailureClass::publish_preflight_failed:
      return "publish_preflight_failed";
    case MacroAleHllcFailureClass::transaction_aborted:
      return "transaction_aborted";
    case MacroAleHllcFailureClass::required_diagnostic_missing:
      return "required_diagnostic_missing";
  }
  return "required_diagnostic_missing";
}

const char* MacroAleHllcRemapOrder() noexcept {
  return "none";
}

std::string BuildMacroAleHllcReportLine(
    const MacroAleHllcDiagnostics& diagnostics) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "mode=direct_moving_face_hllc"
         << "; remap_order=none"
         << "; radial_ale_flux_mode=moving_interface_hllc"
         << "; global_face_begin=" << diagnostics.global_face_begin
         << "; local_face_count=" << diagnostics.local_face_count
         << "; moving_interface_branch_count_left="
         << diagnostics.branch_left_flux_count
         << "; moving_interface_branch_count_left_star="
         << diagnostics.branch_left_star_count
         << "; moving_interface_branch_count_right_star="
         << diagnostics.branch_right_star_count
         << "; moving_interface_branch_count_right="
         << diagnostics.branch_right_flux_count
         << "; moving_interface_branch_count_static_zero_mismatch="
         << diagnostics.branch_static_zero_mismatch_count
         << "; moving_interface_branch_count_invalid="
         << diagnostics.branch_invalid_count
         << "; max_abs_w_face=" << diagnostics.max_abs_w_face
         << "; global_mass_residual=" << diagnostics.global_mass_residual
         << "; global_mom_r_residual=" << diagnostics.global_mom_r_residual
         << "; global_mom_theta_residual=" << diagnostics.global_mom_theta_residual
         << "; global_mom_phi_residual=" << diagnostics.global_mom_phi_residual
         << "; global_e_fluid_total_residual="
         << diagnostics.global_e_fluid_total_residual
         << "; global_e_fluid_total_budget_scale="
         << diagnostics.global_e_fluid_total_budget_scale
         << "; global_e_fluid_total_relative_residual="
         << diagnostics.global_e_fluid_total_relative_residual
         << "; global_chi_e_residual=" << diagnostics.global_chi_e_residual
         << "; global_alpha_chi_residual=" << diagnostics.global_alpha_chi_residual
        << "; global_radiation_chi_residual="
        << diagnostics.global_radiation_chi_residual
        << "; budget_components=rho,mom_r,mom_theta,mom_phi,e_fluid_total,chi_e,alpha_chi,radiation_chi"
        << "; budget_residual_source=computed_full_transport"
        << "; moving_mesh_face_velocity_source=hydro_ale_hllc_face_velocity"
        << "; alpha_hydro_terms_share_hydro_scalar_bundle=true"
        << "; radiation_scalar_bundle_transport="
        << (diagnostics.radiation_group_count > 0u ? "enabled" : "disabled")
        << "; radiation_group_count=" << diagnostics.radiation_group_count
        << "; mpi_rank=" << diagnostics.mpi_rank
        << "; mpi_local_preflight_ok="
        << (diagnostics.mpi_local_preflight_ok ? "true" : "false")
        << "; mpi_preview_ok="
        << (diagnostics.mpi_preview_ok ? "true" : "false")
        << "; mpi_staged_hydro_success="
        << (diagnostics.mpi_staged_hydro_success ? "true" : "false")
        << "; mpi_staged_hydro_executed="
        << (diagnostics.mpi_staged_hydro_executed ? "true" : "false")
        << "; mpi_staged_writeback_available="
        << (diagnostics.mpi_staged_writeback_available ? "true" : "false")
        << "; mpi_staged_geometry_available="
        << (diagnostics.mpi_staged_geometry_available ? "true" : "false")
        << "; mpi_local_face_window_complete="
        << (diagnostics.mpi_local_face_window_complete ? "true" : "false")
        << "; mpi_staged_geometry_matches="
        << (diagnostics.mpi_staged_geometry_matches ? "true" : "false")
        << "; mpi_seam_finite="
        << (diagnostics.mpi_seam_finite ? "true" : "false")
        << "; mpi_residuals_finite="
        << (diagnostics.mpi_residuals_finite ? "true" : "false")
        << "; mpi_seam_within_tolerance="
        << (diagnostics.mpi_seam_within_tolerance ? "true" : "false")
        << "; mpi_residuals_within_tolerance="
        << (diagnostics.mpi_residuals_within_tolerance ? "true" : "false")
        << "; mpi_local_stage_ok="
        << (diagnostics.mpi_local_stage_ok ? "true" : "false")
        << "; mpi_local_publishable="
        << (diagnostics.mpi_local_publishable ? "true" : "false")
        << "; local_max_mass_residual="
         << diagnostics.local_max_mass_residual
         << "; local_max_mom_r_residual="
         << diagnostics.local_max_mom_r_residual
         << "; local_max_mom_theta_residual="
         << diagnostics.local_max_mom_theta_residual
         << "; local_max_mom_phi_residual="
         << diagnostics.local_max_mom_phi_residual
         << "; local_max_e_fluid_total_residual="
         << diagnostics.local_max_e_fluid_total_residual
         << "; local_max_chi_e_residual="
         << diagnostics.local_max_chi_e_residual
         << "; local_max_alpha_chi_residual="
         << diagnostics.local_max_alpha_chi_residual
         << "; local_max_radiation_chi_residual="
         << diagnostics.local_max_radiation_chi_residual
         << "; radiation_hydro_terms_follow_post_H_geometry=true"
         << "; max_seam_w_face_delta=" << diagnostics.max_seam_w_face_delta
         << "; max_seam_proposed_radius_delta="
         << diagnostics.max_seam_proposed_radius_delta
         << "; max_seam_preview_radius_delta="
         << diagnostics.max_seam_preview_radius_delta
         << "; max_seam_flux_mass_delta="
         << diagnostics.max_seam_flux_mass_delta
         << "; max_seam_flux_mass_scale="
         << diagnostics.max_seam_flux_mass_scale
         << "; max_seam_flux_mass_relative_delta="
         << diagnostics.max_seam_flux_mass_relative_delta
         << "; max_seam_flux_mom_r_delta="
         << diagnostics.max_seam_flux_mom_r_delta
         << "; max_seam_flux_mom_r_scale="
         << diagnostics.max_seam_flux_mom_r_scale
         << "; max_seam_flux_mom_r_relative_delta="
         << diagnostics.max_seam_flux_mom_r_relative_delta
         << "; max_seam_flux_mom_theta_delta="
         << diagnostics.max_seam_flux_mom_theta_delta
         << "; max_seam_flux_mom_theta_scale="
         << diagnostics.max_seam_flux_mom_theta_scale
         << "; max_seam_flux_mom_theta_relative_delta="
         << diagnostics.max_seam_flux_mom_theta_relative_delta
         << "; max_seam_flux_mom_phi_delta="
         << diagnostics.max_seam_flux_mom_phi_delta
         << "; max_seam_flux_mom_phi_scale="
         << diagnostics.max_seam_flux_mom_phi_scale
         << "; max_seam_flux_mom_phi_relative_delta="
         << diagnostics.max_seam_flux_mom_phi_relative_delta
         << "; max_seam_flux_e_fluid_total_delta="
         << diagnostics.max_seam_flux_e_fluid_total_delta
         << "; max_seam_flux_e_fluid_total_scale="
         << diagnostics.max_seam_flux_e_fluid_total_scale
         << "; max_seam_flux_e_fluid_total_relative_delta="
         << diagnostics.max_seam_flux_e_fluid_total_relative_delta
         << "; max_seam_flux_chi_e_delta="
         << diagnostics.max_seam_flux_chi_e_delta
         << "; max_seam_flux_chi_e_scale="
         << diagnostics.max_seam_flux_chi_e_scale
         << "; max_seam_flux_chi_e_relative_delta="
         << diagnostics.max_seam_flux_chi_e_relative_delta
         << "; max_seam_flux_alpha_chi_delta="
         << diagnostics.max_seam_flux_alpha_chi_delta
         << "; max_seam_flux_alpha_chi_scale="
         << diagnostics.max_seam_flux_alpha_chi_scale
         << "; max_seam_flux_alpha_chi_relative_delta="
         << diagnostics.max_seam_flux_alpha_chi_relative_delta
         << "; max_seam_flux_radiation_chi_delta="
         << diagnostics.max_seam_flux_radiation_chi_delta
         << "; max_seam_flux_radiation_chi_scale="
         << diagnostics.max_seam_flux_radiation_chi_scale
         << "; max_seam_flux_radiation_chi_relative_delta="
         << diagnostics.max_seam_flux_radiation_chi_relative_delta
         << "; global_stage_ok="
         << (diagnostics.global_stage_ok ? "true" : "false")
         << "; global_publish_ok="
         << (diagnostics.global_publish_ok ? "true" : "false")
         << "; no_partial_canonical_writeback="
         << (diagnostics.no_partial_canonical_writeback ? "true" : "false")
         << "; failure_class="
         << MacroAleHllcFailureClassName(diagnostics.failure_class);
  return report.str();
}

void AppendMacroAleHllcDiagnostics(
    dec3d::core::DiagnosticsPayload& payload,
    const MacroAleHllcDiagnostics& diagnostics) {
  const std::string report =
      diagnostics.report_line.empty()
          ? BuildMacroAleHllcReportLine(diagnostics)
          : diagnostics.report_line;
  AppendDiagnostic(payload, "p1.hydro.macro_ale_hllc.executed", report);
  AppendDiagnostic(
      payload,
      "p1.hydro.macro_ale_hllc.moving_interface_branch",
      report);
  AppendDiagnostic(payload, "p1.hydro.macro_ale_hllc.extensive_budget", report);
  AppendDiagnostic(payload, "p1.hydro.macro_ale_hllc.transaction", report);
}

}  // namespace dec3d::hydro
