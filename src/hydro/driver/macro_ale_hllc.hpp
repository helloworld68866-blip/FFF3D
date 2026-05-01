#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/riemann/hllc_solver.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::hydro {

enum class MacroAleHllcMode {
  direct_moving_face_hllc,
};

enum class MacroAleHllcFailureClass {
  none,
  proposal_not_global_indexed,
  proposal_window_missing,
  proposal_non_monotone,
  geometry_preview_failed,
  coarse_ghost_missing,
  moving_hllc_nonphysical,
  moving_hllc_branch_unclassified,
  extensive_budget_residual_exceeded,
  seam_flux_delta_exceeded,
  physical_recovery_failed,
  publish_preflight_failed,
  transaction_aborted,
  required_diagnostic_missing,
};

struct MacroAleHllcLocalFaceWindow {
  bool available{false};
  std::size_t global_face_begin{0u};
  std::size_t local_face_count{0u};
  std::vector<HydroConservativeState> face_fluxes;
  std::vector<double> face_speeds;
  std::vector<double> proposed_face_radii;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MacroAleHllcDiagnostics {
  bool executed{false};
  std::size_t global_face_begin{0u};
  std::size_t local_face_count{0u};
  std::size_t branch_left_flux_count{0u};
  std::size_t branch_left_star_count{0u};
  std::size_t branch_right_star_count{0u};
  std::size_t branch_right_flux_count{0u};
  std::size_t branch_static_zero_mismatch_count{0u};
  std::size_t branch_invalid_count{0u};
  double max_abs_w_face{0.0};
  double global_mass_residual{0.0};
  double global_mom_r_residual{0.0};
  double global_mom_theta_residual{0.0};
  double global_mom_phi_residual{0.0};
  double global_e_fluid_total_residual{0.0};
  double global_e_fluid_total_budget_scale{1.0};
  double global_e_fluid_total_relative_residual{0.0};
  double global_chi_e_residual{0.0};
  double global_alpha_chi_residual{0.0};
  double global_radiation_chi_residual{0.0};
  double local_max_mass_residual{0.0};
  double local_max_mom_r_residual{0.0};
  double local_max_mom_theta_residual{0.0};
  double local_max_mom_phi_residual{0.0};
  double local_max_e_fluid_total_residual{0.0};
  double local_max_chi_e_residual{0.0};
  double local_max_alpha_chi_residual{0.0};
  double local_max_radiation_chi_residual{0.0};
  double max_seam_w_face_delta{0.0};
  double max_seam_proposed_radius_delta{0.0};
  double max_seam_preview_radius_delta{0.0};
  double max_seam_flux_mass_delta{0.0};
  double max_seam_flux_mass_scale{1.0};
  double max_seam_flux_mass_relative_delta{0.0};
  double max_seam_flux_mom_r_delta{0.0};
  double max_seam_flux_mom_r_scale{1.0};
  double max_seam_flux_mom_r_relative_delta{0.0};
  double max_seam_flux_mom_theta_delta{0.0};
  double max_seam_flux_mom_theta_scale{1.0};
  double max_seam_flux_mom_theta_relative_delta{0.0};
  double max_seam_flux_mom_phi_delta{0.0};
  double max_seam_flux_mom_phi_scale{1.0};
  double max_seam_flux_mom_phi_relative_delta{0.0};
  double max_seam_flux_e_fluid_total_delta{0.0};
  double max_seam_flux_e_fluid_total_scale{1.0};
  double max_seam_flux_e_fluid_total_relative_delta{0.0};
  double max_seam_flux_chi_e_delta{0.0};
  double max_seam_flux_chi_e_scale{1.0};
  double max_seam_flux_chi_e_relative_delta{0.0};
  double max_seam_flux_alpha_chi_delta{0.0};
  double max_seam_flux_alpha_chi_scale{1.0};
  double max_seam_flux_alpha_chi_relative_delta{0.0};
  double max_seam_flux_radiation_chi_delta{0.0};
  double max_seam_flux_radiation_chi_scale{1.0};
  double max_seam_flux_radiation_chi_relative_delta{0.0};
  std::size_t radiation_group_count{0u};
  int mpi_rank{-1};
  bool mpi_local_preflight_ok{true};
  bool mpi_preview_ok{true};
  bool mpi_staged_hydro_success{true};
  bool mpi_staged_hydro_executed{true};
  bool mpi_staged_writeback_available{true};
  bool mpi_staged_geometry_available{true};
  bool mpi_local_face_window_complete{true};
  bool mpi_staged_geometry_matches{true};
  bool mpi_seam_finite{true};
  bool mpi_residuals_finite{true};
  bool mpi_seam_within_tolerance{true};
  bool mpi_residuals_within_tolerance{true};
  bool mpi_local_stage_ok{true};
  bool mpi_local_publishable{true};
  bool global_stage_ok{false};
  bool global_publish_ok{false};
  bool no_partial_canonical_writeback{false};
  MacroAleHllcFailureClass failure_class{MacroAleHllcFailureClass::none};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] const char* MacroAleHllcModeName(MacroAleHllcMode mode) noexcept;
[[nodiscard]] const char* MacroAleHllcFailureClassName(
    MacroAleHllcFailureClass failure_class) noexcept;
[[nodiscard]] const char* MacroAleHllcRemapOrder() noexcept;
[[nodiscard]] std::string BuildMacroAleHllcReportLine(
    const MacroAleHllcDiagnostics& diagnostics);
void AppendMacroAleHllcDiagnostics(
    dec3d::core::DiagnosticsPayload& payload,
    const MacroAleHllcDiagnostics& diagnostics);

}  // namespace dec3d::hydro
