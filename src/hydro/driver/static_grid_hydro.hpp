#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/driver/macro_ale_hllc.hpp"
#include "hydro/driver/radial_overlap_remap.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <string_view>
#include <string>

namespace dec3d::hydro {

struct HydroBudgetResidualSummary {
  double old_mass{0.0};
  double flux_mass_delta{0.0};
  double source_mass_delta{0.0};
  double new_mass{0.0};
  double mass_residual{0.0};
  double old_mom_r{0.0};
  double flux_mom_r_delta{0.0};
  double source_mom_r_delta{0.0};
  double new_mom_r{0.0};
  double mom_r_residual{0.0};
  double old_e_fluid_total{0.0};
  double flux_e_fluid_total_delta{0.0};
  double source_e_fluid_total_delta{0.0};
  double new_e_fluid_total{0.0};
  double e_fluid_total_residual{0.0};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct StaticGridHydroOptions {
  bool apply_geometric_source{true};
  bool apply_radial_sweep{true};
  bool apply_theta_sweep{true};
  bool apply_phi_sweep{true};
  bool use_ppm_reconstruction{false};
  std::size_t reconstruction_ghost_layers{3u};
  bool apply_radial_ale_flux_correction{false};
  bool request_radial_ale_proposal{false};
  std::size_t radial_ale_global_face_begin_index{0u};
  bool use_macro_zoning{false};
  double macro_zoning_coarse_factor{0.5};
  bool macro_zoning_use_radial_ppm{true};
  bool macro_zoning_use_theta_ppm{true};
  bool macro_zoning_use_phi_ppm{true};
  bool defer_macro_ale_compatibility_writeback{false};
  bool use_macro_ale_direct_moving_face_hllc{false};
  bool debug_angular_stage_diagnostics{false};
  bool debug_macro_radial_ppm_face_smoothness{false};
  std::size_t debug_macro_radial_ppm_face{0u};
  bool enable_radiation_hydro_terms{false};
  bool enable_alpha_hydro_terms{false};
};

struct StaticGridHydroResult {
  bool success{false};
  bool geometric_source_executed{false};
  bool radial_executed{false};
  bool theta_executed{false};
  bool phi_executed{false};
  HydroBudgetResidualSummary budget;
  dec3d::core::DiagnosticsPayload diagnostics;
  bool macro_ale_staged_writeback_available{false};
  bool macro_ale_staged_geometry_available{false};
  dec3d::core::Array3D<HydroConservativeState> macro_ale_staged_cells;
  dec3d::mesh::SphericalGeometryMetadata macro_ale_staged_geometry;
  RadialOverlapRemapDiagnostics macro_ale_remap_diagnostics;
  bool macro_ale_hllc_executed{false};
  bool macro_ale_hllc_staged_writeback_available{false};
  bool macro_ale_hllc_staged_geometry_available{false};
  dec3d::core::Array3D<HydroConservativeState> macro_ale_hllc_staged_cells;
  dec3d::mesh::SphericalGeometryMetadata macro_ale_hllc_staged_geometry;
  MacroAleHllcDiagnostics macro_ale_hllc_diagnostics;
  MacroAleHllcLocalFaceWindow macro_ale_hllc_local_face_window;
  double timing_snapshot_wall_s{0.0};
  double timing_scratch_wall_s{0.0};
  double timing_radial_sweep_wall_s{0.0};
  double timing_macro_detect_wall_s{0.0};
  double timing_macro_restrict_wall_s{0.0};
  double timing_macro_update_wall_s{0.0};
  double timing_macro_radial_update_wall_s{0.0};
  double timing_macro_theta_update_wall_s{0.0};
  double timing_macro_phi_update_wall_s{0.0};
  double timing_macro_state_update_wall_s{0.0};
  double timing_macro_prolong_wall_s{0.0};
  double timing_theta_sweep_wall_s{0.0};
  double timing_phi_sweep_wall_s{0.0};
  double timing_commit_wall_s{0.0};
  double timing_source_wall_s{0.0};
  double timing_budget_wall_s{0.0};
  double timing_diagnostics_wall_s{0.0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] StaticGridHydroResult AdvanceStaticGridHydro(
    dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const StaticGridHydroOptions& options = {},
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal = nullptr) noexcept;

[[nodiscard]] StaticGridHydroResult AdvanceStaticGridHydro(
    dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    const StaticGridHydroOptions& options = {},
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal = nullptr) noexcept;

[[nodiscard]] bool ParseHydroBudgetResidualSummary(
    std::string_view report_line,
    HydroBudgetResidualSummary* summary) noexcept;

[[nodiscard]] bool TryExtractHydroBudgetResidualSummary(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    HydroBudgetResidualSummary* summary) noexcept;

[[nodiscard]] bool PublishMacroAleStagedHydroState(
    dec3d::state::HydroStateView& hydro_view,
    const StaticGridHydroResult& staged_result,
    std::string* failure_reason = nullptr) noexcept;

}  // namespace dec3d::hydro
