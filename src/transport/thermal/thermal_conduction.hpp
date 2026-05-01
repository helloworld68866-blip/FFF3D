#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dec3d::transport {

enum class ThermalConductionKappaModel {
  constant_user_supplied,
  spitzer_no_degeneracy,
  lee_more_with_degeneracy,
  unsupported
};

enum class ThermalConductionBackend {
  serial_dense_reference,
  hypre_parcsr_gmres_boomeramg,
  unsupported
};

struct ThermalConductionOptions {
  double dt_s{0.0};
  ThermalConductionKappaModel kappa_model{ThermalConductionKappaModel::constant_user_supplied};
  double kappa_e_cm_inv_s{0.0};
  double kappa_i_cm_inv_s{0.0};
  ThermalConductionBackend backend{ThermalConductionBackend::serial_dense_reference};
  GenericDiffusionBoundaryPolicy boundary_policy;
  double electron_energy_floor{0.0};
  double ion_energy_floor{0.0};
  GenericDiffusionReferenceSolveOptions serial_reference_options{};
};

struct ThermalConductionSolveResult {
  bool success{false};
  std::vector<double> scalar_new;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string backend_executed;
  double residual_l2{0.0};
  double residual_linf{0.0};
  double max_abs_delta{0.0};
};

using ThermalConductionSolveHook =
    std::function<ThermalConductionSolveResult(const GenericDiffusionAssemblyResult&)>;

struct ThermalConductionResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string electron_assembly_report;
  std::string ion_assembly_report;
  std::string electron_solve_report;
  std::string ion_solve_report;
  std::string backend_requested;
  std::string backend_executed;
  std::size_t cell_count{0};
  double dt_s{0.0};
  double kappa_e_cm_inv_s{0.0};
  double kappa_i_cm_inv_s{0.0};
  std::string kappa_model_requested;
  std::string kappa_model_executed;
  std::string provider_model_requested;
  std::string provider_model_executed;
  std::size_t provider_failure_count{0};
  std::string first_provider_failure_cell;
  double min_kappa_e_cm_inv_s{0.0};
  double max_kappa_e_cm_inv_s{0.0};
  double min_kappa_i_cm_inv_s{0.0};
  double max_kappa_i_cm_inv_s{0.0};
  double min_provider_lnLambda{0.0};
  double max_provider_lnLambda{0.0};
  double min_provider_f_LM{0.0};
  double max_provider_f_LM{0.0};
  double max_abs_delta_Te_erg{0.0};
  double max_abs_delta_Ti_erg{0.0};
  double max_electron_energy_residual{0.0};
  double max_ion_energy_residual{0.0};
  double global_thermal_energy_residual{0.0};
  double min_e_electron_after{0.0};
  double min_e_ion_after{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] ThermalConductionResult ApplyThermalConduction(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ThermalConductionOptions& options,
    const ThermalConductionSolveHook& solve_hook = {}) noexcept;

[[nodiscard]] ThermalConductionResult ApplyConstantKappaThermalConduction(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ThermalConductionOptions& options,
    const ThermalConductionSolveHook& solve_hook = {}) noexcept;

[[nodiscard]] bool ValidateThermalConductionDiagnostics(
    const ThermalConductionResult& result) noexcept;

[[nodiscard]] const char* ThermalConductionBackendName(ThermalConductionBackend backend) noexcept;
[[nodiscard]] const char* ThermalConductionKappaModelName(ThermalConductionKappaModel model) noexcept;

}  // namespace dec3d::transport
