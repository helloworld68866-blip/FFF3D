#pragma once

#include "alpha/alpha_coefficients.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cstddef>
#include <string>

namespace dec3d::alpha {

struct OneGroupAlphaTransportOptions {
  double dt_s{0.0};
  AlphaCoefficientProviderOptions provider_options;
  dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy;
  double electron_energy_floor_erg_per_cm3{0.0};
  double ion_energy_floor_erg_per_cm3{0.0};
  double alpha_energy_floor_erg_per_cm3{0.0};
  dec3d::transport::GenericDiffusionReferenceSolveOptions serial_reference_options{};
};

struct OneGroupAlphaTransportResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string provider_report;
  std::string assembly_report;
  std::string solve_report;
  std::string thermodynamic_recovery_report;
  std::string backend_requested{"serial_dense_reference"};
  std::string backend_executed{"none"};
  std::size_t cell_count{0u};
  double dt_s{0.0};
  double delta_alpha_total{0.0};
  double delta_electron_total{0.0};
  double birth_source_total{0.0};
  double drag_deposition_total{0.0};
  double alpha_electron_exchange_residual{0.0};
  double alpha_equation_budget_residual{0.0};
  double global_alpha_plus_electron_budget_residual{0.0};
  double min_epsilon_alpha_before{0.0};
  double max_epsilon_alpha_before{0.0};
  double max_epsilon_alpha_after{0.0};
  double min_alpha_after{0.0};
  double min_e_electron_after{0.0};
  double min_e_ion_after{0.0};
  std::string reactivity_model_requested{"missing"};
  std::string reactivity_model_executed{"missing"};
  std::string bosch_hale_coefficients_source{"none"};
  std::string bosch_hale_temperature_source{"none"};
  double min_dt_reactivity_cm3_s{0.0};
  double max_dt_reactivity_cm3_s{0.0};
  bool thesis_reactivity_claim_allowed{false};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] dec3d::core::AuthoritativeFieldMask AlphaOperatorWriteMask() noexcept;

[[nodiscard]] OneGroupAlphaTransportResult ApplyOneGroupAlphaTransport(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupAlphaTransportOptions& options) noexcept;

[[nodiscard]] bool ValidateAlphaOneGroupOperatorDiagnostics(
    const OneGroupAlphaTransportResult& result) noexcept;

}  // namespace dec3d::alpha
