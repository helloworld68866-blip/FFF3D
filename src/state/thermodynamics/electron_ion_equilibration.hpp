#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/thermodynamics/electron_ion_tau.hpp"

#include <cstddef>
#include <string>

namespace dec3d::state {

struct ElectronIonEquilibrationOptions {
  double dt_s{0.0};
  ElectronIonTauModel tau_model{ElectronIonTauModel::constant_user_supplied};
  double tau_ei_s{0.0};
  double electron_energy_floor{0.0};
  double ion_energy_floor{0.0};
};

struct ElectronIonEquilibrationResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::size_t cell_count{0};
  double max_deltaT_before_erg{0.0};
  double max_deltaT_after_erg{0.0};
  double max_expected_decay_residual_erg{0.0};
  double max_abs_ne_minus_ni{0.0};
  double max_dt_over_tau{0.0};
  double min_tau_ei_s{0.0};
  double max_tau_ei_s{0.0};
  double min_lnLambda_ei_spitzer{0.0};
  double max_lnLambda_ei_spitzer{0.0};
  double max_total_thermal_energy_residual{0.0};
  double min_e_electron_after{0.0};
  double min_e_ion_after{0.0};
  std::size_t invalid_cell_count{0};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] ElectronIonEquilibrationResult ApplyLocalElectronIonEquilibration(
    CanonicalState& state,
    const ElectronIonEquilibrationOptions& options) noexcept;

[[nodiscard]] bool ValidateElectronIonEquilibrationDiagnostics(
    const ElectronIonEquilibrationResult& result) noexcept;

}  // namespace dec3d::state
