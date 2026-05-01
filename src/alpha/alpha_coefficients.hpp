#pragma once

#include "core/array/array3d.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <cstddef>
#include <string>

namespace dec3d::alpha {

enum class AlphaCompositionModel {
  missing,
  equimolar_dt_from_p2_recovery
};

enum class AlphaTauModel {
  missing,
  thesis_spitzer_eq_5_262
};

enum class AlphaReactivityModel {
  missing,
  constant_user_supplied,
  bosch_hale_dt
};

struct AlphaCoefficientProviderOptions {
  AlphaCompositionModel composition_model{AlphaCompositionModel::missing};
  AlphaTauModel tau_model{AlphaTauModel::missing};
  AlphaReactivityModel reactivity_model{AlphaReactivityModel::missing};
  dec3d::state::ThermodynamicRecoveryOptions recovery_options;
  double constant_dt_reactivity_cm3_s{0.0};
  double alpha_energy_floor_erg_cm3{0.0};
};

struct AlphaCoefficientArrays {
  dec3d::core::Array3D<double> tau_alphae_s;
  dec3d::core::Array3D<double> lnLambda_alphae_spitzer;
  dec3d::core::Array3D<double> lambda_drag_cm;
  dec3d::core::Array3D<double> D_alpha_cm2_s;
  dec3d::core::Array3D<double> birth_source_erg_cm3_s;
  dec3d::core::Array3D<double> dt_reactivity_cm3_s;
  dec3d::core::Array3D<double> nD_cm3;
  dec3d::core::Array3D<double> nT_cm3;
};

struct AlphaCoefficientProviderResult {
  bool success{false};
  std::size_t cell_count{0u};
  AlphaCoefficientArrays coefficients;
  dec3d::state::ThermodynamicRecoveryResult recovered;
  double min_tau_alphae_s{0.0};
  double max_tau_alphae_s{0.0};
  double min_D_alpha_cm2_s{0.0};
  double max_D_alpha_cm2_s{0.0};
  double min_birth_source_erg_cm3_s{0.0};
  double max_birth_source_erg_cm3_s{0.0};
  double min_lnLambda_alphae_spitzer{0.0};
  double max_lnLambda_alphae_spitzer{0.0};
  double min_dt_reactivity_cm3_s{0.0};
  double max_dt_reactivity_cm3_s{0.0};
  bool thesis_reactivity_claim_allowed{false};
  bool bosch_hale_reference_locked_locally{false};
  bool double_kB_guard_passed{false};
  std::string reactivity_model_requested{"missing"};
  std::string reactivity_model_executed{"missing"};
  std::string bosch_hale_coefficients_source{"none"};
  std::string bosch_hale_temperature_source{"none"};
  std::string thermodynamic_recovery_report;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] const char* AlphaCompositionModelName(AlphaCompositionModel model) noexcept;
[[nodiscard]] const char* AlphaTauModelName(AlphaTauModel model) noexcept;
[[nodiscard]] const char* AlphaReactivityModelName(AlphaReactivityModel model) noexcept;

[[nodiscard]] AlphaCoefficientProviderResult BuildAlphaCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const AlphaCoefficientProviderOptions& options) noexcept;

[[nodiscard]] double DtReactivityBoschHaleCm3PerS(double Ti_keV) noexcept;

[[nodiscard]] bool ValidateAlphaCoefficientDiagnostics(
    const AlphaCoefficientProviderResult& result) noexcept;

}  // namespace dec3d::alpha
