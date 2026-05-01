#pragma once

#include <string>

namespace dec3d::transport {

enum class ThermalConductivityModel {
  spitzer_no_degeneracy,
  lee_more_with_degeneracy,
  unsupported
};

struct ThermalConductivityInput {
  double Te_erg{0.0};
  double Ti_erg{0.0};
  double ne_cm3{0.0};
  double ni_cm3{0.0};
  double zbar{1.0};
  double mean_ion_mass_g{0.0};
  ThermalConductivityModel model{ThermalConductivityModel::spitzer_no_degeneracy};
};

struct ThermalConductivityResult {
  bool success{false};
  std::string failure_reason;
  std::string report_line;
  std::string failure_diagnostics;
  std::string model_requested;
  std::string model_executed;
  double Te_erg{0.0};
  double Ti_erg{0.0};
  double ne_cm3{0.0};
  double ni_cm3{0.0};
  double zbar{0.0};
  double mean_ion_mass_g{0.0};
  double TF_erg{0.0};
  double electron_thermal_velocity_cm_per_s{0.0};
  double bmax_cm{0.0};
  double bmin_spitzer_cm{0.0};
  double bmin_lee_more_cm{0.0};
  double bmin_used_cm{0.0};
  double lnLambda_raw{0.0};
  double lnLambda{0.0};
  double delta_prefactor{0.0};
  double f_LM{1.0};
  double kappa_e_cm_inv_s{0.0};
  double kappa_i_cm_inv_s{0.0};
  bool lee_more_floor_active{false};
  bool fallback_used{false};

  [[nodiscard]] bool is_complete() const noexcept;
};

struct ThermalConductivityValueResult {
  bool success{false};
  std::string failure_reason;
  double kappa_e_cm_inv_s{0.0};
  double kappa_i_cm_inv_s{0.0};
};

[[nodiscard]] const char* ThermalConductivityModelName(
    ThermalConductivityModel model) noexcept;

[[nodiscard]] ThermalConductivityResult ComputeThermalConductivity(
    const ThermalConductivityInput& input) noexcept;

[[nodiscard]] ThermalConductivityValueResult ComputeThermalConductivityValues(
    const ThermalConductivityInput& input) noexcept;

[[nodiscard]] bool ValidateThermalConductivityDiagnostics(
    const ThermalConductivityResult& result) noexcept;

}  // namespace dec3d::transport
