#pragma once

#include <string>

namespace dec3d::state {

enum class ElectronIonTauModel {
  constant_user_supplied,
  thesis_spitzer_eq_5_241,
  unsupported
};

struct ElectronIonTauInput {
  ElectronIonTauModel model{ElectronIonTauModel::constant_user_supplied};
  double Te_erg_per_particle{0.0};
  double Ti_erg_per_particle{0.0};
  double ne_cm3{0.0};
  double ni_cm3{0.0};
  double zbar{1.0};
  double mean_ion_mass_g{0.0};
  double constant_tau_ei_s{0.0};
};

struct ElectronIonTauResult {
  bool success{false};
  double tau_ei_s{0.0};
  double lnLambda_ei_spitzer{0.0};
  bool double_kB_guard_passed{false};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] const char* ElectronIonTauModelName(ElectronIonTauModel model) noexcept;

[[nodiscard]] ElectronIonTauResult ComputeElectronIonTau(
    const ElectronIonTauInput& input) noexcept;

[[nodiscard]] bool ValidateElectronIonTauDiagnostics(
    const ElectronIonTauResult& result) noexcept;

}  // namespace dec3d::state
