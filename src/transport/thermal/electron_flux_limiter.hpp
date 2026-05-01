#pragma once

#include <string>

namespace dec3d::transport {

enum class ElectronFluxLimiterModel {
  disabled,
  minmax_old_time_face_effective_kappa,
  unsupported
};

struct ElectronFluxLimiterFaceInput {
  ElectronFluxLimiterModel model{ElectronFluxLimiterModel::disabled};
  double alpha_e{0.0};
  double Te_face_erg_per_particle{0.0};
  double ne_face_cm3{0.0};
  double kappa_face_cm_inv_s{0.0};
  double abs_grad_Te_erg_per_cm{0.0};
};

struct ElectronFluxLimiterFaceResult {
  bool success{false};
  bool limited{false};
  double scale{1.0};
  double q_spitzer_abs{0.0};
  double q_max{0.0};
  double flux_ratio_before_limit{0.0};
  double effective_kappa_cm_inv_s{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
};

[[nodiscard]] const char* ElectronFluxLimiterModelName(
    ElectronFluxLimiterModel model) noexcept;

[[nodiscard]] ElectronFluxLimiterFaceResult ApplyElectronFluxLimiter(
    const ElectronFluxLimiterFaceInput& input) noexcept;

[[nodiscard]] ElectronFluxLimiterFaceResult ApplyElectronFluxLimiterValues(
    const ElectronFluxLimiterFaceInput& input) noexcept;

[[nodiscard]] bool ValidateElectronFluxLimiterDiagnostics(
    const ElectronFluxLimiterFaceResult& result) noexcept;

}  // namespace dec3d::transport
