#include "transport/thermal/electron_flux_limiter.hpp"

#include "physics/units/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace dec3d::transport {
namespace {

[[nodiscard]] bool FinitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

}  // namespace

const char* ElectronFluxLimiterModelName(ElectronFluxLimiterModel model) noexcept {
  switch (model) {
    case ElectronFluxLimiterModel::disabled:
      return "disabled";
    case ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa:
      return "minmax_old_time_face_effective_kappa";
    case ElectronFluxLimiterModel::unsupported:
      return "unsupported";
  }
  return "unsupported";
}

namespace {

ElectronFluxLimiterFaceResult ApplyElectronFluxLimiterInternal(
    const ElectronFluxLimiterFaceInput& input,
    bool build_success_report) noexcept {
  ElectronFluxLimiterFaceResult result;
  auto fail = [&](const char* reason) {
    result.success = false;
    result.failure_reason = reason;
    std::ostringstream out;
    out << "diagnostic_id=p2.thermal_flux_limiter.failure"
        << "; failure_reason=" << reason
        << "; canonical_state_mutated=false"
        << "; fallback_used=false";
    result.failure_diagnostics = out.str();
    return result;
  };

  if (input.model == ElectronFluxLimiterModel::unsupported) {
    return fail("unsupported electron flux limiter model");
  }
  if (input.model == ElectronFluxLimiterModel::disabled) {
    result.success = true;
    result.effective_kappa_cm_inv_s = input.kappa_face_cm_inv_s;
    if (build_success_report) {
      result.report_line =
          "diagnostic_id=p2.thermal_flux_limiter; electron_flux_limiter_enabled=false; fallback_used=false";
    }
    return result;
  }
  if (!FinitePositive(input.alpha_e)) {
    return fail("alpha_e must be user-supplied and positive");
  }
  if (!FinitePositive(input.Te_face_erg_per_particle) ||
      !FinitePositive(input.ne_face_cm3) ||
      !FinitePositive(input.kappa_face_cm_inv_s) ||
      !std::isfinite(input.abs_grad_Te_erg_per_cm) ||
      input.abs_grad_Te_erg_per_cm < 0.0) {
    return fail("electron flux limiter inputs are invalid");
  }

  result.q_spitzer_abs = input.kappa_face_cm_inv_s * input.abs_grad_Te_erg_per_cm;
  const double v_te =
      std::sqrt(3.0 * input.Te_face_erg_per_particle /
                dec3d::physics::PhysicsConstantsCGS::electron_mass_g);
  result.q_max =
      input.alpha_e * input.ne_face_cm3 * input.Te_face_erg_per_particle * v_te;
  result.flux_ratio_before_limit =
      result.q_spitzer_abs / std::max(result.q_max, 1.0e-300);
  result.scale = result.q_spitzer_abs > 0.0
      ? std::min(1.0, result.q_max / result.q_spitzer_abs)
      : 1.0;
  result.limited = result.scale < 1.0;
  result.effective_kappa_cm_inv_s = result.scale * input.kappa_face_cm_inv_s;
  result.success = true;

  if (build_success_report) {
    std::ostringstream out;
    out << std::setprecision(17)
        << "diagnostic_id=p2.thermal_flux_limiter"
        << "; electron_flux_limiter_enabled=true"
        << "; electron_flux_limiter_model=" << ElectronFluxLimiterModelName(input.model)
        << "; alpha_e=" << input.alpha_e
        << "; alpha_e_source=user_supplied"
        << "; q_spitzer_abs=" << result.q_spitzer_abs
        << "; q_max=" << result.q_max
        << "; flux_ratio_before_limit=" << result.flux_ratio_before_limit
        << "; limiter_scale=" << result.scale
        << "; limited=" << (result.limited ? "true" : "false")
        << "; effective_kappa_cm_inv_s=" << result.effective_kappa_cm_inv_s
        << "; canonical_state_mutated=false"
        << "; fallback_used=false";
    result.report_line = out.str();
  }
  return result;
}

}  // namespace

ElectronFluxLimiterFaceResult ApplyElectronFluxLimiter(
    const ElectronFluxLimiterFaceInput& input) noexcept {
  return ApplyElectronFluxLimiterInternal(input, true);
}

ElectronFluxLimiterFaceResult ApplyElectronFluxLimiterValues(
    const ElectronFluxLimiterFaceInput& input) noexcept {
  return ApplyElectronFluxLimiterInternal(input, false);
}

bool ValidateElectronFluxLimiterDiagnostics(
    const ElectronFluxLimiterFaceResult& result) noexcept {
  const auto& line = result.report_line;
  if (!result.success) {
    return !result.failure_diagnostics.empty() &&
           result.failure_diagnostics.find("canonical_state_mutated=false") != std::string::npos;
  }
  if (line.find("electron_flux_limiter_enabled=false") != std::string::npos) {
    return line.find("diagnostic_id=p2.thermal_flux_limiter") != std::string::npos &&
           line.find("fallback_used=false") != std::string::npos;
  }
  return line.find("diagnostic_id=p2.thermal_flux_limiter") != std::string::npos &&
         line.find("electron_flux_limiter_model=") != std::string::npos &&
         line.find("alpha_e_source=user_supplied") != std::string::npos &&
         line.find("limiter_scale=") != std::string::npos &&
         line.find("fallback_used=false") != std::string::npos;
}

}  // namespace dec3d::transport
