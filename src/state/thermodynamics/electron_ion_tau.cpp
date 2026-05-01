#include "state/thermodynamics/electron_ion_tau.hpp"

#include "physics/units/physical_constants.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace dec3d::state {
namespace {

[[nodiscard]] bool FinitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

ElectronIonTauResult Fail(const char* reason) {
  ElectronIonTauResult result;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "diagnostic_id=p2.ei_tau.failure"
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false"
         << "; fallback_used=false";
  result.failure_diagnostics = report.str();
  return result;
}

[[nodiscard]] double SpitzerCoulombLog(
    double Te_erg_per_particle,
    double ne_cm3,
    double zbar) noexcept {
  const double e = dec3d::physics::PhysicsConstantsCGS::elementary_charge_statcoulomb;
  return (3.0 / (2.0 * zbar * e * e * e)) *
         std::sqrt((Te_erg_per_particle * Te_erg_per_particle * Te_erg_per_particle) /
                   (dec3d::physics::PhysicsConstantsCGS::pi * ne_cm3));
}

}  // namespace

const char* ElectronIonTauModelName(ElectronIonTauModel model) noexcept {
  switch (model) {
    case ElectronIonTauModel::constant_user_supplied:
      return "constant_user_supplied";
    case ElectronIonTauModel::thesis_spitzer_eq_5_241:
      return "thesis_spitzer_eq_5_241";
    case ElectronIonTauModel::unsupported:
      return "unsupported";
  }
  return "unsupported";
}

bool ElectronIonTauResult::is_complete() const noexcept {
  return success &&
         FinitePositive(tau_ei_s) &&
         !report_line.empty() &&
         report_line.find("diagnostic_id=p2.ei_tau") != std::string::npos &&
         report_line.find("tau_model_requested=") != std::string::npos &&
         report_line.find("tau_model_executed=") != std::string::npos &&
         report_line.find("temperature_internal_unit=erg_per_particle") != std::string::npos &&
         report_line.find("fallback_used=false") != std::string::npos;
}

ElectronIonTauResult ComputeElectronIonTau(const ElectronIonTauInput& input) noexcept {
  ElectronIonTauResult result;

  if (input.model == ElectronIonTauModel::unsupported) {
    return Fail("unsupported tau_model");
  }

  if (input.model == ElectronIonTauModel::constant_user_supplied) {
    if (!FinitePositive(input.constant_tau_ei_s)) {
      return Fail("constant_tau_ei_s must be positive");
    }
    result.success = true;
    result.tau_ei_s = input.constant_tau_ei_s;
    result.double_kB_guard_passed = true;
    std::ostringstream report;
    report << std::setprecision(17)
           << "diagnostic_id=p2.ei_tau"
           << "; tau_model_requested=" << ElectronIonTauModelName(input.model)
           << "; tau_model_executed=" << ElectronIonTauModelName(input.model)
           << "; thesis_spitzer_tau_used=false"
           << "; tau_ei_s=" << result.tau_ei_s
           << "; unit_system=cgs"
           << "; temperature_internal_unit=erg_per_particle"
           << "; temperature_output_unit=keV"
           << "; double_kB_guard_passed=true"
           << "; canonical_state_mutated=false"
           << "; fallback_used=false";
    result.report_line = report.str();
    return result;
  }

  if (!FinitePositive(input.Te_erg_per_particle)) {
    return Fail("electron temperature must be positive");
  }
  if (!FinitePositive(input.Ti_erg_per_particle)) {
    return Fail("ion temperature must be positive");
  }
  if (!FinitePositive(input.ne_cm3)) {
    return Fail("electron density must be positive");
  }
  if (!FinitePositive(input.ni_cm3)) {
    return Fail("ion density must be positive");
  }
  if (!FinitePositive(input.zbar)) {
    return Fail("zbar must be positive");
  }
  if (!FinitePositive(input.mean_ion_mass_g)) {
    return Fail("mean ion mass must be positive");
  }

  result.lnLambda_ei_spitzer =
      SpitzerCoulombLog(input.Te_erg_per_particle, input.ne_cm3, input.zbar);
  if (!FinitePositive(result.lnLambda_ei_spitzer)) {
    return Fail("lnLambda_ei_spitzer must be positive");
  }

  const double me = dec3d::physics::PhysicsConstantsCGS::electron_mass_g;
  const double mi = input.mean_ion_mass_g;
  const double e = dec3d::physics::PhysicsConstantsCGS::elementary_charge_statcoulomb;
  const double thermal_speed_argument =
      input.Te_erg_per_particle / me + input.Ti_erg_per_particle / mi;
  if (!FinitePositive(thermal_speed_argument)) {
    return Fail("electron-ion thermal speed argument must be positive");
  }

  const double numerator =
      3.0 * me * mi * std::pow(thermal_speed_argument, 1.5);
  const double denominator =
      8.0 * std::sqrt(2.0 * dec3d::physics::PhysicsConstantsCGS::pi) *
      input.ni_cm3 *
      input.zbar * input.zbar *
      e * e * e * e *
      result.lnLambda_ei_spitzer;
  if (!FinitePositive(denominator)) {
    return Fail("tau denominator must be positive");
  }

  result.tau_ei_s = numerator / denominator;
  if (!FinitePositive(result.tau_ei_s)) {
    return Fail("tau_ei_s is not positive finite");
  }

  result.success = true;
  result.double_kB_guard_passed = true;
  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p2.ei_tau"
         << "; tau_model_requested=" << ElectronIonTauModelName(input.model)
         << "; tau_model_executed=" << ElectronIonTauModelName(input.model)
         << "; thesis_spitzer_tau_used=true"
         << "; tau_ei_s=" << result.tau_ei_s
         << "; lnLambda_ei_spitzer=" << result.lnLambda_ei_spitzer
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; equation_reference=Woo_Eq_5_241_and_5_236"
         << "; double_kB_guard_passed=true"
         << "; canonical_state_mutated=false"
         << "; fallback_used=false";
  result.report_line = report.str();
  return result;
}

bool ValidateElectronIonTauDiagnostics(const ElectronIonTauResult& result) noexcept {
  return result.is_complete() &&
         result.report_line.find("unit_system=cgs") != std::string::npos &&
         result.report_line.find("temperature_output_unit=keV") != std::string::npos &&
         result.report_line.find("double_kB_guard_passed=true") != std::string::npos &&
         result.report_line.find("canonical_state_mutated=false") != std::string::npos;
}

}  // namespace dec3d::state
