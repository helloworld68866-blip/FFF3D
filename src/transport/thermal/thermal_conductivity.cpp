#include "transport/thermal/thermal_conductivity.hpp"

#include "physics/units/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace dec3d::transport {
namespace {

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

void Fail(ThermalConductivityResult& result, const std::string& reason) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p2.thermal_conductivity.coefficient.failure"
         << "; model_requested=" << result.model_requested
         << "; fallback_used=false"
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false";
  result.failure_diagnostics = report.str();
}

[[nodiscard]] double DeltaPrefactor(double zbar) noexcept {
  return 0.095 * (zbar + 0.24) / (1.0 + 0.24 * zbar);
}

[[nodiscard]] double FermiTemperatureErg(double ne_cm3) noexcept {
  using dec3d::physics::PhysicsConstantsCGS;
  return (PhysicsConstantsCGS::hbar_erg_s * PhysicsConstantsCGS::hbar_erg_s) /
         (2.0 * PhysicsConstantsCGS::electron_mass_g) *
         std::pow(3.0 * PhysicsConstantsCGS::pi * PhysicsConstantsCGS::pi * ne_cm3, 2.0 / 3.0);
}

[[nodiscard]] double ElectronThermalVelocity(double Te_erg) noexcept {
  return std::sqrt(3.0 * Te_erg / dec3d::physics::PhysicsConstantsCGS::electron_mass_g);
}

[[nodiscard]] double DebyeHuckelBmax(
    double Te_erg,
    double Ti_erg,
    double ne_cm3,
    double ni_cm3,
    double zbar,
    double TF_erg) noexcept {
  using dec3d::physics::PhysicsConstantsCGS;
  const double e = PhysicsConstantsCGS::elementary_charge_statcoulomb;
  const double inv_bmax_sq =
      4.0 * PhysicsConstantsCGS::pi * ne_cm3 * e * e /
          std::sqrt(Te_erg * Te_erg + TF_erg * TF_erg) +
      4.0 * PhysicsConstantsCGS::pi * ni_cm3 * zbar * zbar * e * e /
          Ti_erg;
  return 1.0 / std::sqrt(inv_bmax_sq);
}

[[nodiscard]] double SpitzerBmax(
    double Te_erg,
    double Ti_erg,
    double ne_cm3,
    double ni_cm3,
    double zbar) noexcept {
  return DebyeHuckelBmax(Te_erg, Ti_erg, ne_cm3, ni_cm3, zbar, 0.0);
}

[[nodiscard]] double LeeMoreBmax(
    double Te_erg,
    double Ti_erg,
    double ne_cm3,
    double ni_cm3,
    double zbar,
    double TF_erg) noexcept {
  return DebyeHuckelBmax(Te_erg, Ti_erg, ne_cm3, ni_cm3, zbar, TF_erg);
}

[[nodiscard]] double SpitzerBmin(double Te_erg, double zbar) noexcept {
  const double e = dec3d::physics::PhysicsConstantsCGS::elementary_charge_statcoulomb;
  return zbar * e * e / (3.0 * Te_erg);
}

[[nodiscard]] double DebroglieWavelength(double velocity_cm_per_s) noexcept {
  using dec3d::physics::PhysicsConstantsCGS;
  return 2.0 * PhysicsConstantsCGS::pi * PhysicsConstantsCGS::hbar_erg_s /
         (PhysicsConstantsCGS::electron_mass_g * velocity_cm_per_s);
}

[[nodiscard]] double InternalElectronConductivity(
    double Te_erg,
    double zbar,
    double lnLambda,
    double delta,
    double f_LM) noexcept {
  using dec3d::physics::PhysicsConstantsCGS;
  return 20.0 * std::pow(2.0 / PhysicsConstantsCGS::pi, 1.5) *
         std::pow(Te_erg, 2.5) /
         (std::sqrt(PhysicsConstantsCGS::electron_mass_g) *
          zbar *
          std::pow(PhysicsConstantsCGS::elementary_charge_statcoulomb, 4.0) *
          lnLambda) *
         delta *
         f_LM;
}

[[nodiscard]] bool ValidateInput(
    const ThermalConductivityInput& input,
    ThermalConductivityResult& result) {
  if (input.model == ThermalConductivityModel::unsupported) {
    Fail(result, "unsupported conductivity model");
    return false;
  }
  if (!Finite(input.Te_erg) || !Finite(input.Ti_erg) ||
      !Finite(input.ne_cm3) || !Finite(input.ni_cm3) ||
      !Finite(input.zbar) || !Finite(input.mean_ion_mass_g)) {
    Fail(result, "thermal conductivity input must be finite");
    return false;
  }
  if (input.Te_erg <= 0.0) {
    Fail(result, "Te_erg must be positive");
    return false;
  }
  if (input.Ti_erg <= 0.0) {
    Fail(result, "Ti_erg must be positive");
    return false;
  }
  if (input.ne_cm3 <= 0.0) {
    Fail(result, "ne_cm3 must be positive");
    return false;
  }
  if (input.ni_cm3 <= 0.0) {
    Fail(result, "ni_cm3 must be positive");
    return false;
  }
  if (input.zbar <= 0.0) {
    Fail(result, "zbar must be positive");
    return false;
  }
  if (input.mean_ion_mass_g <= 0.0) {
    Fail(result, "mean_ion_mass_g must be positive");
    return false;
  }
  return true;
}

std::string BuildReport(const ThermalConductivityResult& result) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p2.thermal_conductivity.coefficient"
         << "; implementation_id=p2.thermal_conductivity.spitzer_lee_more_v1"
         << "; unit_system=cgs"
         << "; constants_source=PhysicsConstantsCGS"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; kappa_internal_unit=cm^-1_s^-1"
         << "; model_requested=" << result.model_requested
         << "; model_executed=" << result.model_executed
         << "; fallback_used=false"
         << "; zbar=" << result.zbar
         << "; mean_ion_mass_g=" << result.mean_ion_mass_g
         << "; Te_erg=" << result.Te_erg
         << "; Ti_erg=" << result.Ti_erg
         << "; ne_cm3=" << result.ne_cm3
         << "; ni_cm3=" << result.ni_cm3
         << "; TF_erg=" << result.TF_erg
         << "; electron_thermal_velocity_cm_per_s=" << result.electron_thermal_velocity_cm_per_s
         << "; bmax_cm=" << result.bmax_cm
         << "; bmin_spitzer_cm=" << result.bmin_spitzer_cm
         << "; bmin_lee_more_cm=" << result.bmin_lee_more_cm
         << "; bmin_used_cm=" << result.bmin_used_cm
         << "; lnLambda_raw=" << result.lnLambda_raw
         << "; lnLambda=" << result.lnLambda
         << "; lee_more_floor_active=" << (result.lee_more_floor_active ? "true" : "false")
         << "; delta_prefactor=" << result.delta_prefactor
         << "; f_LM=" << result.f_LM
         << "; kappa_e_cm_inv_s=" << result.kappa_e_cm_inv_s
         << "; kappa_i_cm_inv_s=" << result.kappa_i_cm_inv_s
         << "; ion_mass_ratio_rule=kappa_i_eq_kappa_e_sqrt_me_over_mi"
         << "; canonical_state_mutated=false";
  return report.str();
}

}  // namespace

const char* ThermalConductivityModelName(ThermalConductivityModel model) noexcept {
  switch (model) {
    case ThermalConductivityModel::spitzer_no_degeneracy:
      return "spitzer_no_degeneracy";
    case ThermalConductivityModel::lee_more_with_degeneracy:
      return "lee_more_with_degeneracy";
    case ThermalConductivityModel::unsupported:
      return "unsupported";
  }
  return "unsupported";
}

bool ThermalConductivityResult::is_complete() const noexcept {
  return success &&
         !report_line.empty() &&
         report_line.find("diagnostic_id=p2.thermal_conductivity.coefficient") != std::string::npos &&
         report_line.find("implementation_id=p2.thermal_conductivity.spitzer_lee_more_v1") != std::string::npos &&
         report_line.find("fallback_used=false") != std::string::npos &&
         report_line.find("canonical_state_mutated=false") != std::string::npos;
}

ThermalConductivityResult ComputeThermalConductivityInternal(
    const ThermalConductivityInput& input,
    bool build_success_report) noexcept {
  using dec3d::physics::PhysicsConstantsCGS;

  ThermalConductivityResult result;
  result.model_requested = ThermalConductivityModelName(input.model);
  result.Te_erg = input.Te_erg;
  result.Ti_erg = input.Ti_erg;
  result.ne_cm3 = input.ne_cm3;
  result.ni_cm3 = input.ni_cm3;
  result.zbar = input.zbar;
  result.mean_ion_mass_g = input.mean_ion_mass_g;

  if (!ValidateInput(input, result)) {
    return result;
  }

  result.delta_prefactor = DeltaPrefactor(input.zbar);
  result.TF_erg = FermiTemperatureErg(input.ne_cm3);
  result.electron_thermal_velocity_cm_per_s = ElectronThermalVelocity(input.Te_erg);
  result.bmax_cm =
      input.model == ThermalConductivityModel::lee_more_with_degeneracy
          ? LeeMoreBmax(input.Te_erg, input.Ti_erg, input.ne_cm3, input.ni_cm3, input.zbar, result.TF_erg)
          : SpitzerBmax(input.Te_erg, input.Ti_erg, input.ne_cm3, input.ni_cm3, input.zbar);
  result.bmin_spitzer_cm = SpitzerBmin(input.Te_erg, input.zbar);
  result.bmin_lee_more_cm = std::min(
      result.bmin_spitzer_cm, 0.5 * DebroglieWavelength(result.electron_thermal_velocity_cm_per_s));

  result.model_executed = ThermalConductivityModelName(input.model);
  if (!Finite(result.delta_prefactor) || result.delta_prefactor <= 0.0 ||
      !Finite(result.TF_erg) || result.TF_erg < 0.0 ||
      !Finite(result.electron_thermal_velocity_cm_per_s) ||
      result.electron_thermal_velocity_cm_per_s <= 0.0) {
    Fail(result, "thermal conductivity intermediate state is not physical");
    return result;
  }

  if (input.model == ThermalConductivityModel::lee_more_with_degeneracy) {
    result.bmin_used_cm = result.bmin_lee_more_cm;
    if (!Finite(result.bmax_cm) || result.bmax_cm <= 0.0 ||
        !Finite(result.bmin_used_cm) || result.bmin_used_cm <= 0.0) {
      Fail(result, "impact parameters must be positive and finite");
      return result;
    }
    const double ratio = result.bmax_cm / result.bmin_used_cm;
    result.lnLambda_raw = 0.5 * std::log(1.0 + ratio * ratio);
    result.lnLambda = std::max(result.lnLambda_raw, 2.0);
    result.lee_more_floor_active = result.lnLambda_raw < 2.0;
    result.f_LM =
        1.0 +
        (3.0 * std::pow(PhysicsConstantsCGS::pi, 5.0) / 51200.0) *
            std::pow(result.TF_erg / input.Te_erg, 3.0) /
            (result.delta_prefactor * result.delta_prefactor);
  } else {
    result.f_LM = 1.0;
    result.bmin_used_cm = result.bmin_spitzer_cm;
    if (!Finite(result.bmax_cm) || result.bmax_cm <= 0.0 ||
        !Finite(result.bmin_used_cm) || result.bmin_used_cm <= 0.0) {
      Fail(result, "impact parameters must be positive and finite");
      return result;
    }
    if (result.bmax_cm <= result.bmin_used_cm) {
      Fail(result, "Spitzer Coulomb logarithm requires bmax greater than bmin");
      return result;
    }
    result.lnLambda_raw = std::log(result.bmax_cm / result.bmin_used_cm);
    result.lnLambda = result.lnLambda_raw;
  }
  if (!Finite(result.lnLambda) || result.lnLambda <= 0.0) {
    Fail(result, "lnLambda must be positive and finite");
    return result;
  }
  if (!Finite(result.f_LM) || result.f_LM <= 0.0) {
    Fail(result, "f_LM must be positive and finite");
    return result;
  }

  result.kappa_e_cm_inv_s = InternalElectronConductivity(
      input.Te_erg, input.zbar, result.lnLambda, result.delta_prefactor, result.f_LM);
  result.kappa_i_cm_inv_s =
      result.kappa_e_cm_inv_s *
      std::sqrt(PhysicsConstantsCGS::electron_mass_g / input.mean_ion_mass_g);
  if (!Finite(result.kappa_e_cm_inv_s) || result.kappa_e_cm_inv_s < 0.0 ||
      !Finite(result.kappa_i_cm_inv_s) || result.kappa_i_cm_inv_s < 0.0) {
    Fail(result, "conductivity outputs must be non-negative and finite");
    return result;
  }

  result.success = true;
  result.fallback_used = false;
  if (build_success_report) {
    result.report_line = BuildReport(result);
  }
  return result;
}

ThermalConductivityResult ComputeThermalConductivity(
    const ThermalConductivityInput& input) noexcept {
  return ComputeThermalConductivityInternal(input, true);
}

ThermalConductivityValueResult ComputeThermalConductivityValues(
    const ThermalConductivityInput& input) noexcept {
  const auto detailed = ComputeThermalConductivityInternal(input, false);
  ThermalConductivityValueResult result;
  result.success = detailed.success;
  result.failure_reason = detailed.failure_reason;
  result.kappa_e_cm_inv_s = detailed.kappa_e_cm_inv_s;
  result.kappa_i_cm_inv_s = detailed.kappa_i_cm_inv_s;
  return result;
}

bool ValidateThermalConductivityDiagnostics(
    const ThermalConductivityResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }
  const auto& line = result.report_line;
  const char* required[] = {
      "diagnostic_id=p2.thermal_conductivity.coefficient",
      "implementation_id=p2.thermal_conductivity.spitzer_lee_more_v1",
      "unit_system=cgs",
      "constants_source=PhysicsConstantsCGS",
      "temperature_internal_unit=erg_per_particle",
      "temperature_output_unit=keV",
      "kappa_internal_unit=cm^-1_s^-1",
      "model_requested=",
      "model_executed=",
      "fallback_used=false",
      "zbar=",
      "mean_ion_mass_g=",
      "Te_erg=",
      "Ti_erg=",
      "ne_cm3=",
      "ni_cm3=",
      "TF_erg=",
      "electron_thermal_velocity_cm_per_s=",
      "bmax_cm=",
      "bmin_spitzer_cm=",
      "bmin_lee_more_cm=",
      "bmin_used_cm=",
      "lnLambda_raw=",
      "lnLambda=",
      "lee_more_floor_active=",
      "delta_prefactor=",
      "f_LM=",
      "kappa_e_cm_inv_s=",
      "kappa_i_cm_inv_s=",
      "ion_mass_ratio_rule=kappa_i_eq_kappa_e_sqrt_me_over_mi",
      "canonical_state_mutated=false"};
  for (const char* token : required) {
    if (line.find(token) == std::string::npos) {
      return false;
    }
  }
  return true;
}

}  // namespace dec3d::transport
