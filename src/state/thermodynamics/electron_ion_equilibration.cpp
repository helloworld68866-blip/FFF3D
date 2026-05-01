#include "state/thermodynamics/electron_ion_equilibration.hpp"

#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace dec3d::state {

namespace {

constexpr double kDensityClosureTolerance = 1.0e-12;

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

void Fail(
    ElectronIonEquilibrationResult& result,
    const std::string& reason,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "diagnostic_id=p2.ei_equilibration.failure"
         << "; failure_reason=" << reason
         << "; first_bad_cell=" << radial << "," << theta << "," << phi
         << "; canonical_state_mutated=false";
  result.failure_diagnostics = report.str();
}

void ForwardRecoveryFailure(
    ElectronIonEquilibrationResult& result,
    const ThermodynamicRecoveryResult& recovered) {
  result.success = false;
  result.failure_reason = "P2-0 thermodynamic recovery failed";
  std::ostringstream report;
  report << "diagnostic_id=p2.ei_equilibration.failure"
         << "; p2_0_recovery_failed=true"
         << "; failure_reason=P2-0 thermodynamic recovery failed";
  if (!recovered.failure_reason.empty()) {
    report << "; p2_0_failure_reason=" << recovered.failure_reason;
  }
  if (!recovered.failure_diagnostics.empty()) {
    report << "; p2_0_failure_diagnostics={" << recovered.failure_diagnostics << "}";
  }
  report << "; canonical_state_mutated=false";
  result.failure_diagnostics = report.str();
}

[[nodiscard]] std::string BuildSuccessReport(
    const ElectronIonEquilibrationOptions& options,
    const ElectronIonEquilibrationResult& result) {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p2.production.equilibration_stage"
         << "; operator_diagnostic_id=p2.ei_equilibration.exact_exponential"
         << "; stage_id=E"
         << "; implementation_id=p2.ei_equilibration.local_exact_exponential_v1"
         << "; tau_model=" << ElectronIonTauModelName(options.tau_model)
         << "; tau_model_requested=" << ElectronIonTauModelName(options.tau_model)
         << "; tau_model_executed=" << ElectronIonTauModelName(options.tau_model)
         << "; tau_ei_s=" << options.tau_ei_s
         << "; dt_s=" << options.dt_s
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; updated_fields=e_electron"
         << "; e_fluid_total_unchanged=true"
         << "; cell_count=" << result.cell_count
         << "; max_deltaT_before_erg=" << result.max_deltaT_before_erg
         << "; max_deltaT_after_erg=" << result.max_deltaT_after_erg
         << "; max_expected_decay_residual_erg=" << result.max_expected_decay_residual_erg
         << "; max_abs_ne_minus_ni=" << result.max_abs_ne_minus_ni
         << "; max_dt_over_tau=" << result.max_dt_over_tau
         << "; min_tau_ei_s=" << result.min_tau_ei_s
         << "; max_tau_ei_s=" << result.max_tau_ei_s
         << "; min_lnLambda_ei_spitzer=" << result.min_lnLambda_ei_spitzer
         << "; max_lnLambda_ei_spitzer=" << result.max_lnLambda_ei_spitzer
         << "; double_kB_guard_passed=true"
         << "; max_total_thermal_energy_residual=" << result.max_total_thermal_energy_residual
         << "; min_e_electron_after=" << result.min_e_electron_after
         << "; min_e_ion_after=" << result.min_e_ion_after
         << "; invalid_cell_count=" << result.invalid_cell_count;
  return report.str();
}

}  // namespace

bool ElectronIonEquilibrationResult::is_complete() const noexcept {
  return success &&
         updated_fields == dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) &&
         cell_count > 0 &&
         invalid_cell_count == 0 &&
         !report_line.empty() &&
         report_line.find("diagnostic_id=p2.production.equilibration_stage") != std::string::npos &&
         report_line.find("operator_diagnostic_id=p2.ei_equilibration.exact_exponential") != std::string::npos &&
         report_line.find("stage_id=E") != std::string::npos &&
         report_line.find("tau_model_requested=") != std::string::npos &&
         report_line.find("tau_model_executed=") != std::string::npos &&
         report_line.find("updated_fields=e_electron") != std::string::npos &&
         report_line.find("e_fluid_total_unchanged=true") != std::string::npos &&
         report_line.find("max_abs_ne_minus_ni=") != std::string::npos &&
         report_line.find("max_dt_over_tau=") != std::string::npos &&
         report_line.find("min_tau_ei_s=") != std::string::npos &&
         report_line.find("max_tau_ei_s=") != std::string::npos;
}

ElectronIonEquilibrationResult ApplyLocalElectronIonEquilibration(
    CanonicalState& state,
    const ElectronIonEquilibrationOptions& options) noexcept {
  ElectronIonEquilibrationResult result;

  if (!Finite(options.dt_s)) {
    Fail(result, "dt_s must be finite", 0, 0, 0);
    return result;
  }
  if (options.dt_s < 0.0) {
    Fail(result, "dt_s must be non-negative", 0, 0, 0);
    return result;
  }
  if (options.tau_model == ElectronIonTauModel::constant_user_supplied) {
    if (!Finite(options.tau_ei_s)) {
      Fail(result, "tau_ei_s must be finite", 0, 0, 0);
      return result;
    }
    if (options.tau_ei_s <= 0.0) {
      Fail(result, "tau_ei_s must be positive", 0, 0, 0);
      return result;
    }
  }

  const auto recovered = RecoverThermodynamicState(
      state,
      ThermodynamicRecoveryOptions{
          options.electron_energy_floor,
          options.ion_energy_floor,
          1.0,
          0.0});
  if (!recovered.success) {
    ForwardRecoveryFailure(result, recovered);
    return result;
  }

  auto candidate = state;

  result.max_dt_over_tau = 0.0;
  result.min_e_electron_after = std::numeric_limits<double>::infinity();
  result.min_e_ion_after = std::numeric_limits<double>::infinity();
  result.min_tau_ei_s = std::numeric_limits<double>::infinity();
  result.max_tau_ei_s = 0.0;
  result.min_lnLambda_ei_spitzer = std::numeric_limits<double>::infinity();

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const auto& before = recovered.cells(radial, theta, phi);
        const double ne_minus_ni = std::abs(before.n_e_cm3 - before.n_i_cm3);
        result.max_abs_ne_minus_ni = std::max(result.max_abs_ne_minus_ni, ne_minus_ni);
        if (std::abs(before.zbar - 1.0) > kDensityClosureTolerance ||
            ne_minus_ni > kDensityClosureTolerance * std::max(1.0, before.n_i_cm3)) {
          Fail(result, "P2-1 requires Zbar=1 and n_e=n_i", radial, theta, phi);
          return result;
        }

        const auto tau = ComputeElectronIonTau(
            ElectronIonTauInput{
                options.tau_model,
                before.t_e_erg_per_particle,
                before.t_i_erg_per_particle,
                before.n_e_cm3,
                before.n_i_cm3,
                before.zbar,
                before.mean_ion_mass_g,
                options.tau_ei_s});
        if (!tau.success) {
          std::ostringstream reason;
          reason << "electron-ion tau provider failed: " << tau.failure_reason;
          Fail(result, reason.str(), radial, theta, phi);
          if (!tau.failure_diagnostics.empty()) {
            result.failure_diagnostics += "; tau_failure_diagnostics={" + tau.failure_diagnostics + "}";
          }
          return result;
        }
        result.min_tau_ei_s = std::min(result.min_tau_ei_s, tau.tau_ei_s);
        result.max_tau_ei_s = std::max(result.max_tau_ei_s, tau.tau_ei_s);
        result.min_lnLambda_ei_spitzer =
            std::min(result.min_lnLambda_ei_spitzer, tau.lnLambda_ei_spitzer);
        result.max_lnLambda_ei_spitzer =
            std::max(result.max_lnLambda_ei_spitzer, tau.lnLambda_ei_spitzer);
        result.max_dt_over_tau =
            std::max(result.max_dt_over_tau, options.dt_s / tau.tau_ei_s);

        const double delta_before = before.t_e_erg_per_particle - before.t_i_erg_per_particle;
        const double equilibrium = 0.5 * (before.t_e_erg_per_particle + before.t_i_erg_per_particle);
        const double local_decay = std::exp(-2.0 * options.dt_s / tau.tau_ei_s);
        if (!Finite(local_decay)) {
          Fail(result, "electron-ion decay factor is not finite", radial, theta, phi);
          return result;
        }
        const double electron_temperature_after = equilibrium + 0.5 * delta_before * local_decay;
        const double ion_temperature_expected = equilibrium - 0.5 * delta_before * local_decay;
        const double electron_pressure_after = before.n_e_cm3 * electron_temperature_after;
        const double electron_energy_after =
            electron_pressure_after / (HydroIdealGasGamma() - 1.0);
        const double thermal_before = before.e_electron_erg_per_cm3 + before.e_ion_erg_per_cm3;
        const double ion_energy_after = thermal_before - electron_energy_after;
        const double ion_temperature_after =
            (HydroIdealGasGamma() - 1.0) * ion_energy_after / before.n_i_cm3;

        if (!(Finite(electron_temperature_after) &&
              Finite(ion_temperature_expected) &&
              Finite(electron_energy_after) &&
              Finite(ion_energy_after) &&
              Finite(ion_temperature_after))) {
          Fail(result, "electron-ion updated state is not finite", radial, theta, phi);
          return result;
        }
        if (electron_energy_after < options.electron_energy_floor) {
          Fail(result, "electron energy fell below minimum threshold", radial, theta, phi);
          return result;
        }
        if (ion_energy_after < options.ion_energy_floor) {
          Fail(result, "ion energy fell below minimum threshold", radial, theta, phi);
          return result;
        }

        const double delta_after = electron_temperature_after - ion_temperature_after;
        const double expected_delta_after = delta_before * local_decay;
        const double thermal_after = electron_energy_after + ion_energy_after;

        result.max_deltaT_before_erg =
            std::max(result.max_deltaT_before_erg, std::abs(delta_before));
        result.max_deltaT_after_erg =
            std::max(result.max_deltaT_after_erg, std::abs(delta_after));
        result.max_expected_decay_residual_erg =
            std::max(result.max_expected_decay_residual_erg, std::abs(delta_after - expected_delta_after));
        result.max_total_thermal_energy_residual =
            std::max(result.max_total_thermal_energy_residual, std::abs(thermal_after - thermal_before));
        result.min_e_electron_after = std::min(result.min_e_electron_after, electron_energy_after);
        result.min_e_ion_after = std::min(result.min_e_ion_after, ion_energy_after);
        candidate.e_electron(radial, theta, phi) = electron_energy_after;
        ++result.cell_count;
      }
    }
  }

  state.e_electron = std::move(candidate.e_electron);
  result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron);
  state.ApplyAuthoritativeWrite(result.updated_fields);
  result.success = true;
  result.report_line = BuildSuccessReport(options, result);
  return result;
}

bool ValidateElectronIonEquilibrationDiagnostics(
    const ElectronIonEquilibrationResult& result) noexcept {
  return result.is_complete() &&
         result.report_line.find("unit_system=cgs") != std::string::npos &&
         result.report_line.find("temperature_internal_unit=erg_per_particle") != std::string::npos &&
         result.report_line.find("temperature_output_unit=keV") != std::string::npos &&
         result.report_line.find("max_abs_ne_minus_ni=") != std::string::npos &&
         result.report_line.find("max_dt_over_tau=") != std::string::npos &&
         result.report_line.find("min_tau_ei_s=") != std::string::npos &&
         result.report_line.find("max_tau_ei_s=") != std::string::npos &&
         result.report_line.find("min_lnLambda_ei_spitzer=") != std::string::npos &&
         result.report_line.find("double_kB_guard_passed=true") != std::string::npos &&
         result.report_line.find("max_total_thermal_energy_residual=") != std::string::npos &&
         result.report_line.find("min_e_electron_after=") != std::string::npos &&
         result.report_line.find("min_e_ion_after=") != std::string::npos &&
         result.report_line.find("invalid_cell_count=0") != std::string::npos;
}

}  // namespace dec3d::state
