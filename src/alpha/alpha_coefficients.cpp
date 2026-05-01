#include "alpha/alpha_coefficients.hpp"

#include "physics/units/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::alpha {
namespace {

[[nodiscard]] bool FinitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

AlphaCoefficientProviderResult Fail(const std::string& reason) {
  AlphaCoefficientProviderResult result;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "diagnostic_id=p4.alpha.coefficient_provider.failure"
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false"
         << "; fallback_used=false";
  result.failure_diagnostics = report.str();
  return result;
}

[[nodiscard]] double SpitzerAlphaElectronCoulombLog(
    double Te_erg_per_particle,
    double ne_cm3) noexcept {
  const double e = dec3d::physics::PhysicsConstantsCGS::elementary_charge_statcoulomb;
  const double z_alpha = dec3d::physics::PhysicsConstantsCGS::alpha_charge_number;
  return (3.0 / (2.0 * z_alpha * e * e * e)) *
         std::sqrt((Te_erg_per_particle * Te_erg_per_particle * Te_erg_per_particle) /
                   (dec3d::physics::PhysicsConstantsCGS::pi * ne_cm3));
}

[[nodiscard]] double AlphaTauEq5262(
    double Te_erg_per_particle,
    double ne_cm3,
    double lnLambda_alphae) noexcept {
  const double me = dec3d::physics::PhysicsConstantsCGS::electron_mass_g;
  const double ma = dec3d::physics::PhysicsConstantsCGS::alpha_mass_g;
  const double z_alpha = dec3d::physics::PhysicsConstantsCGS::alpha_charge_number;
  const double e = dec3d::physics::PhysicsConstantsCGS::elementary_charge_statcoulomb;
  const double numerator = 3.0 * ma * std::pow(Te_erg_per_particle, 1.5);
  const double denominator =
      8.0 * std::sqrt(2.0 * dec3d::physics::PhysicsConstantsCGS::pi * me) *
      ne_cm3 *
      z_alpha * z_alpha *
      e * e * e * e *
      lnLambda_alphae;
  if (!FinitePositive(denominator)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return numerator / denominator;
}

[[nodiscard]] double AlphaBirthSpeedCmPerS() noexcept {
  return std::sqrt(2.0 * dec3d::physics::PhysicsConstantsCGS::alpha_birth_energy_erg /
                   dec3d::physics::PhysicsConstantsCGS::alpha_mass_g);
}

}  // namespace

double DtReactivityBoschHaleCm3PerS(double Ti_keV) noexcept {
  if (!(std::isfinite(Ti_keV) && Ti_keV > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  constexpr double b_G = 34.3827;
  constexpr double mrc2 = 1124656.0;
  const double T = Ti_keV;
  const double num = 1.51361e-2 * T +
                     4.60643e-3 * T * T -
                     1.06750e-4 * T * T * T;
  const double den = 1.0 +
                     7.51886e-2 * T +
                     1.35000e-2 * T * T +
                     1.3600e-5 * T * T * T;
  if (!(std::isfinite(den) && den != 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double correction = 1.0 - num / den;
  if (!(std::isfinite(correction) && correction > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double theta = T / correction;
  if (!(std::isfinite(theta) && theta > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double zeta = std::pow((b_G * b_G) / (4.0 * theta), 1.0 / 3.0);
  if (!(std::isfinite(zeta) && zeta > 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double reactivity =
      1.17302e-9 * theta *
      std::sqrt(zeta / (mrc2 * T * T * T)) *
      std::exp(-3.0 * zeta);
  return (std::isfinite(reactivity) && reactivity > 0.0)
             ? reactivity
             : std::numeric_limits<double>::quiet_NaN();
}

const char* AlphaCompositionModelName(AlphaCompositionModel model) noexcept {
  switch (model) {
    case AlphaCompositionModel::missing:
      return "missing";
    case AlphaCompositionModel::equimolar_dt_from_p2_recovery:
      return "equimolar_dt_from_p2_recovery";
  }
  return "missing";
}

const char* AlphaTauModelName(AlphaTauModel model) noexcept {
  switch (model) {
    case AlphaTauModel::missing:
      return "missing";
    case AlphaTauModel::thesis_spitzer_eq_5_262:
      return "thesis_spitzer_eq_5_262";
  }
  return "missing";
}

const char* AlphaReactivityModelName(AlphaReactivityModel model) noexcept {
  switch (model) {
    case AlphaReactivityModel::missing:
      return "missing";
    case AlphaReactivityModel::constant_user_supplied:
      return "constant_user_supplied";
    case AlphaReactivityModel::bosch_hale_dt:
      return "bosch_hale_dt";
  }
  return "missing";
}

bool AlphaCoefficientProviderResult::is_complete() const noexcept {
  return success &&
         cell_count > 0u &&
         !coefficients.tau_alphae_s.empty() &&
         !coefficients.D_alpha_cm2_s.empty() &&
         !coefficients.birth_source_erg_cm3_s.empty() &&
         double_kB_guard_passed &&
         !thermodynamic_recovery_report.empty() &&
         !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p4.alpha.coefficient_provider") &&
         Contains(report_line, "phase_id=P4") &&
         Contains(report_line, "stage_id=A") &&
         Contains(report_line, "composition_model=") &&
         Contains(report_line, "tau_alphae_model=") &&
         Contains(report_line, "reactivity_model=") &&
         Contains(report_line, "reactivity_model_requested=") &&
         Contains(report_line, "reactivity_model_executed=") &&
         Contains(report_line, "dt_reactivity_internal_unit=cm3_s") &&
         Contains(report_line, "min_dt_reactivity_cm3_s=") &&
         Contains(report_line, "max_dt_reactivity_cm3_s=") &&
         Contains(report_line, "canonical_state_mutated=false") &&
         Contains(report_line, "fallback_used=false");
}

AlphaCoefficientProviderResult BuildAlphaCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const AlphaCoefficientProviderOptions& options) noexcept {
  if (options.composition_model != AlphaCompositionModel::equimolar_dt_from_p2_recovery) {
    return Fail("composition_model must be equimolar_dt_from_p2_recovery");
  }
  if (options.tau_model != AlphaTauModel::thesis_spitzer_eq_5_262) {
    return Fail("tau_alphae_model must be thesis_spitzer_eq_5_262");
  }
  if (options.reactivity_model == AlphaReactivityModel::missing) {
    return Fail("reactivity_model is missing");
  }
  if (options.reactivity_model == AlphaReactivityModel::constant_user_supplied &&
      !FinitePositive(options.constant_dt_reactivity_cm3_s)) {
    return Fail("constant_dt_reactivity_cm3_s must be positive");
  }
  if (!(std::isfinite(options.alpha_energy_floor_erg_cm3) &&
        options.alpha_energy_floor_erg_cm3 >= 0.0)) {
    return Fail("alpha_energy_floor_erg_cm3 must be finite and nonnegative");
  }
  if (!state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state)) {
    return Fail("canonical alpha_state storage is missing");
  }

  AlphaCoefficientProviderResult result;
  result.recovered = dec3d::state::RecoverThermodynamicState(state, options.recovery_options);
  result.thermodynamic_recovery_report = result.recovered.recovery_diagnostics;
  if (!result.recovered.success) {
    const std::string recovery_report = result.recovered.recovery_diagnostics;
    result = Fail("thermodynamic recovery failed");
    result.thermodynamic_recovery_report = recovery_report;
    return result;
  }

  result.reactivity_model_requested = AlphaReactivityModelName(options.reactivity_model);
  result.reactivity_model_executed = AlphaReactivityModelName(options.reactivity_model);
  if (options.reactivity_model == AlphaReactivityModel::bosch_hale_dt) {
    result.bosch_hale_reference_locked_locally = true;
    result.thesis_reactivity_claim_allowed = true;
    result.bosch_hale_coefficients_source = "project_owner_2026_04_28";
    result.bosch_hale_temperature_source = "Ti_old";
  }

  const std::size_t nr = state.layout.radial_cells;
  const std::size_t nt = state.layout.theta_cells;
  const std::size_t np = state.layout.phi_cells;

  result.coefficients.tau_alphae_s = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.lnLambda_alphae_spitzer = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.lambda_drag_cm = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.D_alpha_cm2_s = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.birth_source_erg_cm3_s = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.dt_reactivity_cm3_s = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.nD_cm3 = dec3d::core::Array3D<double>(nr, nt, np, 0.0);
  result.coefficients.nT_cm3 = dec3d::core::Array3D<double>(nr, nt, np, 0.0);

  result.min_tau_alphae_s = std::numeric_limits<double>::infinity();
  result.max_tau_alphae_s = -std::numeric_limits<double>::infinity();
  result.min_D_alpha_cm2_s = std::numeric_limits<double>::infinity();
  result.max_D_alpha_cm2_s = -std::numeric_limits<double>::infinity();
  result.min_birth_source_erg_cm3_s = std::numeric_limits<double>::infinity();
  result.max_birth_source_erg_cm3_s = -std::numeric_limits<double>::infinity();
  result.min_lnLambda_alphae_spitzer = std::numeric_limits<double>::infinity();
  result.max_lnLambda_alphae_spitzer = -std::numeric_limits<double>::infinity();
  result.min_dt_reactivity_cm3_s = std::numeric_limits<double>::infinity();
  result.max_dt_reactivity_cm3_s = -std::numeric_limits<double>::infinity();

  const double v_alpha0 = AlphaBirthSpeedCmPerS();

  for (std::size_t r = 0; r < nr; ++r) {
    for (std::size_t t = 0; t < nt; ++t) {
      for (std::size_t p = 0; p < np; ++p) {
        const auto& cell = result.recovered.cells(r, t, p);
        if (!FinitePositive(cell.t_e_erg_per_particle) || !FinitePositive(cell.n_e_cm3) ||
            !FinitePositive(cell.n_i_cm3)) {
          return Fail("recovered alpha provider inputs must be positive");
        }
        if (!(std::isfinite(state.alpha_state.storage(r, t, p)) &&
              state.alpha_state.storage(r, t, p) >= options.alpha_energy_floor_erg_cm3)) {
          return Fail("alpha energy density must be finite and above floor");
        }

        const double nD = 0.5 * cell.n_i_cm3;
        const double nT = 0.5 * cell.n_i_cm3;
        const double lnLambda =
            SpitzerAlphaElectronCoulombLog(cell.t_e_erg_per_particle, cell.n_e_cm3);
        const double tau =
            AlphaTauEq5262(cell.t_e_erg_per_particle, cell.n_e_cm3, lnLambda);
        const double lambda_drag = v_alpha0 * tau / 9.0;
        const double D_alpha = v_alpha0 * lambda_drag;
        double dt_reactivity = options.constant_dt_reactivity_cm3_s;
        if (options.reactivity_model == AlphaReactivityModel::bosch_hale_dt) {
          dt_reactivity =
              DtReactivityBoschHaleCm3PerS(dec3d::physics::KeVFromErg(
                  cell.t_i_erg_per_particle));
        }
        const double birth_source =
            nD * nT * dt_reactivity *
            dec3d::physics::PhysicsConstantsCGS::alpha_birth_energy_erg;

        if (!FinitePositive(nD) || !FinitePositive(nT) || !FinitePositive(lnLambda) ||
            !FinitePositive(tau) || !FinitePositive(lambda_drag) ||
            !FinitePositive(D_alpha) || !FinitePositive(dt_reactivity) ||
            !FinitePositive(birth_source)) {
          return Fail("alpha coefficient output is not positive finite");
        }

        result.coefficients.nD_cm3(r, t, p) = nD;
        result.coefficients.nT_cm3(r, t, p) = nT;
        result.coefficients.lnLambda_alphae_spitzer(r, t, p) = lnLambda;
        result.coefficients.tau_alphae_s(r, t, p) = tau;
        result.coefficients.lambda_drag_cm(r, t, p) = lambda_drag;
        result.coefficients.D_alpha_cm2_s(r, t, p) = D_alpha;
        result.coefficients.birth_source_erg_cm3_s(r, t, p) = birth_source;
        result.coefficients.dt_reactivity_cm3_s(r, t, p) = dt_reactivity;

        result.min_tau_alphae_s = std::min(result.min_tau_alphae_s, tau);
        result.max_tau_alphae_s = std::max(result.max_tau_alphae_s, tau);
        result.min_D_alpha_cm2_s = std::min(result.min_D_alpha_cm2_s, D_alpha);
        result.max_D_alpha_cm2_s = std::max(result.max_D_alpha_cm2_s, D_alpha);
        result.min_birth_source_erg_cm3_s =
            std::min(result.min_birth_source_erg_cm3_s, birth_source);
        result.max_birth_source_erg_cm3_s =
            std::max(result.max_birth_source_erg_cm3_s, birth_source);
        result.min_lnLambda_alphae_spitzer =
            std::min(result.min_lnLambda_alphae_spitzer, lnLambda);
        result.max_lnLambda_alphae_spitzer =
            std::max(result.max_lnLambda_alphae_spitzer, lnLambda);
        result.min_dt_reactivity_cm3_s =
            std::min(result.min_dt_reactivity_cm3_s, dt_reactivity);
        result.max_dt_reactivity_cm3_s =
            std::max(result.max_dt_reactivity_cm3_s, dt_reactivity);
        ++result.cell_count;
      }
    }
  }

  result.success = true;
  result.double_kB_guard_passed = true;

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p4.alpha.coefficient_provider"
         << "; phase_id=P4"
         << "; stage_id=A"
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; composition_model=" << AlphaCompositionModelName(options.composition_model)
         << "; equimolar_dt_assumption=true"
         << "; tau_alphae_model=" << AlphaTauModelName(options.tau_model)
         << "; reactivity_model=" << AlphaReactivityModelName(options.reactivity_model)
         << "; reactivity_model_requested=" << result.reactivity_model_requested
         << "; reactivity_model_executed=" << result.reactivity_model_executed
         << "; bosch_hale_reference_locked_locally="
         << (result.bosch_hale_reference_locked_locally ? "true" : "false")
         << "; bosch_hale_coefficients_source="
         << result.bosch_hale_coefficients_source
         << "; bosch_hale_temperature_source="
         << result.bosch_hale_temperature_source
         << "; dt_reactivity_internal_unit=cm3_s"
         << "; min_dt_reactivity_cm3_s=" << result.min_dt_reactivity_cm3_s
         << "; max_dt_reactivity_cm3_s=" << result.max_dt_reactivity_cm3_s
         << "; thesis_reactivity_claim_allowed="
         << (result.thesis_reactivity_claim_allowed ? "true" : "false")
         << "; alpha_birth_energy_erg="
         << dec3d::physics::PhysicsConstantsCGS::alpha_birth_energy_erg
         << "; alpha_mass_g=" << dec3d::physics::PhysicsConstantsCGS::alpha_mass_g
         << "; alpha_charge_number="
         << dec3d::physics::PhysicsConstantsCGS::alpha_charge_number
         << "; alpha_transport_model=atzeni_one_group"
         << "; v_alpha0_source=nonrelativistic_birth_energy"
         << "; v_alpha0_cm_per_s=" << v_alpha0
         << "; lambda_drag_model=v_alpha0_tau_alphae_over_9"
         << "; D_alpha_model=atzeni_drag_one_group"
         << "; min_tau_alphae_s=" << result.min_tau_alphae_s
         << "; max_tau_alphae_s=" << result.max_tau_alphae_s
         << "; min_D_alpha_cm2_s=" << result.min_D_alpha_cm2_s
         << "; max_D_alpha_cm2_s=" << result.max_D_alpha_cm2_s
         << "; min_birth_source_erg_cm3_s="
         << result.min_birth_source_erg_cm3_s
         << "; max_birth_source_erg_cm3_s="
         << result.max_birth_source_erg_cm3_s
         << "; min_lnLambda_alphae_spitzer="
         << result.min_lnLambda_alphae_spitzer
         << "; max_lnLambda_alphae_spitzer="
         << result.max_lnLambda_alphae_spitzer
         << "; thermodynamic_recovery_report_present=true"
         << "; canonical_state_mutated=false"
         << "; fallback_used=false"
         << "; double_kB_guard_passed=true";
  result.report_line = report.str();
  return result;
}

bool ValidateAlphaCoefficientDiagnostics(const AlphaCoefficientProviderResult& result) noexcept {
  return result.is_complete() &&
         Contains(result.report_line, "unit_system=cgs") &&
         Contains(result.report_line, "temperature_internal_unit=erg_per_particle") &&
         Contains(result.report_line, "composition_model=equimolar_dt_from_p2_recovery") &&
         Contains(result.report_line, "tau_alphae_model=thesis_spitzer_eq_5_262") &&
         Contains(result.report_line, "reactivity_model_requested=") &&
         Contains(result.report_line, "reactivity_model_executed=") &&
         Contains(result.report_line, "dt_reactivity_internal_unit=cm3_s") &&
         Contains(result.report_line, "min_dt_reactivity_cm3_s=") &&
         Contains(result.report_line, "max_dt_reactivity_cm3_s=") &&
         Contains(result.report_line, "thesis_reactivity_claim_allowed=") &&
         Contains(result.report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(result.report_line, "v_alpha0_source=nonrelativistic_birth_energy") &&
         Contains(result.report_line, "D_alpha_model=atzeni_drag_one_group") &&
         Contains(result.report_line, "thermodynamic_recovery_report_present=true") &&
         Contains(result.report_line, "double_kB_guard_passed=true");
}

}  // namespace dec3d::alpha
