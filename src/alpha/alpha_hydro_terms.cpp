#include "alpha/alpha_hydro_terms.hpp"

#include "alpha/alpha_contract.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::alpha {
namespace {

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

AlphaHydroTermBundle Fail(AlphaHydroTermBundle result, const char* reason) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p4.alpha.hydro_terms.failure"
      << "; phase_id=P4"
      << "; stage_id=H"
      << "; alpha_hydro_bundle_built=true"
      << "; failure_reason=" << reason
      << "; invalid_alpha_cell_count=" << result.invalid_alpha_cell_count
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  result.failure_diagnostics = out.str();
  return result;
}

void BuildSuccessReport(AlphaHydroTermBundle& result) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p4.alpha.hydro_terms"
      << "; phase_id=P4"
      << "; stage_id=H"
      << "; alpha_hydro_bundle_built=true"
      << "; alpha_hydro_side_terms_model=pressure_scalar_advection_contract"
      << "; passive_epsilon_alpha_advection=false"
      << "; advected_alpha_scalar=P_alpha_power_3_over_5"
      << "; alpha_pressure_relation=P_alpha_eq_2_over_3_epsilon_alpha"
      << "; epsilon_alpha_recovery=three_over_two_times_chi_alpha_power_five_over_three"
      << "; alpha_energy_unit=erg_per_cm3"
      << "; bundle_scope=pre_hydro_view_derivation"
      << "; cell_count=" << result.cell_count
      << "; min_epsilon_alpha_before=" << result.min_epsilon_alpha_before
      << "; max_epsilon_alpha_before=" << result.max_epsilon_alpha_before
      << "; min_chi_alpha_before=" << result.min_chi_alpha_before
      << "; max_chi_alpha_before=" << result.max_chi_alpha_before
      << "; invalid_alpha_cell_count=" << result.invalid_alpha_cell_count
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  result.report_line = out.str();
}

}  // namespace

double AlphaHydroPressureScalarFromEnergyDensity(
    double epsilon_alpha_erg_cm3) noexcept {
  return AlphaPressureScalarFromEnergyDensity(epsilon_alpha_erg_cm3);
}

double AlphaHydroEnergyDensityFromPressureScalar(double chi_alpha) noexcept {
  return AlphaEnergyDensityFromPressureScalar(chi_alpha);
}

double AlphaEpsilonCompressionRatioFromDensityRatio(double rho_ratio) noexcept {
  return (!Finite(rho_ratio) || rho_ratio <= 0.0)
             ? std::numeric_limits<double>::quiet_NaN()
             : std::pow(rho_ratio, 5.0 / 3.0);
}

double AlphaEpsilonAfterPressureScalarCompression(
    double epsilon_old,
    double rho_old,
    double rho_new) noexcept {
  if (!Finite(epsilon_old) || !Finite(rho_old) || !Finite(rho_new) ||
      epsilon_old < 0.0 || rho_old <= 0.0 || rho_new <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return epsilon_old * AlphaEpsilonCompressionRatioFromDensityRatio(rho_new / rho_old);
}

AlphaHydroTermBundle BuildAlphaHydroTermBundle(
    const dec3d::state::CanonicalState& state,
    AlphaHydroTermsOptions options) noexcept {
  AlphaHydroTermBundle result;
  if (!state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state)) {
    return Fail(result, "canonical alpha_state storage is missing");
  }
  if (!Finite(options.alpha_energy_floor_erg_cm3) ||
      options.alpha_energy_floor_erg_cm3 < 0.0) {
    return Fail(result, "alpha energy floor must be finite and non-negative");
  }

  const auto radial = state.layout.radial_cells;
  const auto theta = state.layout.theta_cells;
  const auto phi = state.layout.phi_cells;
  if (state.alpha_state.storage.extent_r() != radial ||
      state.alpha_state.storage.extent_theta() != theta ||
      state.alpha_state.storage.extent_phi() != phi) {
    return Fail(result, "alpha_state storage shape is inconsistent");
  }

  result.cell_count = state.alpha_state.storage.size();
  result.chi_alpha = dec3d::core::Array3D<double>(radial, theta, phi, 0.0);
  result.min_epsilon_alpha_before = std::numeric_limits<double>::infinity();
  result.max_epsilon_alpha_before = -std::numeric_limits<double>::infinity();
  result.min_chi_alpha_before = std::numeric_limits<double>::infinity();
  result.max_chi_alpha_before = -std::numeric_limits<double>::infinity();

  for (std::size_t r = 0; r < radial; ++r) {
    for (std::size_t t = 0; t < theta; ++t) {
      for (std::size_t p = 0; p < phi; ++p) {
        const double epsilon = state.alpha_state.storage(r, t, p);
        if (!Finite(epsilon) || epsilon < options.alpha_energy_floor_erg_cm3) {
          ++result.invalid_alpha_cell_count;
          continue;
        }
        const double chi = AlphaHydroPressureScalarFromEnergyDensity(epsilon);
        if (!Finite(chi)) {
          ++result.invalid_alpha_cell_count;
          continue;
        }
        result.chi_alpha(r, t, p) = chi;
        result.min_epsilon_alpha_before = std::min(result.min_epsilon_alpha_before, epsilon);
        result.max_epsilon_alpha_before = std::max(result.max_epsilon_alpha_before, epsilon);
        result.min_chi_alpha_before = std::min(result.min_chi_alpha_before, chi);
        result.max_chi_alpha_before = std::max(result.max_chi_alpha_before, chi);
      }
    }
  }

  if (result.invalid_alpha_cell_count != 0u) {
    return Fail(result, "alpha_state contains non-finite or below-floor values");
  }

  result.success = true;
  BuildSuccessReport(result);
  return result;
}

bool ValidateAlphaHydroTermsDiagnostics(const std::string& report_line) noexcept {
  return Contains(report_line, "diagnostic_id=p4.alpha.hydro_terms") &&
         Contains(report_line, "phase_id=P4") &&
         Contains(report_line, "stage_id=H") &&
         Contains(report_line, "alpha_hydro_bundle_built=true") &&
         Contains(report_line, "alpha_hydro_side_terms_model=pressure_scalar_advection_contract") &&
         Contains(report_line, "passive_epsilon_alpha_advection=false") &&
         Contains(report_line, "advected_alpha_scalar=P_alpha_power_3_over_5") &&
         Contains(report_line, "alpha_pressure_relation=P_alpha_eq_2_over_3_epsilon_alpha") &&
         Contains(report_line, "epsilon_alpha_recovery=three_over_two_times_chi_alpha_power_five_over_three") &&
         Contains(report_line, "bundle_scope=pre_hydro_view_derivation") &&
         Contains(report_line, "invalid_alpha_cell_count=") &&
         Contains(report_line, "fallback_used=false");
}

}  // namespace dec3d::alpha
