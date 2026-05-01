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

AlphaStateContractResult Fail(const std::string& reason) {
  AlphaStateContractResult result;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "diagnostic_id=p4.alpha.contract.failure"
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false"
         << "; fallback_used=false";
  result.failure_diagnostics = report.str();
  return result;
}

}  // namespace

bool AlphaStateContractResult::is_complete() const noexcept {
  return success &&
         cell_count > 0u &&
         validated_fields == dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state) &&
         updated_fields == 0u &&
         !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p4.alpha.contract") &&
         Contains(report_line, "phase_id=P4") &&
         Contains(report_line, "stage_id=A") &&
         Contains(report_line, "validated_fields=alpha_state") &&
         Contains(report_line, "updated_fields=none") &&
         Contains(report_line, "canonical_state_mutated=false");
}

double AlphaPressureFromEnergyDensity(double epsilon_alpha_erg_cm3) noexcept {
  return (2.0 / 3.0) * epsilon_alpha_erg_cm3;
}

double AlphaEnergyDensityFromPressure(double p_alpha_erg_cm3) noexcept {
  return 1.5 * p_alpha_erg_cm3;
}

double AlphaPressureScalarFromEnergyDensity(double epsilon_alpha_erg_cm3) noexcept {
  const double p_alpha = AlphaPressureFromEnergyDensity(epsilon_alpha_erg_cm3);
  if (!(Finite(p_alpha) && p_alpha >= 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return std::pow(p_alpha, 3.0 / 5.0);
}

double AlphaEnergyDensityFromPressureScalar(double chi_alpha) noexcept {
  if (!(Finite(chi_alpha) && chi_alpha >= 0.0)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return AlphaEnergyDensityFromPressure(std::pow(chi_alpha, 5.0 / 3.0));
}

AlphaStateContractResult ValidateAlphaStateContract(
    const dec3d::state::CanonicalState& state,
    AlphaStateContractOptions options) noexcept {
  if (!state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state)) {
    return Fail("canonical alpha_state storage is missing");
  }

  if (!(Finite(options.alpha_energy_floor_erg_cm3) &&
        options.alpha_energy_floor_erg_cm3 >= 0.0)) {
    return Fail("alpha energy floor must be finite and nonnegative");
  }

  if (state.alpha_state.storage.extent_r() != state.layout.radial_cells ||
      state.alpha_state.storage.extent_theta() != state.layout.theta_cells ||
      state.alpha_state.storage.extent_phi() != state.layout.phi_cells) {
    return Fail("canonical alpha_state shape mismatch");
  }

  AlphaStateContractResult result;
  result.success = true;
  result.cell_count = state.alpha_state.storage.size();
  result.validated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
  result.updated_fields = 0u;
  result.min_alpha_energy_density_erg_cm3 = std::numeric_limits<double>::infinity();
  result.max_alpha_energy_density_erg_cm3 = -std::numeric_limits<double>::infinity();

  for (const double value : state.alpha_state.storage.storage()) {
    if (!(Finite(value) && value >= options.alpha_energy_floor_erg_cm3)) {
      return Fail("alpha energy density must be finite and above floor");
    }
    result.min_alpha_energy_density_erg_cm3 =
        std::min(result.min_alpha_energy_density_erg_cm3, value);
    result.max_alpha_energy_density_erg_cm3 =
        std::max(result.max_alpha_energy_density_erg_cm3, value);
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p4.alpha.contract"
         << "; phase_id=P4"
         << "; stage_id=A"
         << "; alpha_state_authoritative=alpha_energy_density"
         << "; alpha_energy_unit=erg_per_cm3"
         << "; alpha_transport_model=atzeni_one_group"
         << "; alpha_pressure_relation=P_alpha_eq_2_over_3_epsilon_alpha"
         << "; alpha_pressure_scalar=chi_alpha_eq_P_alpha_pow_3_over_5"
         << "; validated_fields=alpha_state"
         << "; updated_fields=none"
         << "; cell_count=" << result.cell_count
         << "; min_alpha_energy_density_erg_cm3="
         << result.min_alpha_energy_density_erg_cm3
         << "; max_alpha_energy_density_erg_cm3="
         << result.max_alpha_energy_density_erg_cm3
         << "; canonical_state_mutated=false"
         << "; fallback_used=false";
  result.report_line = report.str();
  return result;
}

bool ValidateAlphaContractDiagnostics(const AlphaStateContractResult& result) noexcept {
  return result.is_complete() &&
         Contains(result.report_line, "alpha_state_authoritative=alpha_energy_density") &&
         Contains(result.report_line, "alpha_energy_unit=erg_per_cm3") &&
         Contains(result.report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(result.report_line, "alpha_pressure_relation=P_alpha_eq_2_over_3_epsilon_alpha") &&
         Contains(result.report_line, "alpha_pressure_scalar=chi_alpha_eq_P_alpha_pow_3_over_5") &&
         Contains(result.report_line, "fallback_used=false");
}

}  // namespace dec3d::alpha
