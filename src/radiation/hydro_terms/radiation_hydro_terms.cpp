#include "radiation/hydro_terms/radiation_hydro_terms.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::radiation {
namespace {

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

void BuildSuccessReport(RadiationHydroTermBundle& result) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.radiation.hydro_terms"
      << "; stage_id=H"
      << "; radiation_hydro_terms_enabled=true"
      << "; radiation_advection=enabled"
      << "; radiation_pressure_work=enabled"
      << "; passive_Ug_advection=false"
      << "; advected_radiation_scalar=P_g_power_3_over_4"
      << "; radiation_pressure_closure=P_g_equals_U_g_over_3"
      << "; Ug_recovery=three_times_chi_rad_power_four_over_three"
      << "; radiation_group_count=" << result.radiation_group_count
      << "; radiation_energy_unit=erg_per_cm3"
      << "; min_Ug_before=" << result.min_Ug_before
      << "; max_Ug_before=" << result.max_Ug_before
      << "; min_chi_rad_before=" << result.min_chi_rad_before
      << "; max_chi_rad_before=" << result.max_chi_rad_before
      << "; invalid_radiation_cell_count=" << result.invalid_radiation_cell_count
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  result.report_line = out.str();
}

RadiationHydroTermBundle Fail(
    RadiationHydroTermBundle result,
    const char* reason) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.hydro_terms.failure"
      << "; stage_id=H"
      << "; radiation_hydro_terms_enabled=true"
      << "; failure_reason=" << reason
      << "; invalid_radiation_cell_count=" << result.invalid_radiation_cell_count
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  result.failure_diagnostics = out.str();
  return result;
}

}  // namespace

double RadiationPressureScalarFromUg(double ug_erg_per_cm3) noexcept {
  return (!Finite(ug_erg_per_cm3) || ug_erg_per_cm3 < 0.0)
             ? std::numeric_limits<double>::quiet_NaN()
             : std::pow(ug_erg_per_cm3 / 3.0, 3.0 / 4.0);
}

double UgFromRadiationPressureScalar(double chi_rad) noexcept {
  return (!Finite(chi_rad) || chi_rad < 0.0)
             ? std::numeric_limits<double>::quiet_NaN()
             : 3.0 * std::pow(chi_rad, 4.0 / 3.0);
}

double RadiationUgCompressionRatioFromDensityRatio(double rho_ratio) noexcept {
  return (!Finite(rho_ratio) || rho_ratio <= 0.0)
             ? std::numeric_limits<double>::quiet_NaN()
             : std::pow(rho_ratio, 4.0 / 3.0);
}

double RadiationUgAfterPressureScalarCompression(
    double ug_old,
    double rho_old,
    double rho_new) noexcept {
  if (!Finite(ug_old) || !Finite(rho_old) || !Finite(rho_new) ||
      ug_old < 0.0 || rho_old <= 0.0 || rho_new <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return ug_old * RadiationUgCompressionRatioFromDensityRatio(rho_new / rho_old);
}

RadiationHydroTermBundle BuildRadiationHydroTermBundle(
    const dec3d::state::CanonicalState& state,
    const RadiationHydroTermsOptions& options) noexcept {
  RadiationHydroTermBundle result;
  result.radiation_group_count = state.radiation_groups.size();

  if (result.radiation_group_count == 0u ||
      state.layout.radiation_group_count != result.radiation_group_count) {
    return Fail(result, "radiation group layout is missing or inconsistent");
  }
  if (!Finite(options.radiation_energy_floor_erg_per_cm3) ||
      options.radiation_energy_floor_erg_per_cm3 < 0.0) {
    return Fail(result, "radiation energy floor must be finite and non-negative");
  }

  const auto radial = state.layout.radial_cells;
  const auto theta_cells = state.layout.theta_cells;
  const auto phi_cells = state.layout.phi_cells;
  result.chi_rad_by_group.reserve(result.radiation_group_count);
  result.min_Ug_before = std::numeric_limits<double>::infinity();
  result.max_Ug_before = -std::numeric_limits<double>::infinity();
  result.min_chi_rad_before = std::numeric_limits<double>::infinity();
  result.max_chi_rad_before = -std::numeric_limits<double>::infinity();

  for (std::size_t group = 0u; group < result.radiation_group_count; ++group) {
    const auto& ug = state.radiation_groups[group];
    if (ug.extent_r() != radial || ug.extent_theta() != theta_cells ||
        ug.extent_phi() != phi_cells) {
      return Fail(result, "radiation group storage shape is inconsistent");
    }
    dec3d::core::Array3D<double> chi(radial, theta_cells, phi_cells, 0.0);
    for (std::size_t r = 0u; r < radial; ++r) {
      for (std::size_t t = 0u; t < theta_cells; ++t) {
        for (std::size_t p = 0u; p < phi_cells; ++p) {
          const double value = ug(r, t, p);
          if (!Finite(value) || value < options.radiation_energy_floor_erg_per_cm3) {
            ++result.invalid_radiation_cell_count;
            continue;
          }
          const double scalar = RadiationPressureScalarFromUg(value);
          if (!Finite(scalar)) {
            ++result.invalid_radiation_cell_count;
            continue;
          }
          chi(r, t, p) = scalar;
          result.min_Ug_before = std::min(result.min_Ug_before, value);
          result.max_Ug_before = std::max(result.max_Ug_before, value);
          result.min_chi_rad_before = std::min(result.min_chi_rad_before, scalar);
          result.max_chi_rad_before = std::max(result.max_chi_rad_before, scalar);
        }
      }
    }
    result.chi_rad_by_group.push_back(std::move(chi));
  }

  if (result.invalid_radiation_cell_count != 0u) {
    return Fail(result, "radiation group contains non-finite or below-floor values");
  }

  result.success = true;
  BuildSuccessReport(result);
  return result;
}

bool ValidateRadiationHydroTermsDiagnostics(const std::string& report_line) noexcept {
  return Contains(report_line, "diagnostic_id=p3.radiation.hydro_terms") &&
         Contains(report_line, "stage_id=H") &&
         Contains(report_line, "radiation_hydro_terms_enabled=true") &&
         Contains(report_line, "radiation_advection=enabled") &&
         Contains(report_line, "radiation_pressure_work=enabled") &&
         Contains(report_line, "passive_Ug_advection=false") &&
         Contains(report_line, "advected_radiation_scalar=P_g_power_3_over_4") &&
         Contains(report_line, "radiation_pressure_closure=P_g_equals_U_g_over_3") &&
         Contains(report_line, "Ug_recovery=three_times_chi_rad_power_four_over_three") &&
         Contains(report_line, "radiation_group_count=") &&
         Contains(report_line, "invalid_radiation_cell_count=") &&
         Contains(report_line, "fallback_used=false");
}

}  // namespace dec3d::radiation
