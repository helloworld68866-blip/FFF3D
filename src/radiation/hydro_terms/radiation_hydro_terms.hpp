#pragma once

#include "core/array/array3d.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::radiation {

struct RadiationHydroTermsOptions {
  double radiation_energy_floor_erg_per_cm3{0.0};
};

struct RadiationHydroTermBundle {
  bool success{false};
  std::vector<dec3d::core::Array3D<double>> chi_rad_by_group;
  std::size_t radiation_group_count{0u};
  std::size_t invalid_radiation_cell_count{0u};
  double min_Ug_before{0.0};
  double max_Ug_before{0.0};
  double min_chi_rad_before{0.0};
  double max_chi_rad_before{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] double RadiationPressureScalarFromUg(double ug_erg_per_cm3) noexcept;
[[nodiscard]] double UgFromRadiationPressureScalar(double chi_rad) noexcept;
[[nodiscard]] double RadiationUgCompressionRatioFromDensityRatio(double rho_ratio) noexcept;
[[nodiscard]] double RadiationUgAfterPressureScalarCompression(
    double ug_old,
    double rho_old,
    double rho_new) noexcept;

[[nodiscard]] RadiationHydroTermBundle BuildRadiationHydroTermBundle(
    const dec3d::state::CanonicalState& state,
    const RadiationHydroTermsOptions& options) noexcept;

[[nodiscard]] bool ValidateRadiationHydroTermsDiagnostics(
    const std::string& report_line) noexcept;

}  // namespace dec3d::radiation
