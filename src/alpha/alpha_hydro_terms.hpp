#pragma once

#include "core/array/array3d.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <string>

namespace dec3d::alpha {

struct AlphaHydroTermsOptions {
  double alpha_energy_floor_erg_cm3{0.0};
};

struct AlphaHydroTermBundle {
  bool success{false};
  dec3d::core::Array3D<double> chi_alpha;
  std::size_t cell_count{0u};
  std::size_t invalid_alpha_cell_count{0u};
  double min_epsilon_alpha_before{0.0};
  double max_epsilon_alpha_before{0.0};
  double min_chi_alpha_before{0.0};
  double max_chi_alpha_before{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] double AlphaHydroPressureScalarFromEnergyDensity(
    double epsilon_alpha_erg_cm3) noexcept;

[[nodiscard]] double AlphaHydroEnergyDensityFromPressureScalar(double chi_alpha) noexcept;

[[nodiscard]] double AlphaEpsilonCompressionRatioFromDensityRatio(double rho_ratio) noexcept;

[[nodiscard]] double AlphaEpsilonAfterPressureScalarCompression(
    double epsilon_old,
    double rho_old,
    double rho_new) noexcept;

[[nodiscard]] AlphaHydroTermBundle BuildAlphaHydroTermBundle(
    const dec3d::state::CanonicalState& state,
    AlphaHydroTermsOptions options = {}) noexcept;

[[nodiscard]] bool ValidateAlphaHydroTermsDiagnostics(
    const std::string& report_line) noexcept;

}  // namespace dec3d::alpha
