#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <string>

namespace dec3d::alpha {

struct AlphaStateContractOptions {
  double alpha_energy_floor_erg_cm3{0.0};
};

struct AlphaStateContractResult {
  bool success{false};
  std::size_t cell_count{0u};
  dec3d::core::AuthoritativeFieldMask validated_fields{0u};
  dec3d::core::AuthoritativeFieldMask updated_fields{0u};
  double min_alpha_energy_density_erg_cm3{0.0};
  double max_alpha_energy_density_erg_cm3{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] double AlphaPressureFromEnergyDensity(double epsilon_alpha_erg_cm3) noexcept;
[[nodiscard]] double AlphaEnergyDensityFromPressure(double p_alpha_erg_cm3) noexcept;
[[nodiscard]] double AlphaPressureScalarFromEnergyDensity(double epsilon_alpha_erg_cm3) noexcept;
[[nodiscard]] double AlphaEnergyDensityFromPressureScalar(double chi_alpha) noexcept;

[[nodiscard]] AlphaStateContractResult ValidateAlphaStateContract(
    const dec3d::state::CanonicalState& state,
    AlphaStateContractOptions options = {}) noexcept;

[[nodiscard]] bool ValidateAlphaContractDiagnostics(
    const AlphaStateContractResult& result) noexcept;

}  // namespace dec3d::alpha
