#pragma once

#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <iosfwd>
#include <string>
#include <vector>

namespace dec3d::radiation {

enum class RadiationGroupMode {
  missing,
  gray_full_spectrum,
  explicit_frequency_groups
};

struct RadiationGroupLayout {
  RadiationGroupMode mode{RadiationGroupMode::missing};
  std::size_t group_count{0u};
  std::vector<double> frequency_edges_hz;
};

struct RadiationGroupLayoutValidationResult {
  bool success{false};
  RadiationGroupMode mode{RadiationGroupMode::missing};
  std::size_t group_count{0u};
  bool group_count_explicit{false};
  bool gray_full_spectrum{false};
  bool frequency_edges_used{false};
  bool frequency_edges_strictly_increasing{false};
  bool canonical_group_arrays_present{false};
  std::string report_line;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] RadiationGroupLayout MakeGrayFullSpectrumRadiationGroupLayout();

[[nodiscard]] RadiationGroupLayout MakeExplicitFrequencyRadiationGroupLayout(
    const std::vector<double>& frequency_edges_hz);

[[nodiscard]] const char* RadiationGroupModeName(RadiationGroupMode mode) noexcept;

std::ostream& operator<<(std::ostream& os, RadiationGroupMode mode);

[[nodiscard]] RadiationGroupLayoutValidationResult ValidateRadiationGroupLayout(
    const RadiationGroupLayout& layout) noexcept;

[[nodiscard]] RadiationGroupLayoutValidationResult ValidateRadiationGroupStateStorage(
    const RadiationGroupLayout& layout,
    const dec3d::state::CanonicalState& state) noexcept;

[[nodiscard]] bool ValidateRadiationGroupLayoutDiagnostics(
    const RadiationGroupLayoutValidationResult& result) noexcept;

}  // namespace dec3d::radiation
