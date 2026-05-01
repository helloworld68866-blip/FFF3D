#include "radiation/groups/radiation_group_layout.hpp"

#include <cmath>
#include <ostream>
#include <sstream>

namespace dec3d::radiation {
namespace {

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] std::string BuildReport(const RadiationGroupLayoutValidationResult& result) {
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.group_layout"
      << "; stage_id=R"
      << "; radiation_group_mode=" << RadiationGroupModeName(result.mode)
      << "; radiation_group_count=" << result.group_count
      << "; group_count_explicit=" << BoolToken(result.group_count_explicit)
      << "; gray_full_spectrum=" << BoolToken(result.gray_full_spectrum)
      << "; frequency_edges_used=" << BoolToken(result.frequency_edges_used)
      << "; frequency_edges_strictly_increasing=" << BoolToken(result.frequency_edges_strictly_increasing)
      << "; canonical_group_arrays_present=" << BoolToken(result.canonical_group_arrays_present)
      << "; frequency_edges_unit=Hz"
      << "; authoritative_state=radiation_groups";
  return out.str();
}

[[nodiscard]] RadiationGroupLayoutValidationResult Fail(
    RadiationGroupMode mode,
    std::size_t group_count,
    const char* reason) {
  RadiationGroupLayoutValidationResult result;
  result.mode = mode;
  result.group_count = group_count;
  result.failure_reason = reason;
  return result;
}

}  // namespace

bool RadiationGroupLayoutValidationResult::is_complete() const noexcept {
  return success && group_count > 0u && !report_line.empty();
}

RadiationGroupLayout MakeGrayFullSpectrumRadiationGroupLayout() {
  RadiationGroupLayout layout;
  layout.mode = RadiationGroupMode::gray_full_spectrum;
  layout.group_count = 1u;
  return layout;
}

RadiationGroupLayout MakeExplicitFrequencyRadiationGroupLayout(
    const std::vector<double>& frequency_edges_hz) {
  RadiationGroupLayout layout;
  layout.mode = RadiationGroupMode::explicit_frequency_groups;
  layout.group_count = frequency_edges_hz.size() > 0u ? frequency_edges_hz.size() - 1u : 0u;
  layout.frequency_edges_hz = frequency_edges_hz;
  return layout;
}

const char* RadiationGroupModeName(RadiationGroupMode mode) noexcept {
  switch (mode) {
    case RadiationGroupMode::gray_full_spectrum:
      return "gray_full_spectrum";
    case RadiationGroupMode::explicit_frequency_groups:
      return "explicit_frequency_groups";
    case RadiationGroupMode::missing:
    default:
      return "missing";
  }
}

std::ostream& operator<<(std::ostream& os, RadiationGroupMode mode) {
  os << RadiationGroupModeName(mode);
  return os;
}

RadiationGroupLayoutValidationResult ValidateRadiationGroupLayout(
    const RadiationGroupLayout& layout) noexcept {
  if (layout.mode == RadiationGroupMode::missing) {
    return Fail(layout.mode, layout.group_count, "radiation group mode is missing");
  }

  RadiationGroupLayoutValidationResult result;
  result.success = true;
  result.mode = layout.mode;
  result.group_count = layout.group_count;
  result.group_count_explicit = layout.group_count > 0u;

  if (layout.mode == RadiationGroupMode::gray_full_spectrum) {
    if (layout.group_count != 1u) {
      return Fail(layout.mode, layout.group_count, "gray full spectrum requires exactly one group");
    }
    if (!layout.frequency_edges_hz.empty()) {
      return Fail(layout.mode, layout.group_count,
                  "gray full spectrum must not provide finite frequency edges");
    }
    result.gray_full_spectrum = true;
    result.frequency_edges_used = false;
    result.frequency_edges_strictly_increasing = false;
    result.report_line = BuildReport(result);
    return result;
  }

  if (layout.mode == RadiationGroupMode::explicit_frequency_groups) {
    if (layout.frequency_edges_hz.size() < 2u || layout.group_count + 1u != layout.frequency_edges_hz.size()) {
      return Fail(layout.mode, layout.group_count,
                  "explicit frequency groups require group_count plus one frequency edges");
    }

    for (const double edge : layout.frequency_edges_hz) {
      if (!std::isfinite(edge) || edge <= 0.0) {
        return Fail(layout.mode, layout.group_count, "frequency edges must be finite and positive");
      }
    }
    for (std::size_t i = 1u; i < layout.frequency_edges_hz.size(); ++i) {
      if (!(layout.frequency_edges_hz[i] > layout.frequency_edges_hz[i - 1u])) {
        return Fail(layout.mode, layout.group_count, "frequency edges must be strictly increasing");
      }
    }

    result.frequency_edges_used = true;
    result.frequency_edges_strictly_increasing = true;
    result.report_line = BuildReport(result);
    return result;
  }

  return Fail(layout.mode, layout.group_count, "radiation group mode is unsupported");
}

RadiationGroupLayoutValidationResult ValidateRadiationGroupStateStorage(
    const RadiationGroupLayout& layout,
    const dec3d::state::CanonicalState& state) noexcept {
  auto result = ValidateRadiationGroupLayout(layout);
  if (!result.success) {
    return result;
  }

  if (state.radiation_groups.size() != result.group_count ||
      state.layout.radiation_group_count != result.group_count) {
    return Fail(layout.mode, result.group_count, "canonical radiation group count mismatch");
  }

  for (const auto& group : state.radiation_groups) {
    if (group.extent_r() != state.layout.radial_cells ||
        group.extent_theta() != state.layout.theta_cells ||
        group.extent_phi() != state.layout.phi_cells) {
      return Fail(layout.mode, result.group_count, "canonical radiation group shape mismatch");
    }
  }

  result.canonical_group_arrays_present = true;
  result.report_line = BuildReport(result);
  return result;
}

bool ValidateRadiationGroupLayoutDiagnostics(
    const RadiationGroupLayoutValidationResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }

  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.group_layout") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "radiation_group_mode=") &&
         Contains(line, "radiation_group_count=") &&
         Contains(line, "group_count_explicit=") &&
         Contains(line, "gray_full_spectrum=") &&
         Contains(line, "frequency_edges_used=") &&
         Contains(line, "frequency_edges_strictly_increasing=") &&
         Contains(line, "canonical_group_arrays_present=") &&
         Contains(line, "frequency_edges_unit=Hz") &&
         Contains(line, "authoritative_state=radiation_groups");
}

}  // namespace dec3d::radiation
