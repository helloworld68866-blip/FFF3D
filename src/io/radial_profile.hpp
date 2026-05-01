#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace dec3d::io {

struct RadialProfileRow {
  double r_cm{0.0};
  double rho_g_cm3{0.0};
  double Te_keV{0.0};
  double Ti_keV{0.0};
  double vr_cm_s{0.0};
  double vt_cm_s{0.0};
  double vp_cm_s{0.0};
  double epsilon_alpha_erg_cm3{0.0};
  double radiation_scale{1.0};
};

struct RadialProfile {
  std::vector<RadialProfileRow> rows;
  std::vector<std::string> header_columns;
  std::string radius_unit;
};

struct RadialProfileLoadResult {
  bool success{false};
  RadialProfile profile;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct RadialProfileSampleResult {
  bool success{false};
  RadialProfileRow value;
  std::string failure_reason;
};

[[nodiscard]] RadialProfileLoadResult LoadRadialProfile(
    const std::filesystem::path& profile_path,
    const std::string& radius_unit) noexcept;

[[nodiscard]] RadialProfileSampleResult SampleRadialProfile(
    const RadialProfile& profile,
    double radius_cm) noexcept;

[[nodiscard]] bool ValidateRadialProfileDiagnostics(
    const RadialProfileLoadResult& result) noexcept;

}  // namespace dec3d::io
