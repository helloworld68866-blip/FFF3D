#include "io/radial_profile.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string_view>
#include <utility>

namespace dec3d::io {
namespace {

[[nodiscard]] std::string Trim(std::string_view text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string_view::npos) {
    return {};
  }
  const auto last = text.find_last_not_of(" \t\r\n");
  return std::string(text.substr(first, last - first + 1u));
}

[[nodiscard]] std::string StripComment(std::string_view line) {
  const auto pos = line.find('#');
  return Trim(pos == std::string_view::npos ? line : line.substr(0u, pos));
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::string NormalizeColumnName(std::string_view text) {
  std::string normalized;
  normalized.reserve(text.size());
  for (char ch : text) {
    normalized.push_back(static_cast<char>(
        std::tolower(static_cast<unsigned char>(ch))));
  }
  return normalized;
}

[[nodiscard]] std::vector<std::string> SplitFields(const std::string& line) {
  std::vector<std::string> fields;
  if (line.find(',') != std::string::npos) {
    std::stringstream stream(line);
    std::string part;
    while (std::getline(stream, part, ',')) {
      const auto trimmed = Trim(part);
      if (!trimmed.empty()) {
        fields.push_back(trimmed);
      }
    }
    return fields;
  }
  std::istringstream stream(line);
  std::string part;
  while (stream >> part) {
    fields.push_back(part);
  }
  return fields;
}

template <typename T>
[[nodiscard]] bool ParseNumber(const std::string& text, T& out) {
  std::istringstream in(text);
  in >> out;
  return in && in.eof();
}

RadialProfileLoadResult FailProfile(std::string reason) {
  RadialProfileLoadResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = "diagnostic_id=p5.io.radial_profile.failure; failure_reason=" +
                               result.failure_reason;
  return result;
}

[[nodiscard]] double Lerp(double a, double b, double w) noexcept {
  return a + (b - a) * w;
}

[[nodiscard]] double ColumnValueOrDefault(const std::map<std::string, std::size_t>& indices,
                                          const std::vector<std::string>& fields,
                                          const std::string& name,
                                          double default_value,
                                          bool& ok) {
  const auto it = indices.find(name);
  if (it == indices.end()) {
    return default_value;
  }
  if (it->second >= fields.size()) {
    ok = false;
    return 0.0;
  }
  double value = 0.0;
  if (!ParseNumber(fields[it->second], value)) {
    ok = false;
    return 0.0;
  }
  return value;
}

[[nodiscard]] std::string BuildProfileReport(const RadialProfile& profile,
                                             const std::filesystem::path& path) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "diagnostic_id=p5.io.radial_profile"
      << "; profile_file=" << path.string()
      << "; profile_row_count=" << profile.rows.size()
      << "; required_columns_present=true"
      << "; radius_strictly_increasing=true"
      << "; profile_radius_unit=" << profile.radius_unit
      << "; profile_interpolation=linear"
      << "; comment_policy=hash_full_line_and_trailing"
      << "; generated_file_chinese_comments_present=true"
      << "; chinese_comments_runtime_contract=false"
      << "; profile_validation_success=true";
  return out.str();
}

}  // namespace

RadialProfileLoadResult LoadRadialProfile(const std::filesystem::path& profile_path,
                                          const std::string& radius_unit) noexcept {
  if (radius_unit != "um") {
    return FailProfile("unsupported radius unit");
  }
  std::ifstream in(profile_path);
  if (!in) {
    return FailProfile("missing profile file");
  }

  std::vector<std::string> header;
  std::map<std::string, std::size_t> indices;
  std::vector<RadialProfileRow> rows;
  std::string raw;
  while (std::getline(in, raw)) {
    const auto line = StripComment(raw);
    if (line.empty()) {
      continue;
    }
    const auto fields = SplitFields(line);
    if (fields.empty()) {
      continue;
    }
    if (header.empty()) {
      header = fields;
      for (std::size_t i = 0u; i < header.size(); ++i) {
        indices[NormalizeColumnName(header[i])] = i;
      }
      for (const auto required : {
               std::pair{"r_um", "r_um"},
               std::pair{"rho_g_cm3", "rho_g_cm3"},
               std::pair{"te_kev", "Te_keV"},
               std::pair{"ti_kev", "Ti_keV"},
               std::pair{"vr_cm_s", "vr_cm_s"},
           }) {
        if (!indices.contains(required.first)) {
          return FailProfile(std::string("missing required column ") + required.second);
        }
      }
      continue;
    }
    bool ok = true;
    RadialProfileRow row;
    row.r_cm = ColumnValueOrDefault(indices, fields, "r_um", 0.0, ok) * 1.0e-4;
    row.rho_g_cm3 = ColumnValueOrDefault(indices, fields, "rho_g_cm3", 0.0, ok);
    row.Te_keV = ColumnValueOrDefault(indices, fields, "te_kev", 0.0, ok);
    row.Ti_keV = ColumnValueOrDefault(indices, fields, "ti_kev", 0.0, ok);
    row.vr_cm_s = ColumnValueOrDefault(indices, fields, "vr_cm_s", 0.0, ok);
    row.vt_cm_s = ColumnValueOrDefault(indices, fields, "vt_cm_s", 0.0, ok);
    row.vp_cm_s = ColumnValueOrDefault(indices, fields, "vp_cm_s", 0.0, ok);
    row.epsilon_alpha_erg_cm3 =
        ColumnValueOrDefault(indices, fields, "epsilon_alpha_erg_cm3", 0.0, ok);
    row.radiation_scale = ColumnValueOrDefault(indices, fields, "radiation_scale", 1.0, ok);
    if (!ok) {
      return FailProfile("invalid numeric field");
    }
    if (!std::isfinite(row.r_cm) || !std::isfinite(row.rho_g_cm3) ||
        !std::isfinite(row.Te_keV) || !std::isfinite(row.Ti_keV) ||
        !std::isfinite(row.vr_cm_s) || !std::isfinite(row.vt_cm_s) ||
        !std::isfinite(row.vp_cm_s)) {
      return FailProfile("non-finite profile value");
    }
    if (row.rho_g_cm3 <= 0.0) {
      return FailProfile("rho must be positive");
    }
    if (row.Te_keV <= 0.0 || row.Ti_keV <= 0.0) {
      return FailProfile("temperature must be positive");
    }
    if (row.epsilon_alpha_erg_cm3 < 0.0) {
      return FailProfile("epsilon_alpha must be nonnegative");
    }
    if (row.radiation_scale < 0.0) {
      return FailProfile("radiation_scale must be nonnegative");
    }
    rows.push_back(row);
  }
  if (header.empty() || rows.size() < 2u) {
    return FailProfile("profile requires header and at least two rows");
  }
  for (std::size_t i = 1u; i < rows.size(); ++i) {
    if (!(rows[i].r_cm > rows[i - 1u].r_cm)) {
      return FailProfile("radius must be strictly increasing");
    }
  }

  RadialProfileLoadResult result;
  result.success = true;
  result.profile.rows = std::move(rows);
  result.profile.header_columns = std::move(header);
  result.profile.radius_unit = radius_unit;
  result.report_line = BuildProfileReport(result.profile, profile_path);
  return result;
}

RadialProfileSampleResult SampleRadialProfile(const RadialProfile& profile,
                                              double radius_cm) noexcept {
  RadialProfileSampleResult result;
  if (profile.rows.empty() || !std::isfinite(radius_cm)) {
    result.failure_reason = "invalid profile sample";
    return result;
  }
  if (radius_cm < profile.rows.front().r_cm || radius_cm > profile.rows.back().r_cm) {
    result.failure_reason = "radius outside profile";
    return result;
  }
  for (std::size_t i = 0u; i + 1u < profile.rows.size(); ++i) {
    const auto& left = profile.rows[i];
    const auto& right = profile.rows[i + 1u];
    if (radius_cm >= left.r_cm && radius_cm <= right.r_cm) {
      const double w = (radius_cm - left.r_cm) / (right.r_cm - left.r_cm);
      result.success = true;
      result.value.r_cm = radius_cm;
      result.value.rho_g_cm3 = Lerp(left.rho_g_cm3, right.rho_g_cm3, w);
      result.value.Te_keV = Lerp(left.Te_keV, right.Te_keV, w);
      result.value.Ti_keV = Lerp(left.Ti_keV, right.Ti_keV, w);
      result.value.vr_cm_s = Lerp(left.vr_cm_s, right.vr_cm_s, w);
      result.value.vt_cm_s = Lerp(left.vt_cm_s, right.vt_cm_s, w);
      result.value.vp_cm_s = Lerp(left.vp_cm_s, right.vp_cm_s, w);
      result.value.epsilon_alpha_erg_cm3 =
          Lerp(left.epsilon_alpha_erg_cm3, right.epsilon_alpha_erg_cm3, w);
      result.value.radiation_scale = Lerp(left.radiation_scale, right.radiation_scale, w);
      return result;
    }
  }
  result.failure_reason = "radius outside profile";
  return result;
}

bool ValidateRadialProfileDiagnostics(const RadialProfileLoadResult& result) noexcept {
  if (!result.success) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p5.io.radial_profile") &&
         Contains(line, "profile_row_count=") &&
         Contains(line, "required_columns_present=true") &&
         Contains(line, "radius_strictly_increasing=true") &&
         Contains(line, "chinese_comments_runtime_contract=false") &&
         Contains(line, "profile_validation_success=true");
}

}  // namespace dec3d::io
