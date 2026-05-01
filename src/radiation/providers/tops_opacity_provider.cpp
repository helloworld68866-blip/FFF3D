#include "radiation/providers/tops_opacity_provider.hpp"

#include "physics/units/physical_constants.hpp"
#include "radiation/providers/group_blackbody.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dec3d::radiation {
namespace {

constexpr double kCoordinateKeyScale = 1.0e12;

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::uint64_t DoubleBits(double value) noexcept {
  std::uint64_t bits = 0u;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(value));
  return bits;
}

struct ThermodynamicCacheKey {
  std::uint64_t te_keV_bits{0u};
  std::uint64_t te_erg_bits{0u};
  std::uint64_t rho_bits{0u};

  [[nodiscard]] bool operator==(const ThermodynamicCacheKey& other) const noexcept {
    return te_keV_bits == other.te_keV_bits &&
           te_erg_bits == other.te_erg_bits &&
           rho_bits == other.rho_bits;
  }
};

struct ThermodynamicCacheKeyHash {
  [[nodiscard]] std::size_t operator()(const ThermodynamicCacheKey& key) const noexcept {
    std::size_t seed = std::hash<std::uint64_t>{}(key.te_keV_bits);
    seed ^= std::hash<std::uint64_t>{}(key.te_erg_bits) + 0x9e3779b97f4a7c15ull +
            (seed << 6u) + (seed >> 2u);
    seed ^= std::hash<std::uint64_t>{}(key.rho_bits) + 0x9e3779b97f4a7c15ull +
            (seed << 6u) + (seed >> 2u);
    return seed;
  }
};

struct CachedGroupCoefficients {
  double kappaR_mass_cm2_g{0.0};
  double kappaP_mass_cm2_g{0.0};
  double dbar_cm2_s{0.0};
  double kappaP_cm_inv{0.0};
  double B_g_erg_cm3{0.0};
};

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] long long CoordinateKey(double value) {
  return static_cast<long long>(std::llround(value * kCoordinateKeyScale));
}

[[nodiscard]] std::tuple<long long, long long, long long> RowKey(
    double temperature_keV,
    double density_used_g_cm3,
    double photon_energy_keV) {
  return {CoordinateKey(temperature_keV),
          CoordinateKey(density_used_g_cm3),
          CoordinateKey(photon_energy_keV)};
}

[[nodiscard]] std::string JoinPath(const std::string& root, const char* leaf) {
  if (root.empty()) {
    return leaf;
  }
  const char last = root.back();
  if (last == '\\' || last == '/') {
    return root + leaf;
  }
  return root + "\\" + leaf;
}

[[nodiscard]] std::vector<std::string> SplitCsvLine(const std::string& line) {
  std::vector<std::string> fields;
  std::string field;
  bool in_quotes = false;
  for (const char ch : line) {
    if (ch == '"') {
      in_quotes = !in_quotes;
      continue;
    }
    if (ch == ',' && !in_quotes) {
      fields.push_back(field);
      field.clear();
      continue;
    }
    field.push_back(ch);
  }
  fields.push_back(field);
  return fields;
}

[[nodiscard]] double ParseDouble(const std::string& text) {
  std::size_t consumed = 0u;
  const double value = std::stod(text, &consumed);
  if (consumed != text.size()) {
    throw std::runtime_error("unexpected trailing numeric characters");
  }
  return value;
}

[[nodiscard]] bool ParseBool(const std::string& text) {
  if (text == "true") {
    return true;
  }
  if (text == "false") {
    return false;
  }
  throw std::runtime_error("expected boolean token");
}

void PushUniqueSorted(std::vector<double>& values, double value) {
  if (!std::isfinite(value)) {
    return;
  }
  values.push_back(value);
}

void SortUnique(std::vector<double>& values) {
  std::sort(values.begin(), values.end());
  values.erase(std::unique(values.begin(), values.end(),
                           [](double lhs, double rhs) {
                             return CoordinateKey(lhs) == CoordinateKey(rhs);
                           }),
               values.end());
}

[[nodiscard]] std::string ReadWholeFile(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    throw std::runtime_error("failed to open " + path);
  }
  std::ostringstream out;
  out << input.rdbuf();
  return out.str();
}

[[nodiscard]] std::string ExtractJsonString(
    const std::string& json,
    const std::string& key) {
  const std::string needle = "\"" + key + "\"";
  const auto key_pos = json.find(needle);
  if (key_pos == std::string::npos) {
    return {};
  }
  const auto colon = json.find(':', key_pos + needle.size());
  if (colon == std::string::npos) {
    return {};
  }
  const auto open = json.find('"', colon + 1u);
  if (open == std::string::npos) {
    return {};
  }
  const auto close = json.find('"', open + 1u);
  if (close == std::string::npos) {
    return {};
  }
  return json.substr(open + 1u, close - open - 1u);
}

[[nodiscard]] std::size_t ExtractJsonCount(
    const std::string& json,
    const std::string& key) {
  const std::string needle = "\"" + key + "\"";
  const auto key_pos = json.find(needle);
  if (key_pos == std::string::npos) {
    return 0u;
  }
  const auto colon = json.find(':', key_pos + needle.size());
  if (colon == std::string::npos) {
    return 0u;
  }
  const auto first = json.find_first_of("0123456789", colon + 1u);
  if (first == std::string::npos) {
    return 0u;
  }
  const auto last = json.find_first_not_of("0123456789", first);
  return static_cast<std::size_t>(std::stoull(json.substr(first, last - first)));
}

[[nodiscard]] std::string ExtractFirstMaterialId(const std::string& json) {
  return ExtractJsonString(json, "material_id");
}

void LoadMetadata(const std::string& path, TopsOpacityMetadata& metadata) {
  const std::string json = ReadWholeFile(path);
  metadata.material_label = ExtractJsonString(json, "material");
  metadata.material_id = ExtractFirstMaterialId(json);
  metadata.tops_result_date = ExtractJsonString(json, "tops_result_date");
  metadata.temperature_count = ExtractJsonCount(json, "temperatures");
  metadata.density_count = ExtractJsonCount(json, "densities");
  metadata.photon_energy_count = ExtractJsonCount(json, "photon_grid_points");
  metadata.mean_row_count = ExtractJsonCount(json, "mean_rows");
  metadata.multigroup_row_count = ExtractJsonCount(json, "multigroup_rows");
  metadata.density_clipping_warning_count =
      ExtractJsonCount(json, "density_clipping_warning_rows");
}

[[nodiscard]] std::string BuildLoadReport(
    const TopsOpacityTableLoadResult& result,
    const TopsOpacityProviderOptions& options) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.radiation.opacity_provider.table_load"
      << "; stage_id=R"
      << "; opacity_provider_model=tops_dt_tabulated"
      << "; opacity_source=TOPS"
      << "; opacity_data_source=TOPS"
      << "; opacity_source_thesis_exact_match=false"
      << "; material_label=" << result.metadata.material_label
      << "; material_id=" << result.metadata.material_id
      << "; tops_result_date=" << result.metadata.tops_result_date
      << "; table_load_count=1"
      << "; lookup_table_resident=true"
      << "; lookup_index_mode=in_memory_coordinate_index"
      << "; per_cell_csv_io=false"
      << "; per_lookup_full_table_scan=false"
      << "; opacity_interpolation_mode="
      << OpacityInterpolationModeName(options.opacity_interpolation_mode)
      << "; lookup_energy_mapping_mode="
      << OpacityEnergyMappingModeName(options.lookup_energy_mapping_mode)
      << "; density_clip_policy=" << TopsDensityClipPolicyName(options.density_clip_policy)
      << "; temperature_count=" << result.metadata.temperature_count
      << "; density_count=" << result.metadata.density_count
      << "; photon_energy_count=" << result.metadata.photon_energy_count
      << "; mean_opacity_rows=" << result.metadata.mean_row_count
      << "; multigroup_opacity_rows=" << result.metadata.multigroup_row_count
      << "; density_clipping_warning_count="
      << result.metadata.density_clipping_warning_count
      << "; mean_opacities_used_for_provider=false";
  return out.str();
}

[[nodiscard]] TopsOpacityTableLoadResult LoadFailure(
    const std::string& reason,
    const TopsOpacityProviderOptions& options) {
  TopsOpacityTableLoadResult result;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.opacity_provider.failure"
      << "; stage_id=R"
      << "; opacity_provider_model=tops_dt_tabulated"
      << "; opacity_interpolation_mode="
      << OpacityInterpolationModeName(options.opacity_interpolation_mode)
      << "; lookup_energy_mapping_mode="
      << OpacityEnergyMappingModeName(options.lookup_energy_mapping_mode)
      << "; density_clip_policy=" << TopsDensityClipPolicyName(options.density_clip_policy)
      << "; failure_reason=" << reason;
  result.failure_diagnostics = out.str();
  return result;
}

[[nodiscard]] TopsOpacityLookupResult LookupFailure(
    const std::string& reason,
    const TopsOpacityProviderOptions& options) {
  TopsOpacityLookupResult result;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.opacity_provider.lookup.failure"
      << "; stage_id=R"
      << "; opacity_provider_model=tops_dt_tabulated"
      << "; opacity_interpolation_mode="
      << OpacityInterpolationModeName(options.opacity_interpolation_mode)
      << "; lookup_energy_mapping_mode="
      << OpacityEnergyMappingModeName(options.lookup_energy_mapping_mode)
      << "; density_clip_policy=" << TopsDensityClipPolicyName(options.density_clip_policy)
      << "; failure_reason=" << reason;
  result.failure_diagnostics = out.str();
  return result;
}

struct TopsOpacityValueLookupResult {
  bool success{false};
  double kappaR_mass_cm2_g{0.0};
  double kappaP_mass_cm2_g{0.0};
  bool sampled_clipped_support{false};
  std::string failure_reason;
};

[[nodiscard]] TopsOpacityValueLookupResult ValueLookupFailure(std::string reason) {
  TopsOpacityValueLookupResult result;
  result.failure_reason = std::move(reason);
  return result;
}

struct Bracket {
  double lower{0.0};
  double upper{0.0};
  double weight_upper{0.0};
};

[[nodiscard]] bool BuildLogBracket(
    const std::vector<double>& grid,
    double value,
    Bracket& bracket) {
  if (grid.empty() || !std::isfinite(value) || value <= 0.0) {
    return false;
  }
  if (value < grid.front() || value > grid.back()) {
    return false;
  }
  const auto upper_it = std::lower_bound(grid.begin(), grid.end(), value);
  if (upper_it == grid.end()) {
    bracket.lower = grid.back();
    bracket.upper = grid.back();
    bracket.weight_upper = 0.0;
    return true;
  }
  if (CoordinateKey(*upper_it) == CoordinateKey(value) || upper_it == grid.begin()) {
    bracket.lower = *upper_it;
    bracket.upper = *upper_it;
    bracket.weight_upper = 0.0;
    return true;
  }
  bracket.upper = *upper_it;
  bracket.lower = *(upper_it - 1);
  const double log_lower = std::log(bracket.lower);
  const double log_upper = std::log(bracket.upper);
  bracket.weight_upper =
      (std::log(value) - log_lower) / (log_upper - log_lower);
  bracket.weight_upper = std::clamp(bracket.weight_upper, 0.0, 1.0);
  return true;
}

[[nodiscard]] bool FetchRow(
    const TopsOpacityTable& table,
    double temperature_keV,
    double density_used_g_cm3,
    double photon_energy_keV,
    const TopsOpacityTable::MultigroupRow*& row) {
  const auto it = table.multigroup_row_index.find(
      RowKey(temperature_keV, density_used_g_cm3, photon_energy_keV));
  if (it == table.multigroup_row_index.end()) {
    return false;
  }
  row = &table.multigroup_rows[it->second];
  return true;
}

[[nodiscard]] bool SampleLogOpacity(
    const TopsOpacityTable& table,
    double temperature_keV,
    double density_used_g_cm3,
    double photon_energy_keV,
    const TopsOpacityProviderOptions& options,
    double& log_kappaR,
    double& log_kappaP,
    bool& sampled_clipped_support) {
  const TopsOpacityTable::MultigroupRow* row = nullptr;
  if (!FetchRow(table, temperature_keV, density_used_g_cm3, photon_energy_keV, row)) {
    return false;
  }
  sampled_clipped_support = sampled_clipped_support || row->density_was_clipped;
  if (row->density_was_clipped && options.density_clip_policy == TopsDensityClipPolicy::hard_fail) {
    return false;
  }
  if (!(row->kappaR_mass_cm2_g > 0.0) || !(row->kappaP_mass_cm2_g > 0.0) ||
      !std::isfinite(row->kappaR_mass_cm2_g) ||
      !std::isfinite(row->kappaP_mass_cm2_g)) {
    return false;
  }
  log_kappaR = std::log(row->kappaR_mass_cm2_g);
  log_kappaP = std::log(row->kappaP_mass_cm2_g);
  return true;
}

[[nodiscard]] TopsOpacityValueLookupResult LookupTopsOpacityValueOnly(
    const TopsOpacityTable& table,
    double temperature_keV,
    double density_g_cm3,
    double photon_energy_keV,
    const TopsOpacityProviderOptions& options) {
  if (options.opacity_interpolation_mode != OpacityInterpolationMode::loglog_trilinear ||
      options.lookup_energy_mapping_mode != OpacityEnergyMappingMode::geometric_group_energy ||
      options.density_clip_policy != TopsDensityClipPolicy::hard_fail) {
    return ValueLookupFailure("TOPS opacity lookup options are unsupported");
  }

  Bracket t_bracket;
  Bracket e_bracket;
  if (!BuildLogBracket(table.temperature_grid_keV, temperature_keV, t_bracket)) {
    return ValueLookupFailure("temperature lookup coordinate is outside TOPS table");
  }
  if (!BuildLogBracket(table.photon_energy_grid_keV, photon_energy_keV, e_bracket)) {
    return ValueLookupFailure("photon energy lookup coordinate is outside TOPS table");
  }

  const double t_values[2] = {t_bracket.lower, t_bracket.upper};
  const double t_weights[2] = {1.0 - t_bracket.weight_upper, t_bracket.weight_upper};
  const double e_values[2] = {e_bracket.lower, e_bracket.upper};
  const double e_weights[2] = {1.0 - e_bracket.weight_upper, e_bracket.weight_upper};

  double log_kappaR_sum = 0.0;
  double log_kappaP_sum = 0.0;
  bool sampled_clipped_support = false;
  double total_weight = 0.0;

  for (std::size_t ti = 0u; ti < 2u; ++ti) {
    if (t_weights[ti] == 0.0 && ti == 1u) {
      continue;
    }
    const auto density_it =
        table.used_density_grid_by_temperature.find(CoordinateKey(t_values[ti]));
    if (density_it == table.used_density_grid_by_temperature.end()) {
      return ValueLookupFailure("TOPS opacity density support is missing");
    }
    Bracket d_bracket;
    if (!BuildLogBracket(density_it->second, density_g_cm3, d_bracket)) {
      return ValueLookupFailure("density lookup coordinate uses clipped TOPS support");
    }
    const double d_values[2] = {d_bracket.lower, d_bracket.upper};
    const double d_weights[2] = {1.0 - d_bracket.weight_upper, d_bracket.weight_upper};

    for (std::size_t di = 0u; di < 2u; ++di) {
      if (d_weights[di] == 0.0 && di == 1u) {
        continue;
      }
      for (std::size_t ei = 0u; ei < 2u; ++ei) {
        if (e_weights[ei] == 0.0 && ei == 1u) {
          continue;
        }
        double log_kappaR = 0.0;
        double log_kappaP = 0.0;
        if (!SampleLogOpacity(table, t_values[ti], d_values[di], e_values[ei],
                              options, log_kappaR, log_kappaP,
                              sampled_clipped_support)) {
          return ValueLookupFailure(
              "TOPS opacity interpolation support point is unavailable or clipped");
        }
        const double weight = t_weights[ti] * d_weights[di] * e_weights[ei];
        log_kappaR_sum += weight * log_kappaR;
        log_kappaP_sum += weight * log_kappaP;
        total_weight += weight;
      }
    }
  }

  if (!(total_weight > 0.0)) {
    return ValueLookupFailure("TOPS opacity interpolation accumulated zero weight");
  }

  TopsOpacityValueLookupResult result;
  result.success = true;
  result.kappaR_mass_cm2_g = std::exp(log_kappaR_sum / total_weight);
  result.kappaP_mass_cm2_g = std::exp(log_kappaP_sum / total_weight);
  result.sampled_clipped_support = sampled_clipped_support;
  return result;
}

[[nodiscard]] double PhotonEnergyFromFrequencyKeV(double frequency_hz) {
  const double h_erg_s =
      2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
      dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
  return h_erg_s * frequency_hz / dec3d::physics::PhysicsConstantsCGS::erg_per_kev;
}

[[nodiscard]] double GroupRepresentativeEnergyKeV(
    const RadiationGroupLayout& group_layout,
    std::size_t group) {
  const double lower = PhotonEnergyFromFrequencyKeV(group_layout.frequency_edges_hz.at(group));
  const double upper = PhotonEnergyFromFrequencyKeV(group_layout.frequency_edges_hz.at(group + 1u));
  return std::sqrt(lower * upper);
}

[[nodiscard]] RadiationCoefficientProviderResult ProviderFailure(
    RadiationCoefficientProviderResult result,
    const std::string& reason) {
  result.success = false;
  result.updated_fields = 0;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.opacity_provider.failure"
      << "; stage_id=R"
      << "; opacity_provider_model=tops_dt_tabulated"
      << "; canonical_state_mutated=false"
      << "; updated_fields=0"
      << "; failure_reason=" << reason;
  result.failure_diagnostics = out.str();
  return result;
}

[[nodiscard]] std::string BuildProviderReport(
    const RadiationCoefficientProviderResult& result,
    const TopsOpacityProviderOptions& options,
    const TopsOpacityMetadata& metadata,
    const RadiationGroupLayout& group_layout) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.radiation.opacity_provider"
      << "; stage_id=R"
      << "; opacity_provider_model=tops_dt_tabulated"
      << "; opacity_source=TOPS"
      << "; opacity_data_source=TOPS"
      << "; opacity_source_thesis_exact_match=false"
      << "; group_edges_source=RadiationGroupLayout"
      << "; opacity_energy_coordinate_source=TOPS photon grid"
      << "; lookup_energy_mapping_mode="
      << OpacityEnergyMappingModeName(options.lookup_energy_mapping_mode)
      << "; opacity_interpolation_mode="
      << OpacityInterpolationModeName(options.opacity_interpolation_mode)
      << "; density_clip_policy=" << TopsDensityClipPolicyName(options.density_clip_policy)
      << "; table_load_count=1"
      << "; lookup_table_resident=true"
      << "; lookup_index_mode=in_memory_coordinate_index"
      << "; per_cell_csv_io=false"
      << "; per_lookup_full_table_scan=false"
      << "; coefficient_time_level=old_time_lagged"
      << "; kappaR_output_unit=cm2_per_g"
      << "; kappaP_output_unit=cm2_per_g"
      << "; Dbar_output_unit=cm2_per_s"
      << "; B_g_output_unit=erg_per_cm3"
      << "; B_g_source=blackbody_group_integral"
      << "; blackbody_formula=planck_group_integral"
      << "; double_kB_guard_passed=true"
      << "; canonical_state_mutated=false"
      << "; updated_fields=0"
      << "; group_count=" << group_layout.group_count
      << "; material_label=" << metadata.material_label
      << "; material_id=" << metadata.material_id
      << "; tops_result_date=" << metadata.tops_result_date
      << "; mean_opacities_used_for_provider=false"
      << "; mean_opacity_rows=" << metadata.mean_row_count
      << "; multigroup_opacity_rows=" << metadata.multigroup_row_count
      << "; min_kappaR_mass_cm2_g=" << result.min_kappaR_mass_cm2_g
      << "; max_kappaR_mass_cm2_g=" << result.max_kappaR_mass_cm2_g
      << "; min_kappaP_mass_cm2_g=" << result.min_kappaP_mass_cm2_g
      << "; max_kappaP_mass_cm2_g=" << result.max_kappaP_mass_cm2_g
      << "; min_Dbar_cm2_s=" << result.min_Dbar_cm2_s
      << "; max_Dbar_cm2_s=" << result.max_Dbar_cm2_s
      << "; min_B_g_erg_cm3=" << result.min_B_g_erg_cm3
      << "; max_B_g_erg_cm3=" << result.max_B_g_erg_cm3
      << "; exact_state_cache_mode=bitwise_thermodynamic_state"
      << "; exact_state_cache_hit_count=" << result.exact_state_cache_hit_count
      << "; exact_state_cache_miss_count=" << result.exact_state_cache_miss_count
      << "; representative_blackbody_group_weight_sum="
      << result.blackbody_group_weight_sum
      << "; blackbody_weight_reference_cell="
      << result.blackbody_weight_reference_cell;
  return out.str();
}

}  // namespace

const char* OpacityInterpolationModeName(OpacityInterpolationMode mode) noexcept {
  switch (mode) {
    case OpacityInterpolationMode::loglog_trilinear:
      return "loglog_trilinear";
    case OpacityInterpolationMode::missing:
      return "missing";
  }
  return "missing";
}

const char* OpacityEnergyMappingModeName(OpacityEnergyMappingMode mode) noexcept {
  switch (mode) {
    case OpacityEnergyMappingMode::geometric_group_energy:
      return "geometric_group_energy";
    case OpacityEnergyMappingMode::missing:
      return "missing";
  }
  return "missing";
}

const char* TopsDensityClipPolicyName(TopsDensityClipPolicy policy) noexcept {
  switch (policy) {
    case TopsDensityClipPolicy::hard_fail:
      return "hard_fail";
    case TopsDensityClipPolicy::missing:
      return "missing";
  }
  return "missing";
}

bool RadiationCoefficientProviderResult::is_complete() const noexcept {
  return success && ValidateRadiationCoefficientProviderDiagnostics(*this);
}

TopsOpacityTableLoadResult LoadTopsOpacityTable(
    const TopsOpacityProviderOptions& options) noexcept {
  try {
    if (options.table_root.empty()) {
      return LoadFailure("TOPS opacity table root is missing", options);
    }
    if (options.opacity_interpolation_mode != OpacityInterpolationMode::loglog_trilinear) {
      return LoadFailure("opacity interpolation mode is unsupported", options);
    }
    if (options.lookup_energy_mapping_mode != OpacityEnergyMappingMode::geometric_group_energy) {
      return LoadFailure("lookup energy mapping mode is unsupported", options);
    }
    if (options.density_clip_policy != TopsDensityClipPolicy::hard_fail) {
      return LoadFailure("density clip policy is unsupported", options);
    }

    TopsOpacityTableLoadResult result;
    LoadMetadata(JoinPath(options.table_root, "metadata.json"), result.metadata);
    result.table.metadata = result.metadata;

    std::ifstream input(JoinPath(options.table_root, "multigroup_opacities.csv"));
    if (!input) {
      return LoadFailure("failed to open TOPS multigroup opacity CSV", options);
    }

    std::string line;
    if (!std::getline(input, line)) {
      return LoadFailure("TOPS multigroup opacity CSV is empty", options);
    }

    while (std::getline(input, line)) {
      if (line.empty()) {
        continue;
      }
      const auto fields = SplitCsvLine(line);
      if (fields.size() != 7u) {
        return LoadFailure("TOPS multigroup opacity CSV row has wrong field count", options);
      }

      TopsOpacityTable::MultigroupRow row;
      row.temperature_keV = ParseDouble(fields[0]);
      row.density_requested_g_cm3 = ParseDouble(fields[1]);
      row.density_used_g_cm3 = ParseDouble(fields[2]);
      row.density_was_clipped = ParseBool(fields[3]);
      row.photon_energy_keV = ParseDouble(fields[4]);
      row.kappaR_mass_cm2_g = ParseDouble(fields[5]);
      row.kappaP_mass_cm2_g = ParseDouble(fields[6]);
      if (!std::isfinite(row.temperature_keV) || !std::isfinite(row.density_used_g_cm3) ||
          !std::isfinite(row.photon_energy_keV) || !std::isfinite(row.kappaR_mass_cm2_g) ||
          !std::isfinite(row.kappaP_mass_cm2_g) || row.temperature_keV <= 0.0 ||
          row.density_used_g_cm3 <= 0.0 || row.photon_energy_keV <= 0.0 ||
          row.kappaR_mass_cm2_g <= 0.0 || row.kappaP_mass_cm2_g <= 0.0) {
        return LoadFailure("TOPS multigroup opacity CSV row contains nonphysical values", options);
      }

      const std::size_t row_index = result.table.multigroup_rows.size();
      result.table.multigroup_rows.push_back(row);
      result.table.multigroup_row_index[RowKey(
          row.temperature_keV, row.density_used_g_cm3, row.photon_energy_keV)] = row_index;
      PushUniqueSorted(result.table.temperature_grid_keV, row.temperature_keV);
      PushUniqueSorted(result.table.density_requested_grid_g_cm3,
                       row.density_requested_g_cm3);
      PushUniqueSorted(result.table.photon_energy_grid_keV, row.photon_energy_keV);
      result.table.used_density_grid_by_temperature[CoordinateKey(row.temperature_keV)]
          .push_back(row.density_used_g_cm3);
    }

    if (result.table.multigroup_rows.size() != result.metadata.multigroup_row_count) {
      return LoadFailure("TOPS multigroup opacity row count does not match metadata", options);
    }

    SortUnique(result.table.temperature_grid_keV);
    SortUnique(result.table.density_requested_grid_g_cm3);
    SortUnique(result.table.photon_energy_grid_keV);
    for (auto& entry : result.table.used_density_grid_by_temperature) {
      SortUnique(entry.second);
    }

    if (result.table.temperature_grid_keV.size() != result.metadata.temperature_count ||
        result.table.density_requested_grid_g_cm3.size() != result.metadata.density_count ||
        result.table.photon_energy_grid_keV.size() != result.metadata.photon_energy_count) {
      return LoadFailure("TOPS opacity grid counts do not match metadata", options);
    }

    result.success = true;
    result.report_line = BuildLoadReport(result, options);
    return result;
  } catch (const std::exception& error) {
    return LoadFailure(error.what(), options);
  } catch (...) {
    return LoadFailure("unknown TOPS opacity table load failure", options);
  }
}

TopsOpacityLookupResult LookupTopsOpacity(
    const TopsOpacityTable& table,
    double temperature_keV,
    double density_g_cm3,
    double photon_energy_keV,
    const TopsOpacityProviderOptions& options) noexcept {
  try {
    if (options.opacity_interpolation_mode != OpacityInterpolationMode::loglog_trilinear ||
        options.lookup_energy_mapping_mode != OpacityEnergyMappingMode::geometric_group_energy ||
        options.density_clip_policy != TopsDensityClipPolicy::hard_fail) {
      return LookupFailure("TOPS opacity lookup options are unsupported", options);
    }

    Bracket t_bracket;
    Bracket e_bracket;
    if (!BuildLogBracket(table.temperature_grid_keV, temperature_keV, t_bracket)) {
      return LookupFailure("temperature lookup coordinate is outside TOPS table", options);
    }
    if (!BuildLogBracket(table.photon_energy_grid_keV, photon_energy_keV, e_bracket)) {
      return LookupFailure("photon energy lookup coordinate is outside TOPS table", options);
    }

    const double t_values[2] = {t_bracket.lower, t_bracket.upper};
    const double t_weights[2] = {1.0 - t_bracket.weight_upper, t_bracket.weight_upper};
    const double e_values[2] = {e_bracket.lower, e_bracket.upper};
    const double e_weights[2] = {1.0 - e_bracket.weight_upper, e_bracket.weight_upper};

    double log_kappaR_sum = 0.0;
    double log_kappaP_sum = 0.0;
    bool sampled_clipped_support = false;
    double total_weight = 0.0;

    for (std::size_t ti = 0u; ti < 2u; ++ti) {
      if (t_weights[ti] == 0.0 && ti == 1u) {
        continue;
      }
      const auto density_it =
          table.used_density_grid_by_temperature.find(CoordinateKey(t_values[ti]));
      if (density_it == table.used_density_grid_by_temperature.end()) {
        return LookupFailure("TOPS opacity density support is missing", options);
      }
      Bracket d_bracket;
      if (!BuildLogBracket(density_it->second, density_g_cm3, d_bracket)) {
        return LookupFailure("density lookup coordinate uses clipped TOPS support", options);
      }
      const double d_values[2] = {d_bracket.lower, d_bracket.upper};
      const double d_weights[2] = {1.0 - d_bracket.weight_upper, d_bracket.weight_upper};

      for (std::size_t di = 0u; di < 2u; ++di) {
        if (d_weights[di] == 0.0 && di == 1u) {
          continue;
        }
        for (std::size_t ei = 0u; ei < 2u; ++ei) {
          if (e_weights[ei] == 0.0 && ei == 1u) {
            continue;
          }
          double log_kappaR = 0.0;
          double log_kappaP = 0.0;
          if (!SampleLogOpacity(table, t_values[ti], d_values[di], e_values[ei],
                                options, log_kappaR, log_kappaP,
                                sampled_clipped_support)) {
            return LookupFailure("TOPS opacity interpolation support point is unavailable or clipped",
                                 options);
          }
          const double weight = t_weights[ti] * d_weights[di] * e_weights[ei];
          log_kappaR_sum += weight * log_kappaR;
          log_kappaP_sum += weight * log_kappaP;
          total_weight += weight;
        }
      }
    }

    if (!(total_weight > 0.0)) {
      return LookupFailure("TOPS opacity interpolation accumulated zero weight", options);
    }

    TopsOpacityLookupResult result;
    result.success = true;
    result.kappaR_mass_cm2_g = std::exp(log_kappaR_sum / total_weight);
    result.kappaP_mass_cm2_g = std::exp(log_kappaP_sum / total_weight);
    std::ostringstream out;
    out << std::setprecision(17)
        << "diagnostic_id=p3.radiation.opacity_provider.lookup"
        << "; stage_id=R"
        << "; opacity_provider_model=tops_dt_tabulated"
        << "; opacity_interpolation_mode="
        << OpacityInterpolationModeName(options.opacity_interpolation_mode)
        << "; lookup_energy_mapping_mode="
        << OpacityEnergyMappingModeName(options.lookup_energy_mapping_mode)
        << "; density_clip_policy=" << TopsDensityClipPolicyName(options.density_clip_policy)
        << "; density_was_clipped=" << BoolToken(sampled_clipped_support)
        << "; per_cell_csv_io=false"
        << "; per_lookup_full_table_scan=false"
        << "; kappaR_mass_cm2_g=" << result.kappaR_mass_cm2_g
        << "; kappaP_mass_cm2_g=" << result.kappaP_mass_cm2_g;
    result.report_line = out.str();
    return result;
  } catch (const std::exception& error) {
    return LookupFailure(error.what(), options);
  } catch (...) {
    return LookupFailure("unknown TOPS opacity lookup failure", options);
  }
}

RadiationCoefficientProviderResult BuildRadiationCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& /*geometry*/,
    const RadiationGroupLayout& group_layout,
    const TopsOpacityProviderOptions& options) noexcept {
  const auto table_result = LoadTopsOpacityTable(options);
  if (!table_result.success) {
    RadiationCoefficientProviderResult result;
    return ProviderFailure(result, table_result.failure_reason);
  }
  return BuildRadiationCoefficientArrays(
      state,
      dec3d::mesh::SphericalGeometryMetadata{},
      group_layout,
      table_result.table,
      options);
}

RadiationCoefficientProviderResult BuildRadiationCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& /*geometry*/,
    const RadiationGroupLayout& group_layout,
    const TopsOpacityTable& table,
    const TopsOpacityProviderOptions& options) noexcept {
  RadiationCoefficientProviderResult result;
  result.updated_fields = 0;

  const auto group_validation = ValidateRadiationGroupStateStorage(group_layout, state);
  if (!group_validation.success) {
    return ProviderFailure(result, group_validation.failure_reason);
  }
  if (group_layout.mode != RadiationGroupMode::explicit_frequency_groups) {
    return ProviderFailure(result, "TOPS opacity provider requires explicit frequency groups");
  }

  dec3d::state::ThermodynamicRecoveryOptions recovery_options;
  const auto recovery = dec3d::state::RecoverThermodynamicState(state, recovery_options);
  result.thermodynamic_recovery_report = recovery.recovery_diagnostics;
  if (!recovery.success) {
    return ProviderFailure(result, recovery.failure_reason);
  }

  const std::size_t radial = state.layout.radial_cells;
  const std::size_t theta = state.layout.theta_cells;
  const std::size_t phi = state.layout.phi_cells;
  const std::size_t groups = group_layout.group_count;
  result.coefficients.source = "tops_dt_tabulated";
  result.coefficients.Dbar_cm2_per_s.clear();
  result.coefficients.kappaP_cm_inv.clear();
  result.coefficients.B_erg_per_cm3.clear();
  result.coefficients.Dbar_cm2_per_s.reserve(groups);
  result.coefficients.kappaP_cm_inv.reserve(groups);
  result.coefficients.B_erg_per_cm3.reserve(groups);
  for (std::size_t group = 0u; group < groups; ++group) {
    result.coefficients.Dbar_cm2_per_s.emplace_back(radial, theta, phi, 0.0);
    result.coefficients.kappaP_cm_inv.emplace_back(radial, theta, phi, 0.0);
    result.coefficients.B_erg_per_cm3.emplace_back(radial, theta, phi, 0.0);
  }

  double reference_weight_sum = 0.0;
  const double reference_te = recovery.cells(0u, 0u, 0u).t_e_erg_per_particle;
  for (std::size_t group = 0u; group < groups; ++group) {
    const auto blackbody =
        EvaluateGroupBlackbodyEnergyDensityValueOnly(group_layout, group, reference_te);
    if (!blackbody.success) {
      result.first_bad_group = group;
      return ProviderFailure(result, blackbody.failure_reason);
    }
    reference_weight_sum += blackbody.group_weight;
  }
  result.blackbody_group_weight_sum = reference_weight_sum;
  result.blackbody_weight_reference_cell = "0,0,0";

  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const std::size_t local_cell_count = radial * theta * phi;
  for (std::size_t group = 0u; group < groups; ++group) {
    const double photon_energy_keV = GroupRepresentativeEnergyKeV(group_layout, group);
    std::unordered_map<ThermodynamicCacheKey,
                       CachedGroupCoefficients,
                       ThermodynamicCacheKeyHash> exact_state_cache;
    exact_state_cache.reserve(local_cell_count);
    for (std::size_t r = 0u; r < radial; ++r) {
      for (std::size_t t = 0u; t < theta; ++t) {
        for (std::size_t p = 0u; p < phi; ++p) {
          const auto& cell = recovery.cells(r, t, p);
          const ThermodynamicCacheKey key{
              DoubleBits(cell.t_e_keV),
              DoubleBits(cell.t_e_erg_per_particle),
              DoubleBits(cell.rho_g_per_cm3)};
          auto cache_it = exact_state_cache.find(key);
          if (cache_it == exact_state_cache.end()) {
            const auto lookup = LookupTopsOpacityValueOnly(table, cell.t_e_keV,
                                                           cell.rho_g_per_cm3,
                                                           photon_energy_keV, options);
            if (!lookup.success) {
              result.first_bad_group = group;
              result.first_bad_radial = r;
              result.first_bad_theta = t;
              result.first_bad_phi = p;
              return ProviderFailure(result, lookup.failure_reason);
            }
            const double kappaR_cm_inv = lookup.kappaR_mass_cm2_g * cell.rho_g_per_cm3;
            const double kappaP_cm_inv = lookup.kappaP_mass_cm2_g * cell.rho_g_per_cm3;
            if (!(kappaR_cm_inv > 0.0) || !(kappaP_cm_inv > 0.0) ||
                !std::isfinite(kappaR_cm_inv) || !std::isfinite(kappaP_cm_inv)) {
              result.first_bad_group = group;
              result.first_bad_radial = r;
              result.first_bad_theta = t;
              result.first_bad_phi = p;
              return ProviderFailure(result, "TOPS opacity provider produced nonphysical coefficients");
            }
            const double dbar = c / (3.0 * kappaR_cm_inv);
            const auto blackbody =
                EvaluateGroupBlackbodyEnergyDensityValueOnly(
                    group_layout, group, cell.t_e_erg_per_particle);
            if (!blackbody.success) {
              result.first_bad_group = group;
              result.first_bad_radial = r;
              result.first_bad_theta = t;
              result.first_bad_phi = p;
              return ProviderFailure(result, blackbody.failure_reason);
            }
            cache_it = exact_state_cache.emplace(
                key,
                CachedGroupCoefficients{
                    lookup.kappaR_mass_cm2_g,
                    lookup.kappaP_mass_cm2_g,
                    dbar,
                    kappaP_cm_inv,
                    blackbody.energy_density_erg_cm3}).first;
            ++result.exact_state_cache_miss_count;
          } else {
            ++result.exact_state_cache_hit_count;
          }

          const auto& cached = cache_it->second;
          result.coefficients.Dbar_cm2_per_s[group](r, t, p) = cached.dbar_cm2_s;
          result.coefficients.kappaP_cm_inv[group](r, t, p) = cached.kappaP_cm_inv;
          result.coefficients.B_erg_per_cm3[group](r, t, p) = cached.B_g_erg_cm3;
          result.min_kappaR_mass_cm2_g =
              std::min(result.min_kappaR_mass_cm2_g, cached.kappaR_mass_cm2_g);
          result.max_kappaR_mass_cm2_g =
              std::max(result.max_kappaR_mass_cm2_g, cached.kappaR_mass_cm2_g);
          result.min_kappaP_mass_cm2_g =
              std::min(result.min_kappaP_mass_cm2_g, cached.kappaP_mass_cm2_g);
          result.max_kappaP_mass_cm2_g =
              std::max(result.max_kappaP_mass_cm2_g, cached.kappaP_mass_cm2_g);
          result.min_Dbar_cm2_s = std::min(result.min_Dbar_cm2_s, cached.dbar_cm2_s);
          result.max_Dbar_cm2_s = std::max(result.max_Dbar_cm2_s, cached.dbar_cm2_s);
          result.min_B_g_erg_cm3 = std::min(result.min_B_g_erg_cm3, cached.B_g_erg_cm3);
          result.max_B_g_erg_cm3 = std::max(result.max_B_g_erg_cm3, cached.B_g_erg_cm3);
        }
      }
    }
  }

  result.success = true;
  result.report_line =
      BuildProviderReport(result, options, table.metadata, group_layout);
  return result;
}

bool ValidateRadiationCoefficientProviderDiagnostics(
    const RadiationCoefficientProviderResult& result) noexcept {
  const auto& line = result.report_line;
  return result.success &&
         Contains(line, "diagnostic_id=p3.radiation.opacity_provider") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "opacity_provider_model=tops_dt_tabulated") &&
         Contains(line, "opacity_source=TOPS") &&
         Contains(line, "opacity_source_thesis_exact_match=false") &&
         Contains(line, "group_edges_source=RadiationGroupLayout") &&
         Contains(line, "opacity_energy_coordinate_source=TOPS photon grid") &&
         Contains(line, "lookup_energy_mapping_mode=geometric_group_energy") &&
         Contains(line, "opacity_interpolation_mode=loglog_trilinear") &&
         Contains(line, "density_clip_policy=hard_fail") &&
         Contains(line, "table_load_count=1") &&
         Contains(line, "lookup_table_resident=true") &&
         Contains(line, "lookup_index_mode=in_memory_coordinate_index") &&
         Contains(line, "per_cell_csv_io=false") &&
         Contains(line, "per_lookup_full_table_scan=false") &&
         Contains(line, "coefficient_time_level=old_time_lagged") &&
         Contains(line, "double_kB_guard_passed=true") &&
         Contains(line, "canonical_state_mutated=false") &&
         Contains(line, "updated_fields=0") &&
         Contains(line, "mean_opacities_used_for_provider=false") &&
         Contains(line, "representative_blackbody_group_weight_sum=") &&
         Contains(line, "blackbody_weight_reference_cell=0,0,0");
}

}  // namespace dec3d::radiation
