#pragma once

#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace dec3d::radiation {

enum class OpacityInterpolationMode {
  missing,
  loglog_trilinear
};

enum class OpacityEnergyMappingMode {
  missing,
  geometric_group_energy
};

enum class TopsDensityClipPolicy {
  missing,
  hard_fail
};

struct TopsOpacityProviderOptions {
  std::string table_root;
  OpacityInterpolationMode opacity_interpolation_mode{OpacityInterpolationMode::missing};
  OpacityEnergyMappingMode lookup_energy_mapping_mode{OpacityEnergyMappingMode::missing};
  TopsDensityClipPolicy density_clip_policy{TopsDensityClipPolicy::missing};
};

struct TopsOpacityMetadata {
  std::string material_label;
  std::string material_id;
  std::string tops_result_date;
  std::size_t temperature_count{0u};
  std::size_t density_count{0u};
  std::size_t photon_energy_count{0u};
  std::size_t mean_row_count{0u};
  std::size_t multigroup_row_count{0u};
  std::size_t density_clipping_warning_count{0u};
};

struct TopsOpacityTable {
  struct MultigroupRow {
    double temperature_keV{0.0};
    double density_requested_g_cm3{0.0};
    double density_used_g_cm3{0.0};
    bool density_was_clipped{false};
    double photon_energy_keV{0.0};
    double kappaR_mass_cm2_g{0.0};
    double kappaP_mass_cm2_g{0.0};
  };

  TopsOpacityMetadata metadata;
  std::vector<double> temperature_grid_keV;
  std::vector<double> density_requested_grid_g_cm3;
  std::vector<double> photon_energy_grid_keV;
  std::vector<MultigroupRow> multigroup_rows;
  std::map<std::tuple<long long, long long, long long>, std::size_t>
      multigroup_row_index;
  std::map<long long, std::vector<double>> used_density_grid_by_temperature;
};

struct TopsOpacityTableLoadResult {
  bool success{false};
  TopsOpacityTable table;
  TopsOpacityMetadata metadata;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct TopsOpacityLookupResult {
  bool success{false};
  double kappaR_mass_cm2_g{0.0};
  double kappaP_mass_cm2_g{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct RadiationCoefficientProviderResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  MultigroupGrayRadiationCoefficients coefficients;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
  std::string thermodynamic_recovery_report;
  std::size_t first_bad_group{std::numeric_limits<std::size_t>::max()};
  std::size_t first_bad_radial{std::numeric_limits<std::size_t>::max()};
  std::size_t first_bad_theta{std::numeric_limits<std::size_t>::max()};
  std::size_t first_bad_phi{std::numeric_limits<std::size_t>::max()};

  double min_kappaR_mass_cm2_g{std::numeric_limits<double>::infinity()};
  double max_kappaR_mass_cm2_g{0.0};
  double min_kappaP_mass_cm2_g{std::numeric_limits<double>::infinity()};
  double max_kappaP_mass_cm2_g{0.0};
  double min_Dbar_cm2_s{std::numeric_limits<double>::infinity()};
  double max_Dbar_cm2_s{0.0};
  double min_B_g_erg_cm3{std::numeric_limits<double>::infinity()};
  double max_B_g_erg_cm3{0.0};
  double blackbody_group_weight_sum{0.0};
  std::string blackbody_weight_reference_cell{"0,0,0"};
  std::size_t exact_state_cache_hit_count{0u};
  std::size_t exact_state_cache_miss_count{0u};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] const char* OpacityInterpolationModeName(OpacityInterpolationMode mode) noexcept;
[[nodiscard]] const char* OpacityEnergyMappingModeName(OpacityEnergyMappingMode mode) noexcept;
[[nodiscard]] const char* TopsDensityClipPolicyName(TopsDensityClipPolicy policy) noexcept;

[[nodiscard]] TopsOpacityTableLoadResult LoadTopsOpacityTable(
    const TopsOpacityProviderOptions& options) noexcept;

[[nodiscard]] TopsOpacityLookupResult LookupTopsOpacity(
    const TopsOpacityTable& table,
    double temperature_keV,
    double density_g_cm3,
    double photon_energy_keV,
    const TopsOpacityProviderOptions& options) noexcept;

[[nodiscard]] RadiationCoefficientProviderResult BuildRadiationCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const RadiationGroupLayout& group_layout,
    const TopsOpacityProviderOptions& options) noexcept;

[[nodiscard]] RadiationCoefficientProviderResult BuildRadiationCoefficientArrays(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const RadiationGroupLayout& group_layout,
    const TopsOpacityTable& table,
    const TopsOpacityProviderOptions& options) noexcept;

[[nodiscard]] bool ValidateRadiationCoefficientProviderDiagnostics(
    const RadiationCoefficientProviderResult& result) noexcept;

}  // namespace dec3d::radiation
