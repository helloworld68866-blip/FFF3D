#pragma once

#include "radiation/groups/radiation_group_layout.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::benchmarks {

enum class P3RadiationBenchmarkLayer {
  b0_operator_budget_oracle,
  b1_thesis_geometry_slab
};

enum class P3RadiationBenchmarkCase {
  b0_one_cell_source_oracle,
  b0_one_cell_marshak_source_oracle,
  b0_zero_flux_multigroup_conservation,
  b0_marshak_multigroup_leak_closure,
  b0_group_count_one_regression,
  b1_woo_5_25_slab_no_limiter
};

struct P3RadiationReferencePolicy {
  bool lilac_reference_available{false};
  bool parity_claim_allowed{false};
  std::string comparison_mode{"no_external_reference"};
  std::string reference_kind{"none"};
  std::string reference_source{"none"};
};

struct P3RadiationBenchmarkDescriptor {
  P3RadiationBenchmarkCase case_kind{
      P3RadiationBenchmarkCase::b0_one_cell_source_oracle};
  P3RadiationBenchmarkLayer layer{
      P3RadiationBenchmarkLayer::b0_operator_budget_oracle};
  std::string case_id;
  std::string claim_level{"oracle"};
  double r_max_cm{1.0e-2};
  double r0_cm{5.0e-3};
  double te_inner_keV{5.0};
  double te_outer_keV{0.5};
  double rho_inner_g_cm3{50.0};
  double rho_outer_g_cm3{100.0};
  std::string initial_ti_policy{"equal_to_initial_te"};
  std::string initial_radiation_policy{"blackbody_from_initial_te"};
  std::size_t group_count{1u};
  std::string group_edges_source{"RadiationGroupLayout"};
  std::string opacity_provider{"tops_dt_tabulated"};
  std::string opacity_data_source{"TOPS"};
  std::string opacity_energy_coordinate_source{"TOPS photon grid"};
  std::string lookup_energy_mapping_mode{"geometric_group_energy"};
  bool opacity_source_thesis_exact_match{false};
  std::string radiation_boundary_model{"contract_zero_flux_or_scalar_remap"};
  bool marshak_enabled{false};
  std::string radiation_flux_limiter{"disabled"};
  bool advection_enabled{false};
  bool radiation_pressure_work_enabled{false};
  bool production_runtime_registration{false};
  bool lilac_reference_available{false};
  bool parity_claim_allowed{false};
  bool no_limiter_pre_parity_slab{false};
  std::vector<double> output_times_s{1.0e-18};
  double comparison_time_s{1.0e-18};
};

struct P3RadiationBenchmarkValidationResult {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct P3RadiationBenchmarkArtifactManifest {
  std::string output_directory;
  std::vector<std::string> required_files;

  [[nodiscard]] bool ContainsRequiredFile(const std::string& file_name) const noexcept;
};

[[nodiscard]] const char* ToString(P3RadiationBenchmarkLayer layer) noexcept;
[[nodiscard]] const char* ToString(P3RadiationBenchmarkCase case_kind) noexcept;

[[nodiscard]] P3RadiationBenchmarkDescriptor MakeP3RadiationBenchmarkDescriptor(
    P3RadiationBenchmarkCase case_kind);

[[nodiscard]] dec3d::radiation::RadiationGroupLayout MakeP3B1TwelveGroupLayout();

[[nodiscard]] P3RadiationReferencePolicy MakeNoP3RadiationReferencePolicy();

[[nodiscard]] P3RadiationBenchmarkValidationResult
ValidateP3RadiationBenchmarkDescriptor(
    const P3RadiationBenchmarkDescriptor& descriptor) noexcept;

[[nodiscard]] std::string BuildP3RadiationBenchmarkDiagnosticsLine(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy);

[[nodiscard]] bool ValidateP3RadiationBenchmarkDiagnostics(
    const std::string& report_line) noexcept;

[[nodiscard]] P3RadiationBenchmarkArtifactManifest
BuildP3RadiationBenchmarkArtifactManifest(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const std::string& output_root);

}  // namespace dec3d::benchmarks
