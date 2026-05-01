#pragma once

#include <string>
#include <vector>

namespace dec3d::benchmarks {

enum class WooThermalBenchmarkCase {
  spitzer_pure_thermal,
  lee_more_pure_thermal,
  lee_more_thermal_plus_ei,
  exact_l1m1_constant_kappa,
  exact_l0_radial_constant_kappa,
  flux_limiter_stress,
  ei_0d_thesis_spitzer,
  distributed_hte_smoke,
};

enum class BenchmarkBoundaryModel {
  thesis_1d_spherical_zero_flux,
  full_sphere_scalar_remap,
};

enum class BenchmarkThermalKappaModel {
  spitzer_no_degeneracy,
  lee_more_with_degeneracy,
  constant_kappa_exact,
};

enum class BenchmarkTauModel {
  disabled,
  thesis_spitzer_eq_5_241,
};

struct WooThermalReferencePolicy {
  bool lilac_reference_available{false};
  bool digitized_thesis_points_available{false};
  bool parity_claim_allowed{false};
  std::string comparison_mode{"no_external_reference"};
  std::string reference_kind{"none"};
  std::string reference_source{"none"};
};

struct WooThermalBenchmarkDescriptor {
  WooThermalBenchmarkCase case_kind{WooThermalBenchmarkCase::spitzer_pure_thermal};
  std::string case_id;
  std::string woo_figure;
  double r_max_cm{1.0e-2};
  double r0_cm{5.0e-3};
  double te_inner_keV{5.0};
  double te_outer_keV{0.5};
  double rho_inner_g_cm3{50.0};
  double rho_outer_g_cm3{100.0};
  std::string initial_ti_policy{"equal_to_initial_te"};
  BenchmarkBoundaryModel boundary_model{BenchmarkBoundaryModel::thesis_1d_spherical_zero_flux};
  BenchmarkThermalKappaModel thermal_kappa_model{
      BenchmarkThermalKappaModel::spitzer_no_degeneracy};
  BenchmarkTauModel tau_model{BenchmarkTauModel::disabled};
  bool thermal_stage_enabled{true};
  bool equilibration_stage_enabled{false};
  bool momentum_update_enabled{false};
  std::vector<double> output_times_s;
  double comparison_time_s{0.0};
};

struct BenchmarkValidationResult {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct WooThermalArtifactManifest {
  std::string output_directory;
  std::vector<std::string> required_files;

  [[nodiscard]] bool ContainsRequiredFile(const std::string& file_name) const noexcept;
};

[[nodiscard]] const char* ToString(BenchmarkBoundaryModel model) noexcept;
[[nodiscard]] const char* ToString(BenchmarkThermalKappaModel model) noexcept;
[[nodiscard]] const char* ToString(BenchmarkTauModel model) noexcept;

[[nodiscard]] WooThermalBenchmarkDescriptor MakeWooThermalBenchmarkDescriptor(
    WooThermalBenchmarkCase case_kind);

[[nodiscard]] WooThermalReferencePolicy MakeNoReferencePolicy();

[[nodiscard]] BenchmarkValidationResult ValidateWooThermalBenchmarkDescriptor(
    const WooThermalBenchmarkDescriptor& descriptor) noexcept;

[[nodiscard]] std::string BuildWooThermalBenchmarkDiagnosticsLine(
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalReferencePolicy& reference_policy);

[[nodiscard]] bool ValidateWooThermalBenchmarkDiagnostics(
    const std::string& report_line) noexcept;

[[nodiscard]] WooThermalArtifactManifest BuildWooThermalArtifactManifest(
    const WooThermalBenchmarkDescriptor& descriptor,
    const std::string& output_root);

}  // namespace dec3d::benchmarks
