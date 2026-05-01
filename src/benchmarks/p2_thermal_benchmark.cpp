#include "benchmarks/p2_thermal_benchmark.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace dec3d::benchmarks {
namespace {

[[nodiscard]] bool IsFinitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] bool HasToken(const std::string& report, const char* token) noexcept {
  return report.find(token) != std::string::npos;
}

[[nodiscard]] std::string SerializeTimes(const std::vector<double>& times) {
  std::ostringstream out;
  out << std::setprecision(17);
  for (std::size_t index = 0; index < times.size(); ++index) {
    if (index != 0u) {
      out << ",";
    }
    out << times[index];
  }
  return out.str();
}

}  // namespace

const char* ToString(BenchmarkBoundaryModel model) noexcept {
  switch (model) {
    case BenchmarkBoundaryModel::thesis_1d_spherical_zero_flux:
      return "thesis_1d_spherical_zero_flux";
    case BenchmarkBoundaryModel::full_sphere_scalar_remap:
      return "full_sphere_scalar_remap";
  }
  return "unknown";
}

const char* ToString(BenchmarkThermalKappaModel model) noexcept {
  switch (model) {
    case BenchmarkThermalKappaModel::spitzer_no_degeneracy:
      return "spitzer_no_degeneracy";
    case BenchmarkThermalKappaModel::lee_more_with_degeneracy:
      return "lee_more_with_degeneracy";
    case BenchmarkThermalKappaModel::constant_kappa_exact:
      return "constant_kappa_exact";
  }
  return "unknown";
}

const char* ToString(BenchmarkTauModel model) noexcept {
  switch (model) {
    case BenchmarkTauModel::disabled:
      return "disabled";
    case BenchmarkTauModel::thesis_spitzer_eq_5_241:
      return "thesis_spitzer_eq_5_241";
  }
  return "unknown";
}

WooThermalBenchmarkDescriptor MakeWooThermalBenchmarkDescriptor(
    WooThermalBenchmarkCase case_kind) {
  WooThermalBenchmarkDescriptor descriptor;
  descriptor.case_kind = case_kind;
  descriptor.output_times_s = {1.0e-12, 2.0e-12, 5.0e-12};
  descriptor.comparison_time_s = 5.0e-12;

  switch (case_kind) {
    case WooThermalBenchmarkCase::spitzer_pure_thermal:
      descriptor.case_id = "p2_8b_spitzer_pure_thermal";
      descriptor.woo_figure = "5.32";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::spitzer_no_degeneracy;
      descriptor.tau_model = BenchmarkTauModel::disabled;
      descriptor.equilibration_stage_enabled = false;
      break;
    case WooThermalBenchmarkCase::lee_more_pure_thermal:
      descriptor.case_id = "p2_8c_lee_more_pure_thermal";
      descriptor.woo_figure = "5.33";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::lee_more_with_degeneracy;
      descriptor.tau_model = BenchmarkTauModel::disabled;
      descriptor.equilibration_stage_enabled = false;
      break;
    case WooThermalBenchmarkCase::lee_more_thermal_plus_ei:
      descriptor.case_id = "p2_8d_lee_more_thermal_plus_ei";
      descriptor.woo_figure = "5.34";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::lee_more_with_degeneracy;
      descriptor.tau_model = BenchmarkTauModel::thesis_spitzer_eq_5_241;
      descriptor.equilibration_stage_enabled = true;
      break;
    case WooThermalBenchmarkCase::exact_l1m1_constant_kappa:
      descriptor.case_id = "p2_8e_exact_l1m1_constant_kappa";
      descriptor.woo_figure = "exact_l1m1_spherical_bessel_neumann";
      descriptor.boundary_model = BenchmarkBoundaryModel::full_sphere_scalar_remap;
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::constant_kappa_exact;
      descriptor.tau_model = BenchmarkTauModel::disabled;
      descriptor.equilibration_stage_enabled = false;
      descriptor.output_times_s = {2.5e9, 5.0e9, 7.5e9, 1.0e10};
      descriptor.comparison_time_s = 1.0e10;
      break;
    case WooThermalBenchmarkCase::exact_l0_radial_constant_kappa:
      descriptor.case_id = "p2_8f_exact_l0_radial_constant_kappa";
      descriptor.woo_figure = "exact_l0_spherical_bessel_neumann";
      descriptor.boundary_model = BenchmarkBoundaryModel::full_sphere_scalar_remap;
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::constant_kappa_exact;
      descriptor.tau_model = BenchmarkTauModel::disabled;
      descriptor.equilibration_stage_enabled = false;
      descriptor.output_times_s = {2.5e9, 5.0e9, 7.5e9, 1.0e10};
      descriptor.comparison_time_s = 1.0e10;
      break;
    case WooThermalBenchmarkCase::flux_limiter_stress:
      descriptor.case_id = "p2_8g_flux_limiter_stress";
      descriptor.woo_figure = "operator_flux_limiter_stress";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::spitzer_no_degeneracy;
      descriptor.tau_model = BenchmarkTauModel::disabled;
      descriptor.equilibration_stage_enabled = false;
      break;
    case WooThermalBenchmarkCase::ei_0d_thesis_spitzer:
      descriptor.case_id = "p2_8h_ei_0d_thesis_spitzer";
      descriptor.woo_figure = "operator_ei_exact_exponential";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::constant_kappa_exact;
      descriptor.tau_model = BenchmarkTauModel::thesis_spitzer_eq_5_241;
      descriptor.initial_ti_policy = "ion_uniform_0p5keV";
      descriptor.thermal_stage_enabled = false;
      descriptor.equilibration_stage_enabled = true;
      break;
    case WooThermalBenchmarkCase::distributed_hte_smoke:
      descriptor.case_id = "p2_8i_distributed_hte_smoke";
      descriptor.woo_figure = "operator_distributed_hte_smoke";
      descriptor.thermal_kappa_model = BenchmarkThermalKappaModel::lee_more_with_degeneracy;
      descriptor.tau_model = BenchmarkTauModel::thesis_spitzer_eq_5_241;
      descriptor.equilibration_stage_enabled = true;
      break;
  }

  return descriptor;
}

WooThermalReferencePolicy MakeNoReferencePolicy() {
  WooThermalReferencePolicy policy;
  policy.lilac_reference_available = false;
  policy.digitized_thesis_points_available = false;
  policy.parity_claim_allowed = false;
  policy.comparison_mode = "no_external_reference";
  policy.reference_kind = "none";
  policy.reference_source = "none";
  return policy;
}

BenchmarkValidationResult ValidateWooThermalBenchmarkDescriptor(
    const WooThermalBenchmarkDescriptor& descriptor) noexcept {
  BenchmarkValidationResult result;
  if (descriptor.case_id.empty()) {
    result.failure_reason = "case_id is required";
  } else if (descriptor.woo_figure.empty()) {
    result.failure_reason = "woo_figure is required";
  } else if (!IsFinitePositive(descriptor.r_max_cm) ||
             !IsFinitePositive(descriptor.r0_cm) ||
             descriptor.r0_cm >= descriptor.r_max_cm) {
    result.failure_reason = "invalid slab radii";
  } else if (!IsFinitePositive(descriptor.te_inner_keV) ||
             !IsFinitePositive(descriptor.te_outer_keV) ||
             !IsFinitePositive(descriptor.rho_inner_g_cm3) ||
             !IsFinitePositive(descriptor.rho_outer_g_cm3)) {
    result.failure_reason = "invalid slab initial state";
  } else if (descriptor.initial_ti_policy != "equal_to_initial_te" &&
             !(descriptor.case_kind == WooThermalBenchmarkCase::ei_0d_thesis_spitzer &&
               descriptor.initial_ti_policy == "ion_uniform_0p5keV")) {
    result.failure_reason = "unsupported initial_ti_policy";
  } else if (descriptor.case_kind != WooThermalBenchmarkCase::exact_l1m1_constant_kappa &&
             descriptor.case_kind != WooThermalBenchmarkCase::exact_l0_radial_constant_kappa &&
             descriptor.boundary_model !=
                 BenchmarkBoundaryModel::thesis_1d_spherical_zero_flux) {
    result.failure_reason = "direct Woo benchmark requires thesis_1d_spherical_zero_flux";
  } else if ((descriptor.case_kind == WooThermalBenchmarkCase::exact_l1m1_constant_kappa ||
              descriptor.case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) &&
             descriptor.boundary_model != BenchmarkBoundaryModel::full_sphere_scalar_remap) {
    result.failure_reason = "exact spherical harmonic benchmark requires full_sphere_scalar_remap";
  } else if (!descriptor.thermal_stage_enabled &&
             descriptor.case_kind != WooThermalBenchmarkCase::ei_0d_thesis_spitzer) {
    result.failure_reason = "thermal stage must be enabled";
  } else if (descriptor.momentum_update_enabled) {
    result.failure_reason = "momentum update is forbidden for Woo 5.5.4 thermal benchmark";
  } else if (descriptor.output_times_s.empty() ||
             !IsFinitePositive(descriptor.comparison_time_s)) {
    result.failure_reason = "explicit output_times_s and comparison_time_s are required";
  } else {
    bool comparison_time_present = false;
    for (const double time_s : descriptor.output_times_s) {
      if (!IsFinitePositive(time_s)) {
        result.failure_reason = "output_times_s must be finite and positive";
        result.failure_diagnostics =
            "diagnostic_id=p2.thermal_benchmark.descriptor_failure; failure_reason=" +
            result.failure_reason;
        return result;
      }
      comparison_time_present =
          comparison_time_present ||
          std::abs(time_s - descriptor.comparison_time_s) <=
              1.0e-15 * std::max(1.0, descriptor.comparison_time_s);
    }
    if (!comparison_time_present) {
      result.failure_reason = "comparison_time_s must appear in output_times_s";
    }
  }

  result.success = result.failure_reason.empty();
  if (!result.success) {
    result.failure_diagnostics =
        "diagnostic_id=p2.thermal_benchmark.descriptor_failure; failure_reason=" +
        result.failure_reason;
  }
  return result;
}

std::string BuildWooThermalBenchmarkDiagnosticsLine(
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalReferencePolicy& reference_policy) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.thermal_benchmark"
      << "; case_id=" << descriptor.case_id
      << "; woo_figure=" << descriptor.woo_figure
      << "; unit_system=cgs"
      << "; temperature_internal_unit=erg_per_particle"
      << "; temperature_output_unit=keV"
      << "; r_max_cm=" << descriptor.r_max_cm
      << "; r0_cm=" << descriptor.r0_cm
      << "; te_inner_keV=" << descriptor.te_inner_keV
      << "; te_outer_keV=" << descriptor.te_outer_keV
      << "; rho_inner_g_cm3=" << descriptor.rho_inner_g_cm3
      << "; rho_outer_g_cm3=" << descriptor.rho_outer_g_cm3
      << "; initial_ti_policy=" << descriptor.initial_ti_policy
      << "; thermal_stage_executed=" << (descriptor.thermal_stage_enabled ? "true" : "false")
      << "; equilibration_stage_executed="
      << (descriptor.equilibration_stage_enabled ? "true" : "false")
      << "; momentum_update_executed=" << (descriptor.momentum_update_enabled ? "true" : "false")
      << "; thermal_kappa_model_requested=" << ToString(descriptor.thermal_kappa_model)
      << "; thermal_kappa_model_executed=" << ToString(descriptor.thermal_kappa_model)
      << "; tau_model_requested=" << ToString(descriptor.tau_model)
      << "; tau_model_executed=" << ToString(descriptor.tau_model)
      << "; boundary_model=" << ToString(descriptor.boundary_model)
      << "; benchmark_compare_time_policy=explicit_output_times_from_case_descriptor"
      << "; output_times_s=" << SerializeTimes(descriptor.output_times_s)
      << "; comparison_time_s=" << descriptor.comparison_time_s
      << "; comparison_profile_source=profile_at_comparison_time"
      << "; lilac_reference_available="
      << (reference_policy.lilac_reference_available ? "true" : "false")
      << "; parity_claim_allowed=" << (reference_policy.parity_claim_allowed ? "true" : "false")
      << "; comparison_mode=" << reference_policy.comparison_mode
      << "; profile_artifacts_written=false"
      << "; energy_budget_written=false"
      << "; missing_diagnostics_fail=true";
  if (descriptor.case_kind == WooThermalBenchmarkCase::exact_l1m1_constant_kappa) {
    out << "; benchmark_family=operator_exact_solution"
        << "; exact_solution=l1m1_spherical_bessel_neumann"
        << "; exact_projection=cell_volume_average"
        << "; exact_beta=2.081575977818"
        << "; exact_mode_l=1"
        << "; exact_mode_m=1";
  } else if (descriptor.case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
    out << "; benchmark_family=operator_exact_solution"
        << "; exact_solution=l0_spherical_bessel_neumann"
        << "; exact_projection=cell_volume_average"
        << "; exact_beta=4.4934094579090642"
        << "; exact_mode_l=0"
        << "; exact_mode_m=0";
  } else if (descriptor.case_kind == WooThermalBenchmarkCase::flux_limiter_stress) {
    out << "; benchmark_family=operator_flux_limiter_stress"
        << "; electron_flux_limiter_required=true";
  } else if (descriptor.case_kind == WooThermalBenchmarkCase::ei_0d_thesis_spitzer) {
    out << "; benchmark_family=operator_ei_0d_exact_exponential"
        << "; thermal_stage_executed=false"
        << "; equilibration_stage_executed=true";
  } else if (descriptor.case_kind == WooThermalBenchmarkCase::distributed_hte_smoke) {
    out << "; benchmark_family=operator_distributed_hte_smoke"
        << "; stage_order=H,T,E";
  }
  return out.str();
}

bool ValidateWooThermalBenchmarkDiagnostics(const std::string& report_line) noexcept {
  const bool exact_case =
      HasToken(report_line, "case_id=p2_8e_exact_l1m1_constant_kappa") ||
      HasToken(report_line, "case_id=p2_8f_exact_l0_radial_constant_kappa");
  const bool boundary_valid =
      exact_case
          ? HasToken(report_line, "boundary_model=full_sphere_scalar_remap")
          : HasToken(report_line, "boundary_model=thesis_1d_spherical_zero_flux");
  const bool exact_projection_valid =
      !exact_case || HasToken(report_line, "exact_projection=cell_volume_average");
  return HasToken(report_line, "diagnostic_id=p2.thermal_benchmark") &&
         HasToken(report_line, "unit_system=cgs") &&
         HasToken(report_line, "temperature_internal_unit=erg_per_particle") &&
         HasToken(report_line, "temperature_output_unit=keV") &&
         boundary_valid &&
         HasToken(
             report_line,
             "benchmark_compare_time_policy=explicit_output_times_from_case_descriptor") &&
         HasToken(report_line, "output_times_s=") &&
         HasToken(report_line, "comparison_time_s=") &&
         HasToken(report_line, "comparison_profile_source=profile_at_comparison_time") &&
         HasToken(report_line, "lilac_reference_available=") &&
         HasToken(report_line, "parity_claim_allowed=") &&
         HasToken(report_line, "momentum_update_executed=false") &&
         HasToken(report_line, "missing_diagnostics_fail=true") &&
         exact_projection_valid &&
         (!exact_case ||
          (HasToken(report_line, "exact_solution=l1m1_spherical_bessel_neumann") ||
           HasToken(report_line, "exact_solution=l0_spherical_bessel_neumann"))) &&
         !HasToken(report_line, "parity_claim_allowed=true");
}

bool WooThermalArtifactManifest::ContainsRequiredFile(
    const std::string& file_name) const noexcept {
  return std::find(required_files.begin(), required_files.end(), file_name) != required_files.end();
}

WooThermalArtifactManifest BuildWooThermalArtifactManifest(
    const WooThermalBenchmarkDescriptor& descriptor,
    const std::string& output_root) {
  WooThermalArtifactManifest manifest;
  manifest.output_directory = output_root + "/" + descriptor.case_id;
  if (descriptor.case_kind == WooThermalBenchmarkCase::exact_l1m1_constant_kappa ||
      descriptor.case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
    manifest.required_files = {
        "dec3d.out",
        "exact_benchmark_summary.md",
        "benchmark_diagnostics.txt",
        "profile_exact_vs_numeric_step000001.txt",
        "profile_exact_vs_numeric_step000002.txt",
        "profile_exact_vs_numeric_step000003.txt",
        "profile_exact_vs_numeric_step000004.txt",
        "slice_xz_exact_vs_numeric_step000001.txt"};
    return manifest;
  }
  manifest.required_files = {
      "dec3d.out",
      "benchmark_summary.md",
      "benchmark_diagnostics.txt",
      "profile_step000001.txt",
      "te_profile_step000001.txt",
      "ti_profile_step000001.txt",
      "rho_profile_step000001.txt",
      "kappa_e_profile_step000001.txt",
      "kappa_i_profile_step000001.txt",
      "tau_ei_profile_step000001.txt",
      "energy_budget_vs_time.txt",
      "comparison_metrics.txt"};
  return manifest;
}

}  // namespace dec3d::benchmarks
