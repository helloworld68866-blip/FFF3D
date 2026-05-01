#include "benchmarks/p3_radiation_benchmark.hpp"

#include "physics/units/physical_constants.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace dec3d::benchmarks {
namespace {

[[nodiscard]] bool HasToken(const std::string& report, const char* token) noexcept {
  return report.find(token) != std::string::npos;
}

[[nodiscard]] bool IsFinitePositive(double value) noexcept {
  return std::isfinite(value) && value > 0.0;
}

[[nodiscard]] double PlanckConstantErgS() noexcept {
  return 2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
         dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
}

[[nodiscard]] double FrequencyFromKeV(double energy_keV) noexcept {
  return dec3d::physics::ErgFromKeV(energy_keV) / PlanckConstantErgS();
}

}  // namespace

const char* ToString(P3RadiationBenchmarkLayer layer) noexcept {
  switch (layer) {
    case P3RadiationBenchmarkLayer::b0_operator_budget_oracle:
      return "P3-B0";
    case P3RadiationBenchmarkLayer::b1_thesis_geometry_slab:
      return "P3-B1";
  }
  return "unknown";
}

const char* ToString(P3RadiationBenchmarkCase case_kind) noexcept {
  switch (case_kind) {
    case P3RadiationBenchmarkCase::b0_one_cell_source_oracle:
      return "p3_b0_one_cell_source_oracle";
    case P3RadiationBenchmarkCase::b0_one_cell_marshak_source_oracle:
      return "p3_b0_one_cell_marshak_source_oracle";
    case P3RadiationBenchmarkCase::b0_zero_flux_multigroup_conservation:
      return "p3_b0_zero_flux_multigroup_energy_conservation";
    case P3RadiationBenchmarkCase::b0_marshak_multigroup_leak_closure:
      return "p3_b0_marshak_multigroup_leak_closure";
    case P3RadiationBenchmarkCase::b0_group_count_one_regression:
      return "p3_b0_group_count_one_regression";
    case P3RadiationBenchmarkCase::b1_woo_5_25_slab_no_limiter:
      return "p3_b1_woo_5_25_slab_no_limiter";
  }
  return "unknown";
}

P3RadiationBenchmarkDescriptor MakeP3RadiationBenchmarkDescriptor(
    P3RadiationBenchmarkCase case_kind) {
  P3RadiationBenchmarkDescriptor descriptor;
  descriptor.case_kind = case_kind;
  descriptor.case_id = ToString(case_kind);

  switch (case_kind) {
    case P3RadiationBenchmarkCase::b0_one_cell_source_oracle:
    case P3RadiationBenchmarkCase::b0_one_cell_marshak_source_oracle:
    case P3RadiationBenchmarkCase::b0_zero_flux_multigroup_conservation:
    case P3RadiationBenchmarkCase::b0_marshak_multigroup_leak_closure:
    case P3RadiationBenchmarkCase::b0_group_count_one_regression:
      descriptor.layer = P3RadiationBenchmarkLayer::b0_operator_budget_oracle;
      descriptor.claim_level = "oracle";
      descriptor.group_count =
          case_kind == P3RadiationBenchmarkCase::b0_group_count_one_regression
              ? 1u
              : 2u;
      descriptor.marshak_enabled =
          case_kind == P3RadiationBenchmarkCase::b0_one_cell_marshak_source_oracle ||
          case_kind == P3RadiationBenchmarkCase::b0_marshak_multigroup_leak_closure;
      descriptor.radiation_boundary_model =
          descriptor.marshak_enabled ? "thesis_marshak_vacuum"
                                     : "contract_zero_flux_or_scalar_remap";
      break;
    case P3RadiationBenchmarkCase::b1_woo_5_25_slab_no_limiter:
      descriptor.layer = P3RadiationBenchmarkLayer::b1_thesis_geometry_slab;
      descriptor.claim_level = "thesis_geometry_pre_parity";
      descriptor.group_count = 12u;
      descriptor.marshak_enabled = true;
      descriptor.radiation_boundary_model = "thesis_marshak_vacuum";
      descriptor.no_limiter_pre_parity_slab = true;
      break;
  }

  return descriptor;
}

dec3d::radiation::RadiationGroupLayout MakeP3B1TwelveGroupLayout() {
  const double energy_edges_keV[] = {
      0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0,
      3.0, 10.0, 30.0, 60.0, 120.0, 200.0};
  std::vector<double> frequency_edges_hz;
  frequency_edges_hz.reserve(sizeof(energy_edges_keV) / sizeof(energy_edges_keV[0]));
  for (const double energy_keV : energy_edges_keV) {
    frequency_edges_hz.push_back(FrequencyFromKeV(energy_keV));
  }
  return dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
      frequency_edges_hz);
}

P3RadiationReferencePolicy MakeNoP3RadiationReferencePolicy() {
  return {};
}

P3RadiationBenchmarkValidationResult ValidateP3RadiationBenchmarkDescriptor(
    const P3RadiationBenchmarkDescriptor& descriptor) noexcept {
  P3RadiationBenchmarkValidationResult result;

  if (descriptor.case_id.empty()) {
    result.failure_reason = "case_id is required";
  } else if (!IsFinitePositive(descriptor.r_max_cm) ||
             !IsFinitePositive(descriptor.r0_cm) ||
             descriptor.r0_cm >= descriptor.r_max_cm) {
    result.failure_reason = "invalid slab radii";
  } else if (!IsFinitePositive(descriptor.te_inner_keV) ||
             !IsFinitePositive(descriptor.te_outer_keV) ||
             !IsFinitePositive(descriptor.rho_inner_g_cm3) ||
             !IsFinitePositive(descriptor.rho_outer_g_cm3)) {
    result.failure_reason = "invalid slab initial state";
  } else if (descriptor.initial_ti_policy != "equal_to_initial_te") {
    result.failure_reason = "unsupported initial_ti_policy";
  } else if (descriptor.initial_radiation_policy != "blackbody_from_initial_te") {
    result.failure_reason = "unsupported initial_radiation_policy";
  } else if (descriptor.radiation_flux_limiter != "disabled") {
    result.failure_reason = "P3-B first slice requires disabled radiation limiter";
  } else if (descriptor.parity_claim_allowed || descriptor.lilac_reference_available) {
    result.failure_reason = "P3-B first slice forbids parity claims";
  } else if (descriptor.production_runtime_registration) {
    result.failure_reason = "P3-B first slice forbids runtime registration";
  } else if (descriptor.advection_enabled || descriptor.radiation_pressure_work_enabled) {
    result.failure_reason = "P3-B first slice requires frozen-fluid radiation";
  } else if (descriptor.layer == P3RadiationBenchmarkLayer::b1_thesis_geometry_slab &&
             (!descriptor.marshak_enabled ||
              descriptor.radiation_boundary_model != "thesis_marshak_vacuum" ||
              descriptor.group_count != 12u ||
              !descriptor.no_limiter_pre_parity_slab)) {
    result.failure_reason =
        "P3-B1 requires 12 groups, Marshak, and no-limiter pre-parity claim";
  }

  result.success = result.failure_reason.empty();
  if (!result.success) {
    result.failure_diagnostics =
        "diagnostic_id=p3.radiation_benchmark.descriptor_failure; failure_reason=" +
        result.failure_reason;
  }
  return result;
}

std::string BuildP3RadiationBenchmarkDiagnosticsLine(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.radiation_benchmark"
      << "; benchmark_id=" << descriptor.case_id
      << "; case_id=" << descriptor.case_id
      << "; benchmark_layer=" << ToString(descriptor.layer)
      << "; claim_level=" << descriptor.claim_level
      << "; parity_claim_allowed=" << (reference_policy.parity_claim_allowed ? "true" : "false")
      << "; lilac_reference_available="
      << (reference_policy.lilac_reference_available ? "true" : "false")
      << "; comparison_mode=" << reference_policy.comparison_mode
      << "; radiation_boundary_model=" << descriptor.radiation_boundary_model
      << "; marshak_enabled=" << (descriptor.marshak_enabled ? "true" : "false")
      << "; group_count=" << descriptor.group_count
      << "; group_edges_source=" << descriptor.group_edges_source
      << "; opacity_provider=" << descriptor.opacity_provider
      << "; opacity_data_source=" << descriptor.opacity_data_source
      << "; opacity_energy_coordinate_source=" << descriptor.opacity_energy_coordinate_source
      << "; lookup_energy_mapping_mode=" << descriptor.lookup_energy_mapping_mode
      << "; opacity_source_thesis_exact_match="
      << (descriptor.opacity_source_thesis_exact_match ? "true" : "false")
      << "; initial_ti_policy=" << descriptor.initial_ti_policy
      << "; initial_radiation_policy=" << descriptor.initial_radiation_policy
      << "; radiation_flux_limiter=" << descriptor.radiation_flux_limiter
      << "; advection_enabled=" << (descriptor.advection_enabled ? "true" : "false")
      << "; radiation_pressure_work_enabled="
      << (descriptor.radiation_pressure_work_enabled ? "true" : "false")
      << "; production_runtime_registration="
      << (descriptor.production_runtime_registration ? "true" : "false")
      << "; no_limiter_pre_parity_slab="
      << (descriptor.no_limiter_pre_parity_slab ? "true" : "false")
      << "; missing_diagnostics_fail=true";
  return out.str();
}

bool ValidateP3RadiationBenchmarkDiagnostics(const std::string& report_line) noexcept {
  return HasToken(report_line, "diagnostic_id=p3.radiation_benchmark") &&
         HasToken(report_line, "benchmark_layer=") &&
         HasToken(report_line, "claim_level=") &&
         HasToken(report_line, "parity_claim_allowed=false") &&
         HasToken(report_line, "lilac_reference_available=false") &&
         HasToken(report_line, "radiation_boundary_model=") &&
         HasToken(report_line, "marshak_enabled=") &&
         HasToken(report_line, "group_count=") &&
         HasToken(report_line, "group_edges_source=RadiationGroupLayout") &&
         HasToken(report_line, "opacity_provider=") &&
         HasToken(report_line, "opacity_data_source=") &&
         HasToken(report_line, "opacity_source_thesis_exact_match=false") &&
         HasToken(report_line, "initial_ti_policy=") &&
         HasToken(report_line, "initial_radiation_policy=") &&
         HasToken(report_line, "radiation_flux_limiter=disabled") &&
         HasToken(report_line, "advection_enabled=false") &&
         HasToken(report_line, "radiation_pressure_work_enabled=false") &&
         HasToken(report_line, "production_runtime_registration=false") &&
         HasToken(report_line, "missing_diagnostics_fail=true") &&
         !HasToken(report_line, "parity_claim_allowed=true");
}

bool P3RadiationBenchmarkArtifactManifest::ContainsRequiredFile(
    const std::string& file_name) const noexcept {
  return std::find(required_files.begin(), required_files.end(), file_name) !=
         required_files.end();
}

P3RadiationBenchmarkArtifactManifest BuildP3RadiationBenchmarkArtifactManifest(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const std::string& output_root) {
  P3RadiationBenchmarkArtifactManifest manifest;
  manifest.output_directory = output_root + "/" + descriptor.case_id;
  manifest.required_files = {
      "metadata.json",
      "budget_summary.txt",
      "diagnostics.txt"};
  if (descriptor.layer == P3RadiationBenchmarkLayer::b1_thesis_geometry_slab) {
    manifest.required_files.push_back("radial_profiles_step000001.csv");
    manifest.required_files.push_back("group_energy_budget.csv");
    manifest.required_files.push_back("Te_profile.png");
    manifest.required_files.push_back("sum_Ug_profile.png");
  }
  return manifest;
}

}  // namespace dec3d::benchmarks
