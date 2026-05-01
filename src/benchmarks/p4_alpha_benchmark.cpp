#include "benchmarks/p4_alpha_benchmark.hpp"

#include <algorithm>
#include <sstream>

namespace dec3d::benchmarks {
namespace {

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] P4AlphaBenchmarkValidationResult Fail(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const char* reason) {
  P4AlphaBenchmarkValidationResult result;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p4.alpha.benchmark_ladder.failure"
      << "; benchmark_id=" << descriptor.benchmark_id
      << "; case_id=" << descriptor.case_id
      << "; failure_reason=" << reason;
  result.failure_diagnostics = out.str();
  return result;
}

}  // namespace

bool P4AlphaBenchmarkArtifactManifest::ContainsRequiredFile(
    const std::string& file_name) const noexcept {
  return std::find(required_files.begin(), required_files.end(), file_name) !=
         required_files.end();
}

const char* ToString(P4AlphaBenchmarkLayer layer) noexcept {
  switch (layer) {
    case P4AlphaBenchmarkLayer::b0_operator_budget_oracle:
      return "P4-B0";
    case P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot:
      return "P4-B1";
    case P4AlphaBenchmarkLayer::b3_distributed_performance:
      return "P4-B3";
  }
  return "unknown";
}

const char* ToString(P4AlphaBenchmarkCase case_kind) noexcept {
  switch (case_kind) {
    case P4AlphaBenchmarkCase::b0_one_cell_birth_drag_oracle:
      return "p4_b0_one_cell_birth_drag_oracle";
    case P4AlphaBenchmarkCase::b0_one_cell_dt_zero_no_change:
      return "p4_b0_one_cell_dt_zero_no_change";
    case P4AlphaBenchmarkCase::b0_two_cell_diffusion_coupling_oracle:
      return "p4_b0_two_cell_diffusion_coupling_oracle";
    case P4AlphaBenchmarkCase::b0_zero_flux_budget_closure:
      return "p4_b0_zero_flux_budget_closure";
    case P4AlphaBenchmarkCase::b0_distributed_matches_single_rank:
      return "p4_b0_distributed_matches_single_rank";
    case P4AlphaBenchmarkCase::b1_frozen_fluid_alpha_hotspot:
      return "p4_b1_frozen_fluid_alpha_hotspot";
    case P4AlphaBenchmarkCase::b3_distributed_performance:
      return "p4_b3_distributed_performance";
  }
  return "unknown";
}

P4AlphaBenchmarkDescriptor MakeP4AlphaBenchmarkDescriptor(
    P4AlphaBenchmarkCase case_kind) {
  P4AlphaBenchmarkDescriptor descriptor;
  descriptor.case_kind = case_kind;
  descriptor.case_id = ToString(case_kind);

  switch (case_kind) {
    case P4AlphaBenchmarkCase::b0_one_cell_birth_drag_oracle:
    case P4AlphaBenchmarkCase::b0_one_cell_dt_zero_no_change:
    case P4AlphaBenchmarkCase::b0_two_cell_diffusion_coupling_oracle:
    case P4AlphaBenchmarkCase::b0_zero_flux_budget_closure:
    case P4AlphaBenchmarkCase::b0_distributed_matches_single_rank:
      descriptor.layer = P4AlphaBenchmarkLayer::b0_operator_budget_oracle;
      descriptor.benchmark_id = "P4-B0";
      descriptor.claim_level = "operator_oracle";
      descriptor.nr =
          case_kind == P4AlphaBenchmarkCase::b0_two_cell_diffusion_coupling_oracle
              ? 2u
              : 1u;
      descriptor.output_times_s = {1.0e-18};
      descriptor.comparison_times_s = {1.0e-18};
      break;
    case P4AlphaBenchmarkCase::b1_frozen_fluid_alpha_hotspot:
      descriptor.layer = P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot;
      descriptor.benchmark_id = "P4-B1";
      descriptor.claim_level =
          "code_internal_serial_distributed_cross_validation";
      descriptor.nr = 64u;
      descriptor.output_times_s = {1.0e-18, 2.0e-18, 5.0e-18};
      descriptor.comparison_times_s = {5.0e-18};
      break;
    case P4AlphaBenchmarkCase::b3_distributed_performance:
      descriptor.layer = P4AlphaBenchmarkLayer::b3_distributed_performance;
      descriptor.benchmark_id = "P4-B3";
      descriptor.claim_level = "performance_evidence";
      descriptor.nr = 64u;
      descriptor.output_times_s = {1.0e-18};
      descriptor.comparison_times_s = {1.0e-18};
      break;
  }
  return descriptor;
}

P4AlphaBenchmarkValidationResult ValidateP4AlphaBenchmarkDescriptor(
    const P4AlphaBenchmarkDescriptor& descriptor) noexcept {
  if (descriptor.fuel_depletion_enabled) {
    return Fail(descriptor, "fuel depletion is out of scope for P4-B");
  }
  if (descriptor.separate_dt_species_authoritative) {
    return Fail(descriptor, "separate D/T authoritative species are out of scope");
  }
  if (descriptor.parity_claim_allowed) {
    return Fail(descriptor, "P4-B cannot claim external alpha parity");
  }
  if (descriptor.output_times_s.empty()) {
    return Fail(descriptor, "output_times_s must be explicit");
  }
  if (descriptor.comparison_times_s.empty()) {
    return Fail(descriptor, "comparison_times_s must be explicit");
  }
  if (descriptor.reactivity_model != "bosch_hale_dt") {
    return Fail(descriptor, "P4-B requires Bosch-Hale reactivity path");
  }
  if (descriptor.composition_model != "equimolar_dt_from_p2_recovery") {
    return Fail(descriptor, "P4-B requires equimolar DT composition closure");
  }
  if (descriptor.alpha_transport_model != "atzeni_one_group") {
    return Fail(descriptor, "P4-B requires Atzeni one-group alpha transport");
  }
  P4AlphaBenchmarkValidationResult result;
  result.success = true;
  return result;
}

std::string BuildP4AlphaBenchmarkDiagnosticsLine(
    const P4AlphaBenchmarkDescriptor& descriptor) {
  std::ostringstream out;
  out << "diagnostic_id=p4.alpha.benchmark_ladder"
      << "; benchmark_id=" << descriptor.benchmark_id
      << "; benchmark_layer=" << ToString(descriptor.layer)
      << "; case_id=" << descriptor.case_id
      << "; claim_level=" << descriptor.claim_level
      << "; alpha_transport_model=" << descriptor.alpha_transport_model
      << "; reactivity_model_executed=" << descriptor.reactivity_model
      << "; bosch_hale_temperature_source=Ti_old"
      << "; composition_model=" << descriptor.composition_model
      << "; fuel_depletion_enabled=false"
      << "; separate_dt_species_authoritative=false"
      << "; boundary_model=" << descriptor.boundary_model
      << "; initial_ti_policy=" << descriptor.initial_ti_policy
      << "; epsilon_alpha_initial_policy="
      << descriptor.epsilon_alpha_initial_policy
      << "; parity_claim_allowed=false";
  return out.str();
}

bool ValidateP4AlphaBenchmarkDiagnostics(
    const std::string& report_line) noexcept {
  return Contains(report_line, "diagnostic_id=p4.alpha.benchmark_ladder") &&
         Contains(report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(report_line, "reactivity_model_executed=bosch_hale_dt") &&
         Contains(report_line, "bosch_hale_temperature_source=Ti_old") &&
         Contains(report_line, "composition_model=equimolar_dt_from_p2_recovery") &&
         Contains(report_line, "fuel_depletion_enabled=false") &&
         Contains(report_line, "separate_dt_species_authoritative=false") &&
         Contains(report_line, "parity_claim_allowed=false");
}

P4AlphaBenchmarkArtifactManifest BuildP4AlphaBenchmarkArtifactManifest(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const std::string& output_root) {
  P4AlphaBenchmarkArtifactManifest manifest;
  manifest.output_directory = output_root;
  switch (descriptor.layer) {
    case P4AlphaBenchmarkLayer::b0_operator_budget_oracle:
      manifest.required_files = {"p4_b0_operator_oracles.json"};
      break;
    case P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot:
      manifest.required_files = {
          "p4_b1_alpha_hotspot_case_descriptor.json",
          "p4_b1_serial_profiles.csv",
          "p4_b1_budget_summary.csv",
          "p4_benchmark_manifest.json"};
      break;
    case P4AlphaBenchmarkLayer::b3_distributed_performance:
      manifest.required_files = {"p4_b3_distributed_performance.csv"};
      break;
  }
  return manifest;
}

}  // namespace dec3d::benchmarks
