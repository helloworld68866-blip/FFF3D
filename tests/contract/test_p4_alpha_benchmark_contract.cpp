#include "benchmarks/p4_alpha_benchmark.hpp"
#include "benchmarks/p4_alpha_benchmark_runner.hpp"
#include "test_assert.hpp"

#include <filesystem>
#include <exception>
#include <iostream>
#include <string>

namespace {

bool HasToken(const std::string& text, const char* token) {
  return text.find(token) != std::string::npos;
}

}  // namespace

int RunTests() {
  using namespace dec3d::benchmarks;

  const auto b0 = MakeP4AlphaBenchmarkDescriptor(
      P4AlphaBenchmarkCase::b0_one_cell_birth_drag_oracle);
  DEC3D_CHECK_EQ(b0.benchmark_id, std::string{"P4-B0"});
  DEC3D_CHECK_EQ(b0.case_id, std::string{"p4_b0_one_cell_birth_drag_oracle"});
  DEC3D_CHECK(b0.layer == P4AlphaBenchmarkLayer::b0_operator_budget_oracle);
  DEC3D_CHECK(!b0.fuel_depletion_enabled);
  DEC3D_CHECK(!b0.separate_dt_species_authoritative);
  DEC3D_CHECK(!b0.parity_claim_allowed);
  DEC3D_CHECK(ValidateP4AlphaBenchmarkDescriptor(b0).success);

  const auto b1 = MakeP4AlphaBenchmarkDescriptor(
      P4AlphaBenchmarkCase::b1_frozen_fluid_alpha_hotspot);
  DEC3D_CHECK_EQ(b1.benchmark_id, std::string{"P4-B1"});
  DEC3D_CHECK_EQ(b1.case_id, std::string{"p4_b1_frozen_fluid_alpha_hotspot"});
  DEC3D_CHECK_EQ(b1.geometry, std::string{"1d_spherical_slab"});
  DEC3D_CHECK_EQ(
      b1.boundary_model,
      std::string{"benchmark_1d_spherical_scalar_origin_remap_outer_zero_flux"});
  DEC3D_CHECK_EQ(b1.initial_ti_policy, std::string{"equal_to_initial_te"});
  DEC3D_CHECK_EQ(b1.epsilon_alpha_initial_policy, std::string{"zero"});
  DEC3D_CHECK(ValidateP4AlphaBenchmarkDescriptor(b1).success);

  const auto b3 = MakeP4AlphaBenchmarkDescriptor(
      P4AlphaBenchmarkCase::b3_distributed_performance);
  DEC3D_CHECK_EQ(b3.benchmark_id, std::string{"P4-B3"});
  DEC3D_CHECK(b3.layer == P4AlphaBenchmarkLayer::b3_distributed_performance);
  DEC3D_CHECK(ValidateP4AlphaBenchmarkDescriptor(b3).success);

  const auto report = BuildP4AlphaBenchmarkDiagnosticsLine(b1);
  DEC3D_CHECK(HasToken(report, "diagnostic_id=p4.alpha.benchmark_ladder"));
  DEC3D_CHECK(HasToken(report, "benchmark_id=P4-B1"));
  DEC3D_CHECK(HasToken(report, "alpha_transport_model=atzeni_one_group"));
  DEC3D_CHECK(HasToken(report, "reactivity_model_executed=bosch_hale_dt"));
  DEC3D_CHECK(HasToken(report, "composition_model=equimolar_dt_from_p2_recovery"));
  DEC3D_CHECK(HasToken(report, "fuel_depletion_enabled=false"));
  DEC3D_CHECK(HasToken(report, "separate_dt_species_authoritative=false"));
  DEC3D_CHECK(HasToken(report, "parity_claim_allowed=false"));
  DEC3D_CHECK(ValidateP4AlphaBenchmarkDiagnostics(report));

  auto bad = b1;
  bad.fuel_depletion_enabled = true;
  DEC3D_CHECK(!ValidateP4AlphaBenchmarkDescriptor(bad).success);

  auto missing_times = b1;
  missing_times.output_times_s.clear();
  DEC3D_CHECK(!ValidateP4AlphaBenchmarkDescriptor(missing_times).success);

  P4AlphaBenchmarkRunOptions options;
  options.radial_cells = 8u;
  options.theta_cells = 1u;
  options.phi_cells = 1u;
  options.step_count = 2u;
  options.dt_s = 1.0e-18;
  options.output_root = "F:/dec3d/analysis/output/p4_alpha_benchmark_contract";
  options.write_artifacts = true;

  const auto b0_run = RunP4AlphaBenchmarkCase(b0, options);
  if (!b0_run.success) {
    std::cerr << b0_run.failure_reason << '\n'
              << b0_run.failure_diagnostics << '\n';
  }
  DEC3D_CHECK(b0_run.success);
  DEC3D_CHECK(b0_run.b0_operator_oracle_executed);
  DEC3D_CHECK(b0_run.single_rank_run_present);
  DEC3D_CHECK(HasToken(b0_run.report_line, "budget_closed=true"));
  DEC3D_CHECK(std::filesystem::exists(
      options.output_root + "/p4_b0_operator_oracles.json"));

  const auto b1_run = RunP4AlphaBenchmarkCase(b1, options);
  if (!b1_run.success) {
    std::cerr << b1_run.failure_reason << '\n'
              << b1_run.failure_diagnostics << '\n';
  }
  DEC3D_CHECK(b1_run.success);
  DEC3D_CHECK(b1_run.b1_hotspot_executed);
  DEC3D_CHECK(b1_run.profile_artifacts_written);
  DEC3D_CHECK(b1_run.budget_artifacts_written);
  DEC3D_CHECK(std::filesystem::exists(
      options.output_root + "/p4_b1_alpha_hotspot_case_descriptor.json"));
  DEC3D_CHECK(std::filesystem::exists(
      options.output_root + "/p4_b1_serial_profiles.csv"));
  DEC3D_CHECK(std::filesystem::exists(
      options.output_root + "/p4_benchmark_manifest.json"));

  return 0;
}

int main() {
  try {
    return RunTests();
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
