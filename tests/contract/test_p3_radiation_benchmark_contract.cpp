#include "benchmarks/p3_radiation_benchmark.hpp"
#include "benchmarks/p3_radiation_benchmark_runner.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

bool HasToken(const std::string& text, const char* token) {
  return text.find(token) != std::string::npos;
}

}  // namespace

int main() {
  using namespace dec3d::benchmarks;

  const auto source = MakeP3RadiationBenchmarkDescriptor(
      P3RadiationBenchmarkCase::b0_one_cell_source_oracle);
  DEC3D_CHECK_EQ(source.case_id, std::string{"p3_b0_one_cell_source_oracle"});
  DEC3D_CHECK(source.layer == P3RadiationBenchmarkLayer::b0_operator_budget_oracle);
  DEC3D_CHECK_EQ(source.claim_level, std::string{"oracle"});
  DEC3D_CHECK(!source.marshak_enabled);
  DEC3D_CHECK(!source.production_runtime_registration);

  const auto slab = MakeP3RadiationBenchmarkDescriptor(
      P3RadiationBenchmarkCase::b1_woo_5_25_slab_no_limiter);
  DEC3D_CHECK_EQ(slab.case_id, std::string{"p3_b1_woo_5_25_slab_no_limiter"});
  DEC3D_CHECK_EQ(slab.r_max_cm, 1.0e-2);
  DEC3D_CHECK_EQ(slab.r0_cm, 5.0e-3);
  DEC3D_CHECK_EQ(slab.te_inner_keV, 5.0);
  DEC3D_CHECK_EQ(slab.te_outer_keV, 0.5);
  DEC3D_CHECK_EQ(slab.rho_inner_g_cm3, 50.0);
  DEC3D_CHECK_EQ(slab.rho_outer_g_cm3, 100.0);
  DEC3D_CHECK_EQ(slab.initial_ti_policy, std::string{"equal_to_initial_te"});
  DEC3D_CHECK_EQ(slab.initial_radiation_policy,
                 std::string{"blackbody_from_initial_te"});
  DEC3D_CHECK_EQ(slab.group_count, std::size_t{12});
  DEC3D_CHECK(slab.marshak_enabled);
  DEC3D_CHECK_EQ(slab.radiation_boundary_model,
                 std::string{"thesis_marshak_vacuum"});
  DEC3D_CHECK_EQ(slab.radiation_flux_limiter, std::string{"disabled"});
  DEC3D_CHECK(!slab.parity_claim_allowed);
  DEC3D_CHECK(!slab.lilac_reference_available);
  DEC3D_CHECK(!slab.advection_enabled);
  DEC3D_CHECK(!slab.radiation_pressure_work_enabled);
  DEC3D_CHECK(!slab.production_runtime_registration);
  DEC3D_CHECK(ValidateP3RadiationBenchmarkDescriptor(slab).success);

  const auto layout = MakeP3B1TwelveGroupLayout();
  DEC3D_CHECK_EQ(layout.mode,
                 dec3d::radiation::RadiationGroupMode::explicit_frequency_groups);
  DEC3D_CHECK_EQ(layout.group_count, std::size_t{12});
  DEC3D_CHECK_EQ(layout.frequency_edges_hz.size(), std::size_t{13});
  DEC3D_CHECK(dec3d::radiation::ValidateRadiationGroupLayout(layout).success);

  const auto policy = MakeNoP3RadiationReferencePolicy();
  DEC3D_CHECK(!policy.lilac_reference_available);
  DEC3D_CHECK(!policy.parity_claim_allowed);
  DEC3D_CHECK_EQ(policy.comparison_mode, std::string{"no_external_reference"});

  const auto report = BuildP3RadiationBenchmarkDiagnosticsLine(slab, policy);
  DEC3D_CHECK(HasToken(report, "diagnostic_id=p3.radiation_benchmark"));
  DEC3D_CHECK(HasToken(report, "benchmark_layer=P3-B1"));
  DEC3D_CHECK(HasToken(report, "claim_level=thesis_geometry_pre_parity"));
  DEC3D_CHECK(HasToken(report, "parity_claim_allowed=false"));
  DEC3D_CHECK(HasToken(report, "lilac_reference_available=false"));
  DEC3D_CHECK(HasToken(report, "radiation_boundary_model=thesis_marshak_vacuum"));
  DEC3D_CHECK(HasToken(report, "marshak_enabled=true"));
  DEC3D_CHECK(HasToken(report, "group_count=12"));
  DEC3D_CHECK(HasToken(report, "group_edges_source=RadiationGroupLayout"));
  DEC3D_CHECK(HasToken(report, "opacity_energy_coordinate_source=TOPS photon grid"));
  DEC3D_CHECK(HasToken(report, "opacity_source_thesis_exact_match=false"));
  DEC3D_CHECK(HasToken(report, "initial_ti_policy=equal_to_initial_te"));
  DEC3D_CHECK(HasToken(report, "initial_radiation_policy=blackbody_from_initial_te"));
  DEC3D_CHECK(HasToken(report, "radiation_flux_limiter=disabled"));
  DEC3D_CHECK(HasToken(report, "production_runtime_registration=false"));
  DEC3D_CHECK(ValidateP3RadiationBenchmarkDiagnostics(report));

  auto bad = slab;
  bad.marshak_enabled = false;
  bad.radiation_boundary_model = "contract_zero_flux_or_scalar_remap";
  DEC3D_CHECK(!ValidateP3RadiationBenchmarkDescriptor(bad).success);

  auto bad_groups = slab;
  bad_groups.group_count = 11;
  DEC3D_CHECK(!ValidateP3RadiationBenchmarkDescriptor(bad_groups).success);

  auto bad_claim = slab;
  bad_claim.parity_claim_allowed = true;
  DEC3D_CHECK(!ValidateP3RadiationBenchmarkDescriptor(bad_claim).success);

  const auto manifest = BuildP3RadiationBenchmarkArtifactManifest(
      slab,
      "F:/dec3d/analysis/output/p3_radiation_benchmarks");
  DEC3D_CHECK(manifest.ContainsRequiredFile("metadata.json"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("budget_summary.txt"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("diagnostics.txt"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("radial_profiles_step000001.csv"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("group_energy_budget.csv"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("Te_profile.png"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("sum_Ug_profile.png"));

  {
    P3RadiationBenchmarkRunOptions options;
    options.write_artifacts = false;

    const auto source_result = RunP3RadiationBenchmarkCase(
        MakeP3RadiationBenchmarkDescriptor(
            P3RadiationBenchmarkCase::b0_one_cell_source_oracle),
        MakeNoP3RadiationReferencePolicy(),
        options);
    DEC3D_CHECK(source_result.success);
    DEC3D_CHECK(source_result.operator_oracle_executed);
    DEC3D_CHECK(source_result.budget_closed);
    DEC3D_CHECK(source_result.report_line.find("benchmark_layer=P3-B0") !=
                std::string::npos);
    DEC3D_CHECK(source_result.report_line.find("source_oracle_matched=true") !=
                std::string::npos);

    const auto marshak_result = RunP3RadiationBenchmarkCase(
        MakeP3RadiationBenchmarkDescriptor(
            P3RadiationBenchmarkCase::b0_one_cell_marshak_source_oracle),
        MakeNoP3RadiationReferencePolicy(),
        options);
    DEC3D_CHECK(marshak_result.success);
    DEC3D_CHECK(marshak_result.boundary_leak_total > 0.0);
    DEC3D_CHECK(marshak_result.budget_closed);
    DEC3D_CHECK(marshak_result.report_line.find("marshak_enabled=true") !=
                std::string::npos);

    const auto zero_flux_result = RunP3RadiationBenchmarkCase(
        MakeP3RadiationBenchmarkDescriptor(
            P3RadiationBenchmarkCase::b0_zero_flux_multigroup_conservation),
        MakeNoP3RadiationReferencePolicy(),
        options);
    DEC3D_CHECK(zero_flux_result.success);
    DEC3D_CHECK(std::abs(zero_flux_result.boundary_leak_total) < 1.0e-12);
    DEC3D_CHECK(std::abs(zero_flux_result.global_budget_residual) < 1.0e-10);

    const auto group_one_result = RunP3RadiationBenchmarkCase(
        MakeP3RadiationBenchmarkDescriptor(
            P3RadiationBenchmarkCase::b0_group_count_one_regression),
        MakeNoP3RadiationReferencePolicy(),
        options);
    DEC3D_CHECK(group_one_result.success);
    DEC3D_CHECK(group_one_result.group_count_one_regression_passed);
  }

  {
    P3RadiationBenchmarkRunOptions options;
    options.radial_cells = 8;
    options.theta_cells = 1;
    options.phi_cells = 1;
    options.dt_s = 1.0e-18;
    options.output_root = "F:/dec3d/analysis/output/p3_radiation_benchmark_contract";
    options.write_artifacts = true;

    const auto slab_result = RunP3RadiationBenchmarkCase(
        MakeP3RadiationBenchmarkDescriptor(
            P3RadiationBenchmarkCase::b1_woo_5_25_slab_no_limiter),
        MakeNoP3RadiationReferencePolicy(),
        options);
    if (!slab_result.success) {
      std::cerr << slab_result.failure_reason << '\n'
                << slab_result.failure_diagnostics << '\n';
    }
    DEC3D_CHECK(slab_result.success);
    DEC3D_CHECK(slab_result.slab_benchmark_executed);
    DEC3D_CHECK(slab_result.budget_closed);
    DEC3D_CHECK(slab_result.boundary_leak_total >= 0.0);
    DEC3D_CHECK(slab_result.profile_artifacts_written);
    DEC3D_CHECK(slab_result.budget_artifacts_written);
    DEC3D_CHECK(slab_result.metadata_written);
    DEC3D_CHECK(slab_result.diagnostics_written);
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/metadata.json"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/budget_summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/diagnostics.txt"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/radial_profiles_step000001.csv"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/group_energy_budget.csv"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/Te_profile.png"));
    DEC3D_CHECK(std::filesystem::exists(
        slab_result.output_directory + "/sum_Ug_profile.png"));
    DEC3D_CHECK(HasToken(slab_result.report_line, "benchmark_layer=P3-B1"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "case_id=p3_b1_woo_5_25_slab_no_limiter"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "radiation_boundary_model=thesis_marshak_vacuum"));
    DEC3D_CHECK(HasToken(slab_result.report_line, "marshak_enabled=true"));
    DEC3D_CHECK(HasToken(slab_result.report_line, "group_count=12"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "provider_report_present=true"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "radiation_report_present=true"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "opacity_source_thesis_exact_match=false"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "radiation_flux_limiter=disabled"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "parity_claim_allowed=false"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "profile_artifacts_written=true"));
    DEC3D_CHECK(HasToken(slab_result.report_line,
                         "energy_budget_written=true"));
  }

  return 0;
}
