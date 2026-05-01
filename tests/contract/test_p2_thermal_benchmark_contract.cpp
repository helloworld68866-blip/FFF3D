#include "benchmarks/p2_thermal_benchmark.hpp"
#include "test_assert.hpp"

#include <string>

int main() {
  using namespace dec3d::benchmarks;

  const auto spitzer = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::spitzer_pure_thermal);
  DEC3D_CHECK_EQ(spitzer.case_id, std::string{"p2_8b_spitzer_pure_thermal"});
  DEC3D_CHECK_EQ(spitzer.woo_figure, std::string{"5.32"});
  DEC3D_CHECK((
      spitzer.boundary_model ==
      BenchmarkBoundaryModel::thesis_1d_spherical_zero_flux));
  DEC3D_CHECK((
      spitzer.thermal_kappa_model ==
      BenchmarkThermalKappaModel::spitzer_no_degeneracy));
  DEC3D_CHECK(!spitzer.equilibration_stage_enabled);
  DEC3D_CHECK(!spitzer.momentum_update_enabled);
  DEC3D_CHECK(!spitzer.output_times_s.empty());
  DEC3D_CHECK(spitzer.comparison_time_s > 0.0);

  const auto lee_more = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::lee_more_pure_thermal);
  DEC3D_CHECK_EQ(lee_more.case_id, std::string{"p2_8c_lee_more_pure_thermal"});
  DEC3D_CHECK_EQ(lee_more.woo_figure, std::string{"5.33"});
  DEC3D_CHECK((
      lee_more.thermal_kappa_model ==
      BenchmarkThermalKappaModel::lee_more_with_degeneracy));
  DEC3D_CHECK(!lee_more.equilibration_stage_enabled);
  DEC3D_CHECK(!lee_more.momentum_update_enabled);

  const auto thermal_ei = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::lee_more_thermal_plus_ei);
  DEC3D_CHECK_EQ(
      thermal_ei.case_id,
      std::string{"p2_8d_lee_more_thermal_plus_ei"});
  DEC3D_CHECK_EQ(thermal_ei.woo_figure, std::string{"5.34"});
  DEC3D_CHECK(thermal_ei.equilibration_stage_enabled);
  DEC3D_CHECK(!thermal_ei.momentum_update_enabled);
  DEC3D_CHECK((
      thermal_ei.tau_model ==
      BenchmarkTauModel::thesis_spitzer_eq_5_241));

  const auto exact = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::exact_l1m1_constant_kappa);
  DEC3D_CHECK_EQ(
      exact.case_id,
      std::string{"p2_8e_exact_l1m1_constant_kappa"});
  DEC3D_CHECK_EQ(
      exact.woo_figure,
      std::string{"exact_l1m1_spherical_bessel_neumann"});
  DEC3D_CHECK((
      exact.boundary_model ==
      BenchmarkBoundaryModel::full_sphere_scalar_remap));
  DEC3D_CHECK((
      exact.thermal_kappa_model ==
      BenchmarkThermalKappaModel::constant_kappa_exact));
  DEC3D_CHECK(!exact.equilibration_stage_enabled);
  DEC3D_CHECK(!exact.momentum_update_enabled);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDescriptor(exact).success);

  const auto exact_radial = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::exact_l0_radial_constant_kappa);
  DEC3D_CHECK_EQ(
      exact_radial.case_id,
      std::string{"p2_8f_exact_l0_radial_constant_kappa"});
  DEC3D_CHECK_EQ(
      exact_radial.woo_figure,
      std::string{"exact_l0_spherical_bessel_neumann"});
  DEC3D_CHECK((
      exact_radial.boundary_model == BenchmarkBoundaryModel::full_sphere_scalar_remap));
  DEC3D_CHECK((
      exact_radial.thermal_kappa_model == BenchmarkThermalKappaModel::constant_kappa_exact));
  DEC3D_CHECK(!exact_radial.equilibration_stage_enabled);
  DEC3D_CHECK(!exact_radial.momentum_update_enabled);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDescriptor(exact_radial).success);

  const auto limiter = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::flux_limiter_stress);
  DEC3D_CHECK_EQ(limiter.case_id, std::string{"p2_8g_flux_limiter_stress"});
  DEC3D_CHECK_EQ(limiter.woo_figure, std::string{"operator_flux_limiter_stress"});
  DEC3D_CHECK(!limiter.equilibration_stage_enabled);
  DEC3D_CHECK(!limiter.momentum_update_enabled);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDescriptor(limiter).success);

  const auto ei0d = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::ei_0d_thesis_spitzer);
  DEC3D_CHECK_EQ(ei0d.case_id, std::string{"p2_8h_ei_0d_thesis_spitzer"});
  DEC3D_CHECK_EQ(ei0d.woo_figure, std::string{"operator_ei_exact_exponential"});
  DEC3D_CHECK(!ei0d.thermal_stage_enabled);
  DEC3D_CHECK(ei0d.equilibration_stage_enabled);
  DEC3D_CHECK(!ei0d.momentum_update_enabled);
  DEC3D_CHECK_EQ(ei0d.initial_ti_policy, std::string{"ion_uniform_0p5keV"});
  DEC3D_CHECK((ei0d.tau_model == BenchmarkTauModel::thesis_spitzer_eq_5_241));
  DEC3D_CHECK(ValidateWooThermalBenchmarkDescriptor(ei0d).success);

  const auto hte = MakeWooThermalBenchmarkDescriptor(
      WooThermalBenchmarkCase::distributed_hte_smoke);
  DEC3D_CHECK_EQ(hte.case_id, std::string{"p2_8i_distributed_hte_smoke"});
  DEC3D_CHECK_EQ(hte.woo_figure, std::string{"operator_distributed_hte_smoke"});
  DEC3D_CHECK(hte.thermal_stage_enabled);
  DEC3D_CHECK(hte.equilibration_stage_enabled);
  DEC3D_CHECK(!hte.momentum_update_enabled);
  DEC3D_CHECK((hte.tau_model == BenchmarkTauModel::thesis_spitzer_eq_5_241));
  DEC3D_CHECK(ValidateWooThermalBenchmarkDescriptor(hte).success);

  const auto policy = MakeNoReferencePolicy();
  DEC3D_CHECK(!policy.lilac_reference_available);
  DEC3D_CHECK(!policy.parity_claim_allowed);
  DEC3D_CHECK_EQ(policy.comparison_mode, std::string{"no_external_reference"});

  const auto report = BuildWooThermalBenchmarkDiagnosticsLine(spitzer, policy);
  DEC3D_CHECK(report.find("diagnostic_id=p2.thermal_benchmark") != std::string::npos);
  DEC3D_CHECK(report.find("case_id=p2_8b_spitzer_pure_thermal") != std::string::npos);
  DEC3D_CHECK(report.find("boundary_model=thesis_1d_spherical_zero_flux") !=
              std::string::npos);
  DEC3D_CHECK(
      report.find("benchmark_compare_time_policy=explicit_output_times_from_case_descriptor") !=
      std::string::npos);
  DEC3D_CHECK(report.find("lilac_reference_available=false") != std::string::npos);
  DEC3D_CHECK(report.find("parity_claim_allowed=false") != std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(report));

  auto bad = spitzer;
  bad.boundary_model = BenchmarkBoundaryModel::full_sphere_scalar_remap;
  const auto bad_report = BuildWooThermalBenchmarkDiagnosticsLine(bad, policy);
  DEC3D_CHECK(!ValidateWooThermalBenchmarkDiagnostics(bad_report));

  const auto exact_report = BuildWooThermalBenchmarkDiagnosticsLine(exact, policy);
  DEC3D_CHECK(
      exact_report.find("exact_solution=l1m1_spherical_bessel_neumann") !=
      std::string::npos);
  DEC3D_CHECK(
      exact_report.find("exact_projection=cell_volume_average") !=
      std::string::npos);
  DEC3D_CHECK(
      exact_report.find("boundary_model=full_sphere_scalar_remap") !=
      std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(exact_report));

  const auto exact_radial_report =
      BuildWooThermalBenchmarkDiagnosticsLine(exact_radial, policy);
  DEC3D_CHECK(
      exact_radial_report.find("exact_solution=l0_spherical_bessel_neumann") !=
      std::string::npos);
  DEC3D_CHECK(
      exact_radial_report.find("exact_projection=cell_volume_average") !=
      std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(exact_radial_report));

  const auto limiter_report = BuildWooThermalBenchmarkDiagnosticsLine(limiter, policy);
  DEC3D_CHECK(
      limiter_report.find("benchmark_family=operator_flux_limiter_stress") !=
      std::string::npos);
  DEC3D_CHECK(
      limiter_report.find("electron_flux_limiter_required=true") !=
      std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(limiter_report));

  const auto ei0d_report = BuildWooThermalBenchmarkDiagnosticsLine(ei0d, policy);
  DEC3D_CHECK(
      ei0d_report.find("benchmark_family=operator_ei_0d_exact_exponential") !=
      std::string::npos);
  DEC3D_CHECK(
      ei0d_report.find("thermal_stage_executed=false") != std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(ei0d_report));

  const auto hte_report = BuildWooThermalBenchmarkDiagnosticsLine(hte, policy);
  DEC3D_CHECK(
      hte_report.find("benchmark_family=operator_distributed_hte_smoke") !=
      std::string::npos);
  DEC3D_CHECK(hte_report.find("stage_order=H,T,E") != std::string::npos);
  DEC3D_CHECK(ValidateWooThermalBenchmarkDiagnostics(hte_report));

  auto missing_time = spitzer;
  missing_time.output_times_s.clear();
  DEC3D_CHECK(!ValidateWooThermalBenchmarkDescriptor(missing_time).success);

  const auto manifest = BuildWooThermalArtifactManifest(
      spitzer,
      "F:/dec3d/analysis/output/p2_thermal_benchmarks");
  DEC3D_CHECK(manifest.output_directory.find("p2_8b_spitzer_pure_thermal") !=
              std::string::npos);
  DEC3D_CHECK(manifest.required_files.size() >= 8u);
  DEC3D_CHECK(manifest.ContainsRequiredFile("dec3d.out"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("benchmark_summary.md"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("benchmark_diagnostics.txt"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("energy_budget_vs_time.txt"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("comparison_metrics.txt"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("profile_step000001.txt"));

  return 0;
}
