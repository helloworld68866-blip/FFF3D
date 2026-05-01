#include "benchmarks/p2_thermal_benchmark.hpp"
#include "benchmarks/p2_thermal_benchmark_runner.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <exception>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

[[nodiscard]] bool HasToken(const std::string& text, const std::string& token) {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] unsigned long long ExtractUnsignedToken(
    const std::string& text,
    const std::string& key) {
  const auto pos = text.find(key);
  if (pos == std::string::npos) {
    dec3d::test::Fail("missing token", __FILE__, __LINE__, key);
  }
  const auto value_begin = pos + key.size();
  const auto value_end = text.find(';', value_begin);
  return std::stoull(text.substr(value_begin, value_end - value_begin));
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
      dec3d::test::Fail("rank count", __FILE__, __LINE__, "requires mpiexec -n 2");
    }

    for (const auto case_kind : {
             dec3d::benchmarks::WooThermalBenchmarkCase::spitzer_pure_thermal,
             dec3d::benchmarks::WooThermalBenchmarkCase::lee_more_pure_thermal,
             dec3d::benchmarks::WooThermalBenchmarkCase::lee_more_thermal_plus_ei,
             dec3d::benchmarks::WooThermalBenchmarkCase::exact_l1m1_constant_kappa,
             dec3d::benchmarks::WooThermalBenchmarkCase::exact_l0_radial_constant_kappa,
             dec3d::benchmarks::WooThermalBenchmarkCase::flux_limiter_stress,
             dec3d::benchmarks::WooThermalBenchmarkCase::ei_0d_thesis_spitzer,
             dec3d::benchmarks::WooThermalBenchmarkCase::distributed_hte_smoke}) {
      auto descriptor = dec3d::benchmarks::MakeWooThermalBenchmarkDescriptor(case_kind);
      descriptor.output_times_s = {1.0e-14};
      descriptor.comparison_time_s = 1.0e-14;

      auto options = dec3d::benchmarks::WooThermalBenchmarkRunOptions{};
      options.global_radial_cells = 4;
      options.theta_cells = 4;
      options.phi_cells = 4;
      options.output_root = "F:/dec3d/analysis/output/p2_thermal_benchmarks_contract";
      options.write_artifacts = false;

      auto result = dec3d::benchmarks::RunWooThermalBenchmarkCase(
          MPI_COMM_WORLD,
          descriptor,
          dec3d::benchmarks::MakeNoReferencePolicy(),
          options);
      if (!result.success) {
        std::cerr << result.failure_diagnostics << '\n';
      }
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(HasToken(result.report_line, "diagnostic_id=p2.thermal_benchmark"));
      DEC3D_CHECK(HasToken(result.report_line, "profile_artifacts_written=false"));
      DEC3D_CHECK(HasToken(result.report_line, "lilac_reference_available=false"));
      DEC3D_CHECK(HasToken(result.report_line, "parity_claim_allowed=false"));
      if (case_kind ==
          dec3d::benchmarks::WooThermalBenchmarkCase::exact_l1m1_constant_kappa ||
          case_kind ==
          dec3d::benchmarks::WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
        if (case_kind ==
            dec3d::benchmarks::WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
          DEC3D_CHECK(HasToken(result.report_line,
                               "exact_solution=l0_spherical_bessel_neumann"));
        } else {
          DEC3D_CHECK(HasToken(result.report_line,
                               "exact_solution=l1m1_spherical_bessel_neumann"));
        }
        DEC3D_CHECK(HasToken(result.report_line, "boundary_model=full_sphere_scalar_remap"));
        DEC3D_CHECK(HasToken(result.report_line, "te_l1_keV="));
        DEC3D_CHECK(HasToken(result.report_line, "te_linf_keV="));
        DEC3D_CHECK(HasToken(result.report_line, "ti_l1_keV="));
        DEC3D_CHECK(HasToken(result.report_line, "ti_linf_keV="));
        DEC3D_CHECK(HasToken(result.report_line, "origin_remap_used=true"));
        DEC3D_CHECK(HasToken(result.report_line, "pole_remap_used=true"));
      } else if (case_kind ==
                 dec3d::benchmarks::WooThermalBenchmarkCase::flux_limiter_stress) {
        DEC3D_CHECK(HasToken(result.report_line, "electron_flux_limiter_required=true"));
        DEC3D_CHECK(HasToken(result.report_line, "electron_flux_limiter_enabled=true"));
        DEC3D_CHECK(HasToken(result.report_line, "limited_face_count_global="));
        DEC3D_CHECK(ExtractUnsignedToken(result.report_line, "limited_face_count_global=") > 0u);
      } else if (case_kind ==
                 dec3d::benchmarks::WooThermalBenchmarkCase::ei_0d_thesis_spitzer) {
        DEC3D_CHECK(HasToken(result.report_line,
                             "benchmark_family=operator_ei_0d_exact_exponential"));
        DEC3D_CHECK(HasToken(result.report_line, "initial_ti_policy=ion_uniform_0p5keV"));
        DEC3D_CHECK(HasToken(result.report_line, "thermal_stage_executed=false"));
        DEC3D_CHECK(HasToken(result.report_line, "equilibration_stage_report_present=true"));
        DEC3D_CHECK(HasToken(result.report_line, "max_deltaT_before_erg="));
        DEC3D_CHECK(HasToken(result.report_line, "max_dt_over_tau="));
      } else if (case_kind ==
                 dec3d::benchmarks::WooThermalBenchmarkCase::distributed_hte_smoke) {
        DEC3D_CHECK(HasToken(result.report_line, "stage_order=H,T,E"));
        DEC3D_CHECK(HasToken(result.report_line, "thermal_stage_report_present=true"));
        DEC3D_CHECK(HasToken(result.report_line, "equilibration_stage_report_present=true"));
      } else if (
          case_kind != dec3d::benchmarks::WooThermalBenchmarkCase::ei_0d_thesis_spitzer) {
        DEC3D_CHECK(HasToken(result.report_line, "r_max_cm=0.01"));
        DEC3D_CHECK(HasToken(result.report_line, "r0_cm=0.005"));
        DEC3D_CHECK(HasToken(result.report_line, "te_inner_keV=5"));
        DEC3D_CHECK(HasToken(result.report_line, "te_outer_keV=0.5"));
        DEC3D_CHECK(HasToken(result.report_line, "rho_inner_g_cm3=50"));
        DEC3D_CHECK(HasToken(result.report_line, "rho_outer_g_cm3=100"));
        DEC3D_CHECK(HasToken(result.report_line, "initial_ti_policy=equal_to_initial_te"));
      }
      DEC3D_CHECK(!result.momentum_update_executed);
      if (descriptor.equilibration_stage_enabled) {
        DEC3D_CHECK(result.equilibration_stage_executed);
      } else {
        DEC3D_CHECK(!result.equilibration_stage_executed);
      }
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::filesystem::path artifact_root =
        "F:/dec3d/analysis/output/p2_thermal_benchmarks_contract_artifacts";
    if (rank == 0) {
      std::filesystem::remove_all(artifact_root);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    auto artifact_descriptor = dec3d::benchmarks::MakeWooThermalBenchmarkDescriptor(
        dec3d::benchmarks::WooThermalBenchmarkCase::spitzer_pure_thermal);
    artifact_descriptor.output_times_s = {1.0e-14};
    artifact_descriptor.comparison_time_s = 1.0e-14;
    auto artifact_options = dec3d::benchmarks::WooThermalBenchmarkRunOptions{};
    artifact_options.global_radial_cells = 4;
    artifact_options.theta_cells = 4;
    artifact_options.phi_cells = 4;
    artifact_options.output_root = artifact_root.generic_string();
    artifact_options.write_artifacts = true;
    const auto artifact_result = dec3d::benchmarks::RunWooThermalBenchmarkCase(
        MPI_COMM_WORLD,
        artifact_descriptor,
        dec3d::benchmarks::MakeNoReferencePolicy(),
        artifact_options);
    DEC3D_CHECK(artifact_result.success);
    DEC3D_CHECK(artifact_result.profile_artifacts_written);
    DEC3D_CHECK(artifact_result.energy_budget_written);
    DEC3D_CHECK(
        artifact_result.report_line.find("artifact_gather_is_not_production_writeback=true") !=
        std::string::npos);
    if (rank == 0) {
      const std::filesystem::path case_dir =
          artifact_root / "p2_8b_spitzer_pure_thermal";
      DEC3D_CHECK(std::filesystem::exists(case_dir / "dec3d.out"));
      DEC3D_CHECK(std::filesystem::exists(case_dir / "benchmark_summary.md"));
      DEC3D_CHECK(std::filesystem::exists(case_dir / "benchmark_diagnostics.txt"));
      DEC3D_CHECK(std::filesystem::exists(case_dir / "profile_step000001.txt"));
      DEC3D_CHECK(std::filesystem::exists(case_dir / "energy_budget_vs_time.txt"));
      DEC3D_CHECK(std::filesystem::exists(case_dir / "comparison_metrics.txt"));
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
