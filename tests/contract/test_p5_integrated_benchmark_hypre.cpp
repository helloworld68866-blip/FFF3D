#include "benchmarks/p5_integrated_benchmark.hpp"
#include "benchmarks/p5_integrated_benchmark_runner_hypre.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <cstdlib>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool HasToken(const std::string& text, const char* token) {
  return text.find(token) != std::string::npos;
}

std::string PlotPythonExecutable() {
#if defined(_MSC_VER)
  char* configured = nullptr;
  std::size_t configured_size = 0u;
  if (_dupenv_s(&configured, &configured_size, "DEC3D_PLOT_PYTHON") == 0 &&
      configured != nullptr) {
    std::string value{configured};
    std::free(configured);
    if (!value.empty()) {
      return value;
    }
  }
#else
  if (const char* configured = std::getenv("DEC3D_PLOT_PYTHON")) {
    return configured;
  }
#endif
  return "C:/Users/Administrator/anaconda3/envs/spyder55-pip/python.exe";
}

struct P5ProfileDensityEvidence {
  std::size_t unique_step0_density_count{0u};
  double min_step0_density{0.0};
  double max_step0_density{0.0};
};

P5ProfileDensityEvidence ReadP5ProfileDensityEvidence(
    const std::string& csv_path) {
  std::ifstream in(csv_path);
  DEC3D_CHECK(static_cast<bool>(in));
  std::string line;
  std::getline(in, line);
  std::vector<double> densities;
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> cells;
    while (std::getline(ss, cell, ',')) {
      cells.push_back(cell);
    }
    if (cells.size() < 6u || cells[2] != "0") {
      continue;
    }
    densities.push_back(std::stod(cells[5]));
  }
  DEC3D_CHECK(!densities.empty());
  std::sort(densities.begin(), densities.end());
  P5ProfileDensityEvidence evidence;
  evidence.min_step0_density = densities.front();
  evidence.max_step0_density = densities.back();
  double previous = densities.front();
  evidence.unique_step0_density_count = 1u;
  for (double value : densities) {
    if (std::abs(value - previous) > 1.0e-6) {
      ++evidence.unique_step0_density_count;
      previous = value;
    }
  }
  return evidence;
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  try {
    dec3d::benchmarks::P5IntegratedBenchmarkRunOptions options;
    options.output_root =
        "F:/dec3d/analysis/output/p5_integrated_benchmark_hypre_contract";
    options.radial_cells = 64u;
    options.theta_cells = 12u;
    options.phi_cells = 16u;
    options.step_count = 3u;
    options.dt_s = 1.0e-14;
    options.write_artifacts = true;
    options.generate_frame_csv = true;

    const auto primary = dec3d::benchmarks::RunP5IntegratedBenchmarkVariant(
        MPI_COMM_WORLD,
        dec3d::benchmarks::MakeP5IntegratedBenchmarkDescriptor(
            dec3d::benchmarks::P5IntegratedBenchmarkCase::legendre_l2_m0,
            dec3d::benchmarks::P5IntegratedBenchmarkVariant::primary_mpi24_ale),
        options);

    if (!primary.success && rank == 0) {
      std::cerr << primary.failure_reason << '\n'
                << primary.failure_diagnostics << '\n';
    }
    DEC3D_CHECK(primary.success);
    DEC3D_CHECK(primary.production_chain_executed);
    DEC3D_CHECK(primary.all_stage_reports_present);
    DEC3D_CHECK(primary.profile_artifacts_written);
    DEC3D_CHECK(primary.mode_artifacts_written);
    DEC3D_CHECK(primary.budget_artifacts_written);
    DEC3D_CHECK(HasToken(primary.report_line, "stage_order=H,T,E,R,A"));
    DEC3D_CHECK(HasToken(primary.report_line, "mode_l=2"));
    DEC3D_CHECK(HasToken(primary.report_line, "mode_m=0"));
    DEC3D_CHECK(primary.target_mode_amplitude_initial != 0.0);
    DEC3D_CHECK(primary.target_mode_amplitude_final ==
                primary.target_mode_amplitude_final);
    DEC3D_CHECK(primary.phi_leakage_final <= 1.0e-6);

    const auto comparison =
        dec3d::benchmarks::RunP5IntegratedBenchmarkComparisonSet(
            MPI_COMM_WORLD,
            dec3d::benchmarks::P5IntegratedBenchmarkCase::legendre_l2_m0,
            options);
    DEC3D_CHECK(comparison.success);
    DEC3D_CHECK(comparison.primary_present);
    DEC3D_CHECK(comparison.mpi12_comparison_present);
    DEC3D_CHECK(comparison.no_ale_comparison_present);
    DEC3D_CHECK(comparison.mpi12_vs_mpi24_profile_l1_relative <= 5.0e-3);
    DEC3D_CHECK(comparison.mpi12_vs_mpi24_mode_amplitude_relative <= 1.0e-2);
    DEC3D_CHECK(HasToken(
        comparison.report_line,
        "variant_matrix=primary,no_ale_comparison,mpi12_comparison"));

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_case_descriptor.json"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_stage_reports.json"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_radial_profiles.csv"));
      const auto density_evidence = ReadP5ProfileDensityEvidence(
          options.output_root + "/p5_radial_profiles.csv");
      DEC3D_CHECK(density_evidence.unique_step0_density_count > 8u);
      DEC3D_CHECK(density_evidence.min_step0_density < 10.0);
      DEC3D_CHECK(density_evidence.max_step0_density > 250.0);
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_mode_amplitudes.csv"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_budget_history.csv"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p5_variant_comparison_summary.json"));

      const auto python = PlotPythonExecutable();
      if (std::filesystem::exists(python)) {
        const std::string command =
            "\"" + python + "\" "
            "F:/dec3d/analysis/scripts/plot_p5_integrated_benchmark.py "
            "--input-dir " +
            options.output_root;
        DEC3D_CHECK(std::system(command.c_str()) == 0);
        DEC3D_CHECK(std::filesystem::exists(
            options.output_root + "/p5_radial_profiles.png"));
        DEC3D_CHECK(std::filesystem::exists(
            options.output_root + "/p5_mode_amplitudes.png"));
        DEC3D_CHECK(std::filesystem::exists(
            options.output_root + "/p5_budget_residuals.png"));
        DEC3D_CHECK(std::filesystem::exists(
            options.output_root +
            "/p5_p5_legendre_l2_m0_primary_rz_evolution.gif"));
      } else {
        std::cerr << "Skipping P5 plot generation because Python was not found: "
                  << python << '\n';
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
