#include "benchmarks/p4_alpha_benchmark.hpp"
#include "benchmarks/p4_alpha_benchmark_runner_hypre.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>

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

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {

  dec3d::benchmarks::P4AlphaBenchmarkRunOptions options;
  options.radial_cells = 16u;
  options.theta_cells = 1u;
  options.phi_cells = 1u;
  options.dt_s = 1.0e-18;
  options.step_count = 2u;
  options.output_root = "F:/dec3d/analysis/output/p4_alpha_benchmark_hypre_contract";
  options.write_artifacts = true;

  const auto b0 = dec3d::benchmarks::RunP4AlphaBenchmarkSerialDistributedComparison(
      MPI_COMM_WORLD,
      dec3d::benchmarks::MakeP4AlphaBenchmarkDescriptor(
          dec3d::benchmarks::P4AlphaBenchmarkCase::b0_distributed_matches_single_rank),
      options);
  if (!b0.success && rank == 0) {
    std::cerr << b0.failure_reason << '\n' << b0.failure_diagnostics << '\n';
  }
  DEC3D_CHECK(b0.success);
  DEC3D_CHECK(b0.single_rank_run_present);
  DEC3D_CHECK(b0.distributed_run_present);
  DEC3D_CHECK(b0.serial_distributed_cross_validation_present);
  DEC3D_CHECK(std::abs(b0.serial_distributed_epsilon_alpha_linf) < 1.0e-6);
  DEC3D_CHECK(HasToken(b0.report_line, "rank0_gather_solve_used=false"));
  DEC3D_CHECK(HasToken(b0.report_line, "serial_dense_fallback_used=false"));

  const auto b1 = dec3d::benchmarks::RunP4AlphaBenchmarkSerialDistributedComparison(
      MPI_COMM_WORLD,
      dec3d::benchmarks::MakeP4AlphaBenchmarkDescriptor(
          dec3d::benchmarks::P4AlphaBenchmarkCase::b1_frozen_fluid_alpha_hotspot),
      options);
  DEC3D_CHECK(b1.success);
  DEC3D_CHECK(b1.b1_hotspot_executed);
  DEC3D_CHECK(b1.profile_artifacts_written);
  DEC3D_CHECK(b1.serial_distributed_cross_validation_present);
  if (b1.serial_distributed_epsilon_alpha_relative_linf >= 1.0e-8 && rank == 0) {
    std::cerr << "epsilon_alpha_linf="
              << b1.serial_distributed_epsilon_alpha_linf
              << " epsilon_alpha_relative_linf="
              << b1.serial_distributed_epsilon_alpha_relative_linf
              << " Te_linf=" << b1.serial_distributed_Te_linf << '\n'
              << b1.report_line << '\n';
  }
  DEC3D_CHECK(b1.serial_distributed_epsilon_alpha_relative_linf < 1.0e-8);
  DEC3D_CHECK(b1.serial_distributed_Te_linf < 1.0e-10);

  const auto b3 = dec3d::benchmarks::RunP4AlphaBenchmarkCaseDistributed(
      MPI_COMM_WORLD,
      dec3d::benchmarks::MakeP4AlphaBenchmarkDescriptor(
          dec3d::benchmarks::P4AlphaBenchmarkCase::b3_distributed_performance),
      options);
  DEC3D_CHECK(b3.success);
  DEC3D_CHECK(b3.b3_performance_executed);
  DEC3D_CHECK(b3.performance_artifacts_written);

  if (rank == 0) {
    DEC3D_CHECK(std::filesystem::exists(
        options.output_root + "/p4_b1_serial_profiles.csv"));
    DEC3D_CHECK(std::filesystem::exists(
        options.output_root + "/p4_b1_distributed_profiles.csv"));
    DEC3D_CHECK(std::filesystem::exists(
        options.output_root + "/p4_b1_budget_summary.csv"));
    DEC3D_CHECK(std::filesystem::exists(
        options.output_root + "/p4_b1_serial_vs_distributed_summary.json"));
    DEC3D_CHECK(std::filesystem::exists(
        options.output_root + "/p4_b3_distributed_performance.csv"));

  }

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    const auto python = PlotPythonExecutable();
    if (std::filesystem::exists(python)) {
      const std::string command =
          "\"" + python + "\" "
          "F:/dec3d/analysis/scripts/plot_p4_alpha_benchmark.py "
          "--input-dir " +
          options.output_root;
      DEC3D_CHECK(std::system(command.c_str()) == 0);
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p4_b1_profiles.png"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p4_b1_budget_closure.png"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p4_b1_serial_distributed_error.png"));
      DEC3D_CHECK(std::filesystem::exists(
          options.output_root + "/p4_b3_timing_breakdown.png"));
    } else {
      std::cerr << "Skipping P4-B plot generation because DEC3D_PLOT_PYTHON "
                   "is not configured and default Python was not found: "
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
