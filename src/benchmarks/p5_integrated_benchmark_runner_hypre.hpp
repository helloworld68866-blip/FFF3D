#pragma once

#include "benchmarks/p5_integrated_benchmark.hpp"

#include <mpi.h>

#include <cstddef>
#include <string>

namespace dec3d::benchmarks {

struct P5IntegratedBenchmarkRunOptions {
  std::size_t radial_cells{64u};
  std::size_t theta_cells{32u};
  std::size_t phi_cells{64u};
  std::size_t step_count{8u};
  double dt_s{1.0e-14};
  std::string output_root{"F:/dec3d/analysis/output/p5_integrated_benchmark"};
  std::string initial_profile_csv_path{
      "F:/dec3d/analysis/output/p5_initial_profile_from_image_profile_100ps_hydro.csv"};
  std::string initial_profile_yaml_path{
      "F:/python_project/cases/image_profile_100ps_hydro.yaml"};
  bool write_artifacts{true};
  bool generate_frame_csv{true};
};

struct P5IntegratedBenchmarkRunResult {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
  std::string report_line;
  bool production_chain_executed{false};
  bool all_stage_reports_present{false};
  bool profile_artifacts_written{false};
  bool mode_artifacts_written{false};
  bool budget_artifacts_written{false};
  bool timing_artifacts_written{false};
  double phi_leakage_final{0.0};
  double target_mode_amplitude_initial{0.0};
  double target_mode_amplitude_final{0.0};
};

struct P5IntegratedBenchmarkComparisonResult {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
  std::string report_line;
  bool primary_present{false};
  bool mpi12_comparison_present{false};
  bool no_ale_comparison_present{false};
  double mpi12_vs_mpi24_profile_l1_relative{0.0};
  double mpi12_vs_mpi24_profile_linf_relative{0.0};
  double mpi12_vs_mpi24_mode_amplitude_relative{0.0};
  double no_ale_vs_ale_profile_l1_relative{0.0};
  double no_ale_vs_ale_mode_amplitude_relative{0.0};
};

[[nodiscard]] P5IntegratedBenchmarkRunResult RunP5IntegratedBenchmarkVariant(
    MPI_Comm communicator,
    const P5IntegratedBenchmarkDescriptor& descriptor,
    const P5IntegratedBenchmarkRunOptions& options) noexcept;

[[nodiscard]] P5IntegratedBenchmarkComparisonResult
RunP5IntegratedBenchmarkComparisonSet(
    MPI_Comm communicator,
    P5IntegratedBenchmarkCase case_kind,
    const P5IntegratedBenchmarkRunOptions& options) noexcept;

}  // namespace dec3d::benchmarks
