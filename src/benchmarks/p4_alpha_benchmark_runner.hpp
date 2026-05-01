#pragma once

#include "benchmarks/p4_alpha_benchmark.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::benchmarks {

struct P4AlphaBenchmarkRunOptions {
  std::size_t radial_cells{16u};
  std::size_t theta_cells{1u};
  std::size_t phi_cells{1u};
  double dt_s{1.0e-18};
  std::size_t step_count{1u};
  std::string output_root{"analysis/output/p4_alpha_benchmarks"};
  bool write_artifacts{true};
};

struct P4AlphaBenchmarkRunResult {
  bool success{false};
  bool b0_operator_oracle_executed{false};
  bool b1_hotspot_executed{false};
  bool b3_performance_executed{false};
  bool single_rank_run_present{false};
  bool distributed_run_present{false};
  bool serial_distributed_cross_validation_present{false};
  bool profile_artifacts_written{false};
  bool budget_artifacts_written{false};
  bool performance_artifacts_written{false};
  bool plot_artifacts_written{false};
  double delta_alpha_total{0.0};
  double delta_electron_total{0.0};
  double birth_source_total{0.0};
  double drag_deposition_total{0.0};
  double global_alpha_plus_electron_budget_residual{0.0};
  double serial_distributed_epsilon_alpha_linf{0.0};
  double serial_distributed_epsilon_alpha_relative_linf{0.0};
  double serial_distributed_Te_linf{0.0};
  double total_wall_seconds{0.0};
  std::vector<double> step_wall_seconds;
  std::string output_directory;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] P4AlphaBenchmarkRunResult RunP4AlphaBenchmarkCase(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept;

}  // namespace dec3d::benchmarks
