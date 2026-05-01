#pragma once

#include "benchmarks/p3_radiation_benchmark.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::benchmarks {

struct P3RadiationBenchmarkRunOptions {
  std::size_t radial_cells{16u};
  std::size_t theta_cells{1u};
  std::size_t phi_cells{1u};
  double dt_s{1.0e-18};
  std::size_t step_count{1u};
  std::string output_root{"analysis/output/p3_radiation_benchmarks"};
  std::string tops_table_root{"F:/dec3d/data/opacities/tops_dt_2026_04_27"};
  bool write_artifacts{true};
};

struct P3RadiationBenchmarkRunResult {
  bool success{false};
  bool operator_oracle_executed{false};
  bool slab_benchmark_executed{false};
  bool budget_closed{false};
  bool source_oracle_matched{false};
  bool group_count_one_regression_passed{false};
  bool profile_artifacts_written{false};
  bool budget_artifacts_written{false};
  bool metadata_written{false};
  bool diagnostics_written{false};
  double delta_radiation_total{0.0};
  double delta_electron_total{0.0};
  double boundary_leak_total{0.0};
  double global_budget_residual{0.0};
  double total_wall_seconds{0.0};
  std::vector<double> step_wall_seconds;
  std::string output_directory;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] P3RadiationBenchmarkRunResult RunP3RadiationBenchmarkCase(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy,
    const P3RadiationBenchmarkRunOptions& options) noexcept;

}  // namespace dec3d::benchmarks
