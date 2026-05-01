#pragma once

#include "benchmarks/p2_thermal_benchmark.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <mpi.h>

#include <cstddef>
#include <string>

namespace dec3d::benchmarks {

struct WooThermalBenchmarkRunOptions {
  std::size_t global_radial_cells{16};
  std::size_t theta_cells{4};
  std::size_t phi_cells{4};
  std::string output_root{"analysis/output/p2_thermal_benchmarks"};
  bool write_artifacts{true};
  dec3d::transport::DiffusionPoleMetricMode pole_metric_mode{
      dec3d::transport::DiffusionPoleMetricMode::axis_regular_polar_phi};
};

struct WooThermalBenchmarkRunResult {
  bool success{false};
  bool thermal_stage_executed{false};
  bool equilibration_stage_executed{false};
  bool momentum_update_executed{false};
  bool profile_artifacts_written{false};
  bool energy_budget_written{false};
  bool benchmark_summary_written{false};
  bool benchmark_diagnostics_written{false};
  bool comparison_metrics_written{false};
  bool native_profiles_written{false};
  std::string output_directory;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] WooThermalBenchmarkRunResult RunWooThermalBenchmarkCase(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalReferencePolicy& reference_policy,
    const WooThermalBenchmarkRunOptions& options) noexcept;

}  // namespace dec3d::benchmarks
