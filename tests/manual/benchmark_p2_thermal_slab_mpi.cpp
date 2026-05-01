#include "benchmarks/p2_thermal_benchmark.hpp"
#include "benchmarks/p2_thermal_benchmark_runner.hpp"

#include <mpi.h>

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

[[nodiscard]] std::vector<double> UniformOutputTimes(double final_time_s, std::size_t steps) {
  std::vector<double> times;
  times.reserve(steps);
  for (std::size_t step = 1; step <= steps; ++step) {
    times.push_back(final_time_s * static_cast<double>(step) / static_cast<double>(steps));
  }
  return times;
}

[[nodiscard]] std::string ExactOutputRoot(
    std::size_t nr,
    std::size_t ntheta,
    std::size_t nphi,
    std::size_t steps,
    const std::string& suffix) {
  std::ostringstream out;
  out << "F:/dec3d/analysis/output/p2_thermal_benchmarks/exact_convergence"
      << "_nr" << nr
      << "_nt" << ntheta
      << "_np" << nphi
      << "_steps" << steps
      << suffix;
  return out.str();
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    using dec3d::benchmarks::WooThermalBenchmarkCase;
    WooThermalBenchmarkCase case_kind =
        WooThermalBenchmarkCase::spitzer_pure_thermal;
    if (argc > 1) {
      const std::string arg = argv[1];
      if (arg == "lee_more_pure") {
        case_kind = WooThermalBenchmarkCase::lee_more_pure_thermal;
      } else if (arg == "lee_more_te_ti") {
        case_kind = WooThermalBenchmarkCase::lee_more_thermal_plus_ei;
      } else if (arg == "exact" || arg == "exact_l1m1") {
        case_kind = WooThermalBenchmarkCase::exact_l1m1_constant_kappa;
      } else if (arg == "exact_l0") {
        case_kind = WooThermalBenchmarkCase::exact_l0_radial_constant_kappa;
      } else if (arg == "flux_limiter_stress") {
        case_kind = WooThermalBenchmarkCase::flux_limiter_stress;
      } else if (arg == "ei_0d") {
        case_kind = WooThermalBenchmarkCase::ei_0d_thesis_spitzer;
      } else if (arg == "hte_smoke") {
        case_kind = WooThermalBenchmarkCase::distributed_hte_smoke;
      } else if (arg != "spitzer") {
        throw std::runtime_error(
            "unknown case; expected spitzer, lee_more_pure, lee_more_te_ti, "
            "exact_l1m1, exact_l0, flux_limiter_stress, ei_0d, or hte_smoke");
      }
    }

    auto descriptor = dec3d::benchmarks::MakeWooThermalBenchmarkDescriptor(case_kind);
    auto options = dec3d::benchmarks::WooThermalBenchmarkRunOptions{};
    if (case_kind == WooThermalBenchmarkCase::exact_l1m1_constant_kappa ||
        case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
      options.global_radial_cells = 64;
      options.theta_cells = 16;
      options.phi_cells = 16;
      std::size_t steps = descriptor.output_times_s.size();
      if (argc > 2) {
        options.global_radial_cells = static_cast<std::size_t>(std::stoull(argv[2]));
      }
      if (argc > 3) {
        options.theta_cells = static_cast<std::size_t>(std::stoull(argv[3]));
      }
      if (argc > 4) {
        options.phi_cells = static_cast<std::size_t>(std::stoull(argv[4]));
      }
      if (argc > 5) {
        steps = static_cast<std::size_t>(std::stoull(argv[5]));
      }
      if (argc > 6) {
        descriptor.comparison_time_s = std::stod(argv[6]);
      }
      std::string output_suffix;
      if (argc > 7) {
        const std::string pole_metric_arg = argv[7];
        if (pole_metric_arg == "axis_regular_pole") {
          options.pole_metric_mode =
              dec3d::transport::DiffusionPoleMetricMode::axis_regular_polar_phi;
        } else if (pole_metric_arg != "cell_centered_pole") {
          throw std::runtime_error(
              "exact_l1m1 pole metric must be cell_centered_pole or axis_regular_pole");
        } else {
          options.pole_metric_mode =
              dec3d::transport::DiffusionPoleMetricMode::cell_centered_spherical;
          output_suffix = "_cell_centered_pole";
        }
      }
      if (steps == 0u) {
        throw std::runtime_error("exact_l1m1 steps must be positive");
      }
      descriptor.output_times_s = UniformOutputTimes(descriptor.comparison_time_s, steps);
      options.output_root = ExactOutputRoot(
          options.global_radial_cells,
          options.theta_cells,
          options.phi_cells,
          steps,
          output_suffix);
    } else {
      options.global_radial_cells = 256;
      options.theta_cells = 4;
      options.phi_cells = 4;
      options.output_root = "F:/dec3d/analysis/output/p2_thermal_benchmarks";
    }
    options.write_artifacts = true;

    const auto result = dec3d::benchmarks::RunWooThermalBenchmarkCase(
        MPI_COMM_WORLD,
        descriptor,
        dec3d::benchmarks::MakeNoReferencePolicy(),
        options);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      std::cout << result.report_line << '\n';
      if (!result.success) {
        std::cerr << result.failure_diagnostics << '\n';
      }
    }
    MPI_Finalize();
    return result.success ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
