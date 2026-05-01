#include "benchmarks/p3_radiation_benchmark.hpp"
#include "benchmarks/p3_radiation_benchmark_runner.hpp"

#include <cstdlib>
#include <iostream>

int main(int argc, char** argv) {
  dec3d::benchmarks::P3RadiationBenchmarkRunOptions options;
  options.radial_cells = 32u;
  options.theta_cells = 1u;
  options.phi_cells = 2u;
  options.dt_s = 1.0e-18;
  options.output_root = "F:/dec3d/analysis/output/p3_radiation_benchmarks";
  options.write_artifacts = true;

  if (argc > 1) {
    options.radial_cells = static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10));
  }
  if (argc > 2) {
    options.dt_s = std::strtod(argv[2], nullptr);
  }

  const auto result = dec3d::benchmarks::RunP3RadiationBenchmarkCase(
      dec3d::benchmarks::MakeP3RadiationBenchmarkDescriptor(
          dec3d::benchmarks::P3RadiationBenchmarkCase::b1_woo_5_25_slab_no_limiter),
      dec3d::benchmarks::MakeNoP3RadiationReferencePolicy(),
      options);

  if (!result.success) {
    std::cerr << result.failure_reason << '\n'
              << result.failure_diagnostics << '\n';
    return 1;
  }

  std::cout << result.report_line << '\n'
            << "output_directory=" << result.output_directory << '\n';
  return 0;
}
