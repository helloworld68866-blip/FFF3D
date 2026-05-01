#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::benchmarks {

enum class P5IntegratedBenchmarkCase {
  clean_spherical,
  legendre_l2_m0,
  legendre_l4_m0
};

enum class P5IntegratedBenchmarkVariant {
  primary_mpi24_ale,
  no_ale_mpi24,
  mpi12_ale
};

struct P5IntegratedBenchmarkDescriptor {
  P5IntegratedBenchmarkCase case_kind{
      P5IntegratedBenchmarkCase::clean_spherical};
  P5IntegratedBenchmarkVariant variant{
      P5IntegratedBenchmarkVariant::primary_mpi24_ale};
  std::string benchmark_id{"P5"};
  std::string case_id{"p5_clean_spherical"};
  std::string variant_id{"primary"};
  std::size_t radial_cells{64u};
  std::size_t theta_cells{32u};
  std::size_t phi_cells{64u};
  double outer_radius_cm{5.0e-3};
  double perturbation_r0_cm{24.070450097847356e-4};
  double perturbation_shape_width_fraction{0.015};
  double perturbation_amplitude{0.07};
  int mode_l{0};
  int mode_m{0};
  bool ppm_enabled{true};
  bool macro_enabled{true};
  bool ale_enabled{true};
  int mpi_ranks{24};
  bool parity_claim_allowed{false};
  bool lilac_reference_available{false};
  bool one_dimensional_reference_available{false};
  std::string initial_profile_yaml_path{
      "F:/python_project/cases/image_profile_100ps_hydro.yaml"};
  std::string initial_profile_csv_path{
      "F:/dec3d/analysis/output/p5_initial_profile_from_image_profile_100ps_hydro.csv"};
  std::string perturbation_r0_source{"p5_case_descriptor"};
  double woo_reference_r0_um{68.0};
  bool r0_matches_woo_77068{false};
};

struct P5IntegratedBenchmarkArtifactManifest {
  std::string output_directory;
  std::vector<std::string> required_files;
  std::vector<std::string> visual_artifacts;
  std::vector<std::string> comparison_only_files;
  std::vector<std::string> missing_by_design;

  [[nodiscard]] bool ContainsRequiredFile(
      const std::string& file_name) const noexcept;
  [[nodiscard]] bool ContainsVisualArtifact(
      const std::string& file_name) const noexcept;
};

[[nodiscard]] const char* ToString(P5IntegratedBenchmarkCase case_kind) noexcept;
[[nodiscard]] const char* ToString(P5IntegratedBenchmarkVariant variant) noexcept;

[[nodiscard]] P5IntegratedBenchmarkDescriptor MakeP5IntegratedBenchmarkDescriptor(
    P5IntegratedBenchmarkCase case_kind,
    P5IntegratedBenchmarkVariant variant) noexcept;

[[nodiscard]] double LegendreP(int ell, double mu) noexcept;
[[nodiscard]] double WooVelocityPerturbationShape(
    int ell,
    double r_cm,
    double r0_cm) noexcept;
[[nodiscard]] double WooVelocityPerturbationDeltaVr(
    int ell,
    double r_cm,
    double theta_rad,
    double r0_cm,
    double amplitude,
    double reference_velocity_cm_s) noexcept;

[[nodiscard]] std::string BuildP5IntegratedBenchmarkDiagnosticsLine(
    const P5IntegratedBenchmarkDescriptor& descriptor);
[[nodiscard]] bool ValidateP5IntegratedBenchmarkDiagnostics(
    const std::string& report_line) noexcept;

[[nodiscard]] P5IntegratedBenchmarkArtifactManifest
BuildP5IntegratedBenchmarkArtifactManifest(
    const P5IntegratedBenchmarkDescriptor& descriptor,
    const std::string& output_root);

}  // namespace dec3d::benchmarks
