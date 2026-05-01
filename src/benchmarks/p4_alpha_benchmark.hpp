#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::benchmarks {

enum class P4AlphaBenchmarkLayer {
  b0_operator_budget_oracle,
  b1_frozen_fluid_alpha_hotspot,
  b3_distributed_performance
};

enum class P4AlphaBenchmarkCase {
  b0_one_cell_birth_drag_oracle,
  b0_one_cell_dt_zero_no_change,
  b0_two_cell_diffusion_coupling_oracle,
  b0_zero_flux_budget_closure,
  b0_distributed_matches_single_rank,
  b1_frozen_fluid_alpha_hotspot,
  b3_distributed_performance
};

struct P4AlphaBenchmarkDescriptor {
  P4AlphaBenchmarkLayer layer{P4AlphaBenchmarkLayer::b0_operator_budget_oracle};
  P4AlphaBenchmarkCase case_kind{
      P4AlphaBenchmarkCase::b0_one_cell_birth_drag_oracle};
  std::string benchmark_id{"P4-B0"};
  std::string case_id{"p4_b0_one_cell_birth_drag_oracle"};
  std::string claim_level{"operator_oracle"};
  std::string geometry{"1d_spherical_slab"};
  std::size_t nr{16u};
  double R_cm{1.0e-2};
  double r0_cm{5.0e-3};
  double rho_inner_g_cm3{50.0};
  double rho_outer_g_cm3{100.0};
  double Te_inner_keV{5.0};
  double Te_outer_keV{0.5};
  std::string initial_ti_policy{"equal_to_initial_te"};
  std::string epsilon_alpha_initial_policy{"zero"};
  std::vector<double> output_times_s{1.0e-18};
  std::vector<double> comparison_times_s{1.0e-18};
  std::string boundary_model{"benchmark_1d_spherical_scalar_origin_remap_outer_zero_flux"};
  std::string composition_model{"equimolar_dt_from_p2_recovery"};
  std::string reactivity_model{"bosch_hale_dt"};
  std::string alpha_transport_model{"atzeni_one_group"};
  bool fuel_depletion_enabled{false};
  bool separate_dt_species_authoritative{false};
  bool parity_claim_allowed{false};
};

struct P4AlphaBenchmarkValidationResult {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct P4AlphaBenchmarkArtifactManifest {
  std::string output_directory;
  std::vector<std::string> required_files;

  [[nodiscard]] bool ContainsRequiredFile(
      const std::string& file_name) const noexcept;
};

[[nodiscard]] const char* ToString(P4AlphaBenchmarkLayer layer) noexcept;
[[nodiscard]] const char* ToString(P4AlphaBenchmarkCase case_kind) noexcept;

[[nodiscard]] P4AlphaBenchmarkDescriptor MakeP4AlphaBenchmarkDescriptor(
    P4AlphaBenchmarkCase case_kind);

[[nodiscard]] P4AlphaBenchmarkValidationResult ValidateP4AlphaBenchmarkDescriptor(
    const P4AlphaBenchmarkDescriptor& descriptor) noexcept;

[[nodiscard]] std::string BuildP4AlphaBenchmarkDiagnosticsLine(
    const P4AlphaBenchmarkDescriptor& descriptor);

[[nodiscard]] bool ValidateP4AlphaBenchmarkDiagnostics(
    const std::string& report_line) noexcept;

[[nodiscard]] P4AlphaBenchmarkArtifactManifest BuildP4AlphaBenchmarkArtifactManifest(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const std::string& output_root);

}  // namespace dec3d::benchmarks
