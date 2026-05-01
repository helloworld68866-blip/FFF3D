#pragma once

#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace dec3d::hydro {

struct SphericalRadialWaveSanityCheck {
  bool success{false};
  bool angular_symmetry_preserved{false};
  bool tangential_momentum_quiet{false};
  bool radial_profile_changed{false};
  double max_rho_angular_spread{0.0};
  double max_mom_r_angular_spread{0.0};
  double max_e_fluid_total_angular_spread{0.0};
  double max_abs_mom_theta{0.0};
  double max_abs_mom_phi{0.0};
  double max_abs_rho_delta{0.0};
  double max_abs_mom_r{0.0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct SedovShockRadiusSample {
  double time_s{0.0};
  double numerical_shock_radius{0.0};
  double reference_shock_radius{0.0};
  double relative_error{0.0};
};

struct SedovRadialProfileSample {
  std::size_t step_index{0};
  double time_s{0.0};
  double numerical_shock_radius{0.0};
  double reference_shock_radius{0.0};
  std::vector<double> radial_centers;
  std::vector<double> rho_avg;
  std::vector<double> hydro_only_te_avg;
  std::vector<double> velocity_magnitude_avg;
  std::vector<double> mom_r_avg;
  std::vector<double> e_fluid_total_avg;
};

[[nodiscard]] SedovRadialProfileSample BuildSedovRadialProfileSample(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t step_index,
    double time_s) noexcept;

[[nodiscard]] bool WriteSedovRadialProfileSampleFile(
    const std::filesystem::path& path,
    const char* label,
    const SedovRadialProfileSample& sample) noexcept;

struct SedovSphericalBlastCheck {
  bool success{false};
  bool shock_radius_within_tolerance{false};
  bool shell_symmetry_preserved{false};
  bool tangential_momentum_quiet{false};
  double final_time_s{0.0};
  double numerical_shock_radius{0.0};
  double reference_shock_radius{0.0};
  double shock_radius_relative_error{0.0};
  double max_rho_angular_relative_spread{0.0};
  double tangential_to_radial_momentum_ratio{0.0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] SedovSphericalBlastCheck EvaluateSedovSphericalBlast(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double blast_energy,
    double ambient_density,
    double final_time_s,
    double shock_radius_relative_tolerance,
    double shell_symmetry_tolerance,
    double tangential_momentum_ratio_tolerance) noexcept;

[[nodiscard]] bool WriteSedovSphericalBlastOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const std::filesystem::path& output_directory,
    const SedovSphericalBlastCheck& summary,
    const std::vector<SedovShockRadiusSample>& shock_radius_history,
    const std::vector<SedovRadialProfileSample>& profile_history,
    const HydroBudgetResidualSummary* budget = nullptr,
    const char* case_name = "case_sedov_spherical") noexcept;

[[nodiscard]] SphericalRadialWaveSanityCheck EvaluateSphericalRadialWaveSanity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    double angular_tolerance,
    double tangential_momentum_tolerance,
    double radial_change_tolerance) noexcept;

struct SphericalRadialWaveMpiParityCheck {
  bool success{false};
  bool rank_decomposition_valid{false};
  bool parity_within_tolerance{false};
  bool seam_jumps_bounded{false};
  bool ppm_executed{false};
  bool macro_zoning_executed{false};
  std::size_t rank_count{0};
  std::size_t reconstruction_ghost_layers{0};
  double max_rho_difference{0.0};
  double max_mom_r_difference{0.0};
  double max_e_fluid_total_difference{0.0};
  double max_seam_rho_jump{0.0};
  double max_seam_velocity_jump{0.0};
  std::string failure_reason;
  std::string reconstruction_mode;
  std::string macro_zoning_mode;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] SphericalRadialWaveMpiParityCheck EvaluateSphericalRadialWaveMpiParity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    std::size_t rank_count,
    double parity_tolerance,
    double seam_tolerance,
    dec3d::state::CanonicalState* multi_rank_after_out = nullptr) noexcept;

[[nodiscard]] bool WriteSphericalRadialWaveMpiParityOutputs(
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::state::CanonicalState& multi_rank_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const std::filesystem::path& output_directory,
    const SphericalRadialWaveMpiParityCheck& summary,
    const HydroBudgetResidualSummary* single_rank_budget = nullptr,
    const HydroBudgetResidualSummary* multi_rank_budget = nullptr,
    const char* case_name = "case2_radial_wave_mpi_parity") noexcept;

struct RadialAleOnOffCheck {
  bool success{false};
  bool proposal_complete{false};
  bool moving_mesh_applied{false};
  bool ale_changes_solution{false};
  bool ppm_executed{false};
  std::size_t reconstruction_ghost_layers{0};
  double max_rho_difference{0.0};
  double max_mom_r_difference{0.0};
  double max_e_fluid_total_difference{0.0};
  double outer_face_displacement{0.0};
  double max_face_speed{0.0};
  std::string failure_reason;
  std::string reconstruction_mode;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] RadialAleOnOffCheck EvaluateRadialAleOnOff(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& ale_off_after,
    const dec3d::state::CanonicalState& ale_on_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::MeshUpdateProposal& proposal,
    double difference_tolerance,
    bool ppm_executed,
    std::size_t reconstruction_ghost_layers) noexcept;

[[nodiscard]] bool WriteRadialAleOnOffOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& ale_off_after,
    const dec3d::state::CanonicalState& ale_on_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::MeshUpdateProposal& proposal,
    const std::filesystem::path& output_directory,
    const RadialAleOnOffCheck& summary,
    const HydroBudgetResidualSummary* ale_off_budget = nullptr,
    const HydroBudgetResidualSummary* ale_on_budget = nullptr,
    const char* case_name = "case2_radial_wave_ale_on_off") noexcept;

struct DirectionalMpiParityCheck {
  bool success{false};
  bool parity_within_tolerance{false};
  bool ppm_executed{false};
  bool macro_zoning_executed{false};
  std::size_t rank_count{0};
  std::size_t reconstruction_ghost_layers{0};
  double max_rho_difference{0.0};
  double max_hydro_only_te_difference{0.0};
  double max_velocity_magnitude_difference{0.0};
  std::string failure_reason;
  std::string reconstruction_mode;
  std::string macro_zoning_mode;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] DirectionalMpiParityCheck EvaluateDirectionalMpiParity(
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::state::CanonicalState& multi_rank_after,
    std::size_t rank_count,
    double parity_tolerance) noexcept;

struct True3DLowModeSanityCheck {
  bool success{false};
  bool finite_states{false};
  bool full_directional_differs_from_radial_only{false};
  bool mode_projection_nontrivial{false};
  std::size_t shell_index{0};
  double rho_l1_difference{0.0};
  double te_l1_difference{0.0};
  double velocity_l1_difference{0.0};
  double mode_amplitude_before{0.0};
  double mode_amplitude_after_full{0.0};
  double mode_amplitude_after_radial_only{0.0};
  double mode_phase_before{0.0};
  double mode_phase_after_full{0.0};
  double mode_phase_after_radial_only{0.0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] True3DLowModeSanityCheck EvaluateTrue3DLowModeSanity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after_full_directional,
    const dec3d::state::CanonicalState& after_radial_only,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    double difference_tolerance,
    double mode_tolerance) noexcept;

[[nodiscard]] bool WriteTrue3DLowModeOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after_full_directional,
    const dec3d::state::CanonicalState& after_radial_only,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    const std::filesystem::path& output_directory,
    const True3DLowModeSanityCheck& summary,
    const HydroBudgetResidualSummary* budget = nullptr,
    const DirectionalMpiParityCheck* mpi_parity = nullptr,
    const HydroBudgetResidualSummary* single_rank_budget = nullptr,
    const char* case_name = "case3_true3d_low_mode") noexcept;

struct PoleAdjacentSeamStressCheck {
  bool success{false};
  bool finite_solution{false};
  bool seam_growth_bounded{false};
  bool polar_caps_stable{false};
  std::size_t shell_index{0};
  double max_seam_rho_jump_before{0.0};
  double max_seam_rho_jump_after{0.0};
  double max_seam_mom_theta_jump_before{0.0};
  double max_seam_mom_theta_jump_after{0.0};
  double max_seam_mom_phi_jump_before{0.0};
  double max_seam_mom_phi_jump_after{0.0};
  double max_seam_velocity_jump_before{0.0};
  double max_seam_velocity_jump_after{0.0};
  double max_north_cap_velocity{0.0};
  double max_south_cap_velocity{0.0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] PoleAdjacentSeamStressCheck EvaluatePoleAdjacentSeamStress(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    double seam_growth_tolerance,
    double polar_velocity_tolerance) noexcept;

[[nodiscard]] bool WritePoleAdjacentSeamStressOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    const std::filesystem::path& output_directory,
    const PoleAdjacentSeamStressCheck& summary,
    const HydroBudgetResidualSummary* budget = nullptr,
    const DirectionalMpiParityCheck* mpi_parity = nullptr,
    const HydroBudgetResidualSummary* single_rank_budget = nullptr,
    const char* case_name = "case4_pole_seam_stress") noexcept;

}  // namespace dec3d::hydro
