#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/radiation_flux_limiter.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/diffusion/generic_diffusion.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace dec3d::radiation {

inline constexpr std::uint32_t kDistributedRadiationWritesNoFields = 0u;
inline constexpr std::uint32_t kDistributedRadiationWritesRadiationGroups = 1u << 0u;
inline constexpr std::uint32_t kDistributedRadiationWritesElectronEnergy = 1u << 1u;
inline constexpr std::uint32_t kDistributedRadiationWritesFluidTotalEnergy = 1u << 2u;

enum class DistributedRadiationBoundaryModel {
  contract_zero_flux_or_scalar_remap,
  thesis_marshak_vacuum
};

struct DistributedMultigroupRadiationOptions {
  TopsOpacityProviderOptions provider_options{};
  dec3d::transport::GenericDiffusionHypreSolveOptions solve_options{};
  dec3d::state::ThermodynamicRecoveryOptions recovery_options{};
  RadiationFluxLimiterOptions radiation_flux_limiter{};
  DistributedRadiationBoundaryModel boundary_model{
      DistributedRadiationBoundaryModel::contract_zero_flux_or_scalar_remap};
  double radiation_energy_floor_erg_per_cm3{0.0};
  double electron_energy_floor_erg_per_cm3{0.0};
  double ion_energy_floor_erg_per_cm3{0.0};
  bool table_reused_across_steps{false};
};

struct DistributedMultigroupRadiationProblem {
  dec3d::transport::DistributedDiffusionRowOwnership ownership{};
  dec3d::state::CanonicalState local_state{};
  dec3d::mesh::SphericalGeometryMetadata global_geometry{};
  dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy{};
  RadiationGroupLayout group_layout{};
  const TopsOpacityTable* opacity_table{nullptr};
  double dt_s{0.0};
};

struct DistributedMultigroupRadiationResult {
  bool success{false};
  bool canonical_state_mutated{false};
  bool local_stage_ok{false};
  bool global_stage_ok{false};
  bool local_publishable{false};
  bool global_publish_ok{false};
  bool owned_slab_writeback_only{true};
  std::uint32_t updated_fields{kDistributedRadiationWritesNoFields};
  std::string updated_fields_label{"none"};

  std::size_t group_count{0u};
  std::size_t per_group_solve_count{0u};
  std::size_t local_owned_cell_count{0u};
  std::size_t global_owned_cell_count{0u};
  std::size_t first_failing_group_index{std::numeric_limits<std::size_t>::max()};
  int first_failing_rank{-1};

  double delta_radiation_total_all_groups{0.0};
  double delta_electron_total{0.0};
  double source_gain_radiation_total_all_groups{0.0};
  double boundary_leak_total_all_groups{0.0};
  double radiation_electron_exchange_residual{0.0};
  double global_radiation_plus_electron_residual{0.0};
  double coefficient_provider_wall_s{0.0};
  double flux_limiter_wall_s{0.0};
  double assembly_wall_s{0.0};
  double hypre_setup_wall_s{0.0};
  double hypre_solve_wall_s{0.0};
  double writeback_wall_s{0.0};
  int solver_iterations{0};

  std::vector<double> delta_radiation_total_by_group;
  std::vector<double> source_gain_radiation_total_by_group;
  std::vector<double> boundary_leak_total_by_group;

  bool seam_coefficient_halo_exchanged{false};
  bool off_rank_face_conductance_uses_neighbor_Dbar{false};
  std::size_t global_off_rank_column_count{0u};

  bool table_handle_supplied{false};
  bool table_reused_across_steps{false};
  bool per_step_csv_io{false};
  bool per_lookup_full_table_scan{false};
  bool marshak_enabled{false};
  bool radiation_flux_limiter_enabled{false};
  std::string radiation_flux_limiter_model{"disabled"};
  std::size_t limited_face_count_local{0u};
  std::size_t limited_face_count_global{0u};
  double min_limiter_scale{1.0};
  double max_flux_ratio_before_limit{0.0};
  bool outer_marshak_uses_face_effective_D{false};

  std::vector<std::string> per_group_provider_reports;
  std::vector<std::string> per_group_limiter_reports;
  std::vector<std::string> per_group_matrix_reports;
  std::vector<std::string> per_group_solve_reports;
  std::string coefficient_report;
  std::string thermodynamic_recovery_report;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] DistributedMultigroupRadiationResult
ApplyDistributedProviderFedMultigroupRadiation(
    DistributedMultigroupRadiationProblem& problem,
    const DistributedMultigroupRadiationOptions& options) noexcept;

[[nodiscard]] bool ValidateDistributedMultigroupRadiationDiagnostics(
    const DistributedMultigroupRadiationResult& result) noexcept;

}  // namespace dec3d::radiation
