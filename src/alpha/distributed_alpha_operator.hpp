#pragma once

#include "alpha/alpha_coefficients.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/diffusion/generic_diffusion.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <cstddef>
#include <cstdint>
#include <string>

namespace dec3d::alpha {

inline constexpr std::uint32_t kDistributedAlphaWritesNoFields = 0u;
inline constexpr std::uint32_t kDistributedAlphaWritesAlphaState = 1u << 0u;
inline constexpr std::uint32_t kDistributedAlphaWritesElectronEnergy = 1u << 1u;
inline constexpr std::uint32_t kDistributedAlphaWritesFluidTotalEnergy = 1u << 2u;

struct DistributedOneGroupAlphaTransportOptions {
  dec3d::transport::GenericDiffusionHypreSolveOptions solve_options{};
  AlphaCoefficientProviderOptions provider_options{};
  dec3d::state::ThermodynamicRecoveryOptions recovery_options{};
  dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy{};
  double electron_energy_floor_erg_per_cm3{0.0};
  double ion_energy_floor_erg_per_cm3{0.0};
  double alpha_energy_floor_erg_per_cm3{0.0};
  bool force_publish_preflight_failure_for_test{false};
};

struct DistributedOneGroupAlphaTransportProblem {
  MPI_Comm communicator{MPI_COMM_WORLD};
  dec3d::transport::DistributedDiffusionRowOwnership ownership{};
  dec3d::state::CanonicalState* local_state{nullptr};
  dec3d::mesh::SphericalGeometryMetadata global_geometry{};
  double dt_s{0.0};
};

struct DistributedOneGroupAlphaTransportResult {
  bool success{false};
  bool canonical_state_mutated{false};
  bool metadata_written{false};
  bool local_stage_ok{false};
  bool global_stage_ok{false};
  bool local_publishable{false};
  bool global_publish_ok{false};
  bool owned_slab_writeback_only{true};
  std::uint32_t updated_fields{kDistributedAlphaWritesNoFields};
  std::string updated_fields_label{"none"};

  std::size_t local_owned_cell_count{0u};
  std::size_t global_owned_cell_count{0u};
  int first_failing_rank{-1};
  std::size_t solve_count{0u};
  bool distributed_alpha_solve_executed{false};
  bool hypre_solve_executed{false};

  double dt_s{0.0};
  double delta_alpha_total_global{0.0};
  double delta_electron_total_global{0.0};
  double birth_source_total_global{0.0};
  double drag_deposition_total_global{0.0};
  double alpha_equation_budget_residual_global{0.0};
  double alpha_electron_exchange_residual_global{0.0};
  double global_alpha_plus_electron_budget_residual{0.0};
  double coefficient_provider_wall_s{0.0};
  double assembly_wall_s{0.0};
  double hypre_setup_wall_s{0.0};
  double hypre_solve_wall_s{0.0};
  double writeback_wall_s{0.0};
  int solver_iterations{0};

  double min_epsilon_alpha_before_local{0.0};
  double max_epsilon_alpha_before_local{0.0};
  double min_epsilon_alpha_after_local{0.0};
  double max_epsilon_alpha_after_local{0.0};
  double min_epsilon_alpha_after_global{0.0};
  double max_epsilon_alpha_after_global{0.0};
  double min_e_electron_after_global{0.0};
  double min_e_ion_after_global{0.0};
  double min_D_alpha_cm2_s{0.0};
  double max_D_alpha_cm2_s{0.0};
  double min_tau_alphae_s{0.0};
  double max_tau_alphae_s{0.0};
  double min_birth_source_erg_cm3_s{0.0};
  double max_birth_source_erg_cm3_s{0.0};
  double min_dt_reactivity_cm3_s{0.0};
  double max_dt_reactivity_cm3_s{0.0};

  bool seam_coefficient_halo_exchanged{false};
  bool off_rank_face_conductance_uses_neighbor_Dalpha{false};
  std::size_t global_off_rank_column_count{0u};

  std::string backend_requested{"hypre_parcsr_gmres_boomeramg"};
  std::string backend_executed{"none"};
  std::string provider_report;
  std::string assembly_report;
  std::string hypre_solve_report;
  std::string thermodynamic_recovery_report;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] DistributedOneGroupAlphaTransportResult
ApplyDistributedOneGroupAlphaTransport(
    DistributedOneGroupAlphaTransportProblem& problem,
    const DistributedOneGroupAlphaTransportOptions& options) noexcept;

[[nodiscard]] bool ValidateDistributedOneGroupAlphaDiagnostics(
    const DistributedOneGroupAlphaTransportResult& result) noexcept;

}  // namespace dec3d::alpha
