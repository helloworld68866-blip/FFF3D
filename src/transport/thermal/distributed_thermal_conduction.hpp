#pragma once

#include "state/canonical_state/canonical_state.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"
#include "transport/thermal/electron_flux_limiter.hpp"
#include "transport/thermal/thermal_conduction.hpp"

#include <mpi.h>

#include <cstdint>
#include <cstddef>
#include <string>

namespace dec3d::transport {

inline constexpr std::uint32_t kDistributedThermalWritesNoFields = 0u;
inline constexpr std::uint32_t kDistributedThermalWritesElectronEnergy = 1u << 0u;
inline constexpr std::uint32_t kDistributedThermalWritesFluidTotalEnergy = 1u << 1u;

struct DistributedThermalConductionOptions {
  double dt_s{0.0};
  ThermalConductionKappaModel kappa_model{ThermalConductionKappaModel::spitzer_no_degeneracy};
  ElectronFluxLimiterModel electron_flux_limiter_model{ElectronFluxLimiterModel::disabled};
  double electron_flux_limiter_alpha_e{0.0};
  double electron_energy_floor{0.0};
  double ion_energy_floor{0.0};
};

struct DistributedThermalConductionProblem {
  DistributedDiffusionRowOwnership ownership;
  dec3d::state::CanonicalState local_state;
  dec3d::mesh::SphericalGeometryMetadata global_geometry;
  GenericDiffusionBoundaryPolicy boundary_policy;
};

struct DistributedThermalConductionResult {
  bool success{false};
  bool distributed_writeback{false};
  bool gathered_writeback_used{false};
  bool global_stage_ok{false};
  bool global_publish_ok{false};
  std::size_t local_owned_write_count{0};
  std::size_t global_owned_write_count{0};
  std::uint32_t updated_fields{kDistributedThermalWritesNoFields};
  std::string updated_fields_label{"none"};
  double global_thermal_energy_residual{0.0};
  double global_max_cellwise_thermal_energy_mismatch{0.0};
  double coefficient_provider_wall_s{0.0};
  double assembly_wall_s{0.0};
  double hypre_setup_wall_s{0.0};
  double hypre_solve_wall_s{0.0};
  double writeback_wall_s{0.0};
  int solver_iterations{0};
  std::string electron_matrix_report;
  std::string ion_matrix_report;
  std::string electron_solve_report;
  std::string ion_solve_report;
  std::string coefficient_report;
  std::string flux_limiter_report;
  std::string writeback_report;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] DistributedThermalConductionProblem BuildDistributedThermalConductionSmokeProblem(
    MPI_Comm communicator,
    std::size_t global_radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells);

[[nodiscard]] DistributedThermalConductionResult ApplyDistributedVariableKappaThermalConduction(
    DistributedThermalConductionProblem& problem,
    const DistributedThermalConductionOptions& options) noexcept;

[[nodiscard]] bool ValidateDistributedThermalConductionDiagnostics(
    const DistributedThermalConductionResult& result) noexcept;

}  // namespace dec3d::transport
