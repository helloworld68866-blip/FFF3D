#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <mpi.h>

#include <cstddef>
#include <string>

namespace dec3d::hydro {

struct MacroAleMpiCompatibilityDiagnostics {
  bool executed{false};
  int rank{0};
  int rank_count{0};
  std::size_t global_radial_cells{0u};
  std::size_t global_face_begin{0u};
  std::size_t local_face_count{0u};
  std::size_t local_radial_cells{0u};
  double local_mass_residual{0.0};
  double local_mom_r_residual{0.0};
  double local_mom_theta_residual{0.0};
  double local_mom_phi_residual{0.0};
  double local_e_fluid_total_residual{0.0};
  double local_chi_e_residual{0.0};
  double global_mass_residual{0.0};
  double global_mom_r_residual{0.0};
  double global_mom_theta_residual{0.0};
  double global_mom_phi_residual{0.0};
  double global_e_fluid_total_residual{0.0};
  double global_chi_e_residual{0.0};
  double max_seam_face_velocity_delta{0.0};
  double max_seam_proposed_face_delta{0.0};
  double max_seam_preview_radius_delta{0.0};
  bool global_stage_ok{false};
  bool global_publish_ok{false};
  bool canonical_writeback_published{false};
  bool mesh_commit_published{false};
  bool published_from_staged_geometry{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MacroAleMpiCompatibilityResult {
  bool success{false};
  bool local_stage_ok{false};
  bool global_stage_ok{false};
  bool local_publishable{false};
  bool global_publish_ok{false};
  bool canonical_writeback_published{false};
  bool mesh_commit_published{false};
  StaticGridHydroResult staged_hydro_result;
  MacroAleMpiCompatibilityDiagnostics diagnostics;
  dec3d::core::DiagnosticsPayload diagnostics_payload;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] MacroAleMpiCompatibilityResult AdvanceMacroAleMpiCompatibilityStep(
    dec3d::state::HydroStateView& local_hydro_view,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const dec3d::core::MeshUpdateProposal& global_radial_ale_proposal,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    StaticGridHydroOptions options,
    MPI_Comm communicator) noexcept;

}  // namespace dec3d::hydro
