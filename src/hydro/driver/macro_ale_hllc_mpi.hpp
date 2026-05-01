#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/driver/macro_ale_hllc.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <mpi.h>

#include <string>

namespace dec3d::hydro {

struct MacroAleHllcMpiResult {
  bool success{false};
  bool local_stage_ok{false};
  bool global_stage_ok{false};
  bool local_publishable{false};
  bool global_publish_ok{false};
  bool canonical_writeback_published{false};
  bool mesh_commit_published{false};
  StaticGridHydroResult staged_hydro_result;
  MacroAleHllcDiagnostics diagnostics;
  dec3d::core::DiagnosticsPayload diagnostics_payload;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] MacroAleHllcMpiResult AdvanceMacroAleHllcMpiStep(
    dec3d::state::HydroStateView& local_hydro_view,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const dec3d::core::MeshUpdateProposal& global_radial_ale_proposal,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    StaticGridHydroOptions options,
    MPI_Comm communicator) noexcept;

}  // namespace dec3d::hydro
