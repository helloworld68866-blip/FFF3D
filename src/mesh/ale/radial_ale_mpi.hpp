#pragma once

#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <mpi.h>

#include <string>
#include <vector>

namespace dec3d::mesh {

struct RadialAleMpiGlobalProposalResult {
  bool success{false};
  int rank{0};
  int rank_count{0};
  std::vector<double> local_shell_velocity;
  std::vector<double> global_shell_velocity;
  dec3d::core::MeshUpdateProposal proposal;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RadialAleMpiCommitResult {
  bool success{false};
  bool local_preflight_passed{false};
  bool global_preflight_passed{false};
  bool local_stage_attempted{false};
  bool local_stage_success{false};
  bool global_stage_success{false};
  bool published{false};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] RadialAleMpiGlobalProposalResult BuildRadialAleMpiGlobalProposal(
    const dec3d::core::Array3D<double>& local_rho,
    const dec3d::core::Array3D<double>& local_mom_r,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const std::vector<double>& global_radial_faces,
    double dt_s,
    double beta,
    MPI_Comm communicator) noexcept;

[[nodiscard]] RadialAleMpiCommitResult ApplyRadialAleMpiGlobalProposalToLocalGeometry(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    MPI_Comm communicator) noexcept;

}  // namespace dec3d::mesh
