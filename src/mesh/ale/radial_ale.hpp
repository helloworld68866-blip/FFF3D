#pragma once

#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::mesh {

struct RadialAleCommitResult {
  bool attempted{false};
  bool success{false};
  bool proposal_present{false};
  bool geometry_changed{false};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RadialAleLocalProposalWindow {
  bool success{false};
  std::size_t global_face_begin{0u};
  std::size_t local_face_count{0u};
  std::vector<double> radial_face_velocities;
  std::vector<double> proposed_radial_faces;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] std::vector<double> BuildRadialAleShellVelocitySlice(
    const dec3d::core::Array3D<double>& rho,
    const dec3d::core::Array3D<double>& mom_r) noexcept;

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildRadialAleProposalFromGlobalShellVelocity(
    const std::vector<double>& global_shell_velocity,
    const std::vector<double>& global_radial_faces,
    double dt_s,
    double beta = 1.0) noexcept;

[[nodiscard]] RadialAleLocalProposalWindow BuildRadialAleLocalProposalWindow(
    const dec3d::core::MeshUpdateProposal& proposal,
    std::size_t global_face_begin,
    std::size_t local_face_count) noexcept;

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildRadialAleMeshUpdateProposal(
    const dec3d::core::Array3D<double>& rho,
    const dec3d::core::Array3D<double>& mom_r,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    double beta = 1.0) noexcept;

[[nodiscard]] RadialAleCommitResult ApplyRadialAleMeshUpdateProposal(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept;

[[nodiscard]] RadialAleCommitResult ApplyRadialAleMeshUpdateProposal(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t global_face_begin) noexcept;

}  // namespace dec3d::mesh
