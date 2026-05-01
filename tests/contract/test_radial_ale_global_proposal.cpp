#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

[[nodiscard]] bool Expect(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << message << '\n';
    return false;
  }
  return true;
}

[[nodiscard]] bool Near(double lhs, double rhs, double tolerance = 1.0e-14) noexcept {
  return std::abs(lhs - rhs) <= tolerance;
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildGeometry() {
  dec3d::mesh::SphericalMeshDescriptor descriptor;
  descriptor.radial_cells = 8u;
  descriptor.theta_cells = 3u;
  descriptor.phi_cells = 4u;
  descriptor.inner_radius = 0.0;
  descriptor.outer_radius = 1.0;
  return dec3d::mesh::BuildSphericalGeometry(descriptor);
}

void FillState(
    dec3d::core::Array3D<double>& rho,
    dec3d::core::Array3D<double>& mom_r) {
  for (std::size_t radial = 0; radial < rho.extent_r(); ++radial) {
    const double velocity = -0.03 - 0.002 * static_cast<double>(radial);
    for (std::size_t theta = 0; theta < rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.01 * static_cast<double>(radial) +
            0.001 * static_cast<double>(theta) +
            0.0001 * static_cast<double>(phi);
        rho(radial, theta, phi) = density;
        mom_r(radial, theta, phi) = density * velocity;
      }
    }
  }
}

[[nodiscard]] bool TestGlobalBuilderMatchesLegacySingleRank() {
  const auto geometry = BuildGeometry();
  dec3d::core::Array3D<double> rho(8u, 3u, 4u);
  dec3d::core::Array3D<double> mom_r(8u, 3u, 4u);
  FillState(rho, mom_r);

  const auto shell_velocity =
      dec3d::mesh::BuildRadialAleShellVelocitySlice(rho, mom_r);
  const auto global_proposal =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          shell_velocity,
          geometry.radial_faces,
          2.5e-4,
          1.0);
  const auto legacy_proposal =
      dec3d::mesh::BuildRadialAleMeshUpdateProposal(
          rho,
          mom_r,
          geometry,
          2.5e-4,
          1.0);

  bool ok = true;
  ok &= Expect(global_proposal.is_complete_global(geometry.radial_faces.size()),
               "global_proposal_parity: global proposal is incomplete");
  ok &= Expect(legacy_proposal.is_complete(geometry.radial_faces.size()),
               "global_proposal_parity: legacy proposal is incomplete");
  ok &= Expect(global_proposal.radial_face_indexing_is_global,
               "global_proposal_parity: global proposal did not declare global face indexing");
  ok &= Expect(global_proposal.global_radial_face_count == geometry.radial_faces.size(),
               "global_proposal_parity: global proposal face count metadata is wrong");

  for (std::size_t face = 0; face < geometry.radial_faces.size(); ++face) {
    ok &= Expect(Near(global_proposal.radial_face_velocities[face],
                      legacy_proposal.radial_face_velocities[face]),
                 "global_proposal_parity: global and legacy face velocity differ");
    ok &= Expect(Near(global_proposal.proposed_radial_faces[face],
                      legacy_proposal.proposed_radial_faces[face]),
                 "global_proposal_parity: global and legacy proposed face differ");
  }
  return ok;
}

[[nodiscard]] bool TestGlobalBuilderNumericalOracle() {
  const std::vector<double> global_radial_faces{0.0, 1.0, 2.0, 4.0, 8.0};
  const std::vector<double> global_shell_velocity{0.1, -0.2, 0.3, -0.5};
  const double dt_s = 0.25;
  const double beta = 0.4;

  const auto proposal =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          global_shell_velocity,
          global_radial_faces,
          dt_s,
          beta);

  const double expected_outer_face_speed = -0.41333333333333333;
  const double expected_scale_factor = 7.8966666666666665 / 8.0;
  const std::vector<double> expected_proposed_faces{
      0.0,
      0.98708333333333331,
      1.9741666666666666,
      3.9483333333333333,
      7.8966666666666665};
  const std::vector<double> expected_face_velocities{
      0.0,
      -0.051666666666666666,
      -0.10333333333333333,
      -0.20666666666666667,
      -0.41333333333333333};

  bool ok = true;
  ok &= Expect(proposal.is_complete_global(global_radial_faces.size()),
               "global_proposal_parity: numerical oracle proposal is incomplete");
  ok &= Expect(Near(proposal.radial_face_velocities.back(), expected_outer_face_speed),
               "global_proposal_parity: numerical oracle outer face speed mismatch");
  ok &= Expect(Near(proposal.proposed_radial_faces[1u] / global_radial_faces[1u],
                    expected_scale_factor),
               "global_proposal_parity: numerical oracle scale factor mismatch");
  for (std::size_t face = 0; face < expected_proposed_faces.size(); ++face) {
    ok &= Expect(Near(proposal.proposed_radial_faces[face], expected_proposed_faces[face]),
                 "global_proposal_parity: numerical oracle proposed face mismatch");
    ok &= Expect(Near(proposal.radial_face_velocities[face], expected_face_velocities[face]),
                 "global_proposal_parity: numerical oracle face velocity mismatch");
  }
  return ok;
}

[[nodiscard]] bool TestHomologousShellVelocityBuildsHomologousFaceMotion() {
  const std::vector<double> global_radial_faces{0.0, 0.25, 0.50, 0.75, 1.0};
  std::vector<double> shell_velocity;
  for (std::size_t radial = 0; radial + 1u < global_radial_faces.size(); ++radial) {
    const double shell_center =
        0.5 * (global_radial_faces[radial] + global_radial_faces[radial + 1u]);
    shell_velocity.push_back(-shell_center);
  }

  const double dt_s = 1.0e-3;
  const auto proposal =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          shell_velocity,
          global_radial_faces,
          dt_s,
          1.0);

  bool ok = true;
  ok &= Expect(proposal.is_complete_global(global_radial_faces.size()),
               "homologous_face_motion: proposal is incomplete");
  for (std::size_t face = 0; face < global_radial_faces.size(); ++face) {
    const double expected_face_speed = -global_radial_faces[face];
    const double expected_face_radius = global_radial_faces[face] * (1.0 - dt_s);
    ok &= Expect(Near(proposal.radial_face_velocities[face], expected_face_speed, 1.0e-12),
                 "homologous_face_motion: face speed is not homologous");
    ok &= Expect(Near(proposal.proposed_radial_faces[face], expected_face_radius, 1.0e-14),
                 "homologous_face_motion: proposed face is not homologous");
  }
  return ok;
}

[[nodiscard]] bool TestGlobalBuilderIsDeterministic() {
  const auto geometry = BuildGeometry();
  std::vector<double> shell_velocity;
  for (std::size_t radial = 0; radial + 1u < geometry.radial_faces.size(); ++radial) {
    shell_velocity.push_back(-0.02 - 0.001 * static_cast<double>(radial));
  }

  const auto first =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          shell_velocity,
          geometry.radial_faces,
          1.0e-4,
          0.75);
  const auto second =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          shell_velocity,
          geometry.radial_faces,
          1.0e-4,
          0.75);

  bool ok = true;
  ok &= Expect(first.is_complete_global(geometry.radial_faces.size()),
               "global_proposal_parity: first deterministic proposal is incomplete");
  ok &= Expect(second.is_complete_global(geometry.radial_faces.size()),
               "global_proposal_parity: second deterministic proposal is incomplete");
  ok &= Expect(first.radial_face_velocities == second.radial_face_velocities,
               "global_proposal_parity: deterministic builder changed face velocities");
  ok &= Expect(first.proposed_radial_faces == second.proposed_radial_faces,
               "global_proposal_parity: deterministic builder changed proposed faces");
  return ok;
}

[[nodiscard]] bool TestLocalWindowExtraction() {
  const auto geometry = BuildGeometry();
  std::vector<double> shell_velocity(8u, -0.01);
  const auto proposal =
      dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
          shell_velocity,
          geometry.radial_faces,
          1.0e-4,
          1.0);
  const auto window =
      dec3d::mesh::BuildRadialAleLocalProposalWindow(
          proposal,
          3u,
          4u);

  bool ok = true;
  ok &= Expect(window.is_complete(), "proposal_positive_widths: valid local proposal window is incomplete");
  ok &= Expect(window.global_face_begin == 3u, "proposal_positive_widths: window begin index is wrong");
  ok &= Expect(window.local_face_count == 4u, "proposal_positive_widths: window face count is wrong");
  for (std::size_t local_face = 0; local_face < 4u; ++local_face) {
    const std::size_t global_face = 3u + local_face;
    ok &= Expect(Near(window.radial_face_velocities[local_face],
                      proposal.radial_face_velocities[global_face]),
                 "proposal_positive_widths: window face velocity is not a direct global slice");
    ok &= Expect(Near(window.proposed_radial_faces[local_face],
                      proposal.proposed_radial_faces[global_face]),
                 "proposal_positive_widths: window proposed face is not a direct global slice");
  }
  return ok;
}

[[nodiscard]] bool TestFailedCommitDoesNotMutateGeometry() {
  auto geometry = BuildGeometry();
  const auto original_faces = geometry.radial_faces;
  const auto original_volumes = geometry.cell_volumes;

  dec3d::core::MeshUpdateProposal proposal;
  proposal.requested = true;
  proposal.radial_ale = true;
  proposal.radial_face_indexing_is_global = true;
  proposal.global_radial_face_count = geometry.radial_faces.size();
  proposal.dt_s = 1.0e-4;
  proposal.radial_face_velocities.assign(geometry.radial_faces.size(), 0.0);
  proposal.proposed_radial_faces = geometry.radial_faces;
  proposal.proposed_radial_faces[4u] = proposal.proposed_radial_faces[3u];
  proposal.implementation_id = "p1.mesh.ale.radial_proposal";
  proposal.summary = "synthetic nonmonotone proposal";

  const auto result =
      dec3d::mesh::ApplyRadialAleMeshUpdateProposal(proposal, geometry, 0u);

  bool ok = true;
  ok &= Expect(result.attempted, "transaction_failure_no_partial_commit: commit did not attempt validation");
  ok &= Expect(!result.success, "transaction_failure_no_partial_commit: nonmonotone proposal committed");
  ok &= Expect(!result.failure_reason.empty(), "transaction_failure_no_partial_commit: failed commit lacked diagnostics");
  ok &= Expect(geometry.radial_faces == original_faces,
               "transaction_failure_no_partial_commit: failed commit mutated radial faces");
  ok &= Expect(geometry.cell_volumes == original_volumes,
               "transaction_failure_no_partial_commit: failed commit mutated cell volumes");
  return ok;
}

}  // namespace

int main() {
  bool ok = true;
  ok &= TestGlobalBuilderNumericalOracle();
  ok &= TestHomologousShellVelocityBuildsHomologousFaceMotion();
  ok &= TestGlobalBuilderMatchesLegacySingleRank();
  ok &= TestGlobalBuilderIsDeterministic();
  ok &= TestLocalWindowExtraction();
  ok &= TestFailedCommitDoesNotMutateGeometry();
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
