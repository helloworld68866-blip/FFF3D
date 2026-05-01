#include "mesh/ale/radial_ale.hpp"

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <sstream>

namespace dec3d::mesh {

namespace {

[[nodiscard]] double ShellCenterRadius(
    const std::vector<double>& radial_faces,
    std::size_t shell) noexcept {
  return 0.5 * (radial_faces[shell] + radial_faces[shell + 1u]);
}

[[nodiscard]] double ExtrapolateOuterFaceVelocityFromShellCenters(
    const std::vector<double>& global_shell_velocity,
    const std::vector<double>& global_radial_faces) noexcept {
  if (global_shell_velocity.empty() ||
      global_radial_faces.size() != global_shell_velocity.size() + 1u) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const std::size_t last_shell = global_shell_velocity.size() - 1u;
  const double last_velocity = global_shell_velocity[last_shell];
  if (global_shell_velocity.size() == 1u) {
    return last_velocity;
  }

  const std::size_t previous_shell = last_shell - 1u;
  const double previous_center =
      ShellCenterRadius(global_radial_faces, previous_shell);
  const double last_center =
      ShellCenterRadius(global_radial_faces, last_shell);
  const double center_delta = last_center - previous_center;
  if (!(center_delta > 0.0) ||
      !std::isfinite(center_delta) ||
      !std::isfinite(previous_center) ||
      !std::isfinite(last_center) ||
      !std::isfinite(global_shell_velocity[previous_shell]) ||
      !std::isfinite(last_velocity)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double slope =
      (last_velocity - global_shell_velocity[previous_shell]) / center_delta;
  const double outer_radius = global_radial_faces.back();
  return last_velocity + slope * (outer_radius - last_center);
}

void RebuildCellVolumes(SphericalGeometryMetadata& geometry) {
  if (!geometry.valid) {
    geometry.cell_volumes.clear();
    geometry.global_volume = 0.0;
    return;
  }

  const std::size_t radial_cells = geometry.radial_faces.size() - 1u;
  const std::size_t theta_cells = geometry.theta_faces.size() - 1u;
  const std::size_t phi_cells = geometry.phi_faces.size() - 1u;

  geometry.cell_volumes.assign(radial_cells * theta_cells * phi_cells, 0.0);
  geometry.global_volume = 0.0;

  std::size_t linear_index = 0u;
  for (std::size_t radial = 0; radial < radial_cells; ++radial) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[radial + 1u], 3.0) -
         std::pow(geometry.radial_faces[radial], 3.0)) /
        3.0;

    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      const double theta_factor =
          std::cos(geometry.theta_faces[theta]) -
          std::cos(geometry.theta_faces[theta + 1u]);

      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const double phi_factor = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
        const double cell_volume = radial_factor * theta_factor * phi_factor;
        geometry.cell_volumes[linear_index++] = cell_volume;
        geometry.global_volume += cell_volume;
      }
    }
  }
}

}  // namespace

bool RadialAleCommitResult::is_complete() const noexcept {
  return attempted && !report_line.empty();
}

bool RadialAleLocalProposalWindow::is_complete() const noexcept {
  return success &&
         local_face_count > 0u &&
         radial_face_velocities.size() == local_face_count &&
         proposed_radial_faces.size() == local_face_count &&
         failure_reason.empty() &&
         !report_line.empty();
}

std::vector<double> BuildRadialAleShellVelocitySlice(
    const dec3d::core::Array3D<double>& rho,
    const dec3d::core::Array3D<double>& mom_r) noexcept {
  if (rho.empty() ||
      mom_r.empty() ||
      rho.extent_r() != mom_r.extent_r() ||
      rho.extent_theta() != mom_r.extent_theta() ||
      rho.extent_phi() != mom_r.extent_phi()) {
    return {};
  }

  std::vector<double> shell_velocities(rho.extent_r(), 0.0);
  for (std::size_t radial = 0; radial < rho.extent_r(); ++radial) {
    double accumulated_velocity = 0.0;
    std::size_t cell_count = 0u;
    for (std::size_t theta = 0; theta < rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < rho.extent_phi(); ++phi) {
        const double rho_value = rho(radial, theta, phi);
        if (rho_value <= 0.0) {
          continue;
        }

        accumulated_velocity += mom_r(radial, theta, phi) / rho_value;
        ++cell_count;
      }
    }

    shell_velocities[radial] =
        cell_count == 0u ? 0.0 : accumulated_velocity / static_cast<double>(cell_count);
  }

  return shell_velocities;
}

dec3d::core::MeshUpdateProposal BuildRadialAleProposalFromGlobalShellVelocity(
    const std::vector<double>& global_shell_velocity,
    const std::vector<double>& global_radial_faces,
    double dt_s,
    double beta) noexcept {
  dec3d::core::MeshUpdateProposal proposal;
  if (global_shell_velocity.empty() ||
      global_radial_faces.size() != global_shell_velocity.size() + 1u ||
      dt_s <= 0.0 ||
      !std::isfinite(dt_s) ||
      !std::isfinite(beta)) {
    proposal.summary = "radial ALE global proposal inputs are incomplete";
    return proposal;
  }

  for (const double face : global_radial_faces) {
    if (!std::isfinite(face)) {
      proposal.summary = "radial ALE global face array contains non-finite radius";
      return proposal;
    }
  }

  for (std::size_t face = 1u; face < global_radial_faces.size(); ++face) {
    if (!(global_radial_faces[face] > global_radial_faces[face - 1u])) {
      proposal.summary = "radial ALE global face array is not strictly increasing";
      return proposal;
    }
  }

  proposal.requested = true;
  proposal.radial_ale = true;
  proposal.radial_face_indexing_is_global = true;
  proposal.global_radial_face_count = global_radial_faces.size();
  proposal.dt_s = dt_s;
  proposal.radial_face_velocities.assign(global_radial_faces.size(), 0.0);
  proposal.proposed_radial_faces = global_radial_faces;
  proposal.implementation_id = "p1.mesh.ale.radial_proposal";

  const double old_outer_radius = global_radial_faces.back();
  const double outer_face_velocity =
      ExtrapolateOuterFaceVelocityFromShellCenters(
          global_shell_velocity,
          global_radial_faces);
  const double outer_face_speed = beta * outer_face_velocity;
  const double new_outer_radius = old_outer_radius + (dt_s * outer_face_speed);
  if (!std::isfinite(outer_face_velocity) ||
      !(old_outer_radius > 0.0) ||
      !(new_outer_radius > 0.0) ||
      !std::isfinite(new_outer_radius)) {
    proposal.summary = "radial ALE proposal would produce a non-positive outer radius";
    return proposal;
  }

  const double scale_factor = new_outer_radius / old_outer_radius;
  double max_face_speed = 0.0;
  for (std::size_t face = 0; face < global_radial_faces.size(); ++face) {
    proposal.proposed_radial_faces[face] = global_radial_faces[face] * scale_factor;
    proposal.radial_face_velocities[face] =
        (proposal.proposed_radial_faces[face] - global_radial_faces[face]) / dt_s;
    max_face_speed =
        std::max(max_face_speed, std::abs(proposal.radial_face_velocities[face]));
  }

  std::ostringstream summary;
  summary << "radial_face_indexing=global"
          << "; global_radial_faces=" << global_radial_faces.size()
          << "; global_shell_velocity_count=" << global_shell_velocity.size()
          << "; dt_s=" << dt_s
          << "; beta=" << beta
          << "; outer_face_speed=" << outer_face_speed
          << "; outer_face_velocity_source=linear_shell_center_extrapolation"
          << "; scale_factor=" << scale_factor
          << "; max_face_speed=" << max_face_speed;
  proposal.summary = summary.str();
  return proposal;
}

RadialAleLocalProposalWindow BuildRadialAleLocalProposalWindow(
    const dec3d::core::MeshUpdateProposal& proposal,
    std::size_t global_face_begin,
    std::size_t local_face_count) noexcept {
  RadialAleLocalProposalWindow window;
  window.global_face_begin = global_face_begin;
  window.local_face_count = local_face_count;

  if (!proposal.has_radial_face_window(global_face_begin, local_face_count)) {
    window.failure_reason =
        "radial ALE global proposal does not contain requested local face window";
  } else {
    window.radial_face_velocities.assign(
        proposal.radial_face_velocities.begin() + static_cast<std::ptrdiff_t>(global_face_begin),
        proposal.radial_face_velocities.begin() +
            static_cast<std::ptrdiff_t>(global_face_begin + local_face_count));
    window.proposed_radial_faces.assign(
        proposal.proposed_radial_faces.begin() + static_cast<std::ptrdiff_t>(global_face_begin),
        proposal.proposed_radial_faces.begin() +
            static_cast<std::ptrdiff_t>(global_face_begin + local_face_count));
    window.success = true;
  }

  std::ostringstream report;
  report << "global_face_begin=" << global_face_begin
         << "; local_face_count=" << local_face_count
         << "; success=" << (window.success ? "true" : "false");
  if (!window.failure_reason.empty()) {
    report << "; failure_reason=" << window.failure_reason;
  }
  window.report_line = report.str();

  return window;
}

dec3d::core::MeshUpdateProposal BuildRadialAleMeshUpdateProposal(
    const dec3d::core::Array3D<double>& rho,
    const dec3d::core::Array3D<double>& mom_r,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    double beta) noexcept {
  dec3d::core::MeshUpdateProposal proposal;
  if (rho.empty() ||
      mom_r.empty() ||
      !geometry.is_valid() ||
      dt_s <= 0.0 ||
      !std::isfinite(dt_s) ||
      !std::isfinite(beta)) {
    proposal.summary = "radial ALE proposal inputs are incomplete";
    return proposal;
  }

  const auto shell_velocities = BuildRadialAleShellVelocitySlice(rho, mom_r);
  if (shell_velocities.empty() ||
      shell_velocities.size() + 1u != geometry.radial_faces.size()) {
    proposal.summary = "radial ALE proposal inputs are incomplete";
    return proposal;
  }

  return BuildRadialAleProposalFromGlobalShellVelocity(
      shell_velocities,
      geometry.radial_faces,
      dt_s,
      beta);
}

RadialAleCommitResult ApplyRadialAleMeshUpdateProposal(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  return ApplyRadialAleMeshUpdateProposal(proposal, geometry, 0u);
}

RadialAleCommitResult ApplyRadialAleMeshUpdateProposal(
    const dec3d::core::MeshUpdateProposal& proposal,
    dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t global_face_begin) noexcept {
  RadialAleCommitResult result;
  result.attempted = true;
  result.proposal_present = proposal.requested;
  const std::size_t local_face_count = geometry.radial_faces.size();
  const auto window = BuildRadialAleLocalProposalWindow(
      proposal,
      global_face_begin,
      local_face_count);

  if (!proposal.requested) {
    result.failure_reason = "mesh update proposal was not requested";
  } else if (!geometry.is_valid()) {
    result.failure_reason = "canonical mesh geometry is invalid";
  } else if (!window.is_complete()) {
    result.failure_reason = window.failure_reason.empty()
                                ? "mesh update proposal local window is incomplete"
                                : window.failure_reason;
  } else {
    bool monotonic = true;
    for (std::size_t radial_face = 1u; radial_face < window.proposed_radial_faces.size(); ++radial_face) {
      if (!(window.proposed_radial_faces[radial_face] >
            window.proposed_radial_faces[radial_face - 1u])) {
        monotonic = false;
        break;
      }
    }

    if (!monotonic) {
      result.failure_reason = "radial ALE proposal would invert radial face ordering";
    } else {
      result.geometry_changed = window.proposed_radial_faces != geometry.radial_faces;
      geometry.radial_faces = window.proposed_radial_faces;
      geometry.valid = true;
      RebuildCellVolumes(geometry);
      result.success = true;
    }
  }

  std::ostringstream report;
  report << "proposal_present=" << (result.proposal_present ? "true" : "false")
         << "; global_face_begin=" << global_face_begin
         << "; local_face_count=" << local_face_count
         << "; geometry_changed=" << (result.geometry_changed ? "true" : "false")
         << "; success=" << (result.success ? "true" : "false");
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

}  // namespace dec3d::mesh
