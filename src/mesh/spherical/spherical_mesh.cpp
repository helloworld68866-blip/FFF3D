#include "mesh/spherical/spherical_mesh.hpp"

#include <cmath>
#include <numbers>

namespace dec3d::mesh {

[[maybe_unused]] static constexpr const char* kMeshLayerRole = "generic.mesh.infrastructure";

SphericalGeometryMetadata BuildSphericalGeometry(const SphericalMeshDescriptor& descriptor) noexcept {
  SphericalGeometryMetadata metadata;

  if (!descriptor.is_valid()) {
    return metadata;
  }

  metadata.valid = true;

  metadata.radial_faces.resize(descriptor.radial_cells + 1, 0.0);
  metadata.theta_faces.resize(descriptor.theta_cells + 1, 0.0);
  metadata.phi_faces.resize(descriptor.phi_cells + 1, 0.0);
  metadata.cell_volumes.resize(descriptor.radial_cells * descriptor.theta_cells * descriptor.phi_cells, 0.0);

  const double radial_spacing =
      (descriptor.outer_radius - descriptor.inner_radius) / static_cast<double>(descriptor.radial_cells);
  const double theta_spacing = std::numbers::pi / static_cast<double>(descriptor.theta_cells);
  const double phi_spacing = (2.0 * std::numbers::pi) / static_cast<double>(descriptor.phi_cells);

  for (std::size_t radial_index = 0; radial_index < metadata.radial_faces.size(); ++radial_index) {
    metadata.radial_faces[radial_index] =
        descriptor.inner_radius + radial_spacing * static_cast<double>(radial_index);
  }

  for (std::size_t theta_index = 0; theta_index < metadata.theta_faces.size(); ++theta_index) {
    metadata.theta_faces[theta_index] = theta_spacing * static_cast<double>(theta_index);
  }

  for (std::size_t phi_index = 0; phi_index < metadata.phi_faces.size(); ++phi_index) {
    metadata.phi_faces[phi_index] = phi_spacing * static_cast<double>(phi_index);
  }

  std::size_t cell_index = 0;
  for (std::size_t radial_index = 0; radial_index < descriptor.radial_cells; ++radial_index) {
    const double radial_inner = metadata.radial_faces[radial_index];
    const double radial_outer = metadata.radial_faces[radial_index + 1];
    const double radial_factor = (std::pow(radial_outer, 3) - std::pow(radial_inner, 3)) / 3.0;

    for (std::size_t theta_index = 0; theta_index < descriptor.theta_cells; ++theta_index) {
      const double theta_lower = metadata.theta_faces[theta_index];
      const double theta_upper = metadata.theta_faces[theta_index + 1];
      const double polar_factor = std::cos(theta_lower) - std::cos(theta_upper);

      for (std::size_t phi_index = 0; phi_index < descriptor.phi_cells; ++phi_index) {
        const double phi_lower = metadata.phi_faces[phi_index];
        const double phi_upper = metadata.phi_faces[phi_index + 1];
        const double azimuthal_factor = phi_upper - phi_lower;

        metadata.cell_volumes[cell_index] = radial_factor * polar_factor * azimuthal_factor;
        ++cell_index;
      }
    }
  }

  metadata.global_volume =
      (4.0 / 3.0) * std::numbers::pi *
      (std::pow(descriptor.outer_radius, 3) - std::pow(descriptor.inner_radius, 3));

  return metadata;
}

MeshBuildEvidence BuildMeshSubstrateEvidence(const SphericalMeshDescriptor& descriptor) noexcept {
  if (!descriptor.is_valid()) {
    return MeshBuildEvidence{};
  }

  MeshBuildEvidence evidence;
  evidence.built_by_new_code = true;
  evidence.implementation_id = "p0.mesh.substrate";
  return evidence;
}

}  // namespace dec3d::mesh
