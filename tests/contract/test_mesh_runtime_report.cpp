#include "mesh/ghost/ghost_topology.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    dec3d::mesh::SphericalMeshDescriptor descriptor;
    descriptor.radial_cells = 6;
    descriptor.theta_cells = 4;
    descriptor.phi_cells = 8;
    descriptor.inner_radius = 0.5;
    descriptor.outer_radius = 2.0;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(descriptor);
    const auto evidence = dec3d::mesh::BuildMeshSubstrateEvidence(descriptor);
    const auto ownership = dec3d::mesh::BuildRadialOwnership(descriptor.radial_cells, 3);
    const auto topology = dec3d::mesh::BuildGhostTopology(ownership.slices[1], 1, ownership.slices.size());

    const auto report = dec3d::runtime::BuildMeshSubstrateRuntimeReport(
        geometry,
        ownership,
        topology,
        evidence);

    DEC3D_CHECK(report.mesh_geometry_built);
    DEC3D_CHECK(report.ownership_built);
    DEC3D_CHECK(report.ghost_topology_initialized);
    DEC3D_CHECK(report.is_complete());
    DEC3D_CHECK_EQ(report.implementation_id, std::string("p0.mesh.substrate"));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
