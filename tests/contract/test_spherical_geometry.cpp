#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <numbers>

namespace {

bool IsStrictlyIncreasing(const std::vector<double>& values) {
  for (std::size_t index = 1; index < values.size(); ++index) {
    if (!(values[index] > values[index - 1])) {
      return false;
    }
  }

  return true;
}

double Sum(const std::vector<double>& values) {
  double total = 0.0;

  for (const double value : values) {
    total += value;
  }

  return total;
}

double OuterRadialArea(const dec3d::mesh::SphericalGeometryMetadata& geometry,
                       std::size_t theta,
                       std::size_t phi) {
  const double radius = geometry.radial_faces.back();
  const double polar_factor =
      std::cos(geometry.theta_faces[theta]) - std::cos(geometry.theta_faces[theta + 1u]);
  const double delta_phi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return radius * radius * polar_factor * delta_phi;
}

}  // namespace

int main() {
  try {
    dec3d::mesh::SphericalMeshDescriptor descriptor;
    descriptor.radial_cells = 4;
    descriptor.theta_cells = 3;
    descriptor.phi_cells = 5;
    descriptor.inner_radius = 0.0;
    descriptor.outer_radius = 2.0;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(descriptor);
    const auto evidence = dec3d::mesh::BuildMeshSubstrateEvidence(descriptor);

    DEC3D_CHECK(geometry.is_valid());
    DEC3D_CHECK_EQ(geometry.radial_faces.size(), descriptor.radial_cells + 1);
    DEC3D_CHECK_EQ(geometry.theta_faces.size(), descriptor.theta_cells + 1);
    DEC3D_CHECK_EQ(geometry.phi_faces.size(), descriptor.phi_cells + 1);
    DEC3D_CHECK_EQ(
        geometry.cell_volumes.size(),
        descriptor.radial_cells * descriptor.theta_cells * descriptor.phi_cells);

    DEC3D_CHECK(IsStrictlyIncreasing(geometry.radial_faces));
    DEC3D_CHECK(IsStrictlyIncreasing(geometry.theta_faces));
    DEC3D_CHECK(IsStrictlyIncreasing(geometry.phi_faces));

    for (const double volume : geometry.cell_volumes) {
      DEC3D_CHECK(volume > 0.0);
    }

    const double analytic_volume =
        (4.0 / 3.0) * std::numbers::pi * std::pow(descriptor.outer_radius, 3);
    const double summed_volume = Sum(geometry.cell_volumes);
    DEC3D_CHECK(std::abs(geometry.global_volume - analytic_volume) < 1.0e-10);
    DEC3D_CHECK(std::abs(summed_volume - geometry.global_volume) < 1.0e-10);

    DEC3D_CHECK(evidence.is_present());
    DEC3D_CHECK_EQ(evidence.implementation_id, std::string("p0.mesh.substrate"));

    dec3d::mesh::SphericalMeshDescriptor shell_descriptor;
    shell_descriptor.radial_cells = 3;
    shell_descriptor.theta_cells = 4;
    shell_descriptor.phi_cells = 6;
    shell_descriptor.inner_radius = 1.0;
    shell_descriptor.outer_radius = 2.5;

    const auto shell_geometry = dec3d::mesh::BuildSphericalGeometry(shell_descriptor);
    DEC3D_CHECK(shell_geometry.is_valid());
    DEC3D_CHECK(IsStrictlyIncreasing(shell_geometry.radial_faces));

    for (const double volume : shell_geometry.cell_volumes) {
      DEC3D_CHECK(volume > 0.0);
    }

    const double shell_analytic_volume =
        (4.0 / 3.0) * std::numbers::pi *
        (std::pow(shell_descriptor.outer_radius, 3) - std::pow(shell_descriptor.inner_radius, 3));
    DEC3D_CHECK(std::abs(shell_geometry.global_volume - shell_analytic_volume) < 1.0e-10);
    DEC3D_CHECK(std::abs(Sum(shell_geometry.cell_volumes) - shell_geometry.global_volume) < 1.0e-10);

    dec3d::mesh::SphericalMeshDescriptor axisym_descriptor;
    axisym_descriptor.radial_cells = 4;
    axisym_descriptor.theta_cells = 8;
    axisym_descriptor.phi_cells = 1;
    axisym_descriptor.inner_radius = 0.2;
    axisym_descriptor.outer_radius = 1.7;

    const auto axisym_geometry = dec3d::mesh::BuildSphericalGeometry(axisym_descriptor);
    DEC3D_CHECK(axisym_geometry.is_valid());
    DEC3D_CHECK_EQ(axisym_geometry.phi_faces.size(), std::size_t{2});
    DEC3D_CHECK(std::abs(axisym_geometry.phi_faces.front() - 0.0) < 1.0e-15);
    DEC3D_CHECK(std::abs(axisym_geometry.phi_faces.back() - 2.0 * std::numbers::pi) < 1.0e-15);

    const double axisym_expected_volume =
        (4.0 / 3.0) * std::numbers::pi *
        (std::pow(axisym_descriptor.outer_radius, 3) -
         std::pow(axisym_descriptor.inner_radius, 3));
    DEC3D_CHECK(std::abs(Sum(axisym_geometry.cell_volumes) - axisym_expected_volume) <
                1.0e-13);

    double axisym_outer_area = 0.0;
    for (std::size_t theta = 0; theta < axisym_descriptor.theta_cells; ++theta) {
      axisym_outer_area += OuterRadialArea(axisym_geometry, theta, 0u);
    }
    const double axisym_expected_outer_area =
        4.0 * std::numbers::pi * axisym_descriptor.outer_radius *
        axisym_descriptor.outer_radius;
    DEC3D_CHECK(std::abs(axisym_outer_area - axisym_expected_outer_area) < 1.0e-13);

    dec3d::mesh::SphericalMeshDescriptor invalid_descriptor;
    invalid_descriptor.radial_cells = 2;
    invalid_descriptor.theta_cells = 2;
    invalid_descriptor.phi_cells = 2;
    invalid_descriptor.inner_radius = 2.0;
    invalid_descriptor.outer_radius = 1.0;

    const auto invalid_geometry = dec3d::mesh::BuildSphericalGeometry(invalid_descriptor);
    const auto invalid_evidence = dec3d::mesh::BuildMeshSubstrateEvidence(invalid_descriptor);
    DEC3D_CHECK(!invalid_geometry.is_valid());
    DEC3D_CHECK(invalid_geometry.radial_faces.empty());
    DEC3D_CHECK(invalid_geometry.theta_faces.empty());
    DEC3D_CHECK(invalid_geometry.phi_faces.empty());
    DEC3D_CHECK(invalid_geometry.cell_volumes.empty());
    DEC3D_CHECK(!invalid_evidence.is_present());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
