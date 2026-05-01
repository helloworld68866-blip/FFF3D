#include "mesh/macro_zoning/macro_zoning.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"

#include <iostream>
#include <string>

int main() {
  try {
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::DetectMacroZones;
    using dec3d::mesh::SphericalMeshDescriptor;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{1u, 8u, 8u, 0.05, 0.15});
    DEC3D_CHECK(geometry.is_valid());

    const auto map = DetectMacroZones(geometry);
    DEC3D_CHECK(map.is_complete());
    DEC3D_CHECK(map.phi_factor_varies_by_theta_band);
    DEC3D_CHECK_EQ(map.theta_bands.size(), static_cast<std::size_t>(4));

    DEC3D_CHECK_EQ(map.theta_bands[0].theta_factor, static_cast<std::size_t>(2));
    DEC3D_CHECK_EQ(map.theta_bands[1].theta_factor, static_cast<std::size_t>(2));
    DEC3D_CHECK(map.theta_bands[0].phi_factor > map.theta_bands[1].phi_factor);
    DEC3D_CHECK(map.theta_bands[0].phi_factor >= static_cast<std::size_t>(2));
    DEC3D_CHECK_EQ(map.coarse_cells.size(), static_cast<std::size_t>(20));

    for (std::size_t theta = 0; theta < 8u; ++theta) {
      for (std::size_t phi = 0; phi < 8u; ++phi) {
        DEC3D_CHECK(
            map.fine_to_coarse(0u, theta, phi) < map.coarse_cells.size());
      }
    }

    const auto full_sphere_origin_geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{1u, 16u, 16u, 0.0, 5.0e-3});
    DEC3D_CHECK(full_sphere_origin_geometry.is_valid());

    const auto origin_map = DetectMacroZones(full_sphere_origin_geometry);
    DEC3D_CHECK(origin_map.is_complete());
    DEC3D_CHECK(origin_map.phi_factor_saturated_count > 0u);
    DEC3D_CHECK(
        origin_map.report_line.find("phi_factor_saturated_count=") !=
        std::string::npos);
    for (const auto& band : origin_map.theta_bands) {
      if (band.min_delta_s_phi * static_cast<double>(band.phi_factor) +
              1.0e-12 <
          0.5 * (full_sphere_origin_geometry.radial_faces[1] -
                 full_sphere_origin_geometry.radial_faces[0])) {
        DEC3D_CHECK_EQ(band.phi_factor, origin_map.fine_phi_cells);
      }
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
