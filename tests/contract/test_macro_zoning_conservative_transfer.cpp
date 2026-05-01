#include "core/array/array3d.hpp"
#include "mesh/macro_zoning/macro_zoning.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"

#include <array>
#include <cmath>
#include <iostream>

namespace {

double SumExtensive(
    const dec3d::core::Array3D<double>& field,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double sum = 0.0;
  std::size_t linear_index = 0;
  for (std::size_t radial = 0; radial < field.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < field.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < field.extent_phi(); ++phi, ++linear_index) {
        sum += field(radial, theta, phi) * geometry.cell_volumes.at(linear_index);
      }
    }
  }
  return sum;
}

}  // namespace

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::DetectMacroZones;
    using dec3d::mesh::FineHydroPackageView;
    using dec3d::mesh::ProlongMacroZoneHydroPackage;
    using dec3d::mesh::RestrictFineHydroPackage;
    using dec3d::mesh::SphericalMeshDescriptor;

    constexpr std::size_t kRadialCells = 2u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.25});
    DEC3D_CHECK(geometry.is_valid());

    Array3D<double> rho(kRadialCells, kThetaCells, kPhiCells, 0.0);
    Array3D<double> mom_r(kRadialCells, kThetaCells, kPhiCells, 0.0);
    Array3D<double> mom_theta(kRadialCells, kThetaCells, kPhiCells, 0.0);
    Array3D<double> mom_phi(kRadialCells, kThetaCells, kPhiCells, 0.0);
    Array3D<double> e_fluid_total(kRadialCells, kThetaCells, kPhiCells, 0.0);
    Array3D<double> chi_e(kRadialCells, kThetaCells, kPhiCells, 0.0);

    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          const double base =
              1.0 + 0.2 * static_cast<double>(radial) +
              0.03 * static_cast<double>(theta) +
              0.005 * static_cast<double>(phi);
          rho(radial, theta, phi) = base;
          mom_r(radial, theta, phi) = 0.5 * base;
          mom_theta(radial, theta, phi) = -0.25 * base;
          mom_phi(radial, theta, phi) = 0.1 * base;
          e_fluid_total(radial, theta, phi) = 2.0 * base;
          chi_e(radial, theta, phi) = 0.75 * base;
        }
      }
    }

    const auto map = DetectMacroZones(geometry);
    DEC3D_CHECK(map.is_complete());

    const FineHydroPackageView view{
        &rho, &mom_r, &mom_theta, &mom_phi, &e_fluid_total, &chi_e};
    DEC3D_CHECK(view.is_complete());

    const auto coarse = RestrictFineHydroPackage(view, geometry, map);
    DEC3D_CHECK(coarse.is_complete());

    const auto prolonged = ProlongMacroZoneHydroPackage(coarse);
    DEC3D_CHECK(prolonged.is_complete());

    DEC3D_CHECK(std::abs(SumExtensive(rho, geometry) - SumExtensive(prolonged.rho, geometry)) < 1.0e-12);
    DEC3D_CHECK(std::abs(SumExtensive(mom_r, geometry) - SumExtensive(prolonged.mom_r, geometry)) < 1.0e-12);
    DEC3D_CHECK(
        std::abs(SumExtensive(mom_theta, geometry) - SumExtensive(prolonged.mom_theta, geometry)) < 1.0e-12);
    DEC3D_CHECK(std::abs(SumExtensive(mom_phi, geometry) - SumExtensive(prolonged.mom_phi, geometry)) < 1.0e-12);
    DEC3D_CHECK(
        std::abs(SumExtensive(e_fluid_total, geometry) - SumExtensive(prolonged.e_fluid_total, geometry)) <
        1.0e-12);
    DEC3D_CHECK(std::abs(SumExtensive(chi_e, geometry) - SumExtensive(prolonged.chi_e, geometry)) < 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
