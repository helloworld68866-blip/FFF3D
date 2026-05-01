#include "core/array/array3d.hpp"
#include "mesh/macro_zoning/macro_zoning.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::DetectMacroZones;
    using dec3d::mesh::FineHydroPackageView;
    using dec3d::mesh::RejectMacroZoneAuthoritativeCommit;
    using dec3d::mesh::RestrictFineHydroPackage;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ChiEFromElectronPressure;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    constexpr std::size_t kRadialCells = 1u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;

    auto state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    Array3D<double> chi_e(kRadialCells, kThetaCells, kPhiCells, 0.0);

    for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
      for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
        const double density = 1.0 + 0.02 * static_cast<double>(theta) + 0.01 * static_cast<double>(phi);
        const double electron_pressure = 0.2 + 0.01 * static_cast<double>(theta);
        state.rho(0u, theta, phi) = density;
        state.mom_r(0u, theta, phi) = 0.1 * density;
        state.mom_theta(0u, theta, phi) = -0.03 * density;
        state.mom_phi(0u, theta, phi) = 0.02 * density;
        state.e_fluid_total(0u, theta, phi) = 2.0 * density;
        state.e_electron(0u, theta, phi) = ElectronEnergyDensityFromPressure(electron_pressure);
        chi_e(0u, theta, phi) = ChiEFromElectronPressure(electron_pressure);
      }
    }

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.15});
    DEC3D_CHECK(geometry.is_valid());

    const auto map = DetectMacroZones(geometry);
    DEC3D_CHECK(map.is_complete());

    const FineHydroPackageView fine_view{
        &state.rho, &state.mom_r, &state.mom_theta, &state.mom_phi, &state.e_fluid_total, &chi_e};
    DEC3D_CHECK(fine_view.is_complete());

    auto coarse = RestrictFineHydroPackage(fine_view, geometry, map);
    DEC3D_CHECK(coarse.is_complete());
    DEC3D_CHECK(coarse.short_lived_work_view);
    DEC3D_CHECK(!coarse.is_authoritative_truth());

    const double original_rho = state.rho(0u, 0u, 0u);
    coarse.rho.front() += 100.0;
    DEC3D_CHECK(std::abs(state.rho(0u, 0u, 0u) - original_rho) < 1.0e-14);

    const auto reject = RejectMacroZoneAuthoritativeCommit(coarse);
    DEC3D_CHECK(!reject.success);
    DEC3D_CHECK(!reject.failure_reason.empty());
    DEC3D_CHECK(std::abs(state.rho(0u, 0u, 0u) - original_rho) < 1.0e-14);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
