#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::HydroConservativeState;
    using dec3d::hydro::RecoverPrimitiveState;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroStateView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto state = CanonicalState::Create(CanonicalStateLayout{1, 6, 4, 0});
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(0, theta, phi) = 1.0 + 0.02 * static_cast<double>(theta == 0 ? 1 : 0);
        state.mom_r(0, theta, phi) = 0.01 * static_cast<double>(phi);
        state.mom_theta(0, theta, phi) = (theta == 0 ? 0.08 : 0.01) * (phi % 2 == 0 ? 1.0 : -1.0);
        state.mom_phi(0, theta, phi) = (theta == 0 ? -0.06 : 0.02) * static_cast<double>(phi + 1u);
        state.e_fluid_total(0, theta, phi) = 3.0 + 0.05 * static_cast<double>(theta + phi);
        state.e_electron(0, theta, phi) = ElectronEnergyDensityFromPressure(0.5 + 0.01 * static_cast<double>(theta));
      }
    }

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{1, 6, 4, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto hydro_view = BuildHydroStateView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const auto result = AdvanceStaticGridHydro(hydro_view, geometry, 5.0e-4);
    DEC3D_CHECK(result.is_complete());

    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const HydroConservativeState cell{
            (*hydro_view.rho)(0, theta, phi),
            (*hydro_view.mom_r)(0, theta, phi),
            (*hydro_view.mom_theta)(0, theta, phi),
            (*hydro_view.mom_phi)(0, theta, phi),
            (*hydro_view.e_fluid_total)(0, theta, phi),
            hydro_view.chi_e(0, theta, phi)};
        const auto primitive = RecoverPrimitiveState(cell);
        DEC3D_CHECK(cell.is_finite());
        DEC3D_CHECK(primitive.is_physical());
        DEC3D_CHECK(std::abs(cell.mom_theta) < 1.0);
        DEC3D_CHECK(std::abs(cell.mom_phi) < 1.0);
      }
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
