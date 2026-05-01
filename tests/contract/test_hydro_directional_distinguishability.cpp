#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

void SeedAngularStructure(dec3d::state::CanonicalState& state) {
  for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
    for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
      state.rho(0, theta, phi) = 1.0;
      state.mom_r(0, theta, phi) = 0.0;
      state.mom_theta(0, theta, phi) = 0.05 * static_cast<double>(theta);
      state.mom_phi(0, theta, phi) = -0.04 * static_cast<double>(phi);
      state.e_fluid_total(0, theta, phi) = 2.5 + 0.10 * static_cast<double>(theta + phi);
      state.e_electron(0, theta, phi) = 0.6;
    }
  }
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroStateView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto radial_only_state = CanonicalState::Create(CanonicalStateLayout{1, 3, 4, 0});
    auto full_state = CanonicalState::Create(CanonicalStateLayout{1, 3, 4, 0});
    SeedAngularStructure(radial_only_state);
    SeedAngularStructure(full_state);

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{1, 3, 4, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto radial_view = BuildHydroStateView(radial_only_state);
    auto full_view = BuildHydroStateView(full_state);
    DEC3D_CHECK(radial_view.is_complete());
    DEC3D_CHECK(full_view.is_complete());

    const double initial_radial_rho = (*radial_view.rho)(0, 1, 1);
    const double initial_full_rho = (*full_view.rho)(0, 1, 1);

    const auto radial_result = AdvanceStaticGridHydro(
        radial_view,
        geometry,
        1.0e-3,
        StaticGridHydroOptions{false, true, false, false});
    const auto full_result = AdvanceStaticGridHydro(
        full_view,
        geometry,
        1.0e-3,
        StaticGridHydroOptions{false, true, true, true});

    DEC3D_CHECK(radial_result.is_complete());
    DEC3D_CHECK(full_result.is_complete());
    DEC3D_CHECK(std::abs((*radial_view.rho)(0, 1, 1) - initial_radial_rho) < 1.0e-12);
    DEC3D_CHECK(std::abs((*full_view.rho)(0, 1, 1) - initial_full_rho) > 1.0e-8);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
