#include "core/array/array3d.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState MakeState(double seed) noexcept {
  return {
      10.0 + seed,
      20.0 + seed,
      30.0 + seed,
      40.0 + seed,
      50.0 + seed,
      60.0 + seed};
}

void CheckStateEquals(
    const dec3d::hydro::HydroConservativeState& actual,
    const dec3d::hydro::HydroConservativeState& expected) {
  DEC3D_CHECK_EQ(actual.rho, expected.rho);
  DEC3D_CHECK_EQ(actual.mom_r, expected.mom_r);
  DEC3D_CHECK_EQ(actual.mom_theta, expected.mom_theta);
  DEC3D_CHECK_EQ(actual.mom_phi, expected.mom_phi);
  DEC3D_CHECK_EQ(actual.e_fluid_total, expected.e_fluid_total);
  DEC3D_CHECK_EQ(actual.chi_e, expected.chi_e);
}

}  // namespace

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::hydro::FillHydroPhiPeriodicGhosts;
    using dec3d::hydro::HydroDirection;

    Array3D<dec3d::hydro::HydroConservativeState> interior(1, 1, 5);
    for (std::size_t phi = 0; phi < interior.extent_phi(); ++phi) {
      interior(0, 0, phi) = MakeState(static_cast<double>(phi));
    }

    const auto prepared = FillHydroPhiPeriodicGhosts(interior, 3);
    DEC3D_CHECK(prepared.is_complete(1, 1, 5, 3));
    DEC3D_CHECK(prepared.direction == HydroDirection::phi);
    DEC3D_CHECK_EQ(prepared.ghost_layers, 3u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_r(), 1u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_theta(), 1u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_phi(), 11u);

    CheckStateEquals(prepared.ghosted_states(0, 0, 0), interior(0, 0, 2));
    CheckStateEquals(prepared.ghosted_states(0, 0, 1), interior(0, 0, 3));
    CheckStateEquals(prepared.ghosted_states(0, 0, 2), interior(0, 0, 4));

    for (std::size_t phi = 0; phi < interior.extent_phi(); ++phi) {
      CheckStateEquals(prepared.ghosted_states(0, 0, phi + 3u), interior(0, 0, phi));
    }

    CheckStateEquals(prepared.ghosted_states(0, 0, 8), interior(0, 0, 0));
    CheckStateEquals(prepared.ghosted_states(0, 0, 9), interior(0, 0, 1));
    CheckStateEquals(prepared.ghosted_states(0, 0, 10), interior(0, 0, 2));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
