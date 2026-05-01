#include "core/array/array3d.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/driver/theta_pole_boundary.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState MakeState(
    std::size_t theta,
    std::size_t phi) noexcept {
  const double seed = static_cast<double>((theta + 1u) * 10u + phi);
  return {
      1.0 + seed,
      2.0 + seed,
      3.0 + seed,
      4.0 + seed,
      5.0 + seed,
      6.0 + seed};
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
    using dec3d::hydro::BuildThetaPoleGhostState;
    using dec3d::hydro::FillHydroThetaPoleGhosts;
    using dec3d::hydro::HydroDirection;
    using dec3d::hydro::MapPhiAcrossPole;

    Array3D<dec3d::hydro::HydroConservativeState> interior(1, 3, 4);
    for (std::size_t theta = 0; theta < interior.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < interior.extent_phi(); ++phi) {
        interior(0, theta, phi) = MakeState(theta, phi);
      }
    }

    const auto prepared = FillHydroThetaPoleGhosts(interior, 3);
    DEC3D_CHECK(prepared.is_complete(1, 3, 4, 3));
    DEC3D_CHECK(prepared.direction == HydroDirection::theta);
    DEC3D_CHECK_EQ(prepared.ghost_layers, 3u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_r(), 1u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_theta(), 9u);
    DEC3D_CHECK_EQ(prepared.ghosted_states.extent_phi(), 4u);

    for (std::size_t phi = 0; phi < interior.extent_phi(); ++phi) {
      const auto mapped_phi = MapPhiAcrossPole(phi, interior.extent_phi());
      for (std::size_t g = 1; g <= 3u; ++g) {
        const auto expected_lower =
            BuildThetaPoleGhostState(interior(0, g - 1u, mapped_phi));
        const auto expected_upper =
            BuildThetaPoleGhostState(interior(0, interior.extent_theta() - g, mapped_phi));

        CheckStateEquals(prepared.ghosted_states(0, 3u - g, phi), expected_lower);
        CheckStateEquals(
            prepared.ghosted_states(0, interior.extent_theta() + 2u + g, phi),
            expected_upper);
      }
    }

    Array3D<dec3d::hydro::HydroConservativeState> odd_phi(1, 2, 3);
    odd_phi(0, 0, 0) = MakeState(0, 0);
    const auto invalid_prepared = FillHydroThetaPoleGhosts(odd_phi, 3);
    DEC3D_CHECK(!invalid_prepared.success);
    DEC3D_CHECK(!invalid_prepared.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
