#include "core/array/array3d.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState MakeState(double seed) noexcept {
  const dec3d::hydro::HydroPrimitiveState primitive{
      1.0 + 0.1 * seed,
      0.05 * seed,
      -0.02 * seed,
      0.03 * seed,
      1.0 + 0.05 * seed,
      0.6 + 0.01 * seed};
  return dec3d::hydro::MakeConservativeState(primitive);
}

}  // namespace

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::hydro::HydroDirection;
    using dec3d::hydro::PrepareDirectionalGhosts;

    Array3D<dec3d::hydro::HydroConservativeState> interior(3, 3, 4);
    double seed = 0.0;
    for (std::size_t radial = 0; radial < interior.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < interior.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < interior.extent_phi(); ++phi) {
          interior(radial, theta, phi) = MakeState(seed);
          seed += 1.0;
        }
      }
    }

    const auto radial = PrepareDirectionalGhosts(interior, HydroDirection::radial, 3);
    DEC3D_CHECK(radial.is_complete(3, 3, 4, 3));
    DEC3D_CHECK_EQ(radial.ghosted_states.extent_r(), 9u);
    DEC3D_CHECK_EQ(radial.ghosted_states.extent_theta(), 3u);
    DEC3D_CHECK_EQ(radial.ghosted_states.extent_phi(), 4u);
    DEC3D_CHECK(radial.report_line.find("direction=radial") != std::string::npos);

    const auto theta = PrepareDirectionalGhosts(interior, HydroDirection::theta, 3);
    DEC3D_CHECK(theta.is_complete(3, 3, 4, 3));
    DEC3D_CHECK_EQ(theta.ghosted_states.extent_r(), 3u);
    DEC3D_CHECK_EQ(theta.ghosted_states.extent_theta(), 9u);
    DEC3D_CHECK_EQ(theta.ghosted_states.extent_phi(), 4u);
    DEC3D_CHECK(theta.report_line.find("direction=theta") != std::string::npos);

    const auto phi = PrepareDirectionalGhosts(interior, HydroDirection::phi, 3);
    DEC3D_CHECK(phi.is_complete(3, 3, 4, 3));
    DEC3D_CHECK_EQ(phi.ghosted_states.extent_r(), 3u);
    DEC3D_CHECK_EQ(phi.ghosted_states.extent_theta(), 3u);
    DEC3D_CHECK_EQ(phi.ghosted_states.extent_phi(), 10u);
    DEC3D_CHECK(phi.report_line.find("direction=phi") != std::string::npos);

    const auto invalid = PrepareDirectionalGhosts(interior, HydroDirection::phi, 0);
    DEC3D_CHECK(!invalid.success);
    DEC3D_CHECK(!invalid.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
