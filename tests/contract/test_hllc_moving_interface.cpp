#include "hydro/riemann/hllc_solver.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState State(
    double rho,
    double v_r,
    double pressure) {
  return dec3d::hydro::MakeConservativeState(
      dec3d::hydro::HydroPrimitiveState{
          rho,
          v_r,
          0.02,
          -0.01,
          pressure,
          std::pow(pressure * 0.4, 3.0 / 5.0)});
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::HllcActiveRegion;
    using dec3d::hydro::SolveHllcRiemann;
    using dec3d::hydro::SolveMovingInterfaceHllcFlux;

    const auto left = State(1.0, -0.35, 1.0);
    const auto right = State(0.25, 0.45, 0.2);
    const auto static_hllc = SolveHllcRiemann(left, right);
    DEC3D_CHECK(static_hllc.is_complete());

    const double w_face =
        0.5 * (static_hllc.waves.s_star + static_hllc.waves.s_right);
    const auto moving = SolveMovingInterfaceHllcFlux(left, right, w_face);

    DEC3D_CHECK(moving.is_complete());
    DEC3D_CHECK_EQ(moving.mode, std::string("moving_interface_hllc"));
    DEC3D_CHECK_EQ(moving.remap_order, std::string("none"));
    DEC3D_CHECK(moving.active_region == HllcActiveRegion::right_star);
    DEC3D_CHECK(moving.static_zero_region != moving.active_region);
    DEC3D_CHECK(moving.static_zero_region != HllcActiveRegion::invalid);
    DEC3D_CHECK(std::isfinite(moving.flux.rho));
    DEC3D_CHECK(std::abs(moving.flux.rho - static_hllc.interface_flux.rho) > 1.0e-10);

    const auto zero = SolveMovingInterfaceHllcFlux(left, right, 0.0);
    DEC3D_CHECK(zero.is_complete());
    DEC3D_CHECK(zero.active_region == static_hllc.active_region);
    DEC3D_CHECK(std::abs(zero.flux.rho - static_hllc.interface_flux.rho) < 1.0e-12);
    DEC3D_CHECK(std::abs(zero.flux.mom_r - static_hllc.interface_flux.mom_r) < 1.0e-12);
    DEC3D_CHECK(std::abs(zero.flux.e_fluid_total -
                         static_hllc.interface_flux.e_fluid_total) < 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
