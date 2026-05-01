#include "hydro/driver/radial_boundary_contract.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

int main() {
  try {
    using dec3d::hydro::BoundaryFace;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::hydro::PrepareRadialBoundaryStates;
    using dec3d::hydro::RecoverPrimitiveState;

    const HydroPrimitiveState lower_shell_state{
        1.0,
        0.25,
        0.10,
        -0.15,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const auto lower_prepared =
        PrepareRadialBoundaryStates(BoundaryFace::lower, MakeConservativeState(lower_shell_state));
    DEC3D_CHECK(lower_prepared.is_complete());
    const auto lower_ghost = RecoverPrimitiveState(lower_prepared.left_state);
    DEC3D_CHECK(std::abs(lower_ghost.v_r - lower_shell_state.v_r) < 1.0e-12);
    DEC3D_CHECK(std::abs(lower_ghost.v_theta - lower_shell_state.v_theta) < 1.0e-12);
    DEC3D_CHECK(std::abs(lower_ghost.v_phi - lower_shell_state.v_phi) < 1.0e-12);

    const HydroPrimitiveState outer_outflow_state{
        1.0,
        0.30,
        -0.05,
        0.08,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const auto outflow_prepared =
        PrepareRadialBoundaryStates(BoundaryFace::upper, MakeConservativeState(outer_outflow_state));
    DEC3D_CHECK(outflow_prepared.is_complete());
    const auto outflow_ghost = RecoverPrimitiveState(outflow_prepared.right_state);
    DEC3D_CHECK(std::abs(outflow_ghost.v_r - outer_outflow_state.v_r) < 1.0e-12);

    const HydroPrimitiveState outer_inflow_state{
        1.0,
        -0.30,
        -0.05,
        0.08,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const auto inflow_prepared =
        PrepareRadialBoundaryStates(BoundaryFace::upper, MakeConservativeState(outer_inflow_state));
    DEC3D_CHECK(inflow_prepared.is_complete());
    const auto inflow_ghost = RecoverPrimitiveState(inflow_prepared.right_state);
    DEC3D_CHECK(inflow_ghost.v_r >= -1.0e-12);
    DEC3D_CHECK(std::abs(inflow_ghost.v_theta - outer_inflow_state.v_theta) < 1.0e-12);
    DEC3D_CHECK(std::abs(inflow_ghost.v_phi - outer_inflow_state.v_phi) < 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
