#include "hydro/driver/theta_pole_boundary.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::hydro::BoundaryFace;
    using dec3d::hydro::BuildThetaPoleGhostState;
    using dec3d::hydro::HydroConservativeState;
    using dec3d::hydro::MapPhiAcrossPole;
    using dec3d::hydro::PreparePhiPeriodicBoundaryStates;
    using dec3d::hydro::PrepareThetaPoleBoundaryStates;

    DEC3D_CHECK_EQ(MapPhiAcrossPole(0, 4), 2u);
    DEC3D_CHECK_EQ(MapPhiAcrossPole(1, 4), 3u);
    DEC3D_CHECK_EQ(MapPhiAcrossPole(2, 4), 0u);
    DEC3D_CHECK_EQ(MapPhiAcrossPole(3, 4), 1u);
    DEC3D_CHECK_EQ(MapPhiAcrossPole(0, 1), 0u);

    const HydroConservativeState interior{
        1.0,
        2.0,
        3.0,
        -4.0,
        5.0,
        0.6};
    const auto ghost = BuildThetaPoleGhostState(interior);

    DEC3D_CHECK_EQ(ghost.rho, interior.rho);
    DEC3D_CHECK_EQ(ghost.mom_r, interior.mom_r);
    DEC3D_CHECK_EQ(ghost.mom_theta, -interior.mom_theta);
    DEC3D_CHECK_EQ(ghost.mom_phi, -interior.mom_phi);
    DEC3D_CHECK_EQ(ghost.e_fluid_total, interior.e_fluid_total);
    DEC3D_CHECK_EQ(ghost.chi_e, interior.chi_e);

    const HydroConservativeState axisymmetric_interior{
        1.0,
        2.0,
        3.0,
        0.0,
        5.0,
        0.6};
    const auto axisymmetric_ghost = BuildThetaPoleGhostState(axisymmetric_interior);
    DEC3D_CHECK_EQ(axisymmetric_ghost.rho, axisymmetric_interior.rho);
    DEC3D_CHECK_EQ(axisymmetric_ghost.mom_r, axisymmetric_interior.mom_r);
    DEC3D_CHECK_EQ(axisymmetric_ghost.mom_theta, -axisymmetric_interior.mom_theta);
    DEC3D_CHECK_EQ(axisymmetric_ghost.mom_phi, 0.0);
    DEC3D_CHECK_EQ(axisymmetric_ghost.e_fluid_total,
                   axisymmetric_interior.e_fluid_total);
    DEC3D_CHECK_EQ(axisymmetric_ghost.chi_e, axisymmetric_interior.chi_e);

    const auto theta_lower =
        PrepareThetaPoleBoundaryStates(BoundaryFace::lower, interior, interior);
    DEC3D_CHECK(theta_lower.is_complete());
    DEC3D_CHECK_EQ(theta_lower.left_state.mom_theta, ghost.mom_theta);
    DEC3D_CHECK_EQ(theta_lower.right_state.mom_theta, interior.mom_theta);

    const auto theta_upper =
        PrepareThetaPoleBoundaryStates(BoundaryFace::upper, interior, interior);
    DEC3D_CHECK(theta_upper.is_complete());
    DEC3D_CHECK_EQ(theta_upper.left_state.mom_phi, interior.mom_phi);
    DEC3D_CHECK_EQ(theta_upper.right_state.mom_phi, ghost.mom_phi);

    const auto phi_periodic = PreparePhiPeriodicBoundaryStates(interior, ghost);
    DEC3D_CHECK(phi_periodic.is_complete());
    DEC3D_CHECK_EQ(phi_periodic.left_state.mom_r, ghost.mom_r);
    DEC3D_CHECK_EQ(phi_periodic.right_state.mom_r, interior.mom_r);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
