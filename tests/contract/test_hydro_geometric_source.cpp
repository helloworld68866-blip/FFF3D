#include "hydro/riemann/hllc_solver.hpp"
#include "hydro/source/geometric_source.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }

  return false;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::ApplyGeometricSourceStep;
    using dec3d::hydro::CaptureGeometricSourceSnapshot;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroStateView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ChiEFromElectronPressure;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto state = CanonicalState::Create(CanonicalStateLayout{1, 2, 1, 0});
    const HydroPrimitiveState primitive{
        2.0,
        0.3,
        0.4,
        -0.2,
        1.5,
        ChiEFromElectronPressure(0.6)};
    const auto conservative = MakeConservativeState(primitive);

    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      state.rho(0, theta, 0) = conservative.rho;
      state.mom_r(0, theta, 0) = conservative.mom_r;
      state.mom_theta(0, theta, 0) = conservative.mom_theta;
      state.mom_phi(0, theta, 0) = conservative.mom_phi;
      state.e_fluid_total(0, theta, 0) = conservative.e_fluid_total;
      state.e_electron(0, theta, 0) = ElectronEnergyDensityFromPressure(0.6);
    }

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{1, 2, 1, 1.0, 2.0});
    DEC3D_CHECK(geometry.is_valid());

    auto hydro_view = BuildHydroStateView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const double initial_rho = state.rho(0, 0, 0);
    const double initial_e_total = state.e_fluid_total(0, 0, 0);
    const double initial_chi_e = hydro_view.chi_e(0, 0, 0);
    const double initial_mom_r = state.mom_r(0, 0, 0);
    const double initial_mom_theta = state.mom_theta(0, 0, 0);
    const double initial_mom_phi = state.mom_phi(0, 0, 0);
    const double dt_s = 1.0e-2;

    const auto snapshot = CaptureGeometricSourceSnapshot(hydro_view);
    DEC3D_CHECK(snapshot.is_complete());

    const auto result = ApplyGeometricSourceStep(hydro_view, snapshot, geometry, dt_s);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.source.geometric_step"));

    const double radius_center = 0.5 * (geometry.radial_faces[0] + geometry.radial_faces[1]);
    const double theta_center = 0.5 * (geometry.theta_faces[0] + geometry.theta_faces[1]);
    const double cot_theta = std::cos(theta_center) / std::sin(theta_center);
    const double volume = geometry.cell_volumes[0];
    const double radial_face_plus =
        geometry.radial_faces[1] * geometry.radial_faces[1] *
        (std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1])) *
        (geometry.phi_faces[1] - geometry.phi_faces[0]);
    const double radial_face_minus =
        geometry.radial_faces[0] * geometry.radial_faces[0] *
        (std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1])) *
        (geometry.phi_faces[1] - geometry.phi_faces[0]);
    const double theta_face_plus =
        0.5 *
        (geometry.radial_faces[1] * geometry.radial_faces[1] -
         geometry.radial_faces[0] * geometry.radial_faces[0]) *
        std::sin(geometry.theta_faces[1]) *
        (geometry.phi_faces[1] - geometry.phi_faces[0]);
    const double theta_face_minus =
        0.5 *
        (geometry.radial_faces[1] * geometry.radial_faces[1] -
         geometry.radial_faces[0] * geometry.radial_faces[0]) *
        std::sin(geometry.theta_faces[0]) *
        (geometry.phi_faces[1] - geometry.phi_faces[0]);
    const double source_mom_r =
        (primitive.pressure * (radial_face_plus - radial_face_minus) / volume) +
        ((initial_mom_theta * initial_mom_theta + initial_mom_phi * initial_mom_phi) /
         (initial_rho * radius_center));
    const double source_mom_theta =
        (-(initial_mom_theta * initial_mom_r) / (initial_rho * radius_center)) +
        ((cot_theta / radius_center) *
         ((initial_mom_phi * initial_mom_phi) / initial_rho)) +
        (primitive.pressure * (theta_face_plus - theta_face_minus) / volume);
    const double source_mom_phi =
        (-(initial_mom_r * initial_mom_phi) / (initial_rho * radius_center)) -
        ((cot_theta / radius_center) *
         ((initial_mom_theta * initial_mom_phi) / initial_rho));

    DEC3D_CHECK(std::abs(state.rho(0, 0, 0) - initial_rho) < 1.0e-12);
    DEC3D_CHECK(std::abs(state.e_fluid_total(0, 0, 0) - initial_e_total) < 1.0e-12);
    DEC3D_CHECK(std::abs(hydro_view.chi_e(0, 0, 0) - initial_chi_e) < 1.0e-12);
    DEC3D_CHECK(std::abs(state.mom_r(0, 0, 0) - (initial_mom_r + dt_s * source_mom_r)) < 1.0e-12);
    DEC3D_CHECK(std::abs(state.mom_theta(0, 0, 0) - (initial_mom_theta + dt_s * source_mom_theta)) < 1.0e-12);
    DEC3D_CHECK(std::abs(state.mom_phi(0, 0, 0) - (initial_mom_phi + dt_s * source_mom_phi)) < 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
