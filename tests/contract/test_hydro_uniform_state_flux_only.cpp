#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
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
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroStateView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 4, 4, 0});
    const HydroPrimitiveState uniform{
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const auto conservative = MakeConservativeState(uniform);
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          state.rho(radial, theta, phi) = conservative.rho;
          state.mom_r(radial, theta, phi) = conservative.mom_r;
          state.mom_theta(radial, theta, phi) = conservative.mom_theta;
          state.mom_phi(radial, theta, phi) = conservative.mom_phi;
          state.e_fluid_total(radial, theta, phi) = conservative.e_fluid_total;
          state.e_electron(radial, theta, phi) = ElectronEnergyDensityFromPressure(0.4);
        }
      }
    }

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{2, 4, 4, 1.0, 2.0});
    DEC3D_CHECK(geometry.is_valid());

    auto hydro_view = BuildHydroStateView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const auto initial_rho = hydro_view.rho->storage();
    const auto initial_e_total = hydro_view.e_fluid_total->storage();
    const auto initial_chi_e = hydro_view.chi_e.storage();

    const auto result = AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        1.0e-3,
        StaticGridHydroOptions{false, true, true, true});
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.phi"));

    for (std::size_t index = 0; index < initial_rho.size(); ++index) {
      DEC3D_CHECK(std::abs(hydro_view.rho->storage()[index] - initial_rho[index]) < 1.0e-12);
      DEC3D_CHECK(std::abs(hydro_view.e_fluid_total->storage()[index] - initial_e_total[index]) < 1.0e-12);
      DEC3D_CHECK(std::abs(hydro_view.chi_e.storage()[index] - initial_chi_e[index]) < 1.0e-12);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
