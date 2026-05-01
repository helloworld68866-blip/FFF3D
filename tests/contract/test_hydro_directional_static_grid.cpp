#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
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
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::HydroOperator;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto state = CanonicalState::Create(CanonicalStateLayout{1, 3, 4, 0});
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(0, theta, phi) = 1.0;
        state.mom_r(0, theta, phi) = 0.0;
        state.mom_theta(0, theta, phi) = 0.05 * static_cast<double>(theta);
        state.mom_phi(0, theta, phi) = -0.04 * static_cast<double>(phi);
        state.e_fluid_total(0, theta, phi) = 2.5 + 0.15 * static_cast<double>(theta + phi);
        state.e_electron(0, theta, phi) = ElectronEnergyDensityFromPressure(0.4 + 0.05 * static_cast<double>(theta));
      }
    }

    const auto initial_rho = state.rho(0, 1, 1);
    const auto initial_mom_theta = state.mom_theta(0, 1, 1);
    const auto initial_mom_phi = state.mom_phi(0, 1, 1);
    const auto initial_e_fluid_total = state.e_fluid_total(0, 1, 1);

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{1, 3, 4, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    const StageContext context{
        0.0,
        1.0e-3,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.directional",
        "ownership:rank0",
        "diagnostics:p1.hydro.directional"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, state));

    const auto result = hydro.advance();
    DEC3D_CHECK(result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.phi"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.ghost.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.ghost.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.ghost.direction.phi"));
    DEC3D_CHECK(std::abs(state.rho(0, 1, 1) - initial_rho) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.mom_theta(0, 1, 1) - initial_mom_theta) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.mom_phi(0, 1, 1) - initial_mom_phi) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.e_fluid_total(0, 1, 1) - initial_e_fluid_total) > 1.0e-8);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
