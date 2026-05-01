#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
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
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::EvaluateSphericalRadialWaveSanity;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    constexpr std::size_t kRadialCells = 16u;
    constexpr std::size_t kThetaCells = 3u;
    constexpr std::size_t kPhiCells = 4u;

    auto state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto before = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});

    const HydroPrimitiveState inner_state{
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const HydroPrimitiveState outer_state{
        0.125,
        0.0,
        0.0,
        0.0,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};

    const auto inner_conservative = MakeConservativeState(inner_state);
    const auto outer_conservative = MakeConservativeState(outer_state);

    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      const bool inside = radial < (kRadialCells / 2u);
      const auto& cell = inside ? inner_conservative : outer_conservative;
      const double electron_pressure = inside ? 0.4 : 0.04;
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          state.rho(radial, theta, phi) = cell.rho;
          state.mom_r(radial, theta, phi) = cell.mom_r;
          state.mom_theta(radial, theta, phi) = cell.mom_theta;
          state.mom_phi(radial, theta, phi) = cell.mom_phi;
          state.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
          state.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(electron_pressure);

          before.rho(radial, theta, phi) = cell.rho;
          before.mom_r(radial, theta, phi) = cell.mom_r;
          before.mom_theta(radial, theta, phi) = cell.mom_theta;
          before.mom_phi(radial, theta, phi) = cell.mom_phi;
          before.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
          before.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(electron_pressure);
        }
      }
    }

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    const StageContext context{
        0.0,
        1.0e-3,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.radial-wave",
        "ownership:rank0",
        "diagnostics:p1.hydro.radial-wave"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, state));
    const auto stage_result = hydro.advance();
    DEC3D_CHECK(stage_result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.ghost.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.budget.mass"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.budget.e_fluid_total"));

    const auto sanity = EvaluateSphericalRadialWaveSanity(
        before,
        state,
        1.0e-10,
        1.0e-10,
        1.0e-8);
    DEC3D_CHECK(sanity.is_complete());
    DEC3D_CHECK(sanity.angular_symmetry_preserved);
    DEC3D_CHECK(sanity.tangential_momentum_quiet);
    DEC3D_CHECK(sanity.radial_profile_changed);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
