#include "hydro/driver/static_grid_hydro.hpp"
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

void SeedMacroZoningRadialWave(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before) {
  const dec3d::hydro::HydroPrimitiveState inner_state{
      1.0,
      0.0,
      0.0,
      0.0,
      1.0,
      std::pow(0.4, 3.0 / 5.0)};
  const dec3d::hydro::HydroPrimitiveState outer_state{
      0.125,
      0.0,
      0.0,
      0.0,
      0.1,
      std::pow(0.04, 3.0 / 5.0)};

  const auto inner_conservative = dec3d::hydro::MakeConservativeState(inner_state);
  const auto outer_conservative = dec3d::hydro::MakeConservativeState(outer_state);
  const std::size_t radial_midpoint = state.layout.radial_cells / 2u;

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const bool inside = radial < radial_midpoint;
    const auto& base = inside ? inner_conservative : outer_conservative;
    const double electron_pressure = inside ? 0.4 : 0.04;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      const double theta_mode = 1.0 + 0.01 * std::cos(0.7 * static_cast<double>(theta));
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const double phi_mode = 1.0 + 0.01 * std::sin(0.5 * static_cast<double>(phi));
        const double modifier = theta_mode * phi_mode;

        state.rho(radial, theta, phi) = base.rho * modifier;
        state.mom_r(radial, theta, phi) = base.mom_r;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = base.e_fluid_total * modifier;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure * modifier);

        before.rho(radial, theta, phi) = state.rho(radial, theta, phi);
        before.mom_r(radial, theta, phi) = state.mom_r(radial, theta, phi);
        before.mom_theta(radial, theta, phi) = state.mom_theta(radial, theta, phi);
        before.mom_phi(radial, theta, phi) = state.mom_phi(radial, theta, phi);
        before.e_fluid_total(radial, theta, phi) = state.e_fluid_total(radial, theta, phi);
        before.e_electron(radial, theta, phi) = state.e_electron(radial, theta, phi);
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MacroZoningActiveOptions() noexcept {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = false;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_macro_zoning = true;
  options.macro_zoning_coarse_factor = 0.5;
  return options;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    constexpr std::size_t kRadialCells = 16u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;

    auto before = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto state = CanonicalState::Create(before.layout);
    SeedMacroZoningRadialWave(state, before);

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
    DEC3D_CHECK(geometry.is_valid());

    auto hydro_view = BuildHydroWorkView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const auto result = AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        1.0e-3,
        MacroZoningActiveOptions());
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(std::abs(result.budget.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(result.budget.mom_r_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(result.budget.e_fluid_total_residual) < 1.0e-10);
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.detected"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.restrict.executed"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.coarse_update.executed"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.prolong.executed"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.budget.mass"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.budget.mom_r"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.budget.e_fluid_total"));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
