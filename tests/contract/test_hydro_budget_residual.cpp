#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <string_view>

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    constexpr std::size_t kRadialCells = 8u;
    constexpr std::size_t kThetaCells = 3u;
    constexpr std::size_t kPhiCells = 4u;

    auto state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});

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
        }
      }
    }

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto hydro_view = BuildHydroWorkView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const auto result = AdvanceStaticGridHydro(hydro_view, geometry, 1.0e-3);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.budget.is_complete());
    DEC3D_CHECK(std::abs(result.budget.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(result.budget.mom_r_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(result.budget.e_fluid_total_residual) < 1.0e-10);

    dec3d::core::DiagnosticsPayload diagnostics;
    diagnostics.entries.push_back({"p1.hydro.budget.mass", result.budget.report_line});

    dec3d::hydro::HydroBudgetResidualSummary parsed;
    DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(diagnostics, &parsed));
    DEC3D_CHECK(parsed.is_complete());
    DEC3D_CHECK(std::abs(parsed.mass_residual - result.budget.mass_residual) < 1.0e-14);
    DEC3D_CHECK(std::abs(parsed.old_mom_r - result.budget.old_mom_r) < 1.0e-14);
    DEC3D_CHECK(std::abs(parsed.flux_mom_r_delta - result.budget.flux_mom_r_delta) < 1.0e-14);
    DEC3D_CHECK(std::abs(parsed.source_mom_r_delta - result.budget.source_mom_r_delta) < 1.0e-14);
    DEC3D_CHECK(std::abs(parsed.new_mom_r - result.budget.new_mom_r) < 1.0e-14);
    DEC3D_CHECK(std::abs(parsed.mom_r_residual - result.budget.mom_r_residual) < 1.0e-14);
    DEC3D_CHECK(
        std::abs(parsed.e_fluid_total_residual - result.budget.e_fluid_total_residual) < 1.0e-14);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
