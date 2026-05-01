#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] double ComputeTotalMass(
    const dec3d::state::HydroStateView& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double total_mass = 0.0;
  std::size_t linear_index = 0u;
  for (std::size_t radial = 0; radial < state.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho->extent_phi(); ++phi) {
        total_mass += (*state.rho)(radial, theta, phi) * geometry.cell_volumes[linear_index++];
      }
    }
  }
  return total_mass;
}

[[nodiscard]] double ComputeTotalFluidEnergy(
    const dec3d::state::HydroStateView& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double total_energy = 0.0;
  std::size_t linear_index = 0u;
  for (std::size_t radial = 0; radial < state.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho->extent_phi(); ++phi) {
        total_energy +=
            (*state.e_fluid_total)(radial, theta, phi) * geometry.cell_volumes[linear_index++];
      }
    }
  }
  return total_energy;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::ApplyRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto state = CanonicalState::Create(CanonicalStateLayout{4u, 2u, 2u, 0u});
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          const double rho = radial < 2u ? 1.0 : 0.8;
          const double v_r = radial < 2u ? -0.25 : -0.15;
          const double pressure = radial < 2u ? 1.2 : 0.9;
          state.rho(radial, theta, phi) = rho;
          state.mom_r(radial, theta, phi) = rho * v_r;
          state.mom_theta(radial, theta, phi) = 0.0;
          state.mom_phi(radial, theta, phi) = 0.0;
          state.e_fluid_total(radial, theta, phi) =
              pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) +
              0.5 * rho * v_r * v_r;
          state.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(0.4 * pressure);
        }
      }
    }

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{4u, 2u, 2u, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());
    auto updated_geometry = geometry;

    auto hydro_view = BuildHydroWorkView(state);
    DEC3D_CHECK(hydro_view.is_complete());

    const double dt_s = 1.0e-2;
    const auto proposal =
        BuildRadialAleMeshUpdateProposal(state.rho, state.mom_r, geometry, dt_s);
    DEC3D_CHECK(proposal.is_complete(geometry.radial_faces.size()));
    DEC3D_CHECK(ApplyRadialAleMeshUpdateProposal(proposal, updated_geometry).success);

    StaticGridHydroOptions options;
    options.apply_radial_ale_flux_correction = true;
    options.apply_theta_sweep = false;
    options.apply_phi_sweep = false;
    options.apply_geometric_source = true;

    const auto initial_view = dec3d::state::BuildHydroStateView(state);
    DEC3D_CHECK(initial_view.is_complete());
    const double old_mass = ComputeTotalMass(initial_view, geometry);
    const double old_energy = ComputeTotalFluidEnergy(initial_view, geometry);

    const auto result = AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        dt_s,
        options,
        &proposal);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.budget.is_complete());

    const double new_mass = ComputeTotalMass(hydro_view, updated_geometry);
    const double new_energy = ComputeTotalFluidEnergy(hydro_view, updated_geometry);

    DEC3D_CHECK(std::abs(result.budget.old_mass - old_mass) < 1.0e-12);
    DEC3D_CHECK(std::abs(result.budget.old_e_fluid_total - old_energy) < 1.0e-12);
    DEC3D_CHECK(std::abs(result.budget.new_mass - new_mass) < 1.0e-12);
    DEC3D_CHECK(std::abs(result.budget.new_e_fluid_total - new_energy) < 1.0e-12);
    DEC3D_CHECK(std::abs(result.budget.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(result.budget.e_fluid_total_residual) < 1.0e-10);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
