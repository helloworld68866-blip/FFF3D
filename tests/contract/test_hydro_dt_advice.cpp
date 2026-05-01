#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct BoundHydroFixture {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  dec3d::state::CanonicalState state;
};

void FillUniformHydroState(
    dec3d::state::CanonicalState& state,
    const dec3d::hydro::HydroPrimitiveState& primitive,
    double electron_pressure) {
  const auto conservative = dec3d::hydro::MakeConservativeState(primitive);
  const double electron_energy =
      dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(radial, theta, phi) = conservative.rho;
        state.mom_r(radial, theta, phi) = conservative.mom_r;
        state.mom_theta(radial, theta, phi) = conservative.mom_theta;
        state.mom_phi(radial, theta, phi) = conservative.mom_phi;
        state.e_fluid_total(radial, theta, phi) = conservative.e_fluid_total;
        state.e_electron(radial, theta, phi) = electron_energy;
      }
    }
  }
}

[[nodiscard]] BoundHydroFixture BuildBoundHydroFixture(
    const dec3d::mesh::SphericalMeshDescriptor& descriptor,
    const dec3d::hydro::HydroPrimitiveState& primitive,
    double electron_pressure) {
  BoundHydroFixture fixture{
      dec3d::mesh::BuildSphericalGeometry(descriptor),
      dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{
              descriptor.radial_cells,
              descriptor.theta_cells,
              descriptor.phi_cells,
              0u})};

  DEC3D_CHECK(fixture.geometry.is_valid());
  FillUniformHydroState(fixture.state, primitive, electron_pressure);
  return fixture;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::mesh::SphericalMeshDescriptor;

    auto baseline = BuildBoundHydroFixture(
        SphericalMeshDescriptor{4u, 4u, 4u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.10, 0.05, 0.02, 1.0, std::pow(0.4, 3.0 / 5.0)},
        0.4);
    const dec3d::core::StageContext context{
        0.0,
        1.0e-3,
        1u,
        dec3d::core::PhaseId::p1,
        "p1-v0.1-dt",
        "mesh:p1.dt",
        "ownership:rank0",
        "diagnostics:p1.hydro.dt"};
    DEC3D_CHECK(context.is_complete());

    dec3d::hydro::HydroOperator baseline_hydro;
    DEC3D_CHECK(baseline_hydro.bind(context, baseline.geometry, baseline.state));
    const auto baseline_dt = baseline_hydro.estimate_dt();
    DEC3D_CHECK(baseline_dt.is_complete());
    if (!std::isfinite(baseline_dt.hard_cap_dt)) {
      throw std::runtime_error(
          "baseline_dt.hard_cap_dt is non-finite; reason=" + baseline_dt.reason +
          "; evidence=" + baseline_dt.evidence);
    }
    DEC3D_CHECK(std::isfinite(baseline_dt.soft_advice_dt));
    DEC3D_CHECK(baseline_dt.hard_cap_dt < std::numeric_limits<double>::infinity());
    DEC3D_CHECK(baseline_dt.soft_advice_dt < std::numeric_limits<double>::infinity());
    DEC3D_CHECK(baseline_dt.evidence != "p1.hydro.dt.scaffold.placeholder");
    DEC3D_CHECK(baseline_dt.evidence.find("p1.hydro.dt.explicit_cfl") != std::string::npos);

    auto faster_flow = BuildBoundHydroFixture(
        SphericalMeshDescriptor{4u, 4u, 4u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.80, 0.40, 0.20, 1.0, std::pow(0.4, 3.0 / 5.0)},
        0.4);
    dec3d::hydro::HydroOperator faster_flow_hydro;
    DEC3D_CHECK(faster_flow_hydro.bind(context, faster_flow.geometry, faster_flow.state));
    const auto faster_flow_dt = faster_flow_hydro.estimate_dt();
    DEC3D_CHECK(faster_flow_dt.hard_cap_dt < baseline_dt.hard_cap_dt);

    auto higher_pressure = BuildBoundHydroFixture(
        SphericalMeshDescriptor{4u, 4u, 4u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.10, 0.05, 0.02, 4.0, std::pow(1.6, 3.0 / 5.0)},
        1.6);
    dec3d::hydro::HydroOperator higher_pressure_hydro;
    DEC3D_CHECK(higher_pressure_hydro.bind(context, higher_pressure.geometry, higher_pressure.state));
    const auto higher_pressure_dt = higher_pressure_hydro.estimate_dt();
    DEC3D_CHECK(higher_pressure_dt.hard_cap_dt < baseline_dt.hard_cap_dt);

    auto finer_radial = BuildBoundHydroFixture(
        SphericalMeshDescriptor{8u, 4u, 4u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.10, 0.05, 0.02, 1.0, std::pow(0.4, 3.0 / 5.0)},
        0.4);
    auto finer_theta = BuildBoundHydroFixture(
        SphericalMeshDescriptor{4u, 8u, 4u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.10, 0.05, 0.02, 1.0, std::pow(0.4, 3.0 / 5.0)},
        0.4);
    auto finer_phi = BuildBoundHydroFixture(
        SphericalMeshDescriptor{4u, 4u, 8u, 0.5, 1.5},
        HydroPrimitiveState{1.0, 0.10, 0.05, 0.02, 1.0, std::pow(0.4, 3.0 / 5.0)},
        0.4);
    dec3d::hydro::HydroOperator finer_radial_hydro;
    dec3d::hydro::HydroOperator finer_theta_hydro;
    dec3d::hydro::HydroOperator finer_phi_hydro;
    DEC3D_CHECK(finer_radial_hydro.bind(context, finer_radial.geometry, finer_radial.state));
    DEC3D_CHECK(finer_theta_hydro.bind(context, finer_theta.geometry, finer_theta.state));
    DEC3D_CHECK(finer_phi_hydro.bind(context, finer_phi.geometry, finer_phi.state));

    DEC3D_CHECK(finer_radial_hydro.estimate_dt().hard_cap_dt < baseline_dt.hard_cap_dt);
    DEC3D_CHECK(finer_theta_hydro.estimate_dt().hard_cap_dt < baseline_dt.hard_cap_dt);
    DEC3D_CHECK(finer_phi_hydro.estimate_dt().hard_cap_dt < baseline_dt.hard_cap_dt);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
