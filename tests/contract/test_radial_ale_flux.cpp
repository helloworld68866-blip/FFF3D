#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ale/radial_ale.hpp"
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

[[nodiscard]] dec3d::hydro::HydroConservativeState Add(
    const dec3d::hydro::HydroConservativeState& lhs,
    const dec3d::hydro::HydroConservativeState& rhs) noexcept {
  return {
      lhs.rho + rhs.rho,
      lhs.mom_r + rhs.mom_r,
      lhs.mom_theta + rhs.mom_theta,
      lhs.mom_phi + rhs.mom_phi,
      lhs.e_fluid_total + rhs.e_fluid_total,
      lhs.chi_e + rhs.chi_e};
}

[[nodiscard]] dec3d::hydro::HydroConservativeState Subtract(
    const dec3d::hydro::HydroConservativeState& lhs,
    const dec3d::hydro::HydroConservativeState& rhs) noexcept {
  return {
      lhs.rho - rhs.rho,
      lhs.mom_r - rhs.mom_r,
      lhs.mom_theta - rhs.mom_theta,
      lhs.mom_phi - rhs.mom_phi,
      lhs.e_fluid_total - rhs.e_fluid_total,
      lhs.chi_e - rhs.chi_e};
}

[[nodiscard]] dec3d::hydro::HydroConservativeState Scale(
    const dec3d::hydro::HydroConservativeState& state,
    double factor) noexcept {
  return {
      state.rho * factor,
      state.mom_r * factor,
      state.mom_theta * factor,
      state.mom_phi * factor,
      state.e_fluid_total * factor,
      state.chi_e * factor};
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::ComputeAleCorrectedFlux;
    using dec3d::hydro::ComputePhysicalFlux;
    using dec3d::hydro::HllcActiveRegion;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::hydro::SelectMovingInterfaceBranch;
    using dec3d::hydro::SolveHllcRiemann;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    const HydroPrimitiveState left{
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const HydroPrimitiveState right{
        0.125,
        0.0,
        0.0,
        0.0,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};

    const auto left_state = MakeConservativeState(left);
    const auto right_state = MakeConservativeState(right);
    const auto hllc = SolveHllcRiemann(left_state, right_state);
    DEC3D_CHECK(hllc.is_complete());
    DEC3D_CHECK(hllc.active_region == HllcActiveRegion::left_star);
    DEC3D_CHECK(hllc.waves.s_left < 0.0);
    DEC3D_CHECK(hllc.waves.s_star > 0.0);
    DEC3D_CHECK(hllc.waves.s_right > hllc.waves.s_star);

    const double face_speed = 0.5 * (hllc.waves.s_star + hllc.waves.s_right);
    const auto moving_branch =
        SelectMovingInterfaceBranch(hllc, left_state, right_state, face_speed);
    DEC3D_CHECK(moving_branch.is_complete());
    DEC3D_CHECK(moving_branch.active_region == HllcActiveRegion::right_star);

    const auto right_flux = ComputePhysicalFlux(right_state);
    const auto right_star_flux = Add(
        right_flux,
        Scale(Subtract(hllc.right_star_state, right_state), hllc.waves.s_right));
    const auto expected_ale_flux =
        Subtract(right_star_flux, Scale(hllc.right_star_state, face_speed));
    const auto ale_flux =
        ComputeAleCorrectedFlux(hllc, left_state, right_state, face_speed);
    DEC3D_CHECK(std::abs(ale_flux.rho - expected_ale_flux.rho) < 1.0e-12);
    DEC3D_CHECK(std::abs(ale_flux.mom_r - expected_ale_flux.mom_r) < 1.0e-12);
    DEC3D_CHECK(std::abs(ale_flux.e_fluid_total - expected_ale_flux.e_fluid_total) < 1.0e-12);
    DEC3D_CHECK(std::abs(ale_flux.rho - hllc.interface_flux.rho) > 1.0e-12);

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{4u, 2u, 2u, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto proposal_state = CanonicalState::Create(CanonicalStateLayout{4u, 2u, 2u, 0u});
    for (std::size_t radial = 0; radial < proposal_state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < proposal_state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < proposal_state.rho.extent_phi(); ++phi) {
          proposal_state.rho(radial, theta, phi) = 1.0;
          proposal_state.mom_r(radial, theta, phi) = -0.2 - (0.05 * static_cast<double>(radial));
          proposal_state.mom_theta(radial, theta, phi) = 0.0;
          proposal_state.mom_phi(radial, theta, phi) = 0.0;
          proposal_state.e_fluid_total(radial, theta, phi) = 2.5;
          proposal_state.e_electron(radial, theta, phi) = ElectronEnergyDensityFromPressure(0.4);
        }
      }
    }

    const double dt_s = 1.0e-2;
    const auto proposal = BuildRadialAleMeshUpdateProposal(
        proposal_state.rho,
        proposal_state.mom_r,
        geometry,
        dt_s);
    DEC3D_CHECK(proposal.is_complete(geometry.radial_faces.size()));
    DEC3D_CHECK(proposal.radial_face_velocities.back() < 0.0);
    DEC3D_CHECK(proposal.proposed_radial_faces.back() < geometry.radial_faces.back());
    const double scale_factor =
        proposal.proposed_radial_faces.back() / geometry.radial_faces.back();
    for (std::size_t radial_face = 0; radial_face < geometry.radial_faces.size(); ++radial_face) {
      DEC3D_CHECK(std::abs(
          proposal.proposed_radial_faces[radial_face] -
          (geometry.radial_faces[radial_face] * scale_factor)) < 1.0e-12);
    }

    auto static_state = CanonicalState::Create(CanonicalStateLayout{4u, 2u, 2u, 0u});
    auto ale_state = CanonicalState::Create(CanonicalStateLayout{4u, 2u, 2u, 0u});

    for (std::size_t radial = 0; radial < static_state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < static_state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < static_state.rho.extent_phi(); ++phi) {
          const double rho = radial < 2u ? 1.0 : 0.6;
          const double v_r = radial < 2u ? 0.25 : -0.10;
          const double pressure = radial < 2u ? 1.0 : 0.5;
          const double energy =
              pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) +
              0.5 * rho * v_r * v_r;
          const double e_electron = ElectronEnergyDensityFromPressure(0.4 * pressure);

          static_state.rho(radial, theta, phi) = rho;
          static_state.mom_r(radial, theta, phi) = rho * v_r;
          static_state.mom_theta(radial, theta, phi) = 0.0;
          static_state.mom_phi(radial, theta, phi) = 0.0;
          static_state.e_fluid_total(radial, theta, phi) = energy;
          static_state.e_electron(radial, theta, phi) = e_electron;

          ale_state.rho(radial, theta, phi) = rho;
          ale_state.mom_r(radial, theta, phi) = rho * v_r;
          ale_state.mom_theta(radial, theta, phi) = 0.0;
          ale_state.mom_phi(radial, theta, phi) = 0.0;
          ale_state.e_fluid_total(radial, theta, phi) = energy;
          ale_state.e_electron(radial, theta, phi) = e_electron;
        }
      }
    }

    auto static_view = BuildHydroWorkView(static_state);
    auto ale_view = BuildHydroWorkView(ale_state);
    DEC3D_CHECK(static_view.is_complete());
    DEC3D_CHECK(ale_view.is_complete());

    StaticGridHydroOptions static_options;
    static_options.apply_geometric_source = false;
    static_options.apply_theta_sweep = false;
    static_options.apply_phi_sweep = false;

    StaticGridHydroOptions ale_options = static_options;
    ale_options.apply_radial_ale_flux_correction = true;

    const auto static_result = AdvanceStaticGridHydro(
        static_view,
        geometry,
        dt_s,
        static_options);
    DEC3D_CHECK(static_result.success);

    const auto ale_result = AdvanceStaticGridHydro(
        ale_view,
        geometry,
        dt_s,
        dec3d::hydro::RadialGhostOverride{},
        ale_options,
        &proposal);
    DEC3D_CHECK(ale_result.success);
    DEC3D_CHECK(HasDiagnosticCode(
        ale_result.diagnostics,
        "p1.hydro.ale.radial_flux_correction.executed"));
    DEC3D_CHECK(!HasDiagnosticCode(
        static_result.diagnostics,
        "p1.hydro.ale.radial_flux_correction.executed"));
    DEC3D_CHECK(std::abs(
        (*ale_view.rho)(1u, 0u, 0u) - (*static_view.rho)(1u, 0u, 0u)) > 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
