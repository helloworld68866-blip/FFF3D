#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <string_view>

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

[[nodiscard]] bool Contains(std::string_view text, std::string_view token) noexcept {
  return text.find(token) != std::string_view::npos;
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
    const auto initial_mom_r = hydro_view.mom_r->storage();
    const auto initial_mom_theta = hydro_view.mom_theta->storage();
    const auto initial_mom_phi = hydro_view.mom_phi->storage();
    const auto initial_e_total = hydro_view.e_fluid_total->storage();
    const auto initial_chi_e = hydro_view.chi_e.storage();

    const auto result = AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        1.0e-3,
        StaticGridHydroOptions{});
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.direction.phi"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.source.geometric_step"));

    for (std::size_t index = 0; index < initial_rho.size(); ++index) {
      DEC3D_CHECK(std::abs(hydro_view.rho->storage()[index] - initial_rho[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(hydro_view.mom_r->storage()[index] - initial_mom_r[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(hydro_view.mom_theta->storage()[index] - initial_mom_theta[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(hydro_view.mom_phi->storage()[index] - initial_mom_phi[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(hydro_view.e_fluid_total->storage()[index] - initial_e_total[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(hydro_view.chi_e.storage()[index] - initial_chi_e[index]) < 1.0e-10);
    }

    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      const double reference_rho = (*hydro_view.rho)(radial, 0, 0);
      const double reference_mom_r = (*hydro_view.mom_r)(radial, 0, 0);
      const double reference_e_total = (*hydro_view.e_fluid_total)(radial, 0, 0);
      const double reference_chi_e = hydro_view.chi_e(radial, 0, 0);
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        const double reference_mom_theta = (*hydro_view.mom_theta)(radial, theta, 0);
        const double reference_mom_phi = (*hydro_view.mom_phi)(radial, theta, 0);
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          DEC3D_CHECK(std::abs((*hydro_view.rho)(radial, theta, phi) - reference_rho) < 1.0e-10);
          DEC3D_CHECK(std::abs((*hydro_view.mom_r)(radial, theta, phi) - reference_mom_r) < 1.0e-10);
          DEC3D_CHECK(std::abs((*hydro_view.e_fluid_total)(radial, theta, phi) - reference_e_total) < 1.0e-10);
          DEC3D_CHECK(std::abs(hydro_view.chi_e(radial, theta, phi) - reference_chi_e) < 1.0e-10);
          DEC3D_CHECK(std::abs((*hydro_view.mom_theta)(radial, theta, phi) - reference_mom_theta) < 1.0e-10);
          DEC3D_CHECK(std::abs((*hydro_view.mom_phi)(radial, theta, phi) - reference_mom_phi) < 1.0e-10);
        }
      }

      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        const auto mirrored_theta = state.rho.extent_theta() - 1u - theta;
        DEC3D_CHECK(std::abs((*hydro_view.rho)(radial, theta, 0) - (*hydro_view.rho)(radial, mirrored_theta, 0)) < 1.0e-10);
        DEC3D_CHECK(std::abs((*hydro_view.mom_r)(radial, theta, 0) - (*hydro_view.mom_r)(radial, mirrored_theta, 0)) < 1.0e-10);
        DEC3D_CHECK(std::abs((*hydro_view.e_fluid_total)(radial, theta, 0) - (*hydro_view.e_fluid_total)(radial, mirrored_theta, 0)) < 1.0e-10);
        DEC3D_CHECK(std::abs(hydro_view.chi_e(radial, theta, 0) - hydro_view.chi_e(radial, mirrored_theta, 0)) < 1.0e-10);
        DEC3D_CHECK(std::abs((*hydro_view.mom_theta)(radial, theta, 0) + (*hydro_view.mom_theta)(radial, mirrored_theta, 0)) < 1.0e-10);
        DEC3D_CHECK(std::abs((*hydro_view.mom_phi)(radial, theta, 0) - (*hydro_view.mom_phi)(radial, mirrored_theta, 0)) < 1.0e-10);
      }
    }

    auto radiation_state = CanonicalState::Create(CanonicalStateLayout{2, 4, 4, 2});
    for (std::size_t radial = 0; radial < radiation_state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < radiation_state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < radiation_state.rho.extent_phi(); ++phi) {
          radiation_state.rho(radial, theta, phi) = conservative.rho;
          radiation_state.mom_r(radial, theta, phi) = conservative.mom_r;
          radiation_state.mom_theta(radial, theta, phi) = conservative.mom_theta;
          radiation_state.mom_phi(radial, theta, phi) = conservative.mom_phi;
          radiation_state.e_fluid_total(radial, theta, phi) = conservative.e_fluid_total;
          radiation_state.e_electron(radial, theta, phi) = ElectronEnergyDensityFromPressure(0.4);
          radiation_state.radiation_groups[0](radial, theta, phi) = 9.0;
          radiation_state.radiation_groups[1](radial, theta, phi) = 27.0;
        }
      }
    }

    auto radiation_view = BuildHydroStateView(radiation_state);
    DEC3D_CHECK(radiation_view.is_complete());
    DEC3D_CHECK_EQ(radiation_view.radiation_chi.size(), 2u);
    const auto initial_chi_rad0 = radiation_view.radiation_chi[0].storage();
    const auto initial_chi_rad1 = radiation_view.radiation_chi[1].storage();
    StaticGridHydroOptions radiation_options;
    radiation_options.enable_radiation_hydro_terms = true;
    const auto radiation_result = AdvanceStaticGridHydro(
        radiation_view,
        geometry,
        1.0e-3,
        radiation_options);
    DEC3D_CHECK(radiation_result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(
        radiation_result.diagnostics,
        "p3.radiation.hydro_terms.stage_H"));
    DEC3D_CHECK(Contains(radiation_result.report_line, "radiation_hydro_terms_enabled=true"));
    DEC3D_CHECK(Contains(radiation_result.report_line, "advected_radiation_scalar=P_g_power_3_over_4"));
    DEC3D_CHECK(Contains(radiation_result.report_line, "passive_Ug_advection=false"));
    for (std::size_t index = 0; index < initial_chi_rad0.size(); ++index) {
      DEC3D_CHECK(std::abs(radiation_view.radiation_chi[0].storage()[index] - initial_chi_rad0[index]) < 1.0e-10);
      DEC3D_CHECK(std::abs(radiation_view.radiation_chi[1].storage()[index] - initial_chi_rad1[index]) < 1.0e-10);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
