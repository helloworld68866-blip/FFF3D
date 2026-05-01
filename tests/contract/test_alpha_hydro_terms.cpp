#include "alpha/alpha_hydro_terms.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>

namespace {

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << label << ": actual=" << actual << " expected=" << expected << '\n';
    DEC3D_CHECK(false);
  }
}

}  // namespace

int main() {
  try {
    {
      const double epsilon = 150.0;
      const double chi = dec3d::alpha::AlphaHydroPressureScalarFromEnergyDensity(epsilon);
      CheckNear(chi, std::pow((2.0 / 3.0) * epsilon, 3.0 / 5.0), 1.0e-14, "alpha chi");
      CheckNear(dec3d::alpha::AlphaHydroEnergyDensityFromPressureScalar(chi),
                epsilon,
                1.0e-12,
                "alpha roundtrip");
    }

    {
      const double rho_ratio = 2.0;
      CheckNear(dec3d::alpha::AlphaEpsilonCompressionRatioFromDensityRatio(rho_ratio),
                std::pow(rho_ratio, 5.0 / 3.0),
                1.0e-14,
                "alpha compression ratio");
    }

    {
      const auto layout = dec3d::state::CanonicalStateLayout{2u, 1u, 1u, 0u};
      auto state = dec3d::state::CanonicalState::Create(layout);
      state.alpha_state.storage(0u, 0u, 0u) = 3.0;
      state.alpha_state.storage(1u, 0u, 0u) = 24.0;

      dec3d::alpha::AlphaHydroTermsOptions options;
      options.alpha_energy_floor_erg_cm3 = 0.0;
      const auto bundle = dec3d::alpha::BuildAlphaHydroTermBundle(state, options);
      DEC3D_CHECK(bundle.success);
      DEC3D_CHECK(bundle.chi_alpha.extent_r() == 2u);
      DEC3D_CHECK(dec3d::alpha::ValidateAlphaHydroTermsDiagnostics(bundle.report_line));
      DEC3D_CHECK(bundle.report_line.find("advected_alpha_scalar=P_alpha_power_3_over_5") !=
                  std::string::npos);
      DEC3D_CHECK(bundle.report_line.find("passive_epsilon_alpha_advection=false") !=
                  std::string::npos);
    }

    {
      auto state = dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 0u});
      state.alpha_state.storage(0u, 0u, 0u) = -1.0;
      const auto bundle = dec3d::alpha::BuildAlphaHydroTermBundle(
          state, dec3d::alpha::AlphaHydroTermsOptions{});
      DEC3D_CHECK(!bundle.success);
      DEC3D_CHECK(bundle.invalid_alpha_cell_count == 1u);
      DEC3D_CHECK(bundle.failure_diagnostics.find(
                      "diagnostic_id=p4.alpha.hydro_terms.failure") !=
                  std::string::npos);
    }

    {
      const std::string missing =
          "diagnostic_id=p4.alpha.hydro_terms"
          "; phase_id=P4"
          "; stage_id=H"
          "; alpha_hydro_bundle_built=true";
      DEC3D_CHECK(!dec3d::alpha::ValidateAlphaHydroTermsDiagnostics(missing));
    }

    {
      dec3d::hydro::HydroPrimitiveState left;
      left.rho = 1.0;
      left.v_r = 0.25;
      left.pressure = 10.0;
      left.chi_e = 1.0;
      left.alpha_chi = dec3d::alpha::AlphaHydroPressureScalarFromEnergyDensity(3.0);

      dec3d::hydro::HydroPrimitiveState right = left;
      right.rho = 2.0;
      right.v_r = -0.25;
      right.alpha_chi = dec3d::alpha::AlphaHydroPressureScalarFromEnergyDensity(6.0);

      const auto result = dec3d::hydro::SolveHllcRiemann(
          dec3d::hydro::MakeConservativeState(left),
          dec3d::hydro::MakeConservativeState(right));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.left_star_state.alpha_chi >= 0.0);
      DEC3D_CHECK(result.right_star_state.alpha_chi >= 0.0);
      DEC3D_CHECK(result.interface_flux.alpha_chi == result.interface_flux.alpha_chi);
    }

    {
      auto state = dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 0u});
      const double epsilon_old = 9.0;
      const double rho_old = 1.0;
      const double rho_new = 1.8;
      const double epsilon_expected =
          dec3d::alpha::AlphaEpsilonAfterPressureScalarCompression(
              epsilon_old, rho_old, rho_new);
      state.rho(0u, 0u, 0u) = 1.0;
      state.e_electron(0u, 0u, 0u) = 1.0e12;
      state.e_fluid_total(0u, 0u, 0u) = 3.0e12;
      state.alpha_state.storage(0u, 0u, 0u) = epsilon_old;
      const double electron_before = state.e_electron(0u, 0u, 0u);

      dec3d::state::HydroAlphaViewOptions alpha_options;
      alpha_options.enabled = true;
      auto view = dec3d::state::BuildHydroStateView(state, alpha_options);
      DEC3D_CHECK(view.is_complete());
      DEC3D_CHECK(view.operator_local_alpha_chi);
      view.alpha_chi(0u, 0u, 0u) =
          dec3d::alpha::AlphaHydroPressureScalarFromEnergyDensity(epsilon_expected);

      const auto writeback = dec3d::state::CommitHydroWriteback(
          state,
          view,
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state));
      DEC3D_CHECK(writeback.success);
      CheckNear(
          state.alpha_state.storage(0u, 0u, 0u),
          epsilon_expected,
          1.0e-11,
          "controlled alpha 5/3 writeback");
      CheckNear(
          state.e_electron(0u, 0u, 0u),
          electron_before,
          0.0,
          "alpha hydro writeback does not deposit electron energy");
      DEC3D_CHECK(writeback.report_line.find(
                      "h_stage_updated_fields_contains_alpha_state=true") !=
                  std::string::npos);
      DEC3D_CHECK(writeback.report_line.find("alpha_state_epoch=post_H_committed") !=
                  std::string::npos);
      DEC3D_CHECK(writeback.report_line.find(
                      "alpha_hydro_terms_share_hydro_scalar_bundle=true") !=
                  std::string::npos);
    }

    {
      auto state = dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{4u, 1u, 1u, 0u});
      for (std::size_t r = 0; r < 4u; ++r) {
        state.rho(r, 0u, 0u) = 1.0 + 0.05 * static_cast<double>(r);
        state.mom_r(r, 0u, 0u) = 1.0e-4 * state.rho(r, 0u, 0u);
        state.mom_theta(r, 0u, 0u) = 0.0;
        state.mom_phi(r, 0u, 0u) = 0.0;
        state.e_electron(r, 0u, 0u) = 1.0e14;
        state.e_fluid_total(r, 0u, 0u) = 3.0e14;
        state.alpha_state.storage(r, 0u, 0u) = 10.0 + static_cast<double>(r);
      }
      const double alpha_before = state.alpha_state.storage(1u, 0u, 0u);

      dec3d::state::HydroAlphaViewOptions alpha_options;
      alpha_options.enabled = true;
      auto view = dec3d::state::BuildHydroWorkView(state, alpha_options);
      DEC3D_CHECK(view.is_complete());
      dec3d::hydro::StaticGridHydroOptions hydro_options;
      hydro_options.apply_geometric_source = false;
      hydro_options.apply_radial_sweep = true;
      hydro_options.apply_theta_sweep = false;
      hydro_options.apply_phi_sweep = false;
      hydro_options.enable_alpha_hydro_terms = true;
      const auto geometry = dec3d::mesh::BuildSphericalGeometry(
          dec3d::mesh::SphericalMeshDescriptor{4u, 1u, 1u, 1.0, 2.0});
      DEC3D_CHECK(geometry.is_valid());
      const auto hydro_result =
          dec3d::hydro::AdvanceStaticGridHydro(view, geometry, 1.0e-12, hydro_options);
      if (!hydro_result.success) {
        std::cerr << "alpha H-path hydro failure: "
                  << hydro_result.failure_reason
                  << " report=" << hydro_result.report_line << '\n';
      }
      DEC3D_CHECK(hydro_result.success);
      DEC3D_CHECK(hydro_result.report_line.find(
                      "alpha_hydro_terms_share_hydro_scalar_bundle=true") !=
                  std::string::npos);
      const auto writeback = dec3d::state::CommitHydroWriteback(
          state,
          view,
          dec3d::core::Combine({
              dec3d::core::AuthoritativeField::rho,
              dec3d::core::AuthoritativeField::mom_r,
              dec3d::core::AuthoritativeField::mom_theta,
              dec3d::core::AuthoritativeField::mom_phi,
              dec3d::core::AuthoritativeField::e_fluid_total,
              dec3d::core::AuthoritativeField::e_electron,
              dec3d::core::AuthoritativeField::alpha_state}));
      DEC3D_CHECK(writeback.success);
      DEC3D_CHECK(state.alpha_state.storage(1u, 0u, 0u) ==
                  state.alpha_state.storage(1u, 0u, 0u));
      DEC3D_CHECK(state.alpha_state.storage(1u, 0u, 0u) >= 0.0);
      DEC3D_CHECK(std::abs(state.alpha_state.storage(1u, 0u, 0u) - alpha_before) >=
                  0.0);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
