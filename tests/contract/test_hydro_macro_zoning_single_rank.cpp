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

void SeedAngularStructure(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.08 * static_cast<double>(radial) +
            0.03 * std::cos(0.7 * static_cast<double>(theta)) +
            0.02 * std::sin(0.5 * static_cast<double>(phi));
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = 0.02 * density;
        state.mom_theta(radial, theta, phi) =
            0.08 * std::sin(0.6 * static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -0.06 * std::cos(0.4 * static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.5 + 0.12 * density;
        state.e_electron(radial, theta, phi) = 0.4;
      }
    }
  }
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto macro_state = CanonicalState::Create(CanonicalStateLayout{2, 8, 8, 0});
    auto fine_state = CanonicalState::Create(CanonicalStateLayout{2, 8, 8, 0});
    SeedAngularStructure(macro_state);
    SeedAngularStructure(fine_state);

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{2, 8, 8, 0.05, 0.25});
    DEC3D_CHECK(geometry.is_valid());

    auto macro_view = BuildHydroWorkView(macro_state);
    auto fine_view = BuildHydroWorkView(fine_state);
    DEC3D_CHECK(macro_view.is_complete());
    DEC3D_CHECK(fine_view.is_complete());

    StaticGridHydroOptions fine_options{};
    fine_options.apply_geometric_source = false;
    fine_options.apply_radial_sweep = false;
    fine_options.apply_theta_sweep = true;
    fine_options.apply_phi_sweep = true;

    auto macro_options = fine_options;
    macro_options.use_macro_zoning = true;
    macro_options.macro_zoning_coarse_factor = 0.5;

    const double initial_macro_rho = (*macro_view.rho)(0, 0, 0);
    const double initial_fine_rho = (*fine_view.rho)(0, 0, 0);

    const auto fine_result = AdvanceStaticGridHydro(fine_view, geometry, 1.0e-3, fine_options);
    const auto macro_result = AdvanceStaticGridHydro(macro_view, geometry, 1.0e-3, macro_options);

    DEC3D_CHECK(fine_result.is_complete());
    DEC3D_CHECK(macro_result.is_complete());
    DEC3D_CHECK(std::abs(macro_result.budget.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(macro_result.budget.mom_r_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(macro_result.budget.e_fluid_total_residual) < 1.0e-10);
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.macro_zoning.detected"));
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.macro_zoning.restrict.executed"));
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.macro_zoning.coarse_update.executed"));
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.macro_zoning.prolong.executed"));
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(macro_result.diagnostics, "p1.hydro.direction.phi"));
    DEC3D_CHECK(std::abs((*macro_view.rho)(0, 0, 0) - initial_macro_rho) > 1.0e-8);
    DEC3D_CHECK(std::abs((*fine_view.rho)(0, 0, 0) - initial_fine_rho) > 1.0e-8);
    DEC3D_CHECK(std::abs((*macro_view.rho)(0, 0, 0) - (*fine_view.rho)(0, 0, 0)) > 1.0e-8);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
