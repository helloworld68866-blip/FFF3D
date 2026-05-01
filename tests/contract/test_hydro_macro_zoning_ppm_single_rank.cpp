#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
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
            1.0 + 0.04 * static_cast<double>(radial) +
            0.06 * std::cos(0.7 * static_cast<double>(theta)) +
            0.05 * std::sin(0.5 * static_cast<double>(phi));
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = 0.01 * density;
        state.mom_theta(radial, theta, phi) =
            0.09 * std::sin(0.6 * static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -0.07 * std::cos(0.4 * static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.7 + 0.16 * density;
        state.e_electron(radial, theta, phi) = 0.4 + 0.01 * density;
      }
    }
  }
}

[[nodiscard]] double MaxHydroDifference(
    const dec3d::state::HydroStateView& lhs,
    const dec3d::state::HydroStateView& rhs) noexcept {
  double max_difference = 0.0;
  for (std::size_t radial = 0; radial < lhs.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < lhs.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < lhs.rho->extent_phi(); ++phi) {
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.rho)(radial, theta, phi) - (*rhs.rho)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.mom_r)(radial, theta, phi) - (*rhs.mom_r)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.mom_theta)(radial, theta, phi) - (*rhs.mom_theta)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.mom_phi)(radial, theta, phi) - (*rhs.mom_phi)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(
                (*lhs.e_fluid_total)(radial, theta, phi) -
                (*rhs.e_fluid_total)(radial, theta, phi)));
      }
    }
  }
  return max_difference;
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

    constexpr std::size_t kRadialCells = 6u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;

    auto ppm_state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto first_order_state = CanonicalState::Create(ppm_state.layout);
    SeedAngularStructure(ppm_state);
    SeedAngularStructure(first_order_state);

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
    DEC3D_CHECK(geometry.is_valid());

    auto ppm_view = BuildHydroWorkView(ppm_state);
    auto first_order_view = BuildHydroWorkView(first_order_state);
    DEC3D_CHECK(ppm_view.is_complete());
    DEC3D_CHECK(first_order_view.is_complete());

    StaticGridHydroOptions first_order_options{};
    first_order_options.apply_geometric_source = false;
    first_order_options.apply_radial_sweep = false;
    first_order_options.apply_theta_sweep = true;
    first_order_options.apply_phi_sweep = true;
    first_order_options.use_macro_zoning = true;
    first_order_options.macro_zoning_coarse_factor = 0.5;

    auto ppm_options = first_order_options;
    ppm_options.use_ppm_reconstruction = true;
    ppm_options.reconstruction_ghost_layers = 3u;

    const auto first_order_result = AdvanceStaticGridHydro(
        first_order_view,
        geometry,
        1.0e-3,
        first_order_options);
    const auto ppm_result = AdvanceStaticGridHydro(
        ppm_view,
        geometry,
        1.0e-3,
        ppm_options);

    DEC3D_CHECK(first_order_result.is_complete());
    DEC3D_CHECK(ppm_result.is_complete());
    DEC3D_CHECK(std::abs(ppm_result.budget.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(ppm_result.budget.mom_r_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(ppm_result.budget.e_fluid_total_residual) < 1.0e-10);
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.macro_zoning.detected"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.macro_zoning.restrict.executed"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.macro_zoning.coarse_update.executed"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.macro_zoning.prolong.executed"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(ppm_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.phi"));
    DEC3D_CHECK(MaxHydroDifference(ppm_view, first_order_view) > 1.0e-10);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
