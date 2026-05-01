#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

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

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.02 * static_cast<double>(radial) +
            0.003 * std::sin(static_cast<double>(theta + 1u)) +
            0.002 * std::cos(static_cast<double>(phi + 1u));
        const double v_r = 0.01 + 0.001 * static_cast<double>(radial);
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = density * v_r;
        state.mom_theta(radial, theta, phi) =
            0.002 * std::sin(0.5 * static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -0.001 * std::cos(0.25 * static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.0 + 0.05 * density;
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
            std::abs((*lhs.mom_theta)(radial, theta, phi) -
                     (*rhs.mom_theta)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.mom_phi)(radial, theta, phi) -
                     (*rhs.mom_phi)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.e_fluid_total)(radial, theta, phi) -
                     (*rhs.e_fluid_total)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.chi_e(radial, theta, phi) - rhs.chi_e(radial, theta, phi)));
      }
    }
  }
  return max_difference;
}

[[nodiscard]] double MaxStagedDifference(
    const dec3d::state::HydroStateView& view,
    const dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>& staged) noexcept {
  double max_difference = 0.0;
  for (std::size_t radial = 0; radial < staged.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < staged.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < staged.extent_phi(); ++phi) {
        const auto& cell = staged(radial, theta, phi);
        max_difference = std::max(
            max_difference,
            std::abs((*view.rho)(radial, theta, phi) - cell.rho));
        max_difference = std::max(
            max_difference,
            std::abs((*view.mom_r)(radial, theta, phi) - cell.mom_r));
        max_difference = std::max(
            max_difference,
            std::abs((*view.mom_theta)(radial, theta, phi) - cell.mom_theta));
        max_difference = std::max(
            max_difference,
            std::abs((*view.mom_phi)(radial, theta, phi) - cell.mom_phi));
        max_difference = std::max(
            max_difference,
            std::abs((*view.e_fluid_total)(radial, theta, phi) - cell.e_fluid_total));
        max_difference = std::max(
            max_difference,
            std::abs(view.chi_e(radial, theta, phi) - cell.chi_e));
      }
    }
  }
  return max_difference;
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BuildOptions() {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = true;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  options.use_macro_zoning = true;
  options.macro_zoning_coarse_factor = 0.5;
  options.apply_radial_ale_flux_correction = true;
  options.defer_macro_ale_compatibility_writeback = true;
  return options;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::PublishMacroAleStagedHydroState;
    using dec3d::mesh::ApplyRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    constexpr std::size_t kRadialCells = 8u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr double kDt = 2.0e-5;

    auto state =
        CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedState(state);
    const auto before_state = state;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.65});
    DEC3D_CHECK(geometry.is_valid());

    auto view = BuildHydroWorkView(state);
    auto before_view = BuildHydroWorkView(before_state);
    DEC3D_CHECK(view.is_complete());
    DEC3D_CHECK(before_view.is_complete());

    const auto proposal =
        BuildRadialAleMeshUpdateProposal(state.rho, state.mom_r, geometry, kDt);
    DEC3D_CHECK(proposal.is_complete_global(geometry.radial_faces.size()));

    auto expected_preview = geometry;
    const auto preview =
        ApplyRadialAleMeshUpdateProposal(proposal, expected_preview, 0u);
    DEC3D_CHECK(preview.success);

    auto options = BuildOptions();
    const auto result =
        AdvanceStaticGridHydro(view, geometry, kDt, options, &proposal);

    if (!result.success) {
      std::cerr << result.failure_reason << '\n';
    }
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.macro_ale_staged_writeback_available);
    DEC3D_CHECK(result.macro_ale_staged_geometry_available);
    DEC3D_CHECK(result.macro_ale_staged_geometry.radial_faces == expected_preview.radial_faces);
    DEC3D_CHECK(result.macro_ale_staged_cells.extent_r() == kRadialCells);
    DEC3D_CHECK(result.macro_ale_staged_cells.extent_theta() == kThetaCells);
    DEC3D_CHECK(result.macro_ale_staged_cells.extent_phi() == kPhiCells);
    DEC3D_CHECK(result.macro_ale_remap_diagnostics.is_complete());
    DEC3D_CHECK(MaxHydroDifference(view, before_view) == 0.0);
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.transaction.staged"));
    DEC3D_CHECK(!HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.transaction.clean"));

    std::string publish_failure;
    DEC3D_CHECK(PublishMacroAleStagedHydroState(view, result, &publish_failure));
    DEC3D_CHECK(publish_failure.empty());
    DEC3D_CHECK(MaxStagedDifference(view, result.macro_ale_staged_cells) < 1.0e-14);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
