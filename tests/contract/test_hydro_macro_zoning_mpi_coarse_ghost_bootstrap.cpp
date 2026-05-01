#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"
#include "hydro_macro_zoning_mpi_support.hpp"
#include "hydro_mpi_real_case_support.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

[[nodiscard]] const dec3d::core::DiagnosticMessage* FindDiagnostic(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return &entry;
    }
  }
  return nullptr;
}

void ScaleGhostStates(
    dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>& states,
    double scale) noexcept {
  for (std::size_t radial = 0; radial < states.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < states.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < states.extent_phi(); ++phi) {
        auto& state = states(radial, theta, phi);
        state.rho *= scale;
        state.mom_r *= scale;
        state.mom_theta *= scale;
        state.mom_phi *= scale;
        state.e_fluid_total *= scale;
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride MakePoisonedOverride(
    dec3d::hydro::RadialGhostOverride ghost_override) noexcept {
  if (ghost_override.has_inner_neighbor) {
    ScaleGhostStates(ghost_override.inner_ghost_states, 0.70);
  }
  if (ghost_override.has_outer_neighbor) {
    ScaleGhostStates(ghost_override.outer_ghost_states, 1.35);
  }
  ghost_override.report_line = "mpi_radial_halo_exchange_poisoned_for_consumption_check";
  return ghost_override;
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MakeMacroRadialPpmOptions() noexcept {
  auto options = dec3d::testsupport::MacroZoningActiveOptions();
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = false;
  options.apply_phi_sweep = false;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  return options;
}

[[nodiscard]] double MaxHydroDifference(
    const dec3d::state::CanonicalState& lhs,
    const dec3d::state::CanonicalState& rhs) noexcept {
  double max_difference = 0.0;
  for (std::size_t radial = 0; radial < lhs.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < lhs.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < lhs.layout.phi_cells; ++phi) {
        max_difference = std::max(
            max_difference,
            std::abs(lhs.rho(radial, theta, phi) - rhs.rho(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.mom_r(radial, theta, phi) - rhs.mom_r(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.e_fluid_total(radial, theta, phi) - rhs.e_fluid_total(radial, theta, phi)));
      }
    }
  }
  return max_difference;
}

}  // namespace

int main(int argc, char** argv) {
  int mpi_initialized = 0;

  try {
    MPI_Init(&argc, &argv);
    MPI_Initialized(&mpi_initialized);

    int rank = 0;
    int rank_count = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_count);
    DEC3D_CHECK_EQ(rank_count, 3);

    constexpr std::size_t kRadialCells = 12u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kDt = 1.0e-3;

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto seeded = dec3d::state::CanonicalState::Create(before.layout);
    dec3d::testsupport::SeedMacroZoningRadialWave(seeded, before);

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
    DEC3D_CHECK(geometry.is_valid());

    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& local_slice = decomposition.slices[static_cast<std::size_t>(rank)];
    const auto local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        local_slice,
        kThetaCells,
        kPhiCells);
    DEC3D_CHECK(local_geometry.is_valid());

    auto correct_state = dec3d::testsupport::BuildLocalStateSlice(before, local_slice);
    auto poisoned_state = dec3d::testsupport::BuildLocalStateSlice(before, local_slice);

    const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(correct_state.rho, kGhostLayers);
    const auto mom_r_exchange =
        dec3d::mesh::ExchangeRadialHaloField(correct_state.mom_r, kGhostLayers);
    const auto mom_theta_exchange =
        dec3d::mesh::ExchangeRadialHaloField(correct_state.mom_theta, kGhostLayers);
    const auto mom_phi_exchange =
        dec3d::mesh::ExchangeRadialHaloField(correct_state.mom_phi, kGhostLayers);
    const auto e_fluid_total_exchange =
        dec3d::mesh::ExchangeRadialHaloField(correct_state.e_fluid_total, kGhostLayers);
    const auto e_electron_exchange =
        dec3d::mesh::ExchangeRadialHaloField(correct_state.e_electron, kGhostLayers);

    DEC3D_CHECK(rho_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_r_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_theta_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_phi_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(e_fluid_total_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(e_electron_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));

    const auto correct_override = dec3d::testsupport::BuildRadialGhostOverride(
        rho_exchange,
        mom_r_exchange,
        mom_theta_exchange,
        mom_phi_exchange,
        e_fluid_total_exchange,
        e_electron_exchange);
    const auto poisoned_override = MakePoisonedOverride(correct_override);
    const auto options = MakeMacroRadialPpmOptions();

    auto correct_view = dec3d::state::BuildHydroWorkView(correct_state);
    auto poisoned_view = dec3d::state::BuildHydroWorkView(poisoned_state);
    DEC3D_CHECK(correct_view.is_complete());
    DEC3D_CHECK(poisoned_view.is_complete());

    const auto correct_result = dec3d::hydro::AdvanceStaticGridHydro(
        correct_view,
        local_geometry,
        kDt,
        correct_override,
        options);
    const auto poisoned_result = dec3d::hydro::AdvanceStaticGridHydro(
        poisoned_view,
        local_geometry,
        kDt,
        poisoned_override,
        options);

    DEC3D_CHECK(correct_result.is_complete());
    DEC3D_CHECK(poisoned_result.is_complete());
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        correct_result.diagnostics,
        "p1.hydro.macro_zoning.detected"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        correct_result.diagnostics,
        "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        correct_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.radial"));
    DEC3D_CHECK(!dec3d::testsupport::HasDiagnosticCode(
        correct_result.diagnostics,
        "p1.hydro.ghost.direction.radial"));

    const auto* coarse_update = FindDiagnostic(
        correct_result.diagnostics,
        "p1.hydro.macro_zoning.coarse_update.executed");
    DEC3D_CHECK(coarse_update != nullptr);
    DEC3D_CHECK(coarse_update->message.find("radial_executed=true") != std::string::npos);

    const auto correct_writeback = dec3d::state::CommitHydroWriteback(
        correct_state,
        correct_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    const auto poisoned_writeback = dec3d::state::CommitHydroWriteback(
        poisoned_state,
        poisoned_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(correct_writeback.is_complete());
    DEC3D_CHECK(poisoned_writeback.is_complete());

    DEC3D_CHECK(MaxHydroDifference(correct_state, poisoned_state) > 1.0e-12);

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
  }
}
