#include "hydro/driver/macro_ale_mpi_compatibility.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "hydro_mpi_real_case_support.hpp"
#include "mesh/ale/radial_ale_mpi.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr std::size_t kRadialCells = 24u;
constexpr std::size_t kThetaCells = 8u;
constexpr std::size_t kPhiCells = 8u;
constexpr std::size_t kGhostLayers = 3u;
constexpr double kDt = 5.0e-6;

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.01 * static_cast<double>(radial) +
            0.001 * std::sin(static_cast<double>(theta + 1u)) *
                std::cos(static_cast<double>(phi + 1u));
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) =
            density * (0.004 + 0.0002 * static_cast<double>(radial));
        state.mom_theta(radial, theta, phi) =
            1.0e-9 * std::sin(static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -1.0e-9 * std::cos(static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.4 + 0.02 * density;
        state.e_electron(radial, theta, phi) = 0.45 + 0.01 * density;
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BuildOptions() {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = true;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = kGhostLayers;
  options.use_macro_zoning = true;
  options.macro_zoning_coarse_factor = 0.5;
  options.apply_radial_ale_flux_correction = true;
  return options;
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildOverride(
    const dec3d::state::HydroStateView& view) {
  const auto rho_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.rho, kGhostLayers);
  const auto mom_r_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.mom_r, kGhostLayers);
  const auto mom_theta_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.mom_theta, kGhostLayers);
  const auto mom_phi_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.mom_phi, kGhostLayers);
  const auto e_fluid_total_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.e_fluid_total, kGhostLayers);
  const auto e_electron_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*view.e_electron, kGhostLayers);

  DEC3D_CHECK(rho_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
  DEC3D_CHECK(mom_r_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
  DEC3D_CHECK(mom_theta_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
  DEC3D_CHECK(mom_phi_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
  DEC3D_CHECK(e_fluid_total_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
  DEC3D_CHECK(e_electron_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));

  return dec3d::testsupport::BuildRadialGhostOverride(
      rho_exchange,
      mom_r_exchange,
      mom_theta_exchange,
      mom_phi_exchange,
      e_fluid_total_exchange,
      e_electron_exchange);
}

void AssertPhysical(const dec3d::state::HydroStateView& view) {
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const dec3d::hydro::HydroConservativeState state{
            (*view.rho)(radial, theta, phi),
            (*view.mom_r)(radial, theta, phi),
            (*view.mom_theta)(radial, theta, phi),
            (*view.mom_phi)(radial, theta, phi),
            (*view.e_fluid_total)(radial, theta, phi),
            view.chi_e(radial, theta, phi)};
        DEC3D_CHECK(state.is_finite());
        DEC3D_CHECK(dec3d::hydro::RecoverPrimitiveState(state).is_physical());
      }
    }
  }
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

    auto global_state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedState(global_state);
    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
    auto local_state =
        dec3d::testsupport::BuildLocalStateSlice(global_state, slice);
    auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.75});
    DEC3D_CHECK(geometry.is_valid());
    auto local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        slice,
        kThetaCells,
        kPhiCells);
    DEC3D_CHECK(local_geometry.is_valid());
    auto view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(view.is_complete());

    std::vector<double> global_radial_faces = geometry.radial_faces;
    for (int step = 0; step < 3; ++step) {
      const auto proposal = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
          *view.rho,
          *view.mom_r,
          decomposition,
          global_radial_faces,
          kDt,
          1.0,
          MPI_COMM_WORLD);
      DEC3D_CHECK(proposal.is_complete());
      const auto result = dec3d::hydro::AdvanceMacroAleMpiCompatibilityStep(
          view,
          local_geometry,
          decomposition,
          proposal.proposal,
          kDt,
          BuildOverride(view),
          BuildOptions(),
          MPI_COMM_WORLD);
      if (!result.success && rank == 0) {
        std::cerr << result.failure_reason << '\n';
      }
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          result.diagnostics_payload,
          "p1.hydro.macro_ale_mpi.transaction.clean"));
      DEC3D_CHECK(result.diagnostics.published_from_staged_geometry);
      AssertPhysical(view);
      global_radial_faces = proposal.proposal.proposed_radial_faces;
    }

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
