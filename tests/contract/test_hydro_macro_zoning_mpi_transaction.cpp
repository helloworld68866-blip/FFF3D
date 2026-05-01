#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"
#include "hydro_macro_zoning_mpi_support.hpp"

#include <mpi.h>

#include <iostream>
#include <vector>

namespace {

[[nodiscard]] bool SameStorage(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) noexcept {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    if (lhs[index] != rhs[index]) {
      return false;
    }
  }
  return true;
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

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto global_state = dec3d::state::CanonicalState::Create(before.layout);
    dec3d::testsupport::SeedMacroZoningRadialWave(global_state, before);

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
    DEC3D_CHECK(geometry.is_valid());

    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& local_slice = decomposition.slices[static_cast<std::size_t>(rank)];
    auto local_state = dec3d::testsupport::BuildLocalStateSlice(global_state, local_slice);
    const auto local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        local_slice,
        kThetaCells,
        kPhiCells);
    DEC3D_CHECK(local_geometry.is_valid());

    const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.rho, 3u);
    const auto mom_r_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, 3u);
    const auto mom_theta_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, 3u);
    const auto mom_phi_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, 3u);
    const auto e_fluid_total_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, 3u);
    const auto e_electron_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, 3u);

    DEC3D_CHECK(rho_exchange.is_complete(kThetaCells, kPhiCells, 3u));
    DEC3D_CHECK(mom_r_exchange.is_complete(kThetaCells, kPhiCells, 3u));
    DEC3D_CHECK(mom_theta_exchange.is_complete(kThetaCells, kPhiCells, 3u));
    DEC3D_CHECK(mom_phi_exchange.is_complete(kThetaCells, kPhiCells, 3u));
    DEC3D_CHECK(e_fluid_total_exchange.is_complete(kThetaCells, kPhiCells, 3u));
    DEC3D_CHECK(e_electron_exchange.is_complete(kThetaCells, kPhiCells, 3u));

    const auto snapshot_rho = local_state.rho.storage();
    const auto snapshot_mom_r = local_state.mom_r.storage();
    const auto snapshot_mom_theta = local_state.mom_theta.storage();
    const auto snapshot_mom_phi = local_state.mom_phi.storage();
    const auto snapshot_e_fluid_total = local_state.e_fluid_total.storage();
    const auto snapshot_e_electron = local_state.e_electron.storage();

    auto hydro_view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(hydro_view.is_complete());
    const auto local_override = dec3d::testsupport::BuildRadialGhostOverride(
        rho_exchange,
        mom_r_exchange,
        mom_theta_exchange,
        mom_phi_exchange,
        e_fluid_total_exchange,
        e_electron_exchange);

    auto options = dec3d::testsupport::MacroZoningActiveOptions();
    options.apply_radial_sweep = false;
    options.use_ppm_reconstruction = true;
    options.reconstruction_ghost_layers = 3u;

    const auto result = dec3d::hydro::AdvanceStaticGridHydro(
        hydro_view,
        local_geometry,
        1.0e-3,
        local_override,
        options);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics,
        "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics,
        "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(SameStorage(local_state.rho.storage(), snapshot_rho));
    DEC3D_CHECK(SameStorage(local_state.mom_r.storage(), snapshot_mom_r));
    DEC3D_CHECK(SameStorage(local_state.mom_theta.storage(), snapshot_mom_theta));
    DEC3D_CHECK(SameStorage(local_state.mom_phi.storage(), snapshot_mom_phi));
    DEC3D_CHECK(SameStorage(local_state.e_fluid_total.storage(), snapshot_e_fluid_total));
    DEC3D_CHECK(SameStorage(local_state.e_electron.storage(), snapshot_e_electron));

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
