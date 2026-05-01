#include "hydro/driver/macro_ale_hllc_mpi.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro_mpi_real_case_support.hpp"
#include "mesh/ale/radial_ale_mpi.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

constexpr std::size_t kRadialCells = 12u;
constexpr std::size_t kThetaCells = 8u;
constexpr std::size_t kPhiCells = 8u;
constexpr std::size_t kGhostLayers = 3u;
constexpr double kDt = 1.0e-5;

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double theta_mode = std::sin(0.4 * static_cast<double>(theta + 1u));
        const double phi_mode = std::cos(0.3 * static_cast<double>(phi + 1u));
        const double density =
            1.0 + 0.02 * static_cast<double>(radial) +
            0.002 * theta_mode * phi_mode;
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) =
            density * (0.005 + 0.0002 * static_cast<double>(radial));
        state.mom_theta(radial, theta, phi) = 1.0e-8 * density * theta_mode;
        state.mom_phi(radial, theta, phi) = -1.0e-8 * density * phi_mode;
        state.e_fluid_total(radial, theta, phi) = 2.4 + 0.05 * density;
        state.e_electron(radial, theta, phi) = 0.45 + 0.01 * density;
        for (std::size_t group = 0; group < state.radiation_groups.size(); ++group) {
          state.radiation_groups[group](radial, theta, phi) =
              1.0 + 0.1 * static_cast<double>(group + 1u) +
              0.02 * static_cast<double>(radial);
        }
      }
    }
  }
}

[[nodiscard]] double RadiationChiFromUg(double ug) noexcept {
  return std::pow(ug / 3.0, 3.0 / 4.0);
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
  options.use_macro_ale_direct_moving_face_hllc = true;
  return options;
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildOverride(
    const dec3d::state::CanonicalState& local_state) {
  const auto theta_cells = local_state.rho.extent_theta();
  const auto phi_cells = local_state.rho.extent_phi();
  const auto rho_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.rho, kGhostLayers);
  const auto mom_r_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, kGhostLayers);
  const auto mom_theta_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, kGhostLayers);
  const auto mom_phi_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, kGhostLayers);
  const auto e_fluid_total_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, kGhostLayers);
  const auto e_electron_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, kGhostLayers);
  std::vector<dec3d::mesh::RadialHaloFieldExchange> radiation_exchanges;
  radiation_exchanges.reserve(local_state.radiation_groups.size());
  for (const auto& group : local_state.radiation_groups) {
    radiation_exchanges.push_back(
        dec3d::mesh::ExchangeRadialHaloField(group, kGhostLayers));
  }

  DEC3D_CHECK(rho_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_r_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_theta_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_phi_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(e_fluid_total_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(e_electron_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  for (const auto& exchange : radiation_exchanges) {
    DEC3D_CHECK(exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  }

  auto override = dec3d::testsupport::BuildRadialGhostOverride(
      rho_exchange,
      mom_r_exchange,
      mom_theta_exchange,
      mom_phi_exchange,
      e_fluid_total_exchange,
      e_electron_exchange);
  for (std::size_t group = 0; group < radiation_exchanges.size(); ++group) {
    const auto& exchange = radiation_exchanges[group];
    if (override.has_inner_neighbor) {
      for (std::size_t radial = 0; radial < kGhostLayers; ++radial) {
        for (std::size_t theta = 0; theta < theta_cells; ++theta) {
          for (std::size_t phi = 0; phi < phi_cells; ++phi) {
            override.inner_ghost_states(radial, theta, phi).radiation_chi.push_back(
                RadiationChiFromUg(exchange.inner_ghost_values(radial, theta, phi)));
          }
        }
      }
    }
    if (override.has_outer_neighbor) {
      for (std::size_t radial = 0; radial < kGhostLayers; ++radial) {
        for (std::size_t theta = 0; theta < theta_cells; ++theta) {
          for (std::size_t phi = 0; phi < phi_cells; ++phi) {
            override.outer_ghost_states(radial, theta, phi).radiation_chi.push_back(
                RadiationChiFromUg(exchange.outer_ghost_values(radial, theta, phi)));
          }
        }
      }
    }
  }
  DEC3D_CHECK(override.is_complete(theta_cells, phi_cells, kGhostLayers));
  return override;
}

[[nodiscard]] double MaxFieldDifference(
    const dec3d::core::Array3D<double>& lhs,
    const dec3d::core::Array3D<double>& rhs) noexcept {
  double max_difference = 0.0;
  for (std::size_t radial = 0; radial < lhs.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < lhs.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < lhs.extent_phi(); ++phi) {
        max_difference = std::max(
            max_difference,
            std::abs(lhs(radial, theta, phi) - rhs(radial, theta, phi)));
      }
    }
  }
  return max_difference;
}

[[nodiscard]] double SumField(
    const dec3d::core::Array3D<double>& field) noexcept {
  double total = 0.0;
  for (double value : field.storage()) {
    total += value;
  }
  return total;
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
        dec3d::state::CanonicalStateLayout{
            kRadialCells,
            kThetaCells,
            kPhiCells,
            1u});
    SeedState(global_state);

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{
            kRadialCells,
            kThetaCells,
            kPhiCells,
            0.05,
            0.65});
    DEC3D_CHECK(geometry.is_valid());
    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];

    auto local_state =
        dec3d::testsupport::BuildLocalStateSlice(global_state, slice);
    auto local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        slice,
        kThetaCells,
        kPhiCells);
    auto local_view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(local_view.is_complete());
    auto options = BuildOptions();
    options.enable_radiation_hydro_terms = true;

    const auto proposal_result = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
        local_state.rho,
        local_state.mom_r,
        decomposition,
        geometry.radial_faces,
        kDt,
        1.0,
        MPI_COMM_WORLD);
    DEC3D_CHECK(proposal_result.is_complete());

    const auto result = dec3d::hydro::AdvanceMacroAleHllcMpiStep(
        local_view,
        local_geometry,
        decomposition,
        proposal_result.proposal,
        kDt,
        BuildOverride(local_state),
        options,
        MPI_COMM_WORLD);
    if (!result.success && rank == 0) {
      std::cerr << result.failure_reason << '\n';
    }

    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.diagnostics.is_complete());
    DEC3D_CHECK(result.staged_hydro_result.macro_ale_hllc_executed);
    DEC3D_CHECK(result.staged_hydro_result.macro_ale_hllc_local_face_window.is_complete());
    DEC3D_CHECK(result.diagnostics.branch_invalid_count == 0u);
    DEC3D_CHECK(result.diagnostics.max_seam_w_face_delta <= 1.0e-14);
    DEC3D_CHECK(result.diagnostics.max_seam_proposed_radius_delta <= 1.0e-14);
    DEC3D_CHECK(result.diagnostics.max_seam_preview_radius_delta <= 1.0e-14);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_mass_delta <= 1.0e-10);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_mom_r_delta <= 1.0e-10);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_mom_theta_delta <= 1.0e-10);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_mom_phi_delta <= 1.0e-10);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_e_fluid_total_delta <= 1.0e-9);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_chi_e_delta <= 1.0e-10);
    DEC3D_CHECK(result.diagnostics.max_seam_flux_radiation_chi_delta <= 1.0e-10);
    DEC3D_CHECK_EQ(result.diagnostics.radiation_group_count, std::size_t{1});
    DEC3D_CHECK(result.diagnostics.global_stage_ok);
    DEC3D_CHECK(result.diagnostics.global_publish_ok);
    DEC3D_CHECK(result.diagnostics.no_partial_canonical_writeback);
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_hllc.seam_flux_equality"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_hllc.global_transaction.preflight"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_hllc.global_conservation"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_hllc.transaction.clean"));
    DEC3D_CHECK(result.diagnostics.report_line.find(
                    "radiation_scalar_bundle_transport=enabled") != std::string::npos);
    DEC3D_CHECK(result.diagnostics.report_line.find(
                    "global_radiation_chi_residual=") != std::string::npos);
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.staged_hydro_result.diagnostics,
        "p1.hydro.macro_ale_hllc.moving_interface_branch"));
    DEC3D_CHECK(!dec3d::testsupport::HasDiagnosticCode(
        result.staged_hydro_result.diagnostics,
        "p1.hydro.macro_ale.compatibility.executed"));

    dec3d::core::Array3D<double> gathered_rho(kRadialCells, kThetaCells, kPhiCells);
    dec3d::core::Array3D<double> gathered_mom_r(kRadialCells, kThetaCells, kPhiCells);
    dec3d::core::Array3D<double> gathered_mom_theta(kRadialCells, kThetaCells, kPhiCells);
    dec3d::core::Array3D<double> gathered_mom_phi(kRadialCells, kThetaCells, kPhiCells);
    dec3d::core::Array3D<double> gathered_e(kRadialCells, kThetaCells, kPhiCells);
    dec3d::core::Array3D<double> gathered_chi_e(kRadialCells, kThetaCells, kPhiCells);
    dec3d::testsupport::GatherFieldToRoot(
        *local_view.rho,
        decomposition,
        rank,
        gathered_rho);
    dec3d::testsupport::GatherFieldToRoot(
        *local_view.mom_r,
        decomposition,
        rank,
        gathered_mom_r);
    dec3d::testsupport::GatherFieldToRoot(
        *local_view.mom_theta,
        decomposition,
        rank,
        gathered_mom_theta);
    dec3d::testsupport::GatherFieldToRoot(
        *local_view.mom_phi,
        decomposition,
        rank,
        gathered_mom_phi);
    dec3d::testsupport::GatherFieldToRoot(
        *local_view.e_fluid_total,
        decomposition,
        rank,
        gathered_e);
    dec3d::testsupport::GatherFieldToRoot(
        local_view.chi_e,
        decomposition,
        rank,
        gathered_chi_e);

    if (rank == 0) {
      auto single_state = global_state;
      auto single_view = dec3d::state::BuildHydroWorkView(single_state);
      DEC3D_CHECK(single_view.is_complete());
      const auto single_result = dec3d::hydro::AdvanceStaticGridHydro(
          single_view,
          geometry,
          kDt,
          BuildOptions(),
          &proposal_result.proposal);
      if (!single_result.success) {
        std::cerr << single_result.failure_reason << '\n';
      }
      DEC3D_CHECK(single_result.success);
      double max_difference = 0.0;
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_rho, *single_view.rho));
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_mom_r, *single_view.mom_r));
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_mom_theta, *single_view.mom_theta));
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_mom_phi, *single_view.mom_phi));
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_e, *single_view.e_fluid_total));
      max_difference = std::max(
          max_difference,
          MaxFieldDifference(gathered_chi_e, single_view.chi_e));
      DEC3D_CHECK(max_difference < 5.0e-9);
      DEC3D_CHECK(std::abs(SumField(gathered_rho) - SumField(*single_view.rho)) < 1.0e-9);
      DEC3D_CHECK(std::abs(SumField(gathered_e) - SumField(*single_view.e_fluid_total)) < 1.0e-8);
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
