#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

using dec3d::hydro::HydroConservativeState;

[[nodiscard]] std::string ReadTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::string content;
  std::string line;
  while (std::getline(input, line)) {
    content += line;
    content.push_back('\n');
  }

  return content;
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildLocalGeometrySlice(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::RadialOwnershipSlice& slice,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  dec3d::mesh::SphericalGeometryMetadata local_geometry;
  if (!geometry.is_valid() || !slice.initialized || slice.local_cell_count() == 0u) {
    return local_geometry;
  }

  local_geometry.valid = true;
  local_geometry.radial_faces.assign(
      geometry.radial_faces.begin() + static_cast<std::ptrdiff_t>(slice.begin_index),
      geometry.radial_faces.begin() + static_cast<std::ptrdiff_t>(slice.end_index + 1u));
  local_geometry.theta_faces = geometry.theta_faces;
  local_geometry.phi_faces = geometry.phi_faces;
  local_geometry.cell_volumes.reserve(slice.local_cell_count() * theta_cells * phi_cells);

  for (std::size_t radial = slice.begin_index; radial < slice.end_index; ++radial) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const std::size_t global_index =
            ((radial * theta_cells) + theta) * phi_cells + phi;
        const double volume = geometry.cell_volumes[global_index];
        local_geometry.cell_volumes.push_back(volume);
        local_geometry.global_volume += volume;
      }
    }
  }

  return local_geometry;
}

void SeedGlobalRadialWave(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before) {
  const dec3d::hydro::HydroPrimitiveState inner_state{
      1.0,
      0.0,
      0.0,
      0.0,
      1.0,
      std::pow(0.4, 3.0 / 5.0)};
  const dec3d::hydro::HydroPrimitiveState outer_state{
      0.125,
      0.0,
      0.0,
      0.0,
      0.1,
      std::pow(0.04, 3.0 / 5.0)};

  const auto inner_conservative = dec3d::hydro::MakeConservativeState(inner_state);
  const auto outer_conservative = dec3d::hydro::MakeConservativeState(outer_state);

  const std::size_t radial_midpoint = state.layout.radial_cells / 2u;
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const bool inside = radial < radial_midpoint;
    const auto& cell = inside ? inner_conservative : outer_conservative;
    const double electron_pressure = inside ? 0.4 : 0.04;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = cell.rho;
        state.mom_r(radial, theta, phi) = cell.mom_r;
        state.mom_theta(radial, theta, phi) = cell.mom_theta;
        state.mom_phi(radial, theta, phi) = cell.mom_phi;
        state.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);

        before.rho(radial, theta, phi) = cell.rho;
        before.mom_r(radial, theta, phi) = cell.mom_r;
        before.mom_theta(radial, theta, phi) = cell.mom_theta;
        before.mom_phi(radial, theta, phi) = cell.mom_phi;
        before.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
        before.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);
      }
    }
  }
}

[[nodiscard]] dec3d::state::CanonicalState BuildLocalStateSlice(
    const dec3d::state::CanonicalState& global_state,
    const dec3d::mesh::RadialOwnershipSlice& slice) {
  auto local_state = dec3d::state::CanonicalState::Create(
      {slice.local_cell_count(),
       global_state.layout.theta_cells,
       global_state.layout.phi_cells,
       global_state.layout.radiation_group_count});

  for (std::size_t local_radial = 0; local_radial < slice.local_cell_count(); ++local_radial) {
    const std::size_t global_radial = slice.begin_index + local_radial;
    for (std::size_t theta = 0; theta < global_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < global_state.layout.phi_cells; ++phi) {
        local_state.rho(local_radial, theta, phi) = global_state.rho(global_radial, theta, phi);
        local_state.mom_r(local_radial, theta, phi) = global_state.mom_r(global_radial, theta, phi);
        local_state.mom_theta(local_radial, theta, phi) =
            global_state.mom_theta(global_radial, theta, phi);
        local_state.mom_phi(local_radial, theta, phi) =
            global_state.mom_phi(global_radial, theta, phi);
        local_state.e_fluid_total(local_radial, theta, phi) =
            global_state.e_fluid_total(global_radial, theta, phi);
        local_state.e_electron(local_radial, theta, phi) =
            global_state.e_electron(global_radial, theta, phi);
      }
    }
  }

  return local_state;
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildRadialGhostOverride(
    const dec3d::mesh::RadialHaloFieldExchange& rho_exchange,
    const dec3d::mesh::RadialHaloFieldExchange& mom_r_exchange,
    const dec3d::mesh::RadialHaloFieldExchange& mom_theta_exchange,
    const dec3d::mesh::RadialHaloFieldExchange& mom_phi_exchange,
    const dec3d::mesh::RadialHaloFieldExchange& e_fluid_total_exchange,
    const dec3d::mesh::RadialHaloFieldExchange& e_electron_exchange) {
  dec3d::hydro::RadialGhostOverride override;
  override.has_inner_neighbor = rho_exchange.has_inner_neighbor;
  override.has_outer_neighbor = rho_exchange.has_outer_neighbor;
  override.ghost_layers = rho_exchange.ghost_layers;
  override.report_line = "mpi_radial_halo_exchange";

  const std::size_t ghost_layers = rho_exchange.ghost_layers;
  const std::size_t theta_cells = rho_exchange.has_inner_neighbor
                                      ? rho_exchange.inner_ghost_values.extent_theta()
                                      : rho_exchange.outer_ghost_values.extent_theta();
  const std::size_t phi_cells = rho_exchange.has_inner_neighbor
                                    ? rho_exchange.inner_ghost_values.extent_phi()
                                    : rho_exchange.outer_ghost_values.extent_phi();

  if (override.has_inner_neighbor) {
    override.inner_ghost_states = dec3d::core::Array3D<HydroConservativeState>(
        ghost_layers,
        theta_cells,
        phi_cells);
    for (std::size_t radial = 0; radial < ghost_layers; ++radial) {
      for (std::size_t theta = 0; theta < theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < phi_cells; ++phi) {
          override.inner_ghost_states(radial, theta, phi) = {
              rho_exchange.inner_ghost_values(radial, theta, phi),
              mom_r_exchange.inner_ghost_values(radial, theta, phi),
              mom_theta_exchange.inner_ghost_values(radial, theta, phi),
              mom_phi_exchange.inner_ghost_values(radial, theta, phi),
              e_fluid_total_exchange.inner_ghost_values(radial, theta, phi),
              dec3d::state::ChiEFromElectronPressure(
                  dec3d::state::ElectronPressureFromElectronEnergyDensity(
                      e_electron_exchange.inner_ghost_values(radial, theta, phi)))};
        }
      }
    }
  }

  if (override.has_outer_neighbor) {
    override.outer_ghost_states = dec3d::core::Array3D<HydroConservativeState>(
        ghost_layers,
        theta_cells,
        phi_cells);
    for (std::size_t radial = 0; radial < ghost_layers; ++radial) {
      for (std::size_t theta = 0; theta < theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < phi_cells; ++phi) {
          override.outer_ghost_states(radial, theta, phi) = {
              rho_exchange.outer_ghost_values(radial, theta, phi),
              mom_r_exchange.outer_ghost_values(radial, theta, phi),
              mom_theta_exchange.outer_ghost_values(radial, theta, phi),
              mom_phi_exchange.outer_ghost_values(radial, theta, phi),
              e_fluid_total_exchange.outer_ghost_values(radial, theta, phi),
              dec3d::state::ChiEFromElectronPressure(
                  dec3d::state::ElectronPressureFromElectronEnergyDensity(
                      e_electron_exchange.outer_ghost_values(radial, theta, phi)))};
        }
      }
    }
  }

  return override;
}

void GatherFieldToRoot(
    const dec3d::core::Array3D<double>& local_field,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    int rank,
    dec3d::core::Array3D<double>& global_field) {
  std::vector<int> recvcounts;
  std::vector<int> displacements;
  const auto plane_size =
      global_field.extent_theta() * global_field.extent_phi();

  if (rank == 0) {
    recvcounts.resize(decomposition.slices.size(), 0);
    displacements.resize(decomposition.slices.size(), 0);
    for (std::size_t slice_index = 0; slice_index < decomposition.slices.size(); ++slice_index) {
      const auto& slice = decomposition.slices[slice_index];
      recvcounts[slice_index] = static_cast<int>(slice.local_cell_count() * plane_size);
      displacements[slice_index] = static_cast<int>(slice.begin_index * plane_size);
    }
  }

  const int sendcount = static_cast<int>(local_field.size());
  MPI_Gatherv(
      const_cast<double*>(local_field.storage().data()),
      sendcount,
      MPI_DOUBLE,
      rank == 0 ? global_field.storage().data() : nullptr,
      rank == 0 ? recvcounts.data() : nullptr,
      rank == 0 ? displacements.data() : nullptr,
      MPI_DOUBLE,
      0,
      MPI_COMM_WORLD);
}

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

    constexpr std::size_t kRadialCells = 16u;
    constexpr std::size_t kThetaCells = 3u;
    constexpr std::size_t kPhiCells = 4u;
    constexpr double kDt = 1.0e-3;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedGlobalRadialWave(single_rank_after, before);

    const auto decomposition = dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& local_slice = decomposition.slices[static_cast<std::size_t>(rank)];
    const auto local_geometry = BuildLocalGeometrySlice(geometry, local_slice, kThetaCells, kPhiCells);
    DEC3D_CHECK(local_geometry.is_valid());

    auto local_state = BuildLocalStateSlice(before, local_slice);

    const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.rho, 1u);
    const auto mom_r_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, 1u);
    const auto mom_theta_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, 1u);
    const auto mom_phi_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, 1u);
    const auto e_fluid_total_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, 1u);
    const auto e_electron_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, 1u);

    DEC3D_CHECK(rho_exchange.is_complete(kThetaCells, kPhiCells, 1u));
    DEC3D_CHECK(mom_r_exchange.is_complete(kThetaCells, kPhiCells, 1u));
    DEC3D_CHECK(mom_theta_exchange.is_complete(kThetaCells, kPhiCells, 1u));
    DEC3D_CHECK(mom_phi_exchange.is_complete(kThetaCells, kPhiCells, 1u));
    DEC3D_CHECK(e_fluid_total_exchange.is_complete(kThetaCells, kPhiCells, 1u));
    DEC3D_CHECK(e_electron_exchange.is_complete(kThetaCells, kPhiCells, 1u));

    auto local_hydro_view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(local_hydro_view.is_complete());
    const auto local_override = BuildRadialGhostOverride(
        rho_exchange,
        mom_r_exchange,
        mom_theta_exchange,
        mom_phi_exchange,
        e_fluid_total_exchange,
        e_electron_exchange);

    const auto local_result = dec3d::hydro::AdvanceStaticGridHydro(
        local_hydro_view,
        local_geometry,
        kDt,
        local_override);
    DEC3D_CHECK(local_result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.ghost.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.budget.mass"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.budget.e_fluid_total"));

    const auto writeback = dec3d::state::CommitHydroWriteback(
        local_state,
        local_hydro_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(writeback.is_complete());
    DEC3D_CHECK(writeback.success);

    const double local_budget_values[] = {
        local_result.budget.old_mass,
        local_result.budget.flux_mass_delta,
        local_result.budget.source_mass_delta,
        local_result.budget.new_mass,
        local_result.budget.mass_residual,
        local_result.budget.old_e_fluid_total,
        local_result.budget.flux_e_fluid_total_delta,
        local_result.budget.source_e_fluid_total_delta,
        local_result.budget.new_e_fluid_total,
        local_result.budget.e_fluid_total_residual,
    };
    double reduced_budget_values[10] = {};
    MPI_Reduce(
        local_budget_values,
        reduced_budget_values,
        10,
        MPI_DOUBLE,
        MPI_SUM,
        0,
        MPI_COMM_WORLD);

    auto gathered_after = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    GatherFieldToRoot(local_state.rho, decomposition, rank, gathered_after.rho);
    GatherFieldToRoot(local_state.mom_r, decomposition, rank, gathered_after.mom_r);
    GatherFieldToRoot(local_state.mom_theta, decomposition, rank, gathered_after.mom_theta);
    GatherFieldToRoot(local_state.mom_phi, decomposition, rank, gathered_after.mom_phi);
    GatherFieldToRoot(local_state.e_fluid_total, decomposition, rank, gathered_after.e_fluid_total);
    GatherFieldToRoot(local_state.e_electron, decomposition, rank, gathered_after.e_electron);

    int local_ok = 1;
    if (rank == 0) {
      dec3d::hydro::HydroOperator hydro;
      const dec3d::core::StageContext context{
          0.0,
          kDt,
          1,
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          "mesh:p1.radial-wave-mpi-real",
          "ownership:rank0",
          "diagnostics:p1.hydro.radial-wave-mpi-real"};
      DEC3D_CHECK(context.is_complete());
      DEC3D_CHECK(hydro.bind(context, geometry, single_rank_after));
      const auto single_rank_stage_result = hydro.advance();
      DEC3D_CHECK(single_rank_stage_result.is_semantically_complete());
      dec3d::hydro::HydroBudgetResidualSummary single_rank_budget;
      DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
          single_rank_stage_result.diagnostics,
          &single_rank_budget));
      DEC3D_CHECK(single_rank_budget.is_complete());
      DEC3D_CHECK(std::abs(single_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.e_fluid_total_residual) < 1.0e-10);

      dec3d::hydro::HydroBudgetResidualSummary multi_rank_budget{};
      multi_rank_budget.old_mass = reduced_budget_values[0];
      multi_rank_budget.flux_mass_delta = reduced_budget_values[1];
      multi_rank_budget.source_mass_delta = reduced_budget_values[2];
      multi_rank_budget.new_mass = reduced_budget_values[3];
      multi_rank_budget.mass_residual = reduced_budget_values[4];
      multi_rank_budget.old_e_fluid_total = reduced_budget_values[5];
      multi_rank_budget.flux_e_fluid_total_delta = reduced_budget_values[6];
      multi_rank_budget.source_e_fluid_total_delta = reduced_budget_values[7];
      multi_rank_budget.new_e_fluid_total = reduced_budget_values[8];
      multi_rank_budget.e_fluid_total_residual = reduced_budget_values[9];
      multi_rank_budget.report_line =
          "mpi_aggregated_budget_mass_residual=" + std::to_string(multi_rank_budget.mass_residual) +
          "; mpi_aggregated_budget_e_residual=" + std::to_string(multi_rank_budget.e_fluid_total_residual);
      DEC3D_CHECK(std::abs(multi_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.e_fluid_total_residual) < 1.0e-10);

      dec3d::hydro::SphericalRadialWaveMpiParityCheck summary;
      summary.rank_count = static_cast<std::size_t>(rank_count);
      summary.rank_decomposition_valid = decomposition.is_valid();
      for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
        for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
          for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
            summary.max_rho_difference = std::max(
                summary.max_rho_difference,
                std::abs(single_rank_after.rho(radial, theta, phi) -
                         gathered_after.rho(radial, theta, phi)));
            summary.max_mom_r_difference = std::max(
                summary.max_mom_r_difference,
                std::abs(single_rank_after.mom_r(radial, theta, phi) -
                         gathered_after.mom_r(radial, theta, phi)));
            summary.max_e_fluid_total_difference = std::max(
                summary.max_e_fluid_total_difference,
                std::abs(single_rank_after.e_fluid_total(radial, theta, phi) -
                         gathered_after.e_fluid_total(radial, theta, phi)));
          }
        }
      }

      for (std::size_t slice_index = 0; slice_index + 1u < decomposition.slices.size(); ++slice_index) {
        const auto& slice = decomposition.slices[slice_index];
        const std::size_t left_radial = slice.end_index - 1u;
        const std::size_t right_radial = slice.end_index;
        for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
          for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
            const double single_rho_jump =
                std::abs(single_rank_after.rho(left_radial, theta, phi) -
                         single_rank_after.rho(right_radial, theta, phi));
            const double gathered_rho_jump =
                std::abs(gathered_after.rho(left_radial, theta, phi) -
                         gathered_after.rho(right_radial, theta, phi));
            summary.max_seam_rho_jump = std::max(
                summary.max_seam_rho_jump,
                std::abs(gathered_rho_jump - single_rho_jump));

            const double single_velocity_jump =
                std::abs(
                    std::sqrt(std::pow(single_rank_after.mom_r(left_radial, theta, phi) /
                                           single_rank_after.rho(left_radial, theta, phi),
                                       2.0)) -
                    std::sqrt(std::pow(single_rank_after.mom_r(right_radial, theta, phi) /
                                           single_rank_after.rho(right_radial, theta, phi),
                                       2.0)));
            const double gathered_velocity_jump =
                std::abs(
                    std::sqrt(std::pow(gathered_after.mom_r(left_radial, theta, phi) /
                                           gathered_after.rho(left_radial, theta, phi),
                                       2.0)) -
                    std::sqrt(std::pow(gathered_after.mom_r(right_radial, theta, phi) /
                                           gathered_after.rho(right_radial, theta, phi),
                                       2.0)));
            summary.max_seam_velocity_jump = std::max(
                summary.max_seam_velocity_jump,
                std::abs(gathered_velocity_jump - single_velocity_jump));
          }
        }
      }

      summary.parity_within_tolerance =
          summary.max_rho_difference <= 1.0e-10 &&
          summary.max_mom_r_difference <= 1.0e-10 &&
          summary.max_e_fluid_total_difference <= 1.0e-10;
      summary.seam_jumps_bounded =
          summary.max_seam_rho_jump <= 1.0e-10 &&
          summary.max_seam_velocity_jump <= 1.0e-10;
      summary.success =
          summary.rank_decomposition_valid &&
          summary.parity_within_tolerance &&
          summary.seam_jumps_bounded;
      summary.report_line =
          "radial_wave_mpi_real_success=true";

      const std::filesystem::path output_dir =
          "F:\\dec3d\\analysis\\output\\case2_radial_wave_mpi_real";
      std::filesystem::remove_all(output_dir);
      DEC3D_CHECK(dec3d::hydro::WriteSphericalRadialWaveMpiParityOutputs(
          single_rank_after,
          gathered_after,
          geometry,
          output_dir,
          summary,
          &single_rank_budget,
          &multi_rank_budget,
          "case2_radial_wave_mpi_real"));
      DEC3D_CHECK(summary.success);
      DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "seam_diagnostics.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "single_rank_radial_profile_t1.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "multi_rank_radial_profile_t1.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "single_rank_budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "multi_rank_budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual_comparison.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

      const auto summary_text = ReadTextFile(output_dir / "summary.txt");
      DEC3D_CHECK(summary_text.find("case=case2_radial_wave_mpi_real") != std::string::npos);
      DEC3D_CHECK(summary_text.find("rank_count=3") != std::string::npos);
      DEC3D_CHECK(summary_text.find("success=true") != std::string::npos);

      const auto manifest = ReadTextFile(output_dir / "case_manifest.txt");
      DEC3D_CHECK(manifest.find("single_rank_budget=single_rank_budget_residual.txt") != std::string::npos);
      DEC3D_CHECK(manifest.find("multi_rank_budget=multi_rank_budget_residual.txt") != std::string::npos);
      DEC3D_CHECK(manifest.find("budget_comparison=budget_residual_comparison.txt") != std::string::npos);

      const auto budget_comparison = ReadTextFile(output_dir / "budget_residual_comparison.txt");
      DEC3D_CHECK(budget_comparison.find("# kind=budget_residual_comparison") != std::string::npos);
      DEC3D_CHECK(budget_comparison.find("mass_residual_difference=") != std::string::npos);
    }

    MPI_Bcast(&local_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
    return local_ok == 1 ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
  }
}
