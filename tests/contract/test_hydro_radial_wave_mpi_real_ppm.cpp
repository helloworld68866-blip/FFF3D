#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"
#include "hydro_mpi_real_case_support.hpp"

#include <mpi.h>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

namespace {

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

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MakePpmRadialOnlyOptions() noexcept {
  return dec3d::hydro::StaticGridHydroOptions{
      true,
      true,
      false,
      false,
      true,
      3u};
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
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kDt = 1.0e-3;

    const auto ppm_options = MakePpmRadialOnlyOptions();
    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedGlobalRadialWave(single_rank_after, before);

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

    auto local_state = dec3d::testsupport::BuildLocalStateSlice(before, local_slice);

    const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.rho, kGhostLayers);
    const auto mom_r_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, kGhostLayers);
    const auto mom_theta_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, kGhostLayers);
    const auto mom_phi_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, kGhostLayers);
    const auto e_fluid_total_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, kGhostLayers);
    const auto e_electron_exchange =
        dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, kGhostLayers);

    DEC3D_CHECK(rho_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_r_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_theta_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(mom_phi_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(e_fluid_total_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(e_electron_exchange.is_complete(kThetaCells, kPhiCells, kGhostLayers));

    auto local_hydro_view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(local_hydro_view.is_complete());
    const auto local_override = dec3d::testsupport::BuildRadialGhostOverride(
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
        local_override,
        ppm_options);
    DEC3D_CHECK(local_result.is_complete());
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.radial"));

    const auto writeback = dec3d::state::CommitHydroWriteback(
        local_state,
        local_hydro_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(writeback.is_complete());
    DEC3D_CHECK(writeback.success);

    dec3d::hydro::HydroBudgetResidualSummary multi_rank_budget;
    dec3d::testsupport::ReduceBudgetSummaryToRoot(local_result.budget, rank, &multi_rank_budget);

    auto gathered_after = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    dec3d::testsupport::GatherFieldToRoot(local_state.rho, decomposition, rank, gathered_after.rho);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_r, decomposition, rank, gathered_after.mom_r);
    dec3d::testsupport::GatherFieldToRoot(
        local_state.mom_theta,
        decomposition,
        rank,
        gathered_after.mom_theta);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_phi, decomposition, rank, gathered_after.mom_phi);
    dec3d::testsupport::GatherFieldToRoot(
        local_state.e_fluid_total,
        decomposition,
        rank,
        gathered_after.e_fluid_total);
    dec3d::testsupport::GatherFieldToRoot(
        local_state.e_electron,
        decomposition,
        rank,
        gathered_after.e_electron);

    int local_ok = 1;
    if (rank == 0) {
      dec3d::hydro::HydroOperator hydro;
      hydro.SetStaticGridOptions(ppm_options);
      const dec3d::core::StageContext context{
          0.0,
          kDt,
          1,
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          "mesh:p1.radial-wave-mpi-real-ppm",
          "ownership:rank0",
          "diagnostics:p1.hydro.radial-wave-mpi-real-ppm"};
      DEC3D_CHECK(context.is_complete());
      DEC3D_CHECK(hydro.bind(context, geometry, single_rank_after));
      const auto single_rank_stage_result = hydro.advance();
      DEC3D_CHECK(single_rank_stage_result.is_semantically_complete());
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          single_rank_stage_result.diagnostics,
          "p1.hydro.reconstruction.ppm.executed"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          single_rank_stage_result.diagnostics,
          "p1.hydro.reconstruction.ppm.ng3"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          single_rank_stage_result.diagnostics,
          "p1.hydro.reconstruction.ppm.direction.radial"));

      dec3d::hydro::HydroBudgetResidualSummary single_rank_budget;
      DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
          single_rank_stage_result.diagnostics,
          &single_rank_budget));
      DEC3D_CHECK(single_rank_budget.is_complete());
      DEC3D_CHECK(std::abs(single_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.e_fluid_total_residual) < 1.0e-10);
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
      summary.ppm_executed = true;
      summary.reconstruction_ghost_layers = kGhostLayers;
      summary.reconstruction_mode = "radial_only_ppm";
      summary.success =
          summary.rank_decomposition_valid &&
          summary.parity_within_tolerance &&
          summary.seam_jumps_bounded;
      summary.report_line = "radial_wave_mpi_real_ppm_success=true";

      const std::filesystem::path output_dir =
          "F:\\dec3d\\analysis\\output\\case2_radial_wave_mpi_real_ppm";
      std::filesystem::remove_all(output_dir);
      DEC3D_CHECK(dec3d::hydro::WriteSphericalRadialWaveMpiParityOutputs(
          single_rank_after,
          gathered_after,
          geometry,
          output_dir,
          summary,
          &single_rank_budget,
          &multi_rank_budget,
          "case2_radial_wave_mpi_real_ppm"));
      DEC3D_CHECK(summary.success);
      DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "seam_diagnostics.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "single_rank_radial_profile_t1.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "multi_rank_radial_profile_t1.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "single_rank_budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "multi_rank_budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual_comparison.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

      const auto summary_text = dec3d::testsupport::ReadTextFile(output_dir / "summary.txt");
      DEC3D_CHECK(summary_text.find("case=case2_radial_wave_mpi_real_ppm") != std::string::npos);
      DEC3D_CHECK(summary_text.find("rank_count=3") != std::string::npos);
      DEC3D_CHECK(summary_text.find("success=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("ppm_executed=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("reconstruction_ghost_layers=3") != std::string::npos);
      DEC3D_CHECK(summary_text.find("reconstruction_mode=radial_only_ppm") != std::string::npos);

      const auto manifest =
          dec3d::testsupport::ReadTextFile(output_dir / "case_manifest.txt");
      DEC3D_CHECK(
          manifest.find("single_rank_budget=single_rank_budget_residual.txt") != std::string::npos);
      DEC3D_CHECK(
          manifest.find("multi_rank_budget=multi_rank_budget_residual.txt") != std::string::npos);
      DEC3D_CHECK(
          manifest.find("budget_comparison=budget_residual_comparison.txt") != std::string::npos);
      DEC3D_CHECK(manifest.find("reconstruction_mode=radial_only_ppm") != std::string::npos);
      DEC3D_CHECK(manifest.find("reconstruction_ghost_layers=3") != std::string::npos);

      const auto budget_comparison =
          dec3d::testsupport::ReadTextFile(output_dir / "budget_residual_comparison.txt");
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
