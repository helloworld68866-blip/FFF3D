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

namespace {

void SeedTrue3DLowModeCase(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  constexpr double kBaseRho = 1.0;
  constexpr double kBasePressure = 1.0;
  constexpr double kBaseElectronPressure = 0.4;
  constexpr double kAmplitude = 2.0e-2;

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      const double theta_center =
          0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double phi_center =
            0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
        const double low_mode = std::sin(theta_center) * std::cos(phi_center);
        const double v_r = kAmplitude * low_mode;
        const double kinetic = 0.5 * kBaseRho * v_r * v_r;

        state.rho(radial, theta, phi) = kBaseRho;
        state.mom_r(radial, theta, phi) = kBaseRho * v_r;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) =
            kBasePressure / (dec3d::state::HydroIdealGasGamma() - 1.0) + kinetic;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(kBaseElectronPressure);
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MakeFullDirectionalPpmOptions() noexcept {
  return dec3d::hydro::StaticGridHydroOptions{
      true,
      true,
      true,
      true,
      true,
      3u};
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MakeRadialOnlyPpmOptions() noexcept {
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

    constexpr std::size_t kRadialCells = 12u;
    constexpr std::size_t kThetaCells = 6u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kDt = 5.0e-4;
    constexpr std::size_t kShellIndex = 6u;

    const auto full_ppm_options = MakeFullDirectionalPpmOptions();
    const auto radial_only_ppm_options = MakeRadialOnlyPpmOptions();
    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = dec3d::state::CanonicalState::Create(before.layout);
    auto radial_only_after = dec3d::state::CanonicalState::Create(before.layout);
    SeedTrue3DLowModeCase(before, geometry);
    SeedTrue3DLowModeCase(single_rank_after, geometry);
    SeedTrue3DLowModeCase(radial_only_after, geometry);

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
        full_ppm_options);
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
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.theta"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.phi"));

    const auto writeback = dec3d::state::CommitHydroWriteback(
        local_state,
        local_hydro_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(writeback.is_complete());
    DEC3D_CHECK(writeback.success);

    dec3d::hydro::HydroBudgetResidualSummary multi_rank_budget;
    dec3d::testsupport::ReduceBudgetSummaryToRoot(local_result.budget, rank, &multi_rank_budget);

    auto gathered_after = dec3d::state::CanonicalState::Create(before.layout);
    dec3d::testsupport::GatherFieldToRoot(local_state.rho, decomposition, rank, gathered_after.rho);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_r, decomposition, rank, gathered_after.mom_r);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_theta, decomposition, rank, gathered_after.mom_theta);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_phi, decomposition, rank, gathered_after.mom_phi);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_fluid_total, decomposition, rank, gathered_after.e_fluid_total);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_electron, decomposition, rank, gathered_after.e_electron);

    int local_ok = 1;
    if (rank == 0) {
      dec3d::hydro::HydroOperator hydro;
      hydro.SetStaticGridOptions(full_ppm_options);
      const dec3d::core::StageContext context{
          0.0,
          kDt,
          1,
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          "mesh:p1.case3.low-mode-mpi-real-ppm",
          "ownership:rank0",
          "diagnostics:p1.hydro.case3.mpi-real-ppm"};
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
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          single_rank_stage_result.diagnostics,
          "p1.hydro.reconstruction.ppm.direction.theta"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          single_rank_stage_result.diagnostics,
          "p1.hydro.reconstruction.ppm.direction.phi"));

      dec3d::hydro::HydroBudgetResidualSummary single_rank_budget;
      DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
          single_rank_stage_result.diagnostics,
          &single_rank_budget));
      DEC3D_CHECK(single_rank_budget.is_complete());

      auto radial_only_view = dec3d::state::BuildHydroWorkView(radial_only_after);
      DEC3D_CHECK(radial_only_view.is_complete());
      const auto radial_only_result = dec3d::hydro::AdvanceStaticGridHydro(
          radial_only_view,
          geometry,
          kDt,
          radial_only_ppm_options);
      DEC3D_CHECK(radial_only_result.is_complete());
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          radial_only_result.diagnostics,
          "p1.hydro.reconstruction.ppm.executed"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          radial_only_result.diagnostics,
          "p1.hydro.reconstruction.ppm.direction.radial"));
      const auto radial_only_writeback = dec3d::state::CommitHydroWriteback(
          radial_only_after,
          radial_only_view,
          dec3d::state::BuildHydroAuthorizedWriteMask());
      DEC3D_CHECK(radial_only_writeback.is_complete());
      DEC3D_CHECK(radial_only_writeback.success);

      const auto sanity = dec3d::hydro::EvaluateTrue3DLowModeSanity(
          before,
          gathered_after,
          radial_only_after,
          geometry,
          kShellIndex,
          1.0e-8,
          1.0e-8);
      DEC3D_CHECK(sanity.is_complete());
      DEC3D_CHECK(sanity.success);

      auto mpi_parity = dec3d::hydro::EvaluateDirectionalMpiParity(
          single_rank_after,
          gathered_after,
          static_cast<std::size_t>(rank_count),
          1.0e-10);
      DEC3D_CHECK(mpi_parity.is_complete());
      mpi_parity.ppm_executed = true;
      mpi_parity.reconstruction_ghost_layers = kGhostLayers;
      mpi_parity.reconstruction_mode = "full_direction_ppm";
      DEC3D_CHECK(mpi_parity.success);

      DEC3D_CHECK(std::abs(single_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.e_fluid_total_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.e_fluid_total_residual) < 1.0e-10);

      const std::filesystem::path output_dir =
          "F:\\dec3d\\analysis\\output\\case3_true3d_low_mode_mpi_real_ppm";
      std::filesystem::remove_all(output_dir);
      DEC3D_CHECK(dec3d::hydro::WriteTrue3DLowModeOutputs(
          before,
          gathered_after,
          radial_only_after,
          geometry,
          kShellIndex,
          output_dir,
          sanity,
          &multi_rank_budget,
          &mpi_parity,
          &single_rank_budget,
          "case3_true3d_low_mode_mpi_real_ppm"));

      DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "mpi_parity_diagnostics.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual_comparison.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

      const auto summary_text = dec3d::testsupport::ReadTextFile(output_dir / "summary.txt");
      DEC3D_CHECK(summary_text.find("case=case3_true3d_low_mode_mpi_real_ppm") != std::string::npos);
      DEC3D_CHECK(summary_text.find("mpi_parity_within_tolerance=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("ppm_executed=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("reconstruction_ghost_layers=3") != std::string::npos);
      DEC3D_CHECK(summary_text.find("reconstruction_mode=full_direction_ppm") != std::string::npos);

      const auto manifest_text = dec3d::testsupport::ReadTextFile(output_dir / "case_manifest.txt");
      DEC3D_CHECK(manifest_text.find("mpi_parity_diagnostics=mpi_parity_diagnostics.txt") != std::string::npos);
      DEC3D_CHECK(manifest_text.find("budget_comparison=budget_residual_comparison.txt") != std::string::npos);
      DEC3D_CHECK(manifest_text.find("reconstruction_mode=full_direction_ppm") != std::string::npos);
      DEC3D_CHECK(manifest_text.find("reconstruction_ghost_layers=3") != std::string::npos);

      const auto parity_text =
          dec3d::testsupport::ReadTextFile(output_dir / "mpi_parity_diagnostics.txt");
      DEC3D_CHECK(parity_text.find("# kind=mpi_parity_diagnostics") != std::string::npos);
      DEC3D_CHECK(parity_text.find("ppm_executed=true") != std::string::npos);
      DEC3D_CHECK(parity_text.find("reconstruction_ghost_layers=3") != std::string::npos);
      DEC3D_CHECK(parity_text.find("reconstruction_mode=full_direction_ppm") != std::string::npos);
      DEC3D_CHECK(parity_text.find("max_velocity_magnitude_difference=") != std::string::npos);
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
