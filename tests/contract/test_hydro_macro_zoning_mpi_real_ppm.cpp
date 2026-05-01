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
#include "hydro_macro_zoning_mpi_support.hpp"

#include <mpi.h>

#include <filesystem>
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

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions MakeMacroPpmOptions() noexcept {
  auto options = dec3d::testsupport::MacroZoningActiveOptions();
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  return options;
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
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kDt = 1.0e-3;

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = dec3d::state::CanonicalState::Create(before.layout);
    dec3d::testsupport::SeedMacroZoningRadialWave(single_rank_after, before);

    const auto options = MakeMacroPpmOptions();
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
        options);
    if (!local_result.success) {
      throw std::runtime_error(
          local_result.failure_reason.empty()
              ? "local macro-zoning mpi real ppm result failed"
              : local_result.failure_reason);
    }
    DEC3D_CHECK(local_result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.macro_zoning.detected"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.macro_zoning.coarse_update.executed"));
    DEC3D_CHECK(HasDiagnosticCode(
        local_result.diagnostics,
        "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(local_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.phi"));

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
      const dec3d::core::StageContext context{
          0.0,
          kDt,
          1,
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          "mesh:p1.case2.macro-zoning-mpi-real-ppm",
          "ownership:rank0",
          "diagnostics:p1.hydro.case2.macro-zoning.mpi-real.ppm"};
      DEC3D_CHECK(context.is_complete());
      DEC3D_CHECK(hydro.bind(context, geometry, single_rank_after));
      hydro.SetStaticGridOptions(options);
      const auto single_rank_result = hydro.advance();
      DEC3D_CHECK(single_rank_result.is_semantically_complete());
      DEC3D_CHECK(HasDiagnosticCode(single_rank_result.diagnostics, "p1.hydro.stage.macro_zoning"));
      DEC3D_CHECK(HasDiagnosticCode(single_rank_result.diagnostics, "p1.hydro.reconstruction.ppm.executed"));
      DEC3D_CHECK(HasDiagnosticCode(single_rank_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.theta"));
      DEC3D_CHECK(HasDiagnosticCode(single_rank_result.diagnostics, "p1.hydro.reconstruction.ppm.direction.phi"));

      dec3d::hydro::HydroBudgetResidualSummary single_rank_budget;
      DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
          single_rank_result.diagnostics,
          &single_rank_budget));
      DEC3D_CHECK(single_rank_budget.is_complete());
      DEC3D_CHECK(std::abs(single_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.mom_r_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.e_fluid_total_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.mom_r_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.e_fluid_total_residual) < 1.0e-10);

      auto summary = dec3d::testsupport::BuildMacroZoningParitySummary(
          single_rank_after,
          gathered_after,
          decomposition,
          1.0e-10,
          1.0e-10);
      summary.ppm_executed = true;
      summary.reconstruction_ghost_layers = kGhostLayers;
      summary.reconstruction_mode = "macro_zoning_mpi_real_angular_ppm";
      summary.macro_zoning_mode = "mpi_real_angular_macro_ppm";
      DEC3D_CHECK(summary.is_complete());
      DEC3D_CHECK(summary.success);
      DEC3D_CHECK(summary.parity_within_tolerance);
      DEC3D_CHECK(summary.seam_jumps_bounded);

      const std::filesystem::path output_dir =
          "F:\\dec3d\\analysis\\output\\case2_radial_wave_mpi_real_macro_ppm";
      std::filesystem::remove_all(output_dir);
      DEC3D_CHECK(dec3d::hydro::WriteSphericalRadialWaveMpiParityOutputs(
          single_rank_after,
          gathered_after,
          geometry,
          output_dir,
          summary,
          &single_rank_budget,
          &multi_rank_budget,
          "case2_radial_wave_mpi_real_macro_ppm"));

      const auto summary_text = dec3d::testsupport::ReadTextFile(output_dir / "summary.txt");
      DEC3D_CHECK(summary_text.find("case=case2_radial_wave_mpi_real_macro_ppm") != std::string::npos);
      DEC3D_CHECK(summary_text.find("ppm_executed=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("parity_within_tolerance=true") != std::string::npos);
      DEC3D_CHECK(summary_text.find("seam_jumps_bounded=true") != std::string::npos);

      const auto manifest_text = dec3d::testsupport::ReadTextFile(output_dir / "case_manifest.txt");
      DEC3D_CHECK(manifest_text.find("reconstruction_mode=macro_zoning_mpi_real_angular_ppm") != std::string::npos);
      DEC3D_CHECK(manifest_text.find("macro_zoning_mode=mpi_real_angular_macro_ppm") != std::string::npos);
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
