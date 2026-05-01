#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
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
#include <numbers>

namespace {

void SeedPoleSeamStressCase(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  constexpr double kBaseRho = 1.0;
  constexpr double kBasePressure = 1.0;
  constexpr double kBaseElectronPressure = 0.4;
  constexpr double kAmplitude = 3.0e-2;
  constexpr double kSigma = 0.28;
  constexpr double kPolarOffset = 0.25;

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      const double theta_center =
          0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
      const double north = std::exp(
          -std::pow(theta_center - kPolarOffset, 2.0) / (2.0 * kSigma * kSigma));
      const double south = std::exp(
          -std::pow(theta_center - (std::numbers::pi - kPolarOffset), 2.0) /
          (2.0 * kSigma * kSigma));

      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double phi_center =
            0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
        const double phi_mode = std::cos(phi_center) + 0.5 * std::sin(2.0 * phi_center);
        const double delta = kAmplitude * (north - south) * phi_mode;
        const double pressure = kBasePressure * (1.0 + 0.15 * delta);
        const double electron_pressure = kBaseElectronPressure * (1.0 + 0.10 * delta);
        const double v_theta = 2.0e-2 * delta;
        const double v_phi = -2.5e-2 * delta;
        const double kinetic = 0.5 * kBaseRho * (v_theta * v_theta + v_phi * v_phi);

        state.rho(radial, theta, phi) = kBaseRho;
        state.mom_r(radial, theta, phi) = 0.0;
        state.mom_theta(radial, theta, phi) = kBaseRho * v_theta;
        state.mom_phi(radial, theta, phi) = kBaseRho * v_phi;
        state.e_fluid_total(radial, theta, phi) =
            pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) + kinetic;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);
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

    constexpr std::size_t kRadialCells = 3u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr double kDt = 5.0e-4;
    constexpr std::size_t kShellIndex = 1u;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = dec3d::state::CanonicalState::Create(before.layout);
    SeedPoleSeamStressCase(before, geometry);
    SeedPoleSeamStressCase(single_rank_after, geometry);

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
        local_override);
    DEC3D_CHECK(local_result.is_complete());
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(local_result.diagnostics, "p1.hydro.ghost.direction.radial"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(local_result.diagnostics, "p1.hydro.ghost.direction.theta"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(local_result.diagnostics, "p1.hydro.ghost.direction.phi"));

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
          "mesh:p1.case4.pole-stress-mpi-real",
          "ownership:rank0",
          "diagnostics:p1.hydro.case4.mpi-real"};
      DEC3D_CHECK(context.is_complete());
      DEC3D_CHECK(hydro.bind(context, geometry, single_rank_after));
      const auto single_rank_stage_result = hydro.advance();
      DEC3D_CHECK(single_rank_stage_result.is_semantically_complete());

      dec3d::hydro::HydroBudgetResidualSummary single_rank_budget;
      DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
          single_rank_stage_result.diagnostics,
          &single_rank_budget));
      DEC3D_CHECK(single_rank_budget.is_complete());

      const auto stress = dec3d::hydro::EvaluatePoleAdjacentSeamStress(
          before,
          gathered_after,
          geometry,
          kShellIndex,
          1.0e-4,
          2.0e-1);
      DEC3D_CHECK(stress.is_complete());
      DEC3D_CHECK(stress.success);

      const auto mpi_parity = dec3d::hydro::EvaluateDirectionalMpiParity(
          single_rank_after,
          gathered_after,
          static_cast<std::size_t>(rank_count),
          1.0e-10);
      DEC3D_CHECK(mpi_parity.is_complete());
      DEC3D_CHECK(mpi_parity.success);

      DEC3D_CHECK(std::abs(single_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(single_rank_budget.e_fluid_total_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.mass_residual) < 1.0e-10);
      DEC3D_CHECK(std::abs(multi_rank_budget.e_fluid_total_residual) < 1.0e-10);

      const std::filesystem::path output_dir =
          "F:\\dec3d\\analysis\\output\\case4_pole_seam_stress_mpi_real";
      std::filesystem::remove_all(output_dir);
      DEC3D_CHECK(dec3d::hydro::WritePoleAdjacentSeamStressOutputs(
          before,
          gathered_after,
          geometry,
          kShellIndex,
          output_dir,
          stress,
          &multi_rank_budget,
          &mpi_parity,
          &single_rank_budget,
          "case4_pole_seam_stress_mpi_real"));

      DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "mpi_parity_diagnostics.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual_comparison.txt"));
      DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

      const auto summary_text = dec3d::testsupport::ReadTextFile(output_dir / "summary.txt");
      DEC3D_CHECK(summary_text.find("case=case4_pole_seam_stress_mpi_real") != std::string::npos);
      DEC3D_CHECK(summary_text.find("mpi_parity_within_tolerance=true") != std::string::npos);

      const auto manifest_text = dec3d::testsupport::ReadTextFile(output_dir / "case_manifest.txt");
      DEC3D_CHECK(manifest_text.find("mpi_parity_diagnostics=mpi_parity_diagnostics.txt") != std::string::npos);
      DEC3D_CHECK(manifest_text.find("budget_comparison=budget_residual_comparison.txt") != std::string::npos);

      const auto parity_text =
          dec3d::testsupport::ReadTextFile(output_dir / "mpi_parity_diagnostics.txt");
      DEC3D_CHECK(parity_text.find("# kind=mpi_parity_diagnostics") != std::string::npos);
      DEC3D_CHECK(parity_text.find("max_hydro_only_te_difference=") != std::string::npos);
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
