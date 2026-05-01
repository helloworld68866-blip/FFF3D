#include "hydro/driver/macro_ale_mpi_compatibility.hpp"
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
#include <string>
#include <vector>

namespace {

constexpr std::size_t kRadialCells = 12u;
constexpr std::size_t kThetaCells = 8u;
constexpr std::size_t kPhiCells = 8u;
constexpr std::size_t kGhostLayers = 3u;
constexpr double kDt = 1.0e-5;
constexpr double kTolerance = 2.0e-10;

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.015 * static_cast<double>(radial) +
            0.002 * std::sin(0.5 * static_cast<double>(theta + 1u)) +
            0.001 * std::cos(0.25 * static_cast<double>(phi + 1u));
        const double v_r = 0.006 + 0.0005 * static_cast<double>(radial);
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = density * v_r;
        state.mom_theta(radial, theta, phi) =
            1.0e-8 * std::sin(static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -1.0e-8 * std::cos(static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.2 + 0.04 * density;
        state.e_electron(radial, theta, phi) = 0.42 + 0.01 * density;
      }
    }
  }
}

void SeedStressState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double theta_mode = std::sin(static_cast<double>(theta) + 1.0);
        const double phi_mode = std::cos(static_cast<double>(phi) + 0.5);
        const double density =
            1.0 + 0.04 * static_cast<double>(radial) +
            0.003 * theta_mode * phi_mode;
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) =
            density * (-0.02 - 0.001 * static_cast<double>(radial));
        state.mom_theta(radial, theta, phi) =
            1.0e-8 * density * std::sin(static_cast<double>(theta) + 0.25);
        state.mom_phi(radial, theta, phi) =
            -1.0e-8 * density * std::cos(static_cast<double>(phi) + 0.5);
        state.e_fluid_total(radial, theta, phi) = 2.5 + 0.1 * density;
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

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BuildStressOptions() {
  auto options = BuildOptions();
  options.macro_zoning_coarse_factor = 0.5;
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

  DEC3D_CHECK(rho_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_r_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_theta_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(mom_phi_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(e_fluid_total_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));
  DEC3D_CHECK(e_electron_exchange.is_complete(theta_cells, phi_cells, kGhostLayers));

  return dec3d::testsupport::BuildRadialGhostOverride(
      rho_exchange,
      mom_r_exchange,
      mom_theta_exchange,
      mom_phi_exchange,
      e_fluid_total_exchange,
      e_electron_exchange);
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

[[nodiscard]] bool SameStorage(
    const dec3d::core::Array3D<double>& lhs,
    const std::vector<double>& rhs) noexcept {
  return lhs.storage() == rhs;
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

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.65});
    DEC3D_CHECK(geometry.is_valid());
    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    DEC3D_CHECK(decomposition.is_valid());
    const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];

    auto local_state =
        dec3d::testsupport::BuildLocalStateSlice(global_state, slice);
    auto static_state =
        dec3d::testsupport::BuildLocalStateSlice(global_state, slice);
    auto local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        slice,
        kThetaCells,
        kPhiCells);
    auto static_geometry = local_geometry;
    DEC3D_CHECK(local_geometry.is_valid());
    DEC3D_CHECK(static_geometry.is_valid());

    const auto local_override = BuildOverride(local_state);
    const auto proposal_result = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
        local_state.rho,
        local_state.mom_r,
        decomposition,
        geometry.radial_faces,
        kDt,
        1.0,
        MPI_COMM_WORLD);
    DEC3D_CHECK(proposal_result.is_complete());

    auto local_view = dec3d::state::BuildHydroWorkView(local_state);
    DEC3D_CHECK(local_view.is_complete());
    auto result = dec3d::hydro::AdvanceMacroAleMpiCompatibilityStep(
        local_view,
        local_geometry,
        decomposition,
        proposal_result.proposal,
        kDt,
        local_override,
        BuildOptions(),
        MPI_COMM_WORLD);
    if (!result.success && rank == 0) {
      std::cerr << result.failure_reason << '\n';
    }
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.diagnostics.global_stage_ok);
    DEC3D_CHECK(result.diagnostics.global_publish_ok);
    DEC3D_CHECK(result.diagnostics.canonical_writeback_published);
    DEC3D_CHECK(result.diagnostics.mesh_commit_published);
    DEC3D_CHECK(result.diagnostics.published_from_staged_geometry);
    DEC3D_CHECK(result.diagnostics.max_seam_face_velocity_delta <= 1.0e-14);
    DEC3D_CHECK(result.diagnostics.max_seam_proposed_face_delta <= 1.0e-14);
    DEC3D_CHECK(result.diagnostics.max_seam_preview_radius_delta <= 1.0e-14);
    DEC3D_CHECK(std::abs(result.diagnostics.local_mass_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.local_mom_r_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.local_mom_theta_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.local_mom_phi_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.local_e_fluid_total_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.local_chi_e_residual) < kTolerance);
    DEC3D_CHECK(std::abs(result.diagnostics.global_mass_residual) < kTolerance);
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale.local_window"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_mpi.global_transaction.preflight"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_mpi.global_conservation"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.diagnostics_payload,
        "p1.hydro.macro_ale_mpi.transaction.clean"));
    DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
        result.staged_hydro_result.diagnostics,
        "p1.hydro.macro_ale.compatibility.executed"));
    DEC3D_CHECK(
        result.diagnostics.report_line.find(
            "remap_order=first_order_proposal_mapped_overlap") != std::string::npos);

    auto static_view = dec3d::state::BuildHydroWorkView(static_state);
    DEC3D_CHECK(static_view.is_complete());
    auto static_options = BuildOptions();
    static_options.apply_radial_ale_flux_correction = true;
    dec3d::core::MeshUpdateProposal zero_proposal = proposal_result.proposal;
    zero_proposal.radial_face_velocities.assign(geometry.radial_faces.size(), 0.0);
    zero_proposal.proposed_radial_faces = geometry.radial_faces;
    auto zero_result = dec3d::hydro::AdvanceMacroAleMpiCompatibilityStep(
        static_view,
        static_geometry,
        decomposition,
        zero_proposal,
        kDt,
        BuildOverride(static_state),
        static_options,
        MPI_COMM_WORLD);
    DEC3D_CHECK(zero_result.success);

    auto no_ale_view = dec3d::state::BuildHydroWorkView(static_state);
    auto no_ale_options = BuildOptions();
    no_ale_options.apply_radial_ale_flux_correction = false;
    const auto no_ale_result = dec3d::hydro::AdvanceStaticGridHydro(
        no_ale_view,
        dec3d::testsupport::BuildLocalGeometrySlice(geometry, slice, kThetaCells, kPhiCells),
        kDt,
        BuildOverride(static_state),
        no_ale_options);
    DEC3D_CHECK(no_ale_result.success);
    DEC3D_CHECK(MaxHydroDifference(static_view, no_ale_view) < 1.0e-12);

    {
      constexpr std::size_t kStressRadialCells = 15u;
      constexpr std::size_t kStressThetaCells = 16u;
      constexpr std::size_t kStressPhiCells = 16u;
      auto stress_global_state = dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{
              kStressRadialCells,
              kStressThetaCells,
              kStressPhiCells,
              0u});
      SeedStressState(stress_global_state);
      const auto stress_geometry = dec3d::mesh::BuildSphericalGeometry(
          dec3d::mesh::SphericalMeshDescriptor{
              kStressRadialCells,
              kStressThetaCells,
              kStressPhiCells,
              0.04,
              0.70});
      DEC3D_CHECK(stress_geometry.is_valid());
      const auto stress_decomposition = dec3d::mesh::BuildRadialOwnership(
          kStressRadialCells,
          static_cast<std::size_t>(rank_count));
      DEC3D_CHECK(stress_decomposition.is_valid());
      const auto& stress_slice =
          stress_decomposition.slices[static_cast<std::size_t>(rank)];
      auto stress_local_state =
          dec3d::testsupport::BuildLocalStateSlice(stress_global_state, stress_slice);
      auto stress_local_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
          stress_geometry,
          stress_slice,
          kStressThetaCells,
          kStressPhiCells);
      auto stress_view = dec3d::state::BuildHydroWorkView(stress_local_state);
      DEC3D_CHECK(stress_view.is_complete());

      const auto stress_override = BuildOverride(stress_local_state);
      const auto stress_proposal = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
          stress_local_state.rho,
          stress_local_state.mom_r,
          stress_decomposition,
          stress_geometry.radial_faces,
          kDt,
          1.0,
          MPI_COMM_WORLD);
      DEC3D_CHECK(stress_proposal.is_complete());
      const auto stress_result = dec3d::hydro::AdvanceMacroAleMpiCompatibilityStep(
          stress_view,
          stress_local_geometry,
          stress_decomposition,
          stress_proposal.proposal,
          kDt,
          stress_override,
          BuildStressOptions(),
          MPI_COMM_WORLD);
      if (!stress_result.success && rank == 0) {
        std::cerr << stress_result.failure_reason << '\n';
      }
      DEC3D_CHECK(stress_result.success);
      DEC3D_CHECK(stress_result.diagnostics.max_seam_face_velocity_delta <= 1.0e-14);
      DEC3D_CHECK(stress_result.diagnostics.max_seam_proposed_face_delta <= 1.0e-14);
      DEC3D_CHECK(stress_result.diagnostics.max_seam_preview_radius_delta <= 1.0e-14);
      DEC3D_CHECK(std::abs(stress_result.diagnostics.global_mass_residual) <= 1.0e-10);
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          stress_result.diagnostics_payload,
          "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          stress_result.diagnostics_payload,
          "p1.hydro.macro_ale_mpi.seam_preview_geometry_equality"));
      DEC3D_CHECK(dec3d::testsupport::HasDiagnosticCode(
          stress_result.diagnostics_payload,
          "p1.hydro.macro_ale_mpi.transaction.clean"));
    }

    auto bad_state = dec3d::testsupport::BuildLocalStateSlice(global_state, slice);
    auto bad_view = dec3d::state::BuildHydroWorkView(bad_state);
    auto bad_geometry = dec3d::testsupport::BuildLocalGeometrySlice(
        geometry,
        slice,
        kThetaCells,
        kPhiCells);
    const auto before_rho = bad_view.rho->storage();
    const auto before_mom_r = bad_view.mom_r->storage();
    const auto before_faces = bad_geometry.radial_faces;
    auto bad_proposal = proposal_result.proposal;
    if (rank == 1) {
      bad_proposal.proposed_radial_faces.pop_back();
    }
    const auto bad_result = dec3d::hydro::AdvanceMacroAleMpiCompatibilityStep(
        bad_view,
        bad_geometry,
        decomposition,
        bad_proposal,
        kDt,
        BuildOverride(bad_state),
        BuildOptions(),
        MPI_COMM_WORLD);
    DEC3D_CHECK(!bad_result.success);
    DEC3D_CHECK(!bad_result.diagnostics.canonical_writeback_published);
    DEC3D_CHECK(!bad_result.diagnostics.mesh_commit_published);
    DEC3D_CHECK(SameStorage(*bad_view.rho, before_rho));
    DEC3D_CHECK(SameStorage(*bad_view.mom_r, before_mom_r));
    DEC3D_CHECK(bad_geometry.radial_faces == before_faces);

    int local_ok = 1;
    int global_ok = 0;
    MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Finalize();
    return global_ok == 1 ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
  }
}
