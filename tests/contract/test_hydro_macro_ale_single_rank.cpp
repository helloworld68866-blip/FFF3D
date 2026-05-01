#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

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

[[nodiscard]] std::string PayloadForCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return entry.message;
    }
  }
  return {};
}

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.03 * static_cast<double>(radial) +
            0.02 * std::sin(0.5 * static_cast<double>(theta + 1u)) +
            0.01 * std::cos(0.25 * static_cast<double>(phi + 1u));
        const double v_r = 0.012 + 0.002 * static_cast<double>(radial);
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = density * v_r;
        state.mom_theta(radial, theta, phi) =
            0.004 * std::sin(0.4 * static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -0.003 * std::cos(0.3 * static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.4 + 0.08 * density;
        state.e_electron(radial, theta, phi) = 0.45 + 0.02 * density;
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BuildMacroPpmOptions() {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = true;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  options.use_macro_zoning = true;
  options.macro_zoning_coarse_factor = 0.5;
  options.radial_ale_global_face_begin_index = 0u;
  return options;
}

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildZeroMotionProposal(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) {
  dec3d::core::MeshUpdateProposal proposal;
  proposal.requested = true;
  proposal.radial_ale = true;
  proposal.radial_face_indexing_is_global = true;
  proposal.global_radial_face_count = geometry.radial_faces.size();
  proposal.dt_s = dt_s;
  proposal.radial_face_velocities.assign(geometry.radial_faces.size(), 0.0);
  proposal.proposed_radial_faces = geometry.radial_faces;
  proposal.implementation_id = "p1.mesh.ale.radial_proposal";
  proposal.summary =
      "radial_face_indexing=global; test_zero_motion=true; max_face_speed=0";
  return proposal;
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

[[nodiscard]] double MaxCanonicalDifference(
    const dec3d::state::CanonicalState& lhs,
    const dec3d::state::CanonicalState& rhs) noexcept {
  double max_difference = 0.0;
  for (std::size_t radial = 0; radial < lhs.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < lhs.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < lhs.rho.extent_phi(); ++phi) {
        max_difference = std::max(
            max_difference,
            std::abs(lhs.rho(radial, theta, phi) - rhs.rho(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.mom_r(radial, theta, phi) - rhs.mom_r(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.mom_theta(radial, theta, phi) -
                     rhs.mom_theta(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.mom_phi(radial, theta, phi) -
                     rhs.mom_phi(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.e_fluid_total(radial, theta, phi) -
                     rhs.e_fluid_total(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs(lhs.e_electron(radial, theta, phi) -
                     rhs.e_electron(radial, theta, phi)));
      }
    }
  }
  return max_difference;
}

void AssertHydroViewRecoverable(const dec3d::state::HydroStateView& view) {
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const dec3d::hydro::HydroConservativeState conservative{
            (*view.rho)(radial, theta, phi),
            (*view.mom_r)(radial, theta, phi),
            (*view.mom_theta)(radial, theta, phi),
            (*view.mom_phi)(radial, theta, phi),
            (*view.e_fluid_total)(radial, theta, phi),
            view.chi_e(radial, theta, phi)};
        DEC3D_CHECK(conservative.is_finite());
        DEC3D_CHECK(std::isfinite(conservative.chi_e));
        DEC3D_CHECK(conservative.chi_e >= 0.0);
        DEC3D_CHECK(dec3d::hydro::RecoverPrimitiveState(conservative).is_physical());
      }
    }
  }
}

void TestZeroMotionConsistency() {
  using dec3d::hydro::AdvanceStaticGridHydro;
  using dec3d::mesh::BuildSphericalGeometry;
  using dec3d::mesh::SphericalMeshDescriptor;
  using dec3d::state::BuildHydroWorkView;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;

  constexpr std::size_t kRadialCells = 4u;
  constexpr std::size_t kThetaCells = 8u;
  constexpr std::size_t kPhiCells = 8u;
  constexpr double kDt = 2.0e-5;

  auto macro_only_state =
      CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
  auto macro_ale_state = CanonicalState::Create(macro_only_state.layout);
  SeedState(macro_only_state);
  SeedState(macro_ale_state);

  const auto geometry = BuildSphericalGeometry(
      SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
  DEC3D_CHECK(geometry.is_valid());

  auto macro_only_view = BuildHydroWorkView(macro_only_state);
  auto macro_ale_view = BuildHydroWorkView(macro_ale_state);
  DEC3D_CHECK(macro_only_view.is_complete());
  DEC3D_CHECK(macro_ale_view.is_complete());

  auto macro_only_options = BuildMacroPpmOptions();
  auto macro_ale_options = BuildMacroPpmOptions();
  macro_ale_options.apply_radial_ale_flux_correction = true;
  const auto zero_proposal = BuildZeroMotionProposal(geometry, kDt);

  const auto macro_only_result = AdvanceStaticGridHydro(
      macro_only_view,
      geometry,
      kDt,
      macro_only_options);
  const auto macro_ale_result = AdvanceStaticGridHydro(
      macro_ale_view,
      geometry,
      kDt,
      macro_ale_options,
      &zero_proposal);

  DEC3D_CHECK(macro_only_result.success);
  DEC3D_CHECK(macro_ale_result.success);
  DEC3D_CHECK(MaxHydroDifference(macro_only_view, macro_ale_view) < 1.0e-12);
  DEC3D_CHECK(HasDiagnosticCode(
      macro_ale_result.diagnostics,
      "p1.hydro.macro_ale.compatibility.executed"));
  DEC3D_CHECK(HasDiagnosticCode(
      macro_ale_result.diagnostics,
      "p1.hydro.macro_ale.remap.conservation"));
}

void TestPostHydroPathExecuted() {
  using dec3d::hydro::AdvanceStaticGridHydro;
  using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
  using dec3d::mesh::BuildSphericalGeometry;
  using dec3d::mesh::SphericalMeshDescriptor;
  using dec3d::state::BuildHydroAuthorizedWriteMask;
  using dec3d::state::BuildHydroWorkView;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;
  using dec3d::state::CommitHydroWriteback;
  using dec3d::state::ElectronEnergyDensityFromPressure;
  using dec3d::state::ElectronPressureFromChiE;

  constexpr std::size_t kRadialCells = 4u;
  constexpr std::size_t kThetaCells = 8u;
  constexpr std::size_t kPhiCells = 8u;
  constexpr double kDt = 2.0e-5;

  auto state =
      CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
  SeedState(state);
  const auto geometry = BuildSphericalGeometry(
      SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
  DEC3D_CHECK(geometry.is_valid());

  auto view = BuildHydroWorkView(state);
  DEC3D_CHECK(view.is_complete());
  auto options = BuildMacroPpmOptions();
  options.apply_radial_ale_flux_correction = true;
  const auto proposal =
      BuildRadialAleMeshUpdateProposal(state.rho, state.mom_r, geometry, kDt);
  DEC3D_CHECK(proposal.is_complete_global(geometry.radial_faces.size()));

  const auto result = AdvanceStaticGridHydro(
      view,
      geometry,
      kDt,
      options,
      &proposal);

  DEC3D_CHECK(result.success);
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.detected"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.restrict.executed"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.coarse_update.executed"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_zoning.prolong.executed"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.compatibility.executed"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.geometry_preview"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.remap.executed"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.remap.conservation"));
  DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.macro_ale.transaction.clean"));
  DEC3D_CHECK(PayloadForCode(result.diagnostics, "p1.hydro.macro_ale.remap.conservation")
                  .find("remap_order=first_order_proposal_mapped_overlap") != std::string::npos);
  DEC3D_CHECK(PayloadForCode(result.diagnostics, "p1.hydro.macro_ale.geometry_preview")
                  .find("radial_face_indexing=global") != std::string::npos);
  DEC3D_CHECK(PayloadForCode(result.diagnostics, "p1.hydro.macro_ale.geometry_preview")
                  .find("global_face_begin=0") != std::string::npos);

  AssertHydroViewRecoverable(view);
  std::vector<double> expected_electron_energy(state.rho.size(), 0.0);
  bool saw_non_stale_electron_recovery = false;
  std::size_t linear = 0u;
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        expected_electron_energy[linear] =
            ElectronEnergyDensityFromPressure(ElectronPressureFromChiE(
                view.chi_e(radial, theta, phi)));
        if (std::abs(expected_electron_energy[linear] -
                     (*view.e_electron)(radial, theta, phi)) > 1.0e-13) {
          saw_non_stale_electron_recovery = true;
        }
        ++linear;
      }
    }
  }
  DEC3D_CHECK(saw_non_stale_electron_recovery);

  const auto writeback = CommitHydroWriteback(
      state,
      view,
      BuildHydroAuthorizedWriteMask());
  DEC3D_CHECK(writeback.success);
  DEC3D_CHECK(writeback.electron_channel_updated);
  DEC3D_CHECK(HasDiagnosticCode(
      writeback.diagnostics,
      "p1.hydro.writeback.electron_channel"));
  linear = 0u;
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        DEC3D_CHECK(std::abs(
            state.e_electron(radial, theta, phi) -
            expected_electron_energy[linear]) < 1.0e-12);
        ++linear;
      }
    }
  }
}

void TestInvalidProposalDoesNotMutateCanonicalState() {
  using dec3d::hydro::AdvanceStaticGridHydro;
  using dec3d::mesh::BuildSphericalGeometry;
  using dec3d::mesh::SphericalMeshDescriptor;
  using dec3d::state::BuildHydroWorkView;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;

  constexpr std::size_t kRadialCells = 4u;
  constexpr std::size_t kThetaCells = 8u;
  constexpr std::size_t kPhiCells = 8u;
  constexpr double kDt = 2.0e-5;

  auto state =
      CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
  SeedState(state);
  const auto before = state;
  const auto geometry = BuildSphericalGeometry(
      SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.45});
  DEC3D_CHECK(geometry.is_valid());

  auto view = BuildHydroWorkView(state);
  DEC3D_CHECK(view.is_complete());
  auto options = BuildMacroPpmOptions();
  options.apply_radial_ale_flux_correction = true;
  auto invalid_proposal = BuildZeroMotionProposal(geometry, kDt);
  invalid_proposal.radial_face_indexing_is_global = false;

  const auto result = AdvanceStaticGridHydro(
      view,
      geometry,
      kDt,
      options,
      &invalid_proposal);
  DEC3D_CHECK(!result.success);
  DEC3D_CHECK(result.failure_reason.find("global radial face indexing") != std::string::npos);
  DEC3D_CHECK(MaxCanonicalDifference(state, before) == 0.0);
}

}  // namespace

int main() {
  try {
    TestZeroMotionConsistency();
    TestPostHydroPathExecuted();
    TestInvalidProposalDoesNotMutateCanonicalState();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
