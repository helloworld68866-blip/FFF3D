#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

[[nodiscard]] bool HasDiagnostic(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] std::string PayloadForCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& code) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return entry.message;
    }
  }
  return {};
}

[[nodiscard]] std::string ExtractField(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& field) {
  const auto report = PayloadForCode(
      diagnostics,
      "p1.hydro.macro_ale_hllc.executed");
  const std::string prefix = field + "=";
  const auto begin = report.find(prefix);
  if (begin == std::string::npos) {
    return {};
  }
  const auto value_begin = begin + prefix.size();
  const auto value_end = report.find(';', value_begin);
  return report.substr(
      value_begin,
      value_end == std::string::npos ? std::string::npos : value_end - value_begin);
}

void SeedState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double density =
            1.0 + 0.025 * static_cast<double>(radial) +
            0.003 * std::sin(0.5 * static_cast<double>(theta + 1u)) +
            0.002 * std::cos(0.25 * static_cast<double>(phi + 1u));
        const double v_r = 0.01 + 0.001 * static_cast<double>(radial);
        state.rho(radial, theta, phi) = density;
        state.mom_r(radial, theta, phi) = density * v_r;
        state.mom_theta(radial, theta, phi) =
            1.0e-8 * std::sin(static_cast<double>(theta + 1u));
        state.mom_phi(radial, theta, phi) =
            -1.0e-8 * std::cos(static_cast<double>(phi + 1u));
        state.e_fluid_total(radial, theta, phi) = 2.3 + 0.06 * density;
        state.e_electron(radial, theta, phi) = 0.43 + 0.01 * density;
      }
    }
  }
}

void SeedUniformState(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(radial, theta, phi) = 1.0;
        state.mom_r(radial, theta, phi) = 0.0;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = 2.5;
        state.e_electron(radial, theta, phi) = 0.4;
      }
    }
  }
}

void SeedRadiationGroups(dec3d::state::CanonicalState& state) {
  for (std::size_t group = 0; group < state.radiation_groups.size(); ++group) {
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          state.radiation_groups[group](radial, theta, phi) =
              1.0 + 0.2 * static_cast<double>(group + 1u) +
              0.05 * static_cast<double>(radial);
        }
      }
    }
  }
}

void SeedHomologousCompressionState(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  constexpr double kDensity = 1.0;
  constexpr double kSpecificInternalEnergy = 1.0;
  constexpr double kGammaMinusOne = 2.0 / 3.0;
  const double pressure = kGammaMinusOne * kDensity * kSpecificInternalEnergy;
  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    const double radial_center =
        0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
    const double v_r = -radial_center;
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(radial, theta, phi) = kDensity;
        state.mom_r(radial, theta, phi) = kDensity * v_r;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) =
            kDensity * kSpecificInternalEnergy + 0.5 * kDensity * v_r * v_r;
        state.e_electron(radial, theta, phi) = pressure;
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::HydroConservativeState HomologousCompressionStateAtRadius(
    double radius) noexcept {
  constexpr double kDensity = 1.0;
  constexpr double kSpecificInternalEnergy = 1.0;
  constexpr double kGammaMinusOne = 2.0 / 3.0;
  const double pressure = kGammaMinusOne * kDensity * kSpecificInternalEnergy;
  const double v_r = -radius;
  return dec3d::hydro::HydroConservativeState{
      kDensity,
      kDensity * v_r,
      0.0,
      0.0,
      kDensity * kSpecificInternalEnergy + 0.5 * kDensity * v_r * v_r,
      pressure};
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildHomologousOuterGhostOverride(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta_cells,
    std::size_t phi_cells,
    std::size_t ghost_layers) {
  dec3d::hydro::RadialGhostOverride ghost_override;
  ghost_override.has_outer_neighbor = true;
  ghost_override.outer_ghost_reuses_boundary_partition = true;
  ghost_override.ghost_layers = ghost_layers;
  ghost_override.outer_ghost_states =
      dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
          ghost_layers,
          theta_cells,
          phi_cells);
  const double outer_face = geometry.radial_faces.back();
  const double outer_width =
      geometry.radial_faces.back() -
      geometry.radial_faces[geometry.radial_faces.size() - 2u];
  for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
    const double ghost_radius =
        outer_face + (static_cast<double>(ghost) + 0.5) * outer_width;
    const auto ghost_state = HomologousCompressionStateAtRadius(ghost_radius);
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        ghost_override.outer_ghost_states(ghost, theta, phi) = ghost_state;
      }
    }
  }
  ghost_override.report_line =
      "test homologous compression exact outer radial ghosts";
  DEC3D_CHECK(ghost_override.is_complete(theta_cells, phi_cells, ghost_layers));
  return ghost_override;
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BuildStaticMacroOptions() {
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

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildZeroMotionGlobalProposal(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) {
  dec3d::core::MeshUpdateProposal proposal{};
  proposal.requested = true;
  proposal.radial_ale = true;
  proposal.radial_face_indexing_is_global = true;
  proposal.global_radial_face_count = geometry.radial_faces.size();
  proposal.dt_s = dt_s;
  proposal.radial_face_velocities.assign(geometry.radial_faces.size(), 0.0);
  proposal.proposed_radial_faces = geometry.radial_faces;
  proposal.implementation_id = "p1.mesh.ale.radial_proposal";
  proposal.summary =
      "test zero-motion global radial ALE proposal for A4 direct HLLC";
  DEC3D_CHECK(proposal.is_complete_global(geometry.radial_faces.size()));
  return proposal;
}

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildSmoothNonzeroGlobalProposal(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) {
  dec3d::core::MeshUpdateProposal proposal{};
  proposal.requested = true;
  proposal.radial_ale = true;
  proposal.radial_face_indexing_is_global = true;
  proposal.global_radial_face_count = geometry.radial_faces.size();
  proposal.dt_s = dt_s;
  proposal.radial_face_velocities.resize(geometry.radial_faces.size(), 0.0);
  proposal.proposed_radial_faces = geometry.radial_faces;
  for (std::size_t face = 0; face < geometry.radial_faces.size(); ++face) {
    proposal.radial_face_velocities[face] =
        1.0e-2 * static_cast<double>(face) /
        static_cast<double>(geometry.radial_faces.size() - 1u);
    proposal.proposed_radial_faces[face] =
        geometry.radial_faces[face] + dt_s * proposal.radial_face_velocities[face];
  }
  proposal.implementation_id = "p1.mesh.ale.radial_proposal";
  proposal.summary =
      "test smooth nonzero global radial ALE proposal for A4 direct HLLC";
  DEC3D_CHECK(proposal.is_complete_global(geometry.radial_faces.size()));
  return proposal;
}

[[nodiscard]] dec3d::core::MeshUpdateProposal BuildHomologousCompressionProposal(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) {
  std::vector<double> shell_velocities;
  shell_velocities.reserve(geometry.radial_faces.size() - 1u);
  for (std::size_t radial = 0; radial + 1u < geometry.radial_faces.size(); ++radial) {
    const double radial_center =
        0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
    shell_velocities.push_back(-radial_center);
  }
  const auto proposal = dec3d::mesh::BuildRadialAleProposalFromGlobalShellVelocity(
      shell_velocities,
      geometry.radial_faces,
      dt_s,
      1.0);
  DEC3D_CHECK(proposal.is_complete_global(geometry.radial_faces.size()));
  return proposal;
}

[[nodiscard]] double MaxHomologousLabelVelocityError(
    const dec3d::state::HydroStateView& view,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry) noexcept {
  double max_error = 0.0;
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    const double old_radial_center =
        0.5 * (old_geometry.radial_faces[radial] + old_geometry.radial_faces[radial + 1u]);
    const double expected_velocity = -old_radial_center;
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const double rho = (*view.rho)(radial, theta, phi);
        const double velocity = (*view.mom_r)(radial, theta, phi) / rho;
        max_error = std::max(max_error, std::abs(velocity - expected_velocity));
      }
    }
  }
  return max_error;
}

[[nodiscard]] double MaxHomologousDensityRelativeError(
    const dec3d::state::HydroStateView& view,
    double dt_s) noexcept {
  const double exact_density = 1.0 / ((1.0 - dt_s) * (1.0 - dt_s) * (1.0 - dt_s));
  double max_error = 0.0;
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        max_error = std::max(
            max_error,
            std::abs((*view.rho)(radial, theta, phi) - exact_density) / exact_density);
      }
    }
  }
  return max_error;
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
            std::abs((*lhs.rho)(radial, theta, phi) -
                     (*rhs.rho)(radial, theta, phi)));
        max_difference = std::max(
            max_difference,
            std::abs((*lhs.mom_r)(radial, theta, phi) -
                     (*rhs.mom_r)(radial, theta, phi)));
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
            std::abs(lhs.chi_e(radial, theta, phi) -
                     rhs.chi_e(radial, theta, phi)));
      }
    }
  }
  return max_difference;
}

[[nodiscard]] double MaxDensityDeviationFromOne(
    const dec3d::state::HydroStateView& view) noexcept {
  double max_deviation = 0.0;
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        max_deviation = std::max(
            max_deviation,
            std::abs((*view.rho)(radial, theta, phi) - 1.0));
      }
    }
  }
  return max_deviation;
}

[[nodiscard]] double MaxDensityAngularSpread(
    const dec3d::state::HydroStateView& view) noexcept {
  double max_spread = 0.0;
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    double min_density = std::numeric_limits<double>::infinity();
    double max_density = -std::numeric_limits<double>::infinity();
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const double density = (*view.rho)(radial, theta, phi);
        min_density = std::min(min_density, density);
        max_density = std::max(max_density, density);
      }
    }
    max_spread = std::max(max_spread, max_density - min_density);
  }
  return max_spread;
}

[[nodiscard]] double MaxTangentialMomentum(
    const dec3d::state::HydroStateView& view) noexcept {
  double max_momentum = 0.0;
  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        max_momentum = std::max(
            max_momentum,
            std::abs((*view.mom_theta)(radial, theta, phi)));
        max_momentum = std::max(
            max_momentum,
            std::abs((*view.mom_phi)(radial, theta, phi)));
      }
    }
  }
  return max_momentum;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    constexpr std::size_t kRadialCells = 6u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr double kDt = 1.0e-5;

    auto static_state = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto direct_state = CanonicalState::Create(static_state.layout);
    SeedState(static_state);
    SeedState(direct_state);

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.05, 0.55});
    DEC3D_CHECK(geometry.is_valid());

    auto static_view = BuildHydroWorkView(static_state);
    auto direct_view = BuildHydroWorkView(direct_state);
    DEC3D_CHECK(static_view.is_complete());
    DEC3D_CHECK(direct_view.is_complete());

    const auto static_options = BuildStaticMacroOptions();
    auto direct_options = static_options;
    direct_options.apply_radial_ale_flux_correction = true;
    direct_options.use_macro_ale_direct_moving_face_hllc = true;
    const auto zero_proposal = BuildZeroMotionGlobalProposal(geometry, kDt);

    const auto static_result = AdvanceStaticGridHydro(
        static_view,
        geometry,
        kDt,
        static_options);
    const auto direct_result = AdvanceStaticGridHydro(
        direct_view,
        geometry,
        kDt,
        direct_options,
        &zero_proposal);

    if (!static_result.success) {
      std::cerr << static_result.failure_reason << '\n';
    }
    if (!direct_result.success) {
      std::cerr << direct_result.failure_reason << '\n';
    }
    DEC3D_CHECK(static_result.success);
    DEC3D_CHECK(direct_result.success);
    DEC3D_CHECK(direct_result.macro_ale_hllc_executed);
    DEC3D_CHECK(HasDiagnostic(
        direct_result.diagnostics,
        "p1.hydro.macro_ale_hllc.executed"));
    DEC3D_CHECK(HasDiagnostic(
        direct_result.diagnostics,
        "p1.hydro.macro_ale_hllc.coarse_source_budget"));
    DEC3D_CHECK(!HasDiagnostic(
        direct_result.diagnostics,
        "p1.hydro.macro_ale.compatibility.executed"));
    DEC3D_CHECK_EQ(ExtractField(direct_result.diagnostics, "remap_order"), std::string("none"));
    DEC3D_CHECK(
        direct_result.macro_ale_hllc_diagnostics.report_line.find(
            "remap_order=none") != std::string::npos);
    DEC3D_CHECK(
        direct_result.macro_ale_hllc_diagnostics.report_line.find(
            "budget_components=rho,mom_r,mom_theta,mom_phi,e_fluid_total,chi_e") !=
        std::string::npos);
    DEC3D_CHECK(
        direct_result.macro_ale_hllc_diagnostics.report_line.find(
            "budget_residual_source=computed_full_transport") != std::string::npos);
    DEC3D_CHECK_EQ(
        ExtractField(direct_result.diagnostics, "radial_ale_flux_mode"),
        std::string("moving_interface_hllc"));
    DEC3D_CHECK(direct_result.macro_ale_hllc_local_face_window.is_complete());
    DEC3D_CHECK(MaxHydroDifference(static_view, direct_view) < 2.0e-10);

    {
      auto radiation_state = CanonicalState::Create(
          CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 1u});
      SeedState(radiation_state);
      SeedRadiationGroups(radiation_state);
      auto radiation_view = BuildHydroWorkView(radiation_state);
      DEC3D_CHECK(radiation_view.is_complete());
      auto radiation_options = static_options;
      radiation_options.apply_radial_ale_flux_correction = true;
      radiation_options.use_macro_ale_direct_moving_face_hllc = true;
      radiation_options.enable_radiation_hydro_terms = true;
      const auto radiation_result = AdvanceStaticGridHydro(
          radiation_view,
          geometry,
          kDt,
          radiation_options,
          &zero_proposal);
      if (!radiation_result.success) {
        std::cerr << radiation_result.failure_reason << '\n';
      }
      DEC3D_CHECK(radiation_result.success);
      DEC3D_CHECK(radiation_result.macro_ale_hllc_executed);
      DEC3D_CHECK(
          radiation_result.macro_ale_hllc_diagnostics.report_line.find(
              "moving_mesh_face_velocity_source=hydro_ale_hllc_face_velocity") !=
          std::string::npos);
      DEC3D_CHECK(
          radiation_result.macro_ale_hllc_diagnostics.report_line.find(
              "radiation_scalar_bundle_transport=enabled") != std::string::npos);
      DEC3D_CHECK(
          radiation_result.macro_ale_hllc_diagnostics.report_line.find(
              "radiation_group_count=1") != std::string::npos);
      DEC3D_CHECK(
          radiation_result.macro_ale_hllc_diagnostics.report_line.find(
              "global_radiation_chi_residual=") != std::string::npos);
      DEC3D_CHECK(
          radiation_result.macro_ale_hllc_diagnostics.report_line.find(
              "max_seam_flux_radiation_chi_delta=") != std::string::npos);
    }

    auto uniform_state = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedUniformState(uniform_state);
    auto uniform_view = BuildHydroWorkView(uniform_state);
    DEC3D_CHECK(uniform_view.is_complete());

    auto radial_only_direct_options = BuildStaticMacroOptions();
    radial_only_direct_options.apply_geometric_source = false;
    radial_only_direct_options.apply_theta_sweep = false;
    radial_only_direct_options.apply_phi_sweep = false;
    radial_only_direct_options.apply_radial_ale_flux_correction = true;
    radial_only_direct_options.use_macro_ale_direct_moving_face_hllc = true;
    constexpr double kNonzeroDt = 1.0e-5;
    const auto nonzero_proposal = BuildSmoothNonzeroGlobalProposal(geometry, kNonzeroDt);
    const auto nonzero_result = AdvanceStaticGridHydro(
        uniform_view,
        geometry,
        kNonzeroDt,
        radial_only_direct_options,
        &nonzero_proposal);

    if (!nonzero_result.success) {
      std::cerr << nonzero_result.failure_reason << '\n';
    }
    DEC3D_CHECK(nonzero_result.success);
    DEC3D_CHECK(nonzero_result.macro_ale_hllc_diagnostics.is_complete());
    DEC3D_CHECK(std::abs(nonzero_result.macro_ale_hllc_diagnostics.global_mass_residual) < 1.0e-9);
    DEC3D_CHECK(std::abs(nonzero_result.macro_ale_hllc_diagnostics.global_e_fluid_total_residual) < 1.0e-8);
    DEC3D_CHECK(nonzero_result.macro_ale_hllc_diagnostics.max_abs_w_face > 0.0);
    DEC3D_CHECK_EQ(
        ExtractField(nonzero_result.diagnostics, "radial_ale_flux_mode"),
        std::string("moving_interface_hllc"));
    DEC3D_CHECK(MaxDensityDeviationFromOne(uniform_view) < 5.0e-11);
    DEC3D_CHECK(MaxDensityAngularSpread(uniform_view) == 0.0);

    auto full_uniform_state = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedUniformState(full_uniform_state);
    auto full_uniform_view = BuildHydroWorkView(full_uniform_state);
    DEC3D_CHECK(full_uniform_view.is_complete());
    auto full_direct_options = BuildStaticMacroOptions();
    full_direct_options.apply_radial_ale_flux_correction = true;
    full_direct_options.use_macro_ale_direct_moving_face_hllc = true;
    const auto full_nonzero_result = AdvanceStaticGridHydro(
        full_uniform_view,
        geometry,
        kNonzeroDt,
        full_direct_options,
        &nonzero_proposal);
    if (!full_nonzero_result.success) {
      std::cerr << full_nonzero_result.failure_reason << '\n';
    }
    DEC3D_CHECK(full_nonzero_result.success);
    DEC3D_CHECK(MaxDensityAngularSpread(full_uniform_view) == 0.0);
    DEC3D_CHECK(MaxTangentialMomentum(full_uniform_view) == 0.0);

    auto homologous_state = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    const auto homologous_geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 0.6});
    DEC3D_CHECK(homologous_geometry.is_valid());
    SeedHomologousCompressionState(homologous_state, homologous_geometry);
    auto homologous_view = BuildHydroWorkView(homologous_state);
    DEC3D_CHECK(homologous_view.is_complete());
    auto homologous_options = BuildStaticMacroOptions();
    homologous_options.apply_radial_ale_flux_correction = true;
    homologous_options.use_macro_ale_direct_moving_face_hllc = true;
    const auto homologous_proposal =
        BuildHomologousCompressionProposal(homologous_geometry, kNonzeroDt);
    const auto homologous_result = AdvanceStaticGridHydro(
        homologous_view,
        homologous_geometry,
        kNonzeroDt,
        BuildHomologousOuterGhostOverride(
            homologous_geometry,
            kThetaCells,
            kPhiCells,
            homologous_options.reconstruction_ghost_layers),
        homologous_options,
        &homologous_proposal);
    if (!homologous_result.success) {
      std::cerr << homologous_result.failure_reason << '\n';
    }
    DEC3D_CHECK(homologous_result.success);
    const double homologous_density_error =
        MaxHomologousDensityRelativeError(homologous_view, kNonzeroDt);
    const double homologous_velocity_error =
        MaxHomologousLabelVelocityError(homologous_view, homologous_geometry);
    if (homologous_density_error >= 1.0e-12 ||
        homologous_velocity_error >= 1.0e-10) {
      std::cerr << "homologous_density_error=" << homologous_density_error
                << "; homologous_velocity_error=" << homologous_velocity_error << '\n';
    }
    DEC3D_CHECK(homologous_density_error < 1.0e-12);
    DEC3D_CHECK(homologous_velocity_error < 1.0e-10);
    DEC3D_CHECK(MaxDensityAngularSpread(homologous_view) == 0.0);
    DEC3D_CHECK(MaxTangentialMomentum(homologous_view) == 0.0);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
