#include "hydro/driver/hydro_numerical_checks.hpp"

#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace dec3d::hydro {

namespace {

[[nodiscard]] bool LayoutsMatch(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after) noexcept {
  return before.layout.radial_cells == after.layout.radial_cells &&
         before.layout.theta_cells == after.layout.theta_cells &&
         before.layout.phi_cells == after.layout.phi_cells;
}

[[nodiscard]] bool HasAllocatedHydroState(const dec3d::state::CanonicalState& state) noexcept {
  using dec3d::core::AuthoritativeField;
  return state.HasAuthoritativeStorage(AuthoritativeField::rho) &&
         state.HasAuthoritativeStorage(AuthoritativeField::mom_r) &&
         state.HasAuthoritativeStorage(AuthoritativeField::mom_theta) &&
         state.HasAuthoritativeStorage(AuthoritativeField::mom_phi) &&
         state.HasAuthoritativeStorage(AuthoritativeField::e_fluid_total);
}

[[nodiscard]] bool HasAllocatedElectronState(const dec3d::state::CanonicalState& state) noexcept {
  using dec3d::core::AuthoritativeField;
  return state.HasAuthoritativeStorage(AuthoritativeField::e_electron);
}

[[nodiscard]] bool LayoutAndGeometryMatch(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  return geometry.is_valid() &&
         state.layout.radial_cells + 1u == geometry.radial_faces.size() &&
         state.layout.theta_cells + 1u == geometry.theta_faces.size() &&
         state.layout.phi_cells + 1u == geometry.phi_faces.size();
}

[[nodiscard]] double ThetaCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] double PhiCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t phi) noexcept {
  return 0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
}

[[nodiscard]] double HydroOnlyElectronTemperature(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double rho = state.rho(radial, theta, phi);
  if (!(rho > 0.0) || !std::isfinite(rho)) {
    return std::nan("");
  }

  const double electron_pressure =
      dec3d::state::ElectronPressureFromElectronEnergyDensity(
          state.e_electron(radial, theta, phi));
  return electron_pressure / rho;
}

[[nodiscard]] double VelocityMagnitude(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double rho = state.rho(radial, theta, phi);
  if (!(rho > 0.0) || !std::isfinite(rho)) {
    return std::nan("");
  }

  const double v_r = state.mom_r(radial, theta, phi) / rho;
  const double v_theta = state.mom_theta(radial, theta, phi) / rho;
  const double v_phi = state.mom_phi(radial, theta, phi) / rho;
  return std::sqrt(v_r * v_r + v_theta * v_theta + v_phi * v_phi);
}

[[nodiscard]] dec3d::hydro::HydroConservativeState LoadConservativeState(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return {
      state.rho(radial, theta, phi),
      state.mom_r(radial, theta, phi),
      state.mom_theta(radial, theta, phi),
      state.mom_phi(radial, theta, phi),
      state.e_fluid_total(radial, theta, phi),
      dec3d::state::ChiEFromElectronPressure(
          dec3d::state::ElectronPressureFromElectronEnergyDensity(
              state.e_electron(radial, theta, phi)))};
}

[[nodiscard]] bool StateIsFinite(const dec3d::state::CanonicalState& state) noexcept {
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        if (!std::isfinite(state.rho(radial, theta, phi)) ||
            !std::isfinite(state.mom_r(radial, theta, phi)) ||
            !std::isfinite(state.mom_theta(radial, theta, phi)) ||
            !std::isfinite(state.mom_phi(radial, theta, phi)) ||
            !std::isfinite(state.e_fluid_total(radial, theta, phi)) ||
            !std::isfinite(state.e_electron(radial, theta, phi))) {
          return false;
        }
      }
    }
  }

  return true;
}

struct ModeProjection {
  double cosine_coefficient{0.0};
  double sine_coefficient{0.0};
  double amplitude{0.0};
  double phase{0.0};
};

[[nodiscard]] ModeProjection ComputeVrLowModeProjection(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index) noexcept {
  ModeProjection projection;
  const double inverse_count =
      1.0 / static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);

  for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
    const double theta_center = ThetaCenter(geometry, theta);
    const double sin_theta = std::sin(theta_center);
    for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
      const double phi_center = PhiCenter(geometry, phi);
      const double rho = state.rho(shell_index, theta, phi);
      if (!(rho > 0.0) || !std::isfinite(rho)) {
        projection.amplitude = std::nan("");
        projection.phase = std::nan("");
        return projection;
      }

      const double v_r = state.mom_r(shell_index, theta, phi) / rho;
      const double cosine_basis = sin_theta * std::cos(phi_center);
      const double sine_basis = sin_theta * std::sin(phi_center);
      projection.cosine_coefficient += v_r * cosine_basis * inverse_count;
      projection.sine_coefficient += v_r * sine_basis * inverse_count;
    }
  }

  projection.amplitude =
      std::sqrt(projection.cosine_coefficient * projection.cosine_coefficient +
                projection.sine_coefficient * projection.sine_coefficient);
  projection.phase =
      std::atan2(projection.sine_coefficient, projection.cosine_coefficient);
  return projection;
}

[[nodiscard]] bool EnsureOutputDirectory(
    const std::filesystem::path& output_directory) noexcept {
  std::error_code error;
  std::filesystem::create_directories(output_directory, error);
  return !error;
}

[[nodiscard]] bool WriteShellMapFile(
    const std::filesystem::path& path,
    const char* field_name,
    const char* time_label,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    double (*accessor)(const dec3d::state::CanonicalState&, std::size_t, std::size_t, std::size_t) noexcept) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=shell_map\n";
  output << "# field=" << field_name << '\n';
  output << "# time=" << time_label << '\n';
  output << "# shell_index=" << shell_index << '\n';
  output << "# theta_cells=" << state.layout.theta_cells << '\n';
  output << "# phi_cells=" << state.layout.phi_cells << '\n';
  output << "# columns=theta_index phi_index theta_center phi_center value\n";
  output << std::setprecision(17);

  for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
    const double theta_center = ThetaCenter(geometry, theta);
    for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
      const double phi_center = PhiCenter(geometry, phi);
      output << theta << ' '
             << phi << ' '
             << theta_center << ' '
             << phi_center << ' '
             << accessor(state, shell_index, theta, phi) << '\n';
    }
  }

  return true;
}

[[nodiscard]] bool WriteRadialProfileFile(
    const std::filesystem::path& path,
    const char* label,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const auto sample = BuildSedovRadialProfileSample(state, geometry, 0u, 0.0);
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=radial_profile\n";
  output << "# label=" << label << '\n';
  output << "# radial_cells=" << sample.radial_centers.size() << '\n';
  output << "# columns=radial_index radial_center rho_avg hydro_only_te_avg velocity_magnitude_avg mom_r_avg e_fluid_total_avg\n";
  output << std::setprecision(17);

  for (std::size_t radial = 0; radial < sample.radial_centers.size(); ++radial) {
    output << radial << ' '
           << sample.radial_centers[radial] << ' '
           << sample.rho_avg[radial] << ' '
           << sample.hydro_only_te_avg[radial] << ' '
           << sample.velocity_magnitude_avg[radial] << ' '
           << sample.mom_r_avg[radial] << ' '
           << sample.e_fluid_total_avg[radial] << '\n';
  }

  return true;
}

[[nodiscard]] std::string SerializeDoubles(const std::vector<double>& values) {
  std::ostringstream stream;
  stream << std::setprecision(17);
  for (std::size_t index = 0; index < values.size(); ++index) {
    if (index > 0u) {
      stream << ',';
    }
    stream << values[index];
  }
  return stream.str();
}

template <typename Accessor>
[[nodiscard]] bool WriteSedovRadialTimeMapFile(
    const std::filesystem::path& path,
    const char* field_name,
    const std::vector<SedovRadialProfileSample>& profile_history,
    Accessor accessor) {
  if (profile_history.empty()) {
    return false;
  }

  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=radial_time_map\n";
  output << "# field=" << field_name << '\n';
  output << "# sample_count=" << profile_history.size() << '\n';
  output << "# radial_cells=" << profile_history.front().radial_centers.size() << '\n';
  output << "# radial_centers=" << SerializeDoubles(profile_history.front().radial_centers) << '\n';
  output << "# columns=time_s reference_shock_radius numerical_shock_radius values...\n";
  output << std::setprecision(17);

  for (const auto& sample : profile_history) {
    output << sample.time_s << ' '
           << sample.reference_shock_radius << ' '
           << sample.numerical_shock_radius;
    for (double value : accessor(sample)) {
      output << ' ' << value;
    }
    output << '\n';
  }

  return true;
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

void CopyLocalStateBackToGlobal(
    const dec3d::state::CanonicalState& local_state,
    const dec3d::mesh::RadialOwnershipSlice& slice,
    dec3d::state::CanonicalState& global_state) {
  for (std::size_t local_radial = 0; local_radial < slice.local_cell_count(); ++local_radial) {
    const std::size_t global_radial = slice.begin_index + local_radial;
    for (std::size_t theta = 0; theta < global_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < global_state.layout.phi_cells; ++phi) {
        global_state.rho(global_radial, theta, phi) = local_state.rho(local_radial, theta, phi);
        global_state.mom_r(global_radial, theta, phi) = local_state.mom_r(local_radial, theta, phi);
        global_state.mom_theta(global_radial, theta, phi) =
            local_state.mom_theta(local_radial, theta, phi);
        global_state.mom_phi(global_radial, theta, phi) =
            local_state.mom_phi(local_radial, theta, phi);
        global_state.e_fluid_total(global_radial, theta, phi) =
            local_state.e_fluid_total(local_radial, theta, phi);
        global_state.e_electron(global_radial, theta, phi) =
            local_state.e_electron(local_radial, theta, phi);
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildSliceGhostOverride(
    const dec3d::state::CanonicalState& global_state,
    const dec3d::mesh::RadialOwnershipSlice& slice,
    std::size_t ghost_layers) {
  dec3d::hydro::RadialGhostOverride ghost_override;
  ghost_override.ghost_layers = ghost_layers;
  ghost_override.report_line = "simulated_radial_rank_override";

  const std::size_t theta_cells = global_state.layout.theta_cells;
  const std::size_t phi_cells = global_state.layout.phi_cells;

  if (slice.begin_index >= ghost_layers) {
    ghost_override.has_inner_neighbor = true;
    ghost_override.inner_ghost_states = dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
        ghost_layers,
        theta_cells,
        phi_cells);
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      const std::size_t global_radial = slice.begin_index - ghost_layers + ghost;
      for (std::size_t theta = 0; theta < theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < phi_cells; ++phi) {
          ghost_override.inner_ghost_states(ghost, theta, phi) =
              LoadConservativeState(global_state, global_radial, theta, phi);
        }
      }
    }
  }

  if (slice.end_index + ghost_layers <= global_state.layout.radial_cells) {
    ghost_override.has_outer_neighbor = true;
    ghost_override.outer_ghost_states = dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
        ghost_layers,
        theta_cells,
        phi_cells);
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      const std::size_t global_radial = slice.end_index + ghost;
      for (std::size_t theta = 0; theta < theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < phi_cells; ++phi) {
          ghost_override.outer_ghost_states(ghost, theta, phi) =
              LoadConservativeState(global_state, global_radial, theta, phi);
        }
      }
    }
  }

  return ghost_override;
}

struct SimulatedSplitRankHydroResult {
  bool success{false};
  dec3d::state::CanonicalState after;
  dec3d::mesh::RadialOwnershipDecomposition decomposition;
  std::string failure_reason;
};

[[nodiscard]] SimulatedSplitRankHydroResult SimulateSplitRankRadialHydro(
    const dec3d::state::CanonicalState& before,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    std::size_t rank_count) {
  SimulatedSplitRankHydroResult result;
  result.after = dec3d::state::CanonicalState::Create(before.layout);
  result.decomposition = dec3d::mesh::BuildRadialOwnership(before.layout.radial_cells, rank_count);

  if (!result.decomposition.is_valid()) {
    result.failure_reason = result.decomposition.failure_reason.empty()
                                ? "radial ownership decomposition failed"
                                : result.decomposition.failure_reason;
    return result;
  }

  for (const auto& slice : result.decomposition.slices) {
    auto local_state = BuildLocalStateSlice(before, slice);
    const auto local_geometry = BuildLocalGeometrySlice(
        geometry,
        slice,
        before.layout.theta_cells,
        before.layout.phi_cells);
    if (!local_geometry.is_valid()) {
      result.failure_reason = "local spherical geometry slice is invalid";
      return result;
    }

    auto local_hydro_view = dec3d::state::BuildHydroWorkView(local_state);
    if (!local_hydro_view.is_complete()) {
      result.failure_reason = local_hydro_view.failure_reason.empty()
                                  ? "local hydro work view construction failed"
                                  : local_hydro_view.failure_reason;
      return result;
    }

    const auto ghost_override = BuildSliceGhostOverride(before, slice, 1u);
    const auto hydro_result = dec3d::hydro::AdvanceStaticGridHydro(
        local_hydro_view,
        local_geometry,
        dt_s,
        ghost_override);
    if (!hydro_result.success) {
      result.failure_reason = hydro_result.failure_reason.empty()
                                  ? "local split-rank static-grid hydro failed"
                                  : hydro_result.failure_reason;
      return result;
    }

    const auto writeback = dec3d::state::CommitHydroWriteback(
        local_state,
        local_hydro_view,
        dec3d::state::BuildHydroAuthorizedWriteMask());
    if (!writeback.success) {
      result.failure_reason = writeback.failure_reason.empty()
                                  ? "local split-rank hydro writeback failed"
                                  : writeback.failure_reason;
      return result;
    }

    CopyLocalStateBackToGlobal(local_state, slice, result.after);
  }

  result.success = true;
  return result;
}

[[nodiscard]] double AccessRho(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return state.rho(radial, theta, phi);
}

[[nodiscard]] double AccessHydroOnlyTe(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return HydroOnlyElectronTemperature(state, radial, theta, phi);
}

[[nodiscard]] double AccessVelocityMagnitude(
    const dec3d::state::CanonicalState& state,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return VelocityMagnitude(state, radial, theta, phi);
}

[[nodiscard]] bool WritePolarSliceFile(
    const std::filesystem::path& path,
    const char* pole_name,
    const char* time_label,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    std::size_t theta_index) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=polar_slice\n";
  output << "# pole=" << pole_name << '\n';
  output << "# time=" << time_label << '\n';
  output << "# shell_index=" << shell_index << '\n';
  output << "# theta_index=" << theta_index << '\n';
  output << "# columns=phi_index phi_center rho hydro_only_te velocity_magnitude mom_theta mom_phi\n";
  output << std::setprecision(17);

  for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
    output << phi << ' '
           << PhiCenter(geometry, phi) << ' '
           << state.rho(shell_index, theta_index, phi) << ' '
           << HydroOnlyElectronTemperature(state, shell_index, theta_index, phi) << ' '
           << VelocityMagnitude(state, shell_index, theta_index, phi) << ' '
           << state.mom_theta(shell_index, theta_index, phi) << ' '
           << state.mom_phi(shell_index, theta_index, phi) << '\n';
  }

  return true;
}

[[nodiscard]] bool WriteModeProjectionFile(
    const std::filesystem::path& path,
    const ModeProjection& before,
    const ModeProjection& after_full,
    const ModeProjection& after_radial_only) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=mode_projection\n";
  output << "# basis=sin(theta)*cos(phi),sin(theta)*sin(phi)\n";
  output << "# columns=label amplitude phase cosine_coefficient sine_coefficient\n";
  output << std::setprecision(17);
  output << "before "
         << before.amplitude << ' '
         << before.phase << ' '
         << before.cosine_coefficient << ' '
         << before.sine_coefficient << '\n';
  output << "after_full "
         << after_full.amplitude << ' '
         << after_full.phase << ' '
         << after_full.cosine_coefficient << ' '
         << after_full.sine_coefficient << '\n';
  output << "after_radial_only "
         << after_radial_only.amplitude << ' '
         << after_radial_only.phase << ' '
         << after_radial_only.cosine_coefficient << ' '
         << after_radial_only.sine_coefficient << '\n';
  return true;
}

[[nodiscard]] bool WriteKeyValueFile(
    const std::filesystem::path& path,
    const std::vector<std::pair<std::string, std::string>>& entries,
    const char* kind = nullptr) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  if (kind != nullptr) {
    output << "# kind=" << kind << '\n';
  }

  for (const auto& [key, value] : entries) {
    output << key << '=' << value << '\n';
  }

  return true;
}

[[nodiscard]] bool WriteBudgetResidualFile(
    const std::filesystem::path& path,
    const char* label,
    const dec3d::hydro::HydroBudgetResidualSummary& summary) {
  return WriteKeyValueFile(
      path,
      {
          {"label", label},
          {"old_mass", std::to_string(summary.old_mass)},
          {"flux_mass_delta", std::to_string(summary.flux_mass_delta)},
          {"source_mass_delta", std::to_string(summary.source_mass_delta)},
          {"new_mass", std::to_string(summary.new_mass)},
          {"mass_residual", std::to_string(summary.mass_residual)},
          {"old_mom_r", std::to_string(summary.old_mom_r)},
          {"flux_mom_r_delta", std::to_string(summary.flux_mom_r_delta)},
          {"source_mom_r_delta", std::to_string(summary.source_mom_r_delta)},
          {"new_mom_r", std::to_string(summary.new_mom_r)},
          {"mom_r_residual", std::to_string(summary.mom_r_residual)},
          {"old_e_fluid_total", std::to_string(summary.old_e_fluid_total)},
          {"flux_e_fluid_total_delta", std::to_string(summary.flux_e_fluid_total_delta)},
          {"source_e_fluid_total_delta", std::to_string(summary.source_e_fluid_total_delta)},
          {"new_e_fluid_total", std::to_string(summary.new_e_fluid_total)},
          {"e_fluid_total_residual", std::to_string(summary.e_fluid_total_residual)},
      },
      "budget_residual");
}

[[nodiscard]] bool WriteBudgetResidualComparisonFile(
    const std::filesystem::path& path,
    const dec3d::hydro::HydroBudgetResidualSummary& single_rank_budget,
    const dec3d::hydro::HydroBudgetResidualSummary& multi_rank_budget) {
  return WriteKeyValueFile(
      path,
      {
          {"single_rank_mass_residual", std::to_string(single_rank_budget.mass_residual)},
          {"multi_rank_mass_residual", std::to_string(multi_rank_budget.mass_residual)},
          {"mass_residual_difference", std::to_string(std::abs(
               single_rank_budget.mass_residual - multi_rank_budget.mass_residual))},
          {"single_rank_mom_r_residual", std::to_string(single_rank_budget.mom_r_residual)},
          {"multi_rank_mom_r_residual", std::to_string(multi_rank_budget.mom_r_residual)},
          {"mom_r_residual_difference", std::to_string(std::abs(
               single_rank_budget.mom_r_residual - multi_rank_budget.mom_r_residual))},
          {"single_rank_e_fluid_total_residual", std::to_string(single_rank_budget.e_fluid_total_residual)},
          {"multi_rank_e_fluid_total_residual", std::to_string(multi_rank_budget.e_fluid_total_residual)},
          {"e_fluid_total_residual_difference", std::to_string(std::abs(
               single_rank_budget.e_fluid_total_residual - multi_rank_budget.e_fluid_total_residual))},
      },
      "budget_residual_comparison");
}

}  // namespace

bool SphericalRadialWaveSanityCheck::is_complete() const noexcept {
  return success && !report_line.empty();
}

bool SedovSphericalBlastCheck::is_complete() const noexcept {
  return !report_line.empty();
}

SedovRadialProfileSample BuildSedovRadialProfileSample(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t step_index,
    double time_s) noexcept {
  SedovRadialProfileSample sample;
  sample.step_index = step_index;
  sample.time_s = time_s;

  if (!LayoutAndGeometryMatch(state, geometry)) {
    return sample;
  }

  const double inverse_shell_count =
      1.0 / static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);
  sample.radial_centers.reserve(state.layout.radial_cells);
  sample.rho_avg.reserve(state.layout.radial_cells);
  sample.hydro_only_te_avg.reserve(state.layout.radial_cells);
  sample.velocity_magnitude_avg.reserve(state.layout.radial_cells);
  sample.mom_r_avg.reserve(state.layout.radial_cells);
  sample.e_fluid_total_avg.reserve(state.layout.radial_cells);

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    double rho_sum = 0.0;
    double te_sum = 0.0;
    double velocity_sum = 0.0;
    double mom_r_sum = 0.0;
    double e_fluid_total_sum = 0.0;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        rho_sum += state.rho(radial, theta, phi);
        te_sum += HydroOnlyElectronTemperature(state, radial, theta, phi);
        velocity_sum += VelocityMagnitude(state, radial, theta, phi);
        mom_r_sum += state.mom_r(radial, theta, phi);
        e_fluid_total_sum += state.e_fluid_total(radial, theta, phi);
      }
    }

    sample.radial_centers.push_back(
        0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]));
    sample.rho_avg.push_back(rho_sum * inverse_shell_count);
    sample.hydro_only_te_avg.push_back(te_sum * inverse_shell_count);
    sample.velocity_magnitude_avg.push_back(velocity_sum * inverse_shell_count);
    sample.mom_r_avg.push_back(mom_r_sum * inverse_shell_count);
    sample.e_fluid_total_avg.push_back(e_fluid_total_sum * inverse_shell_count);
  }

  return sample;
}

bool WriteSedovRadialProfileSampleFile(
    const std::filesystem::path& path,
    const char* label,
    const SedovRadialProfileSample& sample) noexcept {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=radial_profile\n";
  output << "# label=" << label << '\n';
  output << "# step_index=" << sample.step_index << '\n';
  output << "# time_s=" << sample.time_s << '\n';
  output << "# reference_shock_radius=" << sample.reference_shock_radius << '\n';
  output << "# numerical_shock_radius=" << sample.numerical_shock_radius << '\n';
  output << "# radial_cells=" << sample.radial_centers.size() << '\n';
  output << "# columns=radial_index radial_center rho_avg hydro_only_te_avg velocity_magnitude_avg mom_r_avg e_fluid_total_avg\n";
  output << std::setprecision(17);

  for (std::size_t radial = 0; radial < sample.radial_centers.size(); ++radial) {
    output << radial << ' '
           << sample.radial_centers[radial] << ' '
           << sample.rho_avg[radial] << ' '
           << sample.hydro_only_te_avg[radial] << ' '
           << sample.velocity_magnitude_avg[radial] << ' '
           << sample.mom_r_avg[radial] << ' '
           << sample.e_fluid_total_avg[radial] << '\n';
  }

  return true;
}

namespace {

[[nodiscard]] double SedovShockRadiusCoefficientGammaFiveThirds() noexcept {
  return 1.15;
}

[[nodiscard]] std::vector<double> BuildShellAveragedField(
    const dec3d::state::CanonicalState& state,
    double (*accessor)(const dec3d::state::CanonicalState&, std::size_t, std::size_t, std::size_t) noexcept) {
  std::vector<double> shell_average(state.layout.radial_cells, 0.0);
  const double inverse_shell_count =
      1.0 / static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        sum += accessor(state, radial, theta, phi);
      }
    }
    shell_average[radial] = sum * inverse_shell_count;
  }
  return shell_average;
}

[[nodiscard]] double ShellMeanRho(
    const dec3d::state::CanonicalState& state,
    std::size_t radial) noexcept {
  double sum = 0.0;
  const double inverse_shell_count =
      1.0 / static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);
  for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
    for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
      sum += state.rho(radial, theta, phi);
    }
  }
  return sum * inverse_shell_count;
}

[[nodiscard]] double EstimateSedovShockRadiusFromState(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  if (!LayoutAndGeometryMatch(state, geometry) || state.layout.radial_cells < 2u) {
    return 0.0;
  }

  const auto rho_profile = BuildShellAveragedField(state, &AccessRho);
  std::size_t shock_index = 0u;
  double max_positive_jump = -1.0;
  for (std::size_t radial = 0; radial + 1u < rho_profile.size(); ++radial) {
    const double jump = rho_profile[radial + 1u] - rho_profile[radial];
    if (jump > max_positive_jump) {
      max_positive_jump = jump;
      shock_index = radial;
    }
  }

  return 0.5 * (geometry.radial_faces[shock_index + 1u] + geometry.radial_faces[shock_index + 2u]);
}

[[nodiscard]] bool WriteShockRadiusHistoryFile(
    const std::filesystem::path& path,
    const std::vector<SedovShockRadiusSample>& history) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }
  output << "# kind=shock_radius_history\n";
  output << "# columns=time_s numerical_radius reference_radius relative_error\n";
  output << std::setprecision(17);
  for (const auto& sample : history) {
    output << sample.time_s << ' '
           << sample.numerical_shock_radius << ' '
           << sample.reference_shock_radius << ' '
           << sample.relative_error << '\n';
  }
  return true;
}

[[nodiscard]] bool WriteRunProgressFile(
    const std::filesystem::path& path,
    const std::vector<SedovShockRadiusSample>& history) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }
  output << "# kind=run_progress\n";
  output << "# columns=step time_s status\n";
  output << std::setprecision(17);
  for (std::size_t index = 0; index < history.size(); ++index) {
    output << (index + 1u) << ' '
           << history[index].time_s << ' '
           << "pending\n";
  }
  return true;
}

}  // namespace

SedovSphericalBlastCheck EvaluateSedovSphericalBlast(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double blast_energy,
    double ambient_density,
    double final_time_s,
    double shock_radius_relative_tolerance,
    double shell_symmetry_tolerance,
    double tangential_momentum_ratio_tolerance) noexcept {
  SedovSphericalBlastCheck check;
  check.final_time_s = final_time_s;

  if (!HasAllocatedHydroState(before) ||
      !HasAllocatedHydroState(after) ||
      !HasAllocatedElectronState(before) ||
      !HasAllocatedElectronState(after)) {
    check.failure_reason = "Sedov benchmark requires allocated hydro and electron state";
  } else if (!LayoutsMatch(before, after) || !LayoutAndGeometryMatch(before, geometry)) {
    check.failure_reason = "Sedov benchmark requires matching layouts and geometry";
  } else if (!(blast_energy > 0.0) || !(ambient_density > 0.0) || !(final_time_s > 0.0)) {
    check.failure_reason = "Sedov benchmark requires positive blast energy, ambient density, and final time";
  } else {
    check.numerical_shock_radius = EstimateSedovShockRadiusFromState(after, geometry);
    check.reference_shock_radius =
        SedovShockRadiusCoefficientGammaFiveThirds() *
        std::pow(blast_energy * final_time_s * final_time_s / ambient_density, 0.2);
    check.shock_radius_relative_error =
        check.reference_shock_radius > 0.0
            ? std::abs(check.numerical_shock_radius - check.reference_shock_radius) /
                  check.reference_shock_radius
            : 0.0;
    check.shock_radius_within_tolerance =
        check.shock_radius_relative_error <= shock_radius_relative_tolerance;

    double radial_momentum_l1 = 0.0;
    double tangential_momentum_l1 = 0.0;
    for (std::size_t radial = 0; radial < after.layout.radial_cells; ++radial) {
      const double shell_mean = ShellMeanRho(after, radial);
      if (!(shell_mean > 0.0) || !std::isfinite(shell_mean)) {
        check.failure_reason = "Sedov benchmark encountered non-physical shell-mean density";
        break;
      }
      for (std::size_t theta = 0; theta < after.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < after.layout.phi_cells; ++phi) {
          check.max_rho_angular_relative_spread = std::max(
              check.max_rho_angular_relative_spread,
              std::abs(after.rho(radial, theta, phi) - shell_mean) / shell_mean);
          radial_momentum_l1 += std::abs(after.mom_r(radial, theta, phi));
          tangential_momentum_l1 +=
              std::abs(after.mom_theta(radial, theta, phi)) +
              std::abs(after.mom_phi(radial, theta, phi));
        }
      }
    }

    if (check.failure_reason.empty()) {
      check.shell_symmetry_preserved =
          check.max_rho_angular_relative_spread <= shell_symmetry_tolerance;
      check.tangential_to_radial_momentum_ratio =
          radial_momentum_l1 > 0.0 ? tangential_momentum_l1 / radial_momentum_l1
                                   : std::numeric_limits<double>::infinity();
      check.tangential_momentum_quiet =
          check.tangential_to_radial_momentum_ratio <= tangential_momentum_ratio_tolerance;
      check.success =
          check.shock_radius_within_tolerance &&
          check.shell_symmetry_preserved &&
          check.tangential_momentum_quiet;

      if (!check.success) {
        if (!check.shock_radius_within_tolerance) {
          check.failure_reason = "Sedov shock radius exceeded tolerance";
        } else if (!check.shell_symmetry_preserved) {
          check.failure_reason = "Sedov shell symmetry drift exceeded tolerance";
        } else if (!check.tangential_momentum_quiet) {
          check.failure_reason = "Sedov tangential momentum ratio exceeded tolerance";
        }
      }
    }
  }

  std::ostringstream report;
  report << "sedov_success=" << (check.success ? "true" : "false")
         << "; shock_radius_within_tolerance=" << (check.shock_radius_within_tolerance ? "true" : "false")
         << "; shell_symmetry_preserved=" << (check.shell_symmetry_preserved ? "true" : "false")
         << "; tangential_momentum_quiet=" << (check.tangential_momentum_quiet ? "true" : "false")
         << "; final_time_s=" << check.final_time_s
         << "; numerical_shock_radius=" << check.numerical_shock_radius
         << "; reference_shock_radius=" << check.reference_shock_radius
         << "; shock_radius_relative_error=" << check.shock_radius_relative_error
         << "; max_rho_angular_relative_spread=" << check.max_rho_angular_relative_spread
         << "; tangential_to_radial_momentum_ratio=" << check.tangential_to_radial_momentum_ratio;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();
  return check;
}

bool WriteSedovSphericalBlastOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const std::filesystem::path& output_directory,
    const SedovSphericalBlastCheck& summary,
    const std::vector<SedovShockRadiusSample>& shock_radius_history,
    const std::vector<SedovRadialProfileSample>& profile_history,
    const HydroBudgetResidualSummary* budget,
    const char* case_name) noexcept {
  if (!EnsureOutputDirectory(output_directory) || !geometry.is_valid()) {
    return false;
  }

  const bool have_budget = budget != nullptr && budget->is_complete();
  const std::string case_name_string =
      case_name == nullptr ? "case_sedov_spherical" : case_name;
  const std::size_t shell_index = std::min<std::size_t>(after.layout.radial_cells / 2u, after.layout.radial_cells - 1u);

  bool wrote_profile_history = true;
  for (const auto& sample : profile_history) {
    std::ostringstream file_name;
    file_name << "radial_profile_step"
              << std::setw(6) << std::setfill('0') << sample.step_index
              << ".txt";
    std::ostringstream label;
    label << "step_" << sample.step_index;
    wrote_profile_history =
        wrote_profile_history &&
        WriteSedovRadialProfileSampleFile(
            output_directory / file_name.str(),
            label.str().c_str(),
            sample);
  }

  return
      wrote_profile_history &&
      WriteShockRadiusHistoryFile(output_directory / "shock_radius_vs_time.txt", shock_radius_history) &&
      WriteRunProgressFile(output_directory / "dec3d.out", shock_radius_history) &&
      WriteRadialProfileFile(output_directory / "radial_profile_t0.txt", "t0", before, geometry) &&
      WriteRadialProfileFile(output_directory / "radial_profile_t1.txt", "t1", after, geometry) &&
      WriteSedovRadialTimeMapFile(
          output_directory / "rho_rt_map.txt",
          "rho",
          profile_history,
          [](const SedovRadialProfileSample& sample) -> const std::vector<double>& {
            return sample.rho_avg;
          }) &&
      WriteSedovRadialTimeMapFile(
          output_directory / "hydro_only_te_rt_map.txt",
          "hydro_only_te",
          profile_history,
          [](const SedovRadialProfileSample& sample) -> const std::vector<double>& {
            return sample.hydro_only_te_avg;
          }) &&
      WriteSedovRadialTimeMapFile(
          output_directory / "velocity_magnitude_rt_map.txt",
          "velocity_magnitude",
          profile_history,
          [](const SedovRadialProfileSample& sample) -> const std::vector<double>& {
            return sample.velocity_magnitude_avg;
          }) &&
      WriteSedovRadialTimeMapFile(
          output_directory / "mom_r_rt_map.txt",
          "mom_r",
          profile_history,
          [](const SedovRadialProfileSample& sample) -> const std::vector<double>& {
            return sample.mom_r_avg;
          }) &&
      WriteSedovRadialTimeMapFile(
          output_directory / "e_fluid_total_rt_map.txt",
          "e_fluid_total",
          profile_history,
          [](const SedovRadialProfileSample& sample) -> const std::vector<double>& {
            return sample.e_fluid_total_avg;
          }) &&
      WriteShellMapFile(output_directory / "shell_map_rho_t0.txt", "rho", "t0", before, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "shell_map_rho_t1.txt", "rho", "t1", after, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "shell_map_te_t0.txt", "hydro_only_te", "t0", before, geometry, shell_index, &AccessHydroOnlyTe) &&
      WriteShellMapFile(output_directory / "shell_map_te_t1.txt", "hydro_only_te", "t1", after, geometry, shell_index, &AccessHydroOnlyTe) &&
      WriteShellMapFile(output_directory / "shell_map_velocity_t0.txt", "velocity_magnitude", "t0", before, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteShellMapFile(output_directory / "shell_map_velocity_t1.txt", "velocity_magnitude", "t1", after, geometry, shell_index, &AccessVelocityMagnitude) &&
      (!have_budget || WriteBudgetResidualFile(output_directory / "budget_residual.txt", "sedov", *budget)) &&
      WriteKeyValueFile(
          output_directory / "summary.txt",
          {
              {"case", case_name_string},
              {"success", summary.success ? "true" : "false"},
              {"shock_radius_within_tolerance", summary.shock_radius_within_tolerance ? "true" : "false"},
              {"shell_symmetry_preserved", summary.shell_symmetry_preserved ? "true" : "false"},
              {"tangential_momentum_quiet", summary.tangential_momentum_quiet ? "true" : "false"},
              {"final_time_s", std::to_string(summary.final_time_s)},
              {"numerical_shock_radius", std::to_string(summary.numerical_shock_radius)},
              {"reference_shock_radius", std::to_string(summary.reference_shock_radius)},
              {"shock_radius_relative_error", std::to_string(summary.shock_radius_relative_error)},
              {"report_line", summary.report_line},
          }) &&
      WriteKeyValueFile(
          output_directory / "case_manifest.txt",
          {
              {"case", case_name_string},
              {"summary", "summary.txt"},
              {"run_progress", "dec3d.out"},
              {"shock_radius_history", "shock_radius_vs_time.txt"},
              {"radial_profile_series", "radial_profile_step*.txt"},
              {"radial_profile_t0", "radial_profile_t0.txt"},
              {"radial_profile_t1", "radial_profile_t1.txt"},
              {"rho_rt_map", "rho_rt_map.txt"},
              {"hydro_only_te_rt_map", "hydro_only_te_rt_map.txt"},
              {"velocity_magnitude_rt_map", "velocity_magnitude_rt_map.txt"},
              {"mom_r_rt_map", "mom_r_rt_map.txt"},
              {"e_fluid_total_rt_map", "e_fluid_total_rt_map.txt"},
              {"shell_map_rho_t0", "shell_map_rho_t0.txt"},
              {"shell_map_rho_t1", "shell_map_rho_t1.txt"},
              {"shell_map_te_t0", "shell_map_te_t0.txt"},
              {"shell_map_te_t1", "shell_map_te_t1.txt"},
              {"shell_map_velocity_t0", "shell_map_velocity_t0.txt"},
              {"shell_map_velocity_t1", "shell_map_velocity_t1.txt"},
              {"budget_residual", have_budget ? "budget_residual.txt" : ""},
          },
          "case_manifest");
}

SphericalRadialWaveSanityCheck EvaluateSphericalRadialWaveSanity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    double angular_tolerance,
    double tangential_momentum_tolerance,
    double radial_change_tolerance) noexcept {
  SphericalRadialWaveSanityCheck check;

  if (!HasAllocatedHydroState(before) || !HasAllocatedHydroState(after)) {
    check.failure_reason = "radial wave sanity requires allocated hydro state";
  } else if (!LayoutsMatch(before, after)) {
    check.failure_reason = "radial wave sanity requires matching before/after layouts";
  } else if (before.layout.radial_cells == 0u ||
             before.layout.theta_cells == 0u ||
             before.layout.phi_cells == 0u) {
    check.failure_reason = "radial wave sanity requires non-empty spherical state";
  } else {
    const auto radial_cells = before.layout.radial_cells;
    const auto theta_cells = before.layout.theta_cells;
    const auto phi_cells = before.layout.phi_cells;

    for (std::size_t radial = 0; radial < radial_cells; ++radial) {
      double rho_min = after.rho(radial, 0, 0);
      double rho_max = rho_min;
      double mom_r_min = after.mom_r(radial, 0, 0);
      double mom_r_max = mom_r_min;
      double e_min = after.e_fluid_total(radial, 0, 0);
      double e_max = e_min;

      for (std::size_t theta = 0; theta < theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < phi_cells; ++phi) {
          rho_min = std::min(rho_min, after.rho(radial, theta, phi));
          rho_max = std::max(rho_max, after.rho(radial, theta, phi));
          mom_r_min = std::min(mom_r_min, after.mom_r(radial, theta, phi));
          mom_r_max = std::max(mom_r_max, after.mom_r(radial, theta, phi));
          e_min = std::min(e_min, after.e_fluid_total(radial, theta, phi));
          e_max = std::max(e_max, after.e_fluid_total(radial, theta, phi));

          check.max_abs_mom_theta = std::max(
              check.max_abs_mom_theta,
              std::abs(after.mom_theta(radial, theta, phi)));
          check.max_abs_mom_phi = std::max(
              check.max_abs_mom_phi,
              std::abs(after.mom_phi(radial, theta, phi)));
          check.max_abs_rho_delta = std::max(
              check.max_abs_rho_delta,
              std::abs(after.rho(radial, theta, phi) - before.rho(radial, theta, phi)));
          check.max_abs_mom_r = std::max(
              check.max_abs_mom_r,
              std::abs(after.mom_r(radial, theta, phi)));
        }
      }

      check.max_rho_angular_spread = std::max(check.max_rho_angular_spread, rho_max - rho_min);
      check.max_mom_r_angular_spread = std::max(check.max_mom_r_angular_spread, mom_r_max - mom_r_min);
      check.max_e_fluid_total_angular_spread = std::max(check.max_e_fluid_total_angular_spread, e_max - e_min);
    }

    check.angular_symmetry_preserved =
        check.max_rho_angular_spread <= angular_tolerance &&
        check.max_mom_r_angular_spread <= angular_tolerance &&
        check.max_e_fluid_total_angular_spread <= angular_tolerance;
    check.tangential_momentum_quiet =
        check.max_abs_mom_theta <= tangential_momentum_tolerance &&
        check.max_abs_mom_phi <= tangential_momentum_tolerance;
    check.radial_profile_changed =
        check.max_abs_rho_delta > radial_change_tolerance &&
        check.max_abs_mom_r > radial_change_tolerance;
    check.success =
        check.angular_symmetry_preserved &&
        check.tangential_momentum_quiet &&
        check.radial_profile_changed;

    if (!check.success) {
      if (!check.angular_symmetry_preserved) {
        check.failure_reason = "angular shell symmetry was not preserved";
      } else if (!check.tangential_momentum_quiet) {
        check.failure_reason = "tangential momentum grew beyond tolerance";
      } else if (!check.radial_profile_changed) {
        check.failure_reason = "radial profile did not evolve enough to qualify as a wave/shock sanity case";
      }
    }
  }

  std::ostringstream report;
  report << "spherical_radial_wave_success=" << (check.success ? "true" : "false")
         << "; angular_symmetry_preserved=" << (check.angular_symmetry_preserved ? "true" : "false")
         << "; tangential_momentum_quiet=" << (check.tangential_momentum_quiet ? "true" : "false")
         << "; radial_profile_changed=" << (check.radial_profile_changed ? "true" : "false")
         << "; max_rho_angular_spread=" << check.max_rho_angular_spread
         << "; max_mom_r_angular_spread=" << check.max_mom_r_angular_spread
         << "; max_e_fluid_total_angular_spread=" << check.max_e_fluid_total_angular_spread
         << "; max_abs_mom_theta=" << check.max_abs_mom_theta
         << "; max_abs_mom_phi=" << check.max_abs_mom_phi
         << "; max_abs_rho_delta=" << check.max_abs_rho_delta
         << "; max_abs_mom_r=" << check.max_abs_mom_r;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

bool SphericalRadialWaveMpiParityCheck::is_complete() const noexcept {
  return !report_line.empty();
}

bool RadialAleOnOffCheck::is_complete() const noexcept {
  return !report_line.empty();
}

bool DirectionalMpiParityCheck::is_complete() const noexcept {
  return !report_line.empty();
}

RadialAleOnOffCheck EvaluateRadialAleOnOff(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& ale_off_after,
    const dec3d::state::CanonicalState& ale_on_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::MeshUpdateProposal& proposal,
    double difference_tolerance,
    bool ppm_executed,
    std::size_t reconstruction_ghost_layers) noexcept {
  RadialAleOnOffCheck check;
  check.ppm_executed = ppm_executed;
  check.reconstruction_ghost_layers = reconstruction_ghost_layers;
  check.reconstruction_mode =
      ppm_executed ? "radial_only_ppm_ale" : "radial_only_first_order_ale";

  if (!HasAllocatedHydroState(before) ||
      !HasAllocatedHydroState(ale_off_after) ||
      !HasAllocatedHydroState(ale_on_after) ||
      !HasAllocatedElectronState(before) ||
      !HasAllocatedElectronState(ale_off_after) ||
      !HasAllocatedElectronState(ale_on_after)) {
    check.failure_reason = "radial ALE on/off sanity requires allocated hydro and electron state";
  } else if (!LayoutsMatch(before, ale_off_after) || !LayoutsMatch(before, ale_on_after)) {
    check.failure_reason = "radial ALE on/off sanity requires matching layouts";
  } else if (!LayoutAndGeometryMatch(before, geometry)) {
    check.failure_reason = "radial ALE on/off sanity requires matching spherical geometry";
  } else {
    check.proposal_complete = proposal.is_complete(geometry.radial_faces.size());
    if (check.proposal_complete) {
      check.outer_face_displacement =
          proposal.proposed_radial_faces.back() - geometry.radial_faces.back();
      for (double face_speed : proposal.radial_face_velocities) {
        check.max_face_speed = std::max(check.max_face_speed, std::abs(face_speed));
      }
      check.moving_mesh_applied =
          std::abs(check.outer_face_displacement) > difference_tolerance &&
          check.max_face_speed > difference_tolerance;
    }

    for (std::size_t radial = 0; radial < before.layout.radial_cells; ++radial) {
      for (std::size_t theta = 0; theta < before.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < before.layout.phi_cells; ++phi) {
          check.max_rho_difference = std::max(
              check.max_rho_difference,
              std::abs(ale_on_after.rho(radial, theta, phi) -
                       ale_off_after.rho(radial, theta, phi)));
          check.max_mom_r_difference = std::max(
              check.max_mom_r_difference,
              std::abs(ale_on_after.mom_r(radial, theta, phi) -
                       ale_off_after.mom_r(radial, theta, phi)));
          check.max_e_fluid_total_difference = std::max(
              check.max_e_fluid_total_difference,
              std::abs(ale_on_after.e_fluid_total(radial, theta, phi) -
                       ale_off_after.e_fluid_total(radial, theta, phi)));
        }
      }
    }

    check.ale_changes_solution =
        check.max_rho_difference > difference_tolerance ||
        check.max_mom_r_difference > difference_tolerance ||
        check.max_e_fluid_total_difference > difference_tolerance;
    check.success =
        check.proposal_complete &&
        check.moving_mesh_applied &&
        check.ale_changes_solution;

    if (!check.success) {
      if (!check.proposal_complete) {
        check.failure_reason = "radial ALE proposal is incomplete for on/off sanity";
      } else if (!check.moving_mesh_applied) {
        check.failure_reason = "radial ALE proposal did not move the mesh enough";
      } else if (!check.ale_changes_solution) {
        check.failure_reason = "ALE on/off radial wave produced indistinguishable results";
      }
    }
  }

  std::ostringstream report;
  report << "radial_ale_on_off_success=" << (check.success ? "true" : "false")
         << "; proposal_complete=" << (check.proposal_complete ? "true" : "false")
         << "; moving_mesh_applied=" << (check.moving_mesh_applied ? "true" : "false")
         << "; ale_changes_solution=" << (check.ale_changes_solution ? "true" : "false")
         << "; ppm_executed=" << (check.ppm_executed ? "true" : "false")
         << "; reconstruction_ghost_layers=" << check.reconstruction_ghost_layers
         << "; reconstruction_mode=" << check.reconstruction_mode
         << "; max_rho_difference=" << check.max_rho_difference
         << "; max_mom_r_difference=" << check.max_mom_r_difference
         << "; max_e_fluid_total_difference=" << check.max_e_fluid_total_difference
         << "; outer_face_displacement=" << check.outer_face_displacement
         << "; max_face_speed=" << check.max_face_speed;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

SphericalRadialWaveMpiParityCheck EvaluateSphericalRadialWaveMpiParity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    std::size_t rank_count,
    double parity_tolerance,
    double seam_tolerance,
    dec3d::state::CanonicalState* multi_rank_after_out) noexcept {
  SphericalRadialWaveMpiParityCheck check;
  check.rank_count = rank_count;

  if (!HasAllocatedHydroState(before) || !HasAllocatedHydroState(single_rank_after)) {
    check.failure_reason = "radial wave mpi parity requires allocated hydro state";
  } else if (!HasAllocatedElectronState(before) || !HasAllocatedElectronState(single_rank_after)) {
    check.failure_reason = "radial wave mpi parity requires allocated electron state";
  } else if (!LayoutsMatch(before, single_rank_after) || !LayoutAndGeometryMatch(before, geometry)) {
    check.failure_reason = "radial wave mpi parity requires matching layouts and geometry";
  } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
    check.failure_reason = "radial wave mpi parity requires a finite positive dt";
  } else {
    const auto simulated = SimulateSplitRankRadialHydro(before, geometry, dt_s, rank_count);
    check.rank_decomposition_valid = simulated.decomposition.is_valid();
    if (!simulated.success) {
      check.failure_reason = simulated.failure_reason;
    } else {
      const auto& multi_rank_after = simulated.after;
      if (multi_rank_after_out != nullptr) {
        *multi_rank_after_out = multi_rank_after;
      }

      for (std::size_t radial = 0; radial < before.layout.radial_cells; ++radial) {
        for (std::size_t theta = 0; theta < before.layout.theta_cells; ++theta) {
          for (std::size_t phi = 0; phi < before.layout.phi_cells; ++phi) {
            check.max_rho_difference = std::max(
                check.max_rho_difference,
                std::abs(single_rank_after.rho(radial, theta, phi) -
                         multi_rank_after.rho(radial, theta, phi)));
            check.max_mom_r_difference = std::max(
                check.max_mom_r_difference,
                std::abs(single_rank_after.mom_r(radial, theta, phi) -
                         multi_rank_after.mom_r(radial, theta, phi)));
            check.max_e_fluid_total_difference = std::max(
                check.max_e_fluid_total_difference,
                std::abs(single_rank_after.e_fluid_total(radial, theta, phi) -
                         multi_rank_after.e_fluid_total(radial, theta, phi)));
          }
        }
      }

      for (std::size_t slice_index = 0; slice_index + 1u < simulated.decomposition.slices.size(); ++slice_index) {
        const auto& left_slice = simulated.decomposition.slices[slice_index];
        const std::size_t left_radial = left_slice.end_index - 1u;
        const std::size_t right_radial = left_slice.end_index;
        for (std::size_t theta = 0; theta < before.layout.theta_cells; ++theta) {
          for (std::size_t phi = 0; phi < before.layout.phi_cells; ++phi) {
            const double single_rho_jump =
                std::abs(single_rank_after.rho(left_radial, theta, phi) -
                         single_rank_after.rho(right_radial, theta, phi));
            const double multi_rho_jump =
                std::abs(multi_rank_after.rho(left_radial, theta, phi) -
                         multi_rank_after.rho(right_radial, theta, phi));
            check.max_seam_rho_jump = std::max(
                check.max_seam_rho_jump,
                std::abs(multi_rho_jump - single_rho_jump));

            const double single_velocity_jump =
                std::abs(VelocityMagnitude(single_rank_after, left_radial, theta, phi) -
                         VelocityMagnitude(single_rank_after, right_radial, theta, phi));
            const double multi_velocity_jump =
                std::abs(VelocityMagnitude(multi_rank_after, left_radial, theta, phi) -
                         VelocityMagnitude(multi_rank_after, right_radial, theta, phi));
            check.max_seam_velocity_jump = std::max(
                check.max_seam_velocity_jump,
                std::abs(multi_velocity_jump - single_velocity_jump));
          }
        }
      }

      check.parity_within_tolerance =
          check.max_rho_difference <= parity_tolerance &&
          check.max_mom_r_difference <= parity_tolerance &&
          check.max_e_fluid_total_difference <= parity_tolerance;
      check.seam_jumps_bounded =
          check.max_seam_rho_jump <= seam_tolerance &&
          check.max_seam_velocity_jump <= seam_tolerance;
      check.success =
          check.rank_decomposition_valid &&
          check.parity_within_tolerance &&
          check.seam_jumps_bounded;

      if (!check.success) {
        if (!check.rank_decomposition_valid) {
          check.failure_reason = "radial ownership decomposition is invalid";
        } else if (!check.parity_within_tolerance) {
          check.failure_reason = "multi-rank radial wave parity drift exceeded tolerance";
        } else if (!check.seam_jumps_bounded) {
          check.failure_reason = "radial seam jump mismatch exceeded tolerance";
        }
      }
    }
  }

  std::ostringstream report;
  report << "radial_wave_mpi_parity_success=" << (check.success ? "true" : "false")
         << "; rank_count=" << check.rank_count
         << "; rank_decomposition_valid=" << (check.rank_decomposition_valid ? "true" : "false")
         << "; parity_within_tolerance=" << (check.parity_within_tolerance ? "true" : "false")
         << "; seam_jumps_bounded=" << (check.seam_jumps_bounded ? "true" : "false")
         << "; macro_zoning_executed=" << (check.macro_zoning_executed ? "true" : "false")
         << "; max_rho_difference=" << check.max_rho_difference
         << "; max_mom_r_difference=" << check.max_mom_r_difference
         << "; max_e_fluid_total_difference=" << check.max_e_fluid_total_difference
         << "; max_seam_rho_jump_mismatch=" << check.max_seam_rho_jump
         << "; max_seam_velocity_jump_mismatch=" << check.max_seam_velocity_jump
         << "; macro_zoning_mode=" << check.macro_zoning_mode;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

DirectionalMpiParityCheck EvaluateDirectionalMpiParity(
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::state::CanonicalState& multi_rank_after,
    std::size_t rank_count,
    double parity_tolerance) noexcept {
  DirectionalMpiParityCheck check;
  check.rank_count = rank_count;

  if (!HasAllocatedHydroState(single_rank_after) ||
      !HasAllocatedHydroState(multi_rank_after) ||
      !HasAllocatedElectronState(single_rank_after) ||
      !HasAllocatedElectronState(multi_rank_after)) {
    check.failure_reason = "directional mpi parity requires allocated hydro and electron state";
  } else if (!LayoutsMatch(single_rank_after, multi_rank_after)) {
    check.failure_reason = "directional mpi parity requires matching layouts";
  } else {
    for (std::size_t radial = 0; radial < single_rank_after.layout.radial_cells; ++radial) {
      for (std::size_t theta = 0; theta < single_rank_after.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < single_rank_after.layout.phi_cells; ++phi) {
          check.max_rho_difference = std::max(
              check.max_rho_difference,
              std::abs(single_rank_after.rho(radial, theta, phi) -
                       multi_rank_after.rho(radial, theta, phi)));
          check.max_hydro_only_te_difference = std::max(
              check.max_hydro_only_te_difference,
              std::abs(
                  HydroOnlyElectronTemperature(single_rank_after, radial, theta, phi) -
                  HydroOnlyElectronTemperature(multi_rank_after, radial, theta, phi)));
          check.max_velocity_magnitude_difference = std::max(
              check.max_velocity_magnitude_difference,
              std::abs(
                  VelocityMagnitude(single_rank_after, radial, theta, phi) -
                  VelocityMagnitude(multi_rank_after, radial, theta, phi)));
        }
      }
    }

    check.parity_within_tolerance =
        check.max_rho_difference <= parity_tolerance &&
        check.max_hydro_only_te_difference <= parity_tolerance &&
        check.max_velocity_magnitude_difference <= parity_tolerance;
    check.success = check.parity_within_tolerance;
    if (!check.success) {
      check.failure_reason = "directional mpi parity exceeded tolerance";
    }
  }

  std::ostringstream report;
  report << "directional_mpi_parity_success=" << (check.success ? "true" : "false")
         << "; rank_count=" << check.rank_count
         << "; parity_within_tolerance=" << (check.parity_within_tolerance ? "true" : "false")
         << "; macro_zoning_executed=" << (check.macro_zoning_executed ? "true" : "false")
         << "; max_rho_difference=" << check.max_rho_difference
         << "; max_hydro_only_te_difference=" << check.max_hydro_only_te_difference
         << "; max_velocity_magnitude_difference=" << check.max_velocity_magnitude_difference
         << "; macro_zoning_mode=" << check.macro_zoning_mode;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

bool WriteSphericalRadialWaveMpiParityOutputs(
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::state::CanonicalState& multi_rank_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const std::filesystem::path& output_directory,
    const SphericalRadialWaveMpiParityCheck& summary,
    const HydroBudgetResidualSummary* single_rank_budget,
    const HydroBudgetResidualSummary* multi_rank_budget,
    const char* case_name) noexcept {
  if (!EnsureOutputDirectory(output_directory)) {
    return false;
  }
  if (!geometry.is_valid()) {
    return false;
  }

  const std::size_t shell_index = single_rank_after.layout.radial_cells / 2u;
  const bool have_budgets =
      single_rank_budget != nullptr && single_rank_budget->is_complete() &&
      multi_rank_budget != nullptr && multi_rank_budget->is_complete();
  const std::string case_name_string =
      case_name == nullptr ? "case2_radial_wave_mpi_parity" : case_name;

  return
      WriteRadialProfileFile(
          output_directory / "single_rank_radial_profile_t1.txt",
          "single_rank_t1",
          single_rank_after,
          geometry) &&
      WriteRadialProfileFile(
          output_directory / "multi_rank_radial_profile_t1.txt",
          "multi_rank_t1",
          multi_rank_after,
          geometry) &&
      WriteShellMapFile(output_directory / "rho_shell_map_t1_single.txt", "rho", "t1_single", single_rank_after, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "rho_shell_map_t1_multi.txt", "rho", "t1_multi", multi_rank_after, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "velocity_magnitude_shell_map_t1_single.txt", "velocity_magnitude", "t1_single", single_rank_after, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteShellMapFile(output_directory / "velocity_magnitude_shell_map_t1_multi.txt", "velocity_magnitude", "t1_multi", multi_rank_after, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteKeyValueFile(
          output_directory / "seam_diagnostics.txt",
          {
              {"rank_count", std::to_string(summary.rank_count)},
              {"max_rho_difference", std::to_string(summary.max_rho_difference)},
              {"max_mom_r_difference", std::to_string(summary.max_mom_r_difference)},
              {"max_e_fluid_total_difference", std::to_string(summary.max_e_fluid_total_difference)},
              {"max_seam_rho_jump_mismatch", std::to_string(summary.max_seam_rho_jump)},
              {"max_seam_velocity_jump_mismatch", std::to_string(summary.max_seam_velocity_jump)},
          },
          "seam_diagnostics") &&
      (!have_budgets ||
       (WriteBudgetResidualFile(
            output_directory / "single_rank_budget_residual.txt",
            "single_rank",
            *single_rank_budget) &&
        WriteBudgetResidualFile(
            output_directory / "multi_rank_budget_residual.txt",
            "multi_rank",
            *multi_rank_budget) &&
        WriteBudgetResidualComparisonFile(
            output_directory / "budget_residual_comparison.txt",
            *single_rank_budget,
            *multi_rank_budget))) &&
      WriteKeyValueFile(
          output_directory / "summary.txt",
          {
              {"case", case_name_string},
              {"success", summary.success ? "true" : "false"},
              {"rank_count", std::to_string(summary.rank_count)},
              {"rank_decomposition_valid", summary.rank_decomposition_valid ? "true" : "false"},
              {"parity_within_tolerance", summary.parity_within_tolerance ? "true" : "false"},
              {"seam_jumps_bounded", summary.seam_jumps_bounded ? "true" : "false"},
              {"ppm_executed", summary.ppm_executed ? "true" : "false"},
              {"macro_zoning_executed", summary.macro_zoning_executed ? "true" : "false"},
              {"reconstruction_ghost_layers", std::to_string(summary.reconstruction_ghost_layers)},
              {"reconstruction_mode", summary.reconstruction_mode},
              {"macro_zoning_mode", summary.macro_zoning_mode},
              {"report_line", summary.report_line},
          }) &&
      WriteKeyValueFile(
          output_directory / "case_manifest.txt",
          {
              {"case", case_name_string},
              {"reconstruction_mode", summary.reconstruction_mode},
              {"macro_zoning_mode", summary.macro_zoning_mode},
              {"reconstruction_ghost_layers", std::to_string(summary.reconstruction_ghost_layers)},
              {"summary", "summary.txt"},
              {"single_rank_profile", "single_rank_radial_profile_t1.txt"},
              {"multi_rank_profile", "multi_rank_radial_profile_t1.txt"},
              {"rho_shell_map_single", "rho_shell_map_t1_single.txt"},
              {"rho_shell_map_multi", "rho_shell_map_t1_multi.txt"},
              {"velocity_shell_map_single", "velocity_magnitude_shell_map_t1_single.txt"},
              {"velocity_shell_map_multi", "velocity_magnitude_shell_map_t1_multi.txt"},
              {"seam_diagnostics", "seam_diagnostics.txt"},
              {"single_rank_budget", have_budgets ? "single_rank_budget_residual.txt" : ""},
              {"multi_rank_budget", have_budgets ? "multi_rank_budget_residual.txt" : ""},
              {"budget_comparison", have_budgets ? "budget_residual_comparison.txt" : ""},
          },
          "case_manifest");
}

bool WriteRadialAleOnOffOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& ale_off_after,
    const dec3d::state::CanonicalState& ale_on_after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::MeshUpdateProposal& proposal,
    const std::filesystem::path& output_directory,
    const RadialAleOnOffCheck& summary,
    const HydroBudgetResidualSummary* ale_off_budget,
    const HydroBudgetResidualSummary* ale_on_budget,
    const char* case_name) noexcept {
  if (!EnsureOutputDirectory(output_directory) || !geometry.is_valid()) {
    return false;
  }

  const bool have_budgets =
      ale_off_budget != nullptr && ale_off_budget->is_complete() &&
      ale_on_budget != nullptr && ale_on_budget->is_complete();
  const std::string case_name_string =
      case_name == nullptr ? "case2_radial_wave_ale_on_off" : case_name;

  return
      WriteRadialProfileFile(
          output_directory / "before_radial_profile_t0.txt",
          "before_t0",
          before,
          geometry) &&
      WriteRadialProfileFile(
          output_directory / "ale_off_radial_profile_t1.txt",
          "ale_off_t1",
          ale_off_after,
          geometry) &&
      WriteRadialProfileFile(
          output_directory / "ale_on_radial_profile_t1.txt",
          "ale_on_t1",
          ale_on_after,
          geometry) &&
      WriteKeyValueFile(
          output_directory / "moving_mesh_diagnostics.txt",
          {
              {"proposal_complete", summary.proposal_complete ? "true" : "false"},
              {"moving_mesh_applied", summary.moving_mesh_applied ? "true" : "false"},
              {"outer_face_old", geometry.radial_faces.empty() ? "" : std::to_string(geometry.radial_faces.back())},
              {"outer_face_new", proposal.proposed_radial_faces.empty() ? "" : std::to_string(proposal.proposed_radial_faces.back())},
              {"outer_face_displacement", std::to_string(summary.outer_face_displacement)},
              {"max_face_speed", std::to_string(summary.max_face_speed)},
              {"proposal_summary", proposal.summary},
          },
          "moving_mesh_diagnostics") &&
      (!have_budgets ||
       (WriteBudgetResidualFile(
            output_directory / "budget_residual_off.txt",
            "ale_off",
            *ale_off_budget) &&
        WriteBudgetResidualFile(
            output_directory / "budget_residual_on.txt",
            "ale_on",
            *ale_on_budget) &&
        WriteBudgetResidualComparisonFile(
            output_directory / "budget_residual_comparison.txt",
            *ale_off_budget,
            *ale_on_budget))) &&
      WriteKeyValueFile(
          output_directory / "summary.txt",
          {
              {"case", case_name_string},
              {"success", summary.success ? "true" : "false"},
              {"proposal_complete", summary.proposal_complete ? "true" : "false"},
              {"moving_mesh_applied", summary.moving_mesh_applied ? "true" : "false"},
              {"ale_changes_solution", summary.ale_changes_solution ? "true" : "false"},
              {"ppm_executed", summary.ppm_executed ? "true" : "false"},
              {"reconstruction_ghost_layers", std::to_string(summary.reconstruction_ghost_layers)},
              {"reconstruction_mode", summary.reconstruction_mode},
              {"max_rho_difference", std::to_string(summary.max_rho_difference)},
              {"max_mom_r_difference", std::to_string(summary.max_mom_r_difference)},
              {"max_e_fluid_total_difference", std::to_string(summary.max_e_fluid_total_difference)},
              {"outer_face_displacement", std::to_string(summary.outer_face_displacement)},
              {"max_face_speed", std::to_string(summary.max_face_speed)},
              {"report_line", summary.report_line},
          }) &&
      WriteKeyValueFile(
          output_directory / "case_manifest.txt",
          {
              {"case", case_name_string},
              {"reconstruction_mode", summary.reconstruction_mode},
              {"reconstruction_ghost_layers", std::to_string(summary.reconstruction_ghost_layers)},
              {"summary", "summary.txt"},
              {"moving_mesh_diagnostics", "moving_mesh_diagnostics.txt"},
              {"before_profile", "before_radial_profile_t0.txt"},
              {"ale_off_profile", "ale_off_radial_profile_t1.txt"},
              {"ale_on_profile", "ale_on_radial_profile_t1.txt"},
              {"budget_residual_off", have_budgets ? "budget_residual_off.txt" : ""},
              {"budget_residual_on", have_budgets ? "budget_residual_on.txt" : ""},
              {"budget_comparison", have_budgets ? "budget_residual_comparison.txt" : ""},
          },
          "case_manifest");
}

bool True3DLowModeSanityCheck::is_complete() const noexcept {
  return !report_line.empty();
}

True3DLowModeSanityCheck EvaluateTrue3DLowModeSanity(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after_full_directional,
    const dec3d::state::CanonicalState& after_radial_only,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    double difference_tolerance,
    double mode_tolerance) noexcept {
  True3DLowModeSanityCheck check;
  check.shell_index = shell_index;

  if (!HasAllocatedHydroState(before) ||
      !HasAllocatedHydroState(after_full_directional) ||
      !HasAllocatedHydroState(after_radial_only) ||
      !HasAllocatedElectronState(before) ||
      !HasAllocatedElectronState(after_full_directional) ||
      !HasAllocatedElectronState(after_radial_only)) {
    check.failure_reason = "true-3d low-mode sanity requires allocated hydro and electron state";
  } else if (!LayoutsMatch(before, after_full_directional) ||
             !LayoutsMatch(before, after_radial_only) ||
             !LayoutAndGeometryMatch(before, geometry)) {
    check.failure_reason = "true-3d low-mode sanity requires matching layouts and geometry";
  } else if (shell_index >= before.layout.radial_cells) {
    check.failure_reason = "true-3d low-mode sanity shell index is out of range";
  } else {
    check.finite_states =
        StateIsFinite(before) &&
        StateIsFinite(after_full_directional) &&
        StateIsFinite(after_radial_only);

    for (std::size_t theta = 0; theta < before.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < before.layout.phi_cells; ++phi) {
        check.rho_l1_difference += std::abs(
            after_full_directional.rho(shell_index, theta, phi) -
            after_radial_only.rho(shell_index, theta, phi));
        check.te_l1_difference += std::abs(
            HydroOnlyElectronTemperature(after_full_directional, shell_index, theta, phi) -
            HydroOnlyElectronTemperature(after_radial_only, shell_index, theta, phi));
        check.velocity_l1_difference += std::abs(
            VelocityMagnitude(after_full_directional, shell_index, theta, phi) -
            VelocityMagnitude(after_radial_only, shell_index, theta, phi));
      }
    }

    const auto before_projection =
        ComputeVrLowModeProjection(before, geometry, shell_index);
    const auto full_projection =
        ComputeVrLowModeProjection(after_full_directional, geometry, shell_index);
    const auto radial_projection =
        ComputeVrLowModeProjection(after_radial_only, geometry, shell_index);

    check.mode_amplitude_before = before_projection.amplitude;
    check.mode_amplitude_after_full = full_projection.amplitude;
    check.mode_amplitude_after_radial_only = radial_projection.amplitude;
    check.mode_phase_before = before_projection.phase;
    check.mode_phase_after_full = full_projection.phase;
    check.mode_phase_after_radial_only = radial_projection.phase;

    check.full_directional_differs_from_radial_only =
        check.rho_l1_difference > difference_tolerance ||
        check.te_l1_difference > difference_tolerance ||
        check.velocity_l1_difference > difference_tolerance;
    check.mode_projection_nontrivial =
        std::abs(check.mode_amplitude_after_full - check.mode_amplitude_after_radial_only) > mode_tolerance ||
        check.mode_amplitude_after_full > mode_tolerance;
    check.success =
        check.finite_states &&
        check.full_directional_differs_from_radial_only &&
        check.mode_projection_nontrivial;

    if (!check.success) {
      if (!check.finite_states) {
        check.failure_reason = "true-3d low-mode sanity encountered non-finite state values";
      } else if (!check.full_directional_differs_from_radial_only) {
        check.failure_reason = "full directional low-mode case does not differ enough from radial-only reference";
      } else if (!check.mode_projection_nontrivial) {
        check.failure_reason = "low-mode projection remained trivial after the true-3d update";
      }
    }
  }

  std::ostringstream report;
  report << "true3d_low_mode_success=" << (check.success ? "true" : "false")
         << "; finite_states=" << (check.finite_states ? "true" : "false")
         << "; full_directional_differs_from_radial_only="
         << (check.full_directional_differs_from_radial_only ? "true" : "false")
         << "; mode_projection_nontrivial="
         << (check.mode_projection_nontrivial ? "true" : "false")
         << "; shell_index=" << check.shell_index
         << "; rho_l1_difference=" << check.rho_l1_difference
         << "; te_l1_difference=" << check.te_l1_difference
         << "; velocity_l1_difference=" << check.velocity_l1_difference
         << "; mode_amplitude_before=" << check.mode_amplitude_before
         << "; mode_amplitude_after_full=" << check.mode_amplitude_after_full
         << "; mode_amplitude_after_radial_only=" << check.mode_amplitude_after_radial_only
         << "; mode_phase_before=" << check.mode_phase_before
         << "; mode_phase_after_full=" << check.mode_phase_after_full
         << "; mode_phase_after_radial_only=" << check.mode_phase_after_radial_only;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

bool WriteTrue3DLowModeOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after_full_directional,
    const dec3d::state::CanonicalState& after_radial_only,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    const std::filesystem::path& output_directory,
    const True3DLowModeSanityCheck& summary,
    const HydroBudgetResidualSummary* budget,
    const DirectionalMpiParityCheck* mpi_parity,
    const HydroBudgetResidualSummary* single_rank_budget,
    const char* case_name) noexcept {
  if (!EnsureOutputDirectory(output_directory)) {
    return false;
  }

  const auto before_projection = ComputeVrLowModeProjection(before, geometry, shell_index);
  const auto full_projection = ComputeVrLowModeProjection(after_full_directional, geometry, shell_index);
  const auto radial_projection = ComputeVrLowModeProjection(after_radial_only, geometry, shell_index);
  const bool have_budget = budget != nullptr && budget->is_complete();
  const bool have_mpi_parity = mpi_parity != nullptr && mpi_parity->is_complete();
  const bool have_single_rank_budget =
      single_rank_budget != nullptr && single_rank_budget->is_complete();
  const std::string case_name_string =
      case_name == nullptr ? "case3_true3d_low_mode" : case_name;

  return
      WriteShellMapFile(output_directory / "rho_shell_map_t0.txt", "rho", "t0", before, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "rho_shell_map_t1_full.txt", "rho", "t1_full", after_full_directional, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "hydro_only_te_shell_map_t0.txt", "hydro_only_te", "t0", before, geometry, shell_index, &AccessHydroOnlyTe) &&
      WriteShellMapFile(output_directory / "hydro_only_te_shell_map_t1_full.txt", "hydro_only_te", "t1_full", after_full_directional, geometry, shell_index, &AccessHydroOnlyTe) &&
      WriteShellMapFile(output_directory / "velocity_magnitude_shell_map_t0.txt", "velocity_magnitude", "t0", before, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteShellMapFile(output_directory / "velocity_magnitude_shell_map_t1_full.txt", "velocity_magnitude", "t1_full", after_full_directional, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteModeProjectionFile(output_directory / "mode_projection.txt", before_projection, full_projection, radial_projection) &&
      WriteKeyValueFile(
          output_directory / "radial_only_comparison.txt",
          {
              {"rho_l1_difference", std::to_string(summary.rho_l1_difference)},
              {"te_l1_difference", std::to_string(summary.te_l1_difference)},
              {"velocity_l1_difference", std::to_string(summary.velocity_l1_difference)},
              {"mode_amplitude_after_full", std::to_string(summary.mode_amplitude_after_full)},
              {"mode_amplitude_after_radial_only", std::to_string(summary.mode_amplitude_after_radial_only)},
          },
          "radial_only_comparison") &&
      (!have_budget ||
       WriteBudgetResidualFile(
           output_directory / "budget_residual.txt",
           "full_directional",
           *budget)) &&
      (!have_mpi_parity ||
       WriteKeyValueFile(
           output_directory / "mpi_parity_diagnostics.txt",
           {
                {"rank_count", std::to_string(mpi_parity->rank_count)},
                {"parity_within_tolerance", mpi_parity->parity_within_tolerance ? "true" : "false"},
                {"ppm_executed", mpi_parity->ppm_executed ? "true" : "false"},
                {"macro_zoning_executed", mpi_parity->macro_zoning_executed ? "true" : "false"},
                {"reconstruction_ghost_layers", std::to_string(mpi_parity->reconstruction_ghost_layers)},
                {"reconstruction_mode", mpi_parity->reconstruction_mode},
                {"macro_zoning_mode", mpi_parity->macro_zoning_mode},
                {"max_rho_difference", std::to_string(mpi_parity->max_rho_difference)},
                {"max_hydro_only_te_difference", std::to_string(mpi_parity->max_hydro_only_te_difference)},
                {"max_velocity_magnitude_difference", std::to_string(mpi_parity->max_velocity_magnitude_difference)},
           },
           "mpi_parity_diagnostics")) &&
      (!(have_budget && have_single_rank_budget) ||
       WriteBudgetResidualComparisonFile(
           output_directory / "budget_residual_comparison.txt",
           *single_rank_budget,
           *budget)) &&
      WriteKeyValueFile(
          output_directory / "summary.txt",
          {
              {"case", case_name_string},
              {"success", summary.success ? "true" : "false"},
              {"finite_states", summary.finite_states ? "true" : "false"},
              {"full_directional_differs_from_radial_only", summary.full_directional_differs_from_radial_only ? "true" : "false"},
              {"mode_projection_nontrivial", summary.mode_projection_nontrivial ? "true" : "false"},
              {"shell_index", std::to_string(summary.shell_index)},
              {"rho_l1_difference", std::to_string(summary.rho_l1_difference)},
              {"te_l1_difference", std::to_string(summary.te_l1_difference)},
               {"velocity_l1_difference", std::to_string(summary.velocity_l1_difference)},
               {"mpi_parity_within_tolerance", have_mpi_parity ? (mpi_parity->parity_within_tolerance ? "true" : "false") : ""},
               {"ppm_executed", have_mpi_parity ? (mpi_parity->ppm_executed ? "true" : "false") : ""},
               {"macro_zoning_executed", have_mpi_parity ? (mpi_parity->macro_zoning_executed ? "true" : "false") : ""},
               {"reconstruction_ghost_layers", have_mpi_parity ? std::to_string(mpi_parity->reconstruction_ghost_layers) : ""},
               {"reconstruction_mode", have_mpi_parity ? mpi_parity->reconstruction_mode : ""},
               {"macro_zoning_mode", have_mpi_parity ? mpi_parity->macro_zoning_mode : ""},
               {"report_line", summary.report_line},
           }) &&
      WriteKeyValueFile(
          output_directory / "case_manifest.txt",
          {
               {"case", case_name_string},
               {"reconstruction_mode", have_mpi_parity ? mpi_parity->reconstruction_mode : ""},
               {"macro_zoning_mode", have_mpi_parity ? mpi_parity->macro_zoning_mode : ""},
               {"reconstruction_ghost_layers", have_mpi_parity ? std::to_string(mpi_parity->reconstruction_ghost_layers) : ""},
               {"summary", "summary.txt"},
              {"rho_t0", "rho_shell_map_t0.txt"},
              {"rho_t1_full", "rho_shell_map_t1_full.txt"},
              {"te_t0", "hydro_only_te_shell_map_t0.txt"},
              {"te_t1_full", "hydro_only_te_shell_map_t1_full.txt"},
              {"velocity_t0", "velocity_magnitude_shell_map_t0.txt"},
              {"velocity_t1_full", "velocity_magnitude_shell_map_t1_full.txt"},
              {"mode_projection", "mode_projection.txt"},
              {"radial_only_comparison", "radial_only_comparison.txt"},
              {"budget_residual", have_budget ? "budget_residual.txt" : ""},
              {"mpi_parity_diagnostics", have_mpi_parity ? "mpi_parity_diagnostics.txt" : ""},
              {"budget_comparison", (have_budget && have_single_rank_budget) ? "budget_residual_comparison.txt" : ""},
          },
          "case_manifest");
}

bool PoleAdjacentSeamStressCheck::is_complete() const noexcept {
  return !report_line.empty();
}

PoleAdjacentSeamStressCheck EvaluatePoleAdjacentSeamStress(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    double seam_growth_tolerance,
    double polar_velocity_tolerance) noexcept {
  PoleAdjacentSeamStressCheck check;
  check.shell_index = shell_index;

  if (!HasAllocatedHydroState(before) ||
      !HasAllocatedHydroState(after) ||
      !HasAllocatedElectronState(before) ||
      !HasAllocatedElectronState(after)) {
    check.failure_reason = "pole-adjacent seam stress requires allocated hydro and electron state";
  } else if (!LayoutsMatch(before, after) || !LayoutAndGeometryMatch(before, geometry)) {
    check.failure_reason = "pole-adjacent seam stress requires matching layouts and geometry";
  } else if (shell_index >= before.layout.radial_cells) {
    check.failure_reason = "pole-adjacent seam stress shell index is out of range";
  } else {
    check.finite_solution = StateIsFinite(after);

    for (std::size_t theta = 0; theta < before.layout.theta_cells; ++theta) {
      const double rho_jump_before =
          std::abs(before.rho(shell_index, theta, before.layout.phi_cells - 1u) -
                   before.rho(shell_index, theta, 0u));
      const double rho_jump_after =
          std::abs(after.rho(shell_index, theta, after.layout.phi_cells - 1u) -
                   after.rho(shell_index, theta, 0u));
      const double mom_theta_jump_before =
          std::abs(before.mom_theta(shell_index, theta, before.layout.phi_cells - 1u) -
                   before.mom_theta(shell_index, theta, 0u));
      const double mom_theta_jump_after =
          std::abs(after.mom_theta(shell_index, theta, after.layout.phi_cells - 1u) -
                   after.mom_theta(shell_index, theta, 0u));
      const double mom_phi_jump_before =
          std::abs(before.mom_phi(shell_index, theta, before.layout.phi_cells - 1u) -
                   before.mom_phi(shell_index, theta, 0u));
      const double mom_phi_jump_after =
          std::abs(after.mom_phi(shell_index, theta, after.layout.phi_cells - 1u) -
                   after.mom_phi(shell_index, theta, 0u));
      const double velocity_jump_before =
          std::abs(VelocityMagnitude(before, shell_index, theta, before.layout.phi_cells - 1u) -
                   VelocityMagnitude(before, shell_index, theta, 0u));
      const double velocity_jump_after =
          std::abs(VelocityMagnitude(after, shell_index, theta, after.layout.phi_cells - 1u) -
                   VelocityMagnitude(after, shell_index, theta, 0u));

      check.max_seam_rho_jump_before = std::max(check.max_seam_rho_jump_before, rho_jump_before);
      check.max_seam_rho_jump_after = std::max(check.max_seam_rho_jump_after, rho_jump_after);
      check.max_seam_mom_theta_jump_before =
          std::max(check.max_seam_mom_theta_jump_before, mom_theta_jump_before);
      check.max_seam_mom_theta_jump_after =
          std::max(check.max_seam_mom_theta_jump_after, mom_theta_jump_after);
      check.max_seam_mom_phi_jump_before =
          std::max(check.max_seam_mom_phi_jump_before, mom_phi_jump_before);
      check.max_seam_mom_phi_jump_after =
          std::max(check.max_seam_mom_phi_jump_after, mom_phi_jump_after);
      check.max_seam_velocity_jump_before = std::max(check.max_seam_velocity_jump_before, velocity_jump_before);
      check.max_seam_velocity_jump_after = std::max(check.max_seam_velocity_jump_after, velocity_jump_after);
    }

    const std::size_t cap_depth = std::min<std::size_t>(2u, before.layout.theta_cells / 2u);
    for (std::size_t theta = 0; theta < cap_depth; ++theta) {
      for (std::size_t phi = 0; phi < before.layout.phi_cells; ++phi) {
        check.max_north_cap_velocity = std::max(
            check.max_north_cap_velocity,
            VelocityMagnitude(after, shell_index, theta, phi));
        check.max_south_cap_velocity = std::max(
            check.max_south_cap_velocity,
            VelocityMagnitude(after, shell_index, after.layout.theta_cells - 1u - theta, phi));
      }
    }

    check.seam_growth_bounded =
        check.max_seam_rho_jump_after <= check.max_seam_rho_jump_before + seam_growth_tolerance &&
        check.max_seam_velocity_jump_after <= check.max_seam_velocity_jump_before + seam_growth_tolerance;
    check.polar_caps_stable =
        check.max_north_cap_velocity <= polar_velocity_tolerance &&
        check.max_south_cap_velocity <= polar_velocity_tolerance;
    check.success =
        check.finite_solution &&
        check.seam_growth_bounded &&
        check.polar_caps_stable;

    if (!check.success) {
      if (!check.finite_solution) {
        check.failure_reason = "pole-adjacent seam stress produced non-finite state values";
      } else if (!check.seam_growth_bounded) {
        check.failure_reason = "seam jump growth exceeded tolerance";
      } else if (!check.polar_caps_stable) {
        check.failure_reason = "polar cap velocity exceeded tolerance";
      }
    }
  }

  std::ostringstream report;
  report << "pole_adjacent_seam_stress_success=" << (check.success ? "true" : "false")
         << "; finite_solution=" << (check.finite_solution ? "true" : "false")
         << "; seam_growth_bounded=" << (check.seam_growth_bounded ? "true" : "false")
         << "; polar_caps_stable=" << (check.polar_caps_stable ? "true" : "false")
         << "; shell_index=" << check.shell_index
         << "; max_seam_rho_jump_before=" << check.max_seam_rho_jump_before
         << "; max_seam_rho_jump_after=" << check.max_seam_rho_jump_after
         << "; max_seam_mom_theta_jump_before=" << check.max_seam_mom_theta_jump_before
         << "; max_seam_mom_theta_jump_after=" << check.max_seam_mom_theta_jump_after
         << "; max_seam_mom_phi_jump_before=" << check.max_seam_mom_phi_jump_before
         << "; max_seam_mom_phi_jump_after=" << check.max_seam_mom_phi_jump_after
         << "; max_seam_velocity_jump_before=" << check.max_seam_velocity_jump_before
         << "; max_seam_velocity_jump_after=" << check.max_seam_velocity_jump_after
         << "; max_north_cap_velocity=" << check.max_north_cap_velocity
         << "; max_south_cap_velocity=" << check.max_south_cap_velocity;
  if (!check.failure_reason.empty()) {
    report << "; failure_reason=" << check.failure_reason;
  }
  check.report_line = report.str();

  return check;
}

bool WritePoleAdjacentSeamStressOutputs(
    const dec3d::state::CanonicalState& before,
    const dec3d::state::CanonicalState& after,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t shell_index,
    const std::filesystem::path& output_directory,
    const PoleAdjacentSeamStressCheck& summary,
    const HydroBudgetResidualSummary* budget,
    const DirectionalMpiParityCheck* mpi_parity,
    const HydroBudgetResidualSummary* single_rank_budget,
    const char* case_name) noexcept {
  if (!EnsureOutputDirectory(output_directory)) {
    return false;
  }
  const bool have_budget = budget != nullptr && budget->is_complete();
  const bool have_mpi_parity = mpi_parity != nullptr && mpi_parity->is_complete();
  const bool have_single_rank_budget =
      single_rank_budget != nullptr && single_rank_budget->is_complete();
  const std::string case_name_string =
      case_name == nullptr ? "case4_pole_seam_stress" : case_name;

  return
      WritePolarSliceFile(output_directory / "north_polar_slice_t0.txt", "north", "t0", before, geometry, shell_index, 0u) &&
      WritePolarSliceFile(output_directory / "north_polar_slice_t1.txt", "north", "t1", after, geometry, shell_index, 0u) &&
      WritePolarSliceFile(output_directory / "south_polar_slice_t0.txt", "south", "t0", before, geometry, shell_index, before.layout.theta_cells - 1u) &&
      WritePolarSliceFile(output_directory / "south_polar_slice_t1.txt", "south", "t1", after, geometry, shell_index, after.layout.theta_cells - 1u) &&
      WriteShellMapFile(output_directory / "rho_shell_map_t1.txt", "rho", "t1", after, geometry, shell_index, &AccessRho) &&
      WriteShellMapFile(output_directory / "velocity_magnitude_shell_map_t1.txt", "velocity_magnitude", "t1", after, geometry, shell_index, &AccessVelocityMagnitude) &&
      WriteKeyValueFile(
          output_directory / "seam_polar_diagnostics.txt",
          {
              {"max_seam_rho_jump_before", std::to_string(summary.max_seam_rho_jump_before)},
              {"max_seam_rho_jump_after", std::to_string(summary.max_seam_rho_jump_after)},
              {"max_seam_mom_theta_jump_before", std::to_string(summary.max_seam_mom_theta_jump_before)},
              {"max_seam_mom_theta_jump_after", std::to_string(summary.max_seam_mom_theta_jump_after)},
              {"max_seam_mom_phi_jump_before", std::to_string(summary.max_seam_mom_phi_jump_before)},
              {"max_seam_mom_phi_jump_after", std::to_string(summary.max_seam_mom_phi_jump_after)},
              {"max_seam_velocity_jump_before", std::to_string(summary.max_seam_velocity_jump_before)},
              {"max_seam_velocity_jump_after", std::to_string(summary.max_seam_velocity_jump_after)},
              {"max_north_cap_velocity", std::to_string(summary.max_north_cap_velocity)},
              {"max_south_cap_velocity", std::to_string(summary.max_south_cap_velocity)},
          },
          "seam_diagnostics") &&
      (!have_budget ||
       WriteBudgetResidualFile(
           output_directory / "budget_residual.txt",
           "pole_seam_stress",
           *budget)) &&
      (!have_mpi_parity ||
       WriteKeyValueFile(
           output_directory / "mpi_parity_diagnostics.txt",
           {
                {"rank_count", std::to_string(mpi_parity->rank_count)},
                {"parity_within_tolerance", mpi_parity->parity_within_tolerance ? "true" : "false"},
                {"ppm_executed", mpi_parity->ppm_executed ? "true" : "false"},
                {"macro_zoning_executed", mpi_parity->macro_zoning_executed ? "true" : "false"},
                {"reconstruction_ghost_layers", std::to_string(mpi_parity->reconstruction_ghost_layers)},
                {"reconstruction_mode", mpi_parity->reconstruction_mode},
                {"macro_zoning_mode", mpi_parity->macro_zoning_mode},
                {"max_rho_difference", std::to_string(mpi_parity->max_rho_difference)},
                {"max_hydro_only_te_difference", std::to_string(mpi_parity->max_hydro_only_te_difference)},
                {"max_velocity_magnitude_difference", std::to_string(mpi_parity->max_velocity_magnitude_difference)},
           },
           "mpi_parity_diagnostics")) &&
      (!(have_budget && have_single_rank_budget) ||
       WriteBudgetResidualComparisonFile(
           output_directory / "budget_residual_comparison.txt",
           *single_rank_budget,
           *budget)) &&
      WriteKeyValueFile(
          output_directory / "summary.txt",
          {
              {"case", case_name_string},
              {"success", summary.success ? "true" : "false"},
              {"finite_solution", summary.finite_solution ? "true" : "false"},
              {"seam_growth_bounded", summary.seam_growth_bounded ? "true" : "false"},
              {"polar_caps_stable", summary.polar_caps_stable ? "true" : "false"},
               {"shell_index", std::to_string(summary.shell_index)},
               {"mpi_parity_within_tolerance", have_mpi_parity ? (mpi_parity->parity_within_tolerance ? "true" : "false") : ""},
               {"ppm_executed", have_mpi_parity ? (mpi_parity->ppm_executed ? "true" : "false") : ""},
               {"macro_zoning_executed", have_mpi_parity ? (mpi_parity->macro_zoning_executed ? "true" : "false") : ""},
               {"reconstruction_ghost_layers", have_mpi_parity ? std::to_string(mpi_parity->reconstruction_ghost_layers) : ""},
               {"reconstruction_mode", have_mpi_parity ? mpi_parity->reconstruction_mode : ""},
               {"macro_zoning_mode", have_mpi_parity ? mpi_parity->macro_zoning_mode : ""},
               {"report_line", summary.report_line},
           }) &&
      WriteKeyValueFile(
          output_directory / "case_manifest.txt",
          {
               {"case", case_name_string},
               {"reconstruction_mode", have_mpi_parity ? mpi_parity->reconstruction_mode : ""},
               {"macro_zoning_mode", have_mpi_parity ? mpi_parity->macro_zoning_mode : ""},
               {"reconstruction_ghost_layers", have_mpi_parity ? std::to_string(mpi_parity->reconstruction_ghost_layers) : ""},
              {"summary", "summary.txt"},
              {"north_t0", "north_polar_slice_t0.txt"},
              {"north_t1", "north_polar_slice_t1.txt"},
              {"south_t0", "south_polar_slice_t0.txt"},
              {"south_t1", "south_polar_slice_t1.txt"},
              {"rho_shell_t1", "rho_shell_map_t1.txt"},
              {"velocity_shell_t1", "velocity_magnitude_shell_map_t1.txt"},
              {"seam_polar_diagnostics", "seam_polar_diagnostics.txt"},
              {"budget_residual", have_budget ? "budget_residual.txt" : ""},
              {"mpi_parity_diagnostics", have_mpi_parity ? "mpi_parity_diagnostics.txt" : ""},
              {"budget_comparison", (have_budget && have_single_rank_budget) ? "budget_residual_comparison.txt" : ""},
          },
          "case_manifest");
}

}  // namespace dec3d::hydro
