#pragma once

#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace dec3d::testsupport {

using dec3d::hydro::HydroConservativeState;

[[nodiscard]] inline std::string ReadTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::string content;
  std::string line;
  while (std::getline(input, line)) {
    content += line;
    content.push_back('\n');
  }

  return content;
}

[[nodiscard]] inline bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }

  return false;
}

[[nodiscard]] inline dec3d::mesh::SphericalGeometryMetadata BuildLocalGeometrySlice(
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

[[nodiscard]] inline dec3d::state::CanonicalState BuildLocalStateSlice(
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
        for (std::size_t group = 0; group < global_state.radiation_groups.size(); ++group) {
          local_state.radiation_groups[group](local_radial, theta, phi) =
              global_state.radiation_groups[group](global_radial, theta, phi);
        }
      }
    }
  }

  return local_state;
}

[[nodiscard]] inline dec3d::hydro::RadialGhostOverride BuildRadialGhostOverride(
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

inline void GatherFieldToRoot(
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

inline void ReduceBudgetSummaryToRoot(
    const dec3d::hydro::HydroBudgetResidualSummary& local_budget,
    int rank,
    dec3d::hydro::HydroBudgetResidualSummary* reduced_budget) {
  const double local_budget_values[] = {
      local_budget.old_mass,
      local_budget.flux_mass_delta,
      local_budget.source_mass_delta,
      local_budget.new_mass,
      local_budget.mass_residual,
      local_budget.old_mom_r,
      local_budget.flux_mom_r_delta,
      local_budget.source_mom_r_delta,
      local_budget.new_mom_r,
      local_budget.mom_r_residual,
      local_budget.old_e_fluid_total,
      local_budget.flux_e_fluid_total_delta,
      local_budget.source_e_fluid_total_delta,
      local_budget.new_e_fluid_total,
      local_budget.e_fluid_total_residual,
  };
  double reduced_budget_values[15] = {};
  MPI_Reduce(
      local_budget_values,
      reduced_budget_values,
      15,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      MPI_COMM_WORLD);

  if (rank != 0 || reduced_budget == nullptr) {
    return;
  }

  reduced_budget->old_mass = reduced_budget_values[0];
  reduced_budget->flux_mass_delta = reduced_budget_values[1];
  reduced_budget->source_mass_delta = reduced_budget_values[2];
  reduced_budget->new_mass = reduced_budget_values[3];
  reduced_budget->mass_residual = reduced_budget_values[4];
  reduced_budget->old_mom_r = reduced_budget_values[5];
  reduced_budget->flux_mom_r_delta = reduced_budget_values[6];
  reduced_budget->source_mom_r_delta = reduced_budget_values[7];
  reduced_budget->new_mom_r = reduced_budget_values[8];
  reduced_budget->mom_r_residual = reduced_budget_values[9];
  reduced_budget->old_e_fluid_total = reduced_budget_values[10];
  reduced_budget->flux_e_fluid_total_delta = reduced_budget_values[11];
  reduced_budget->source_e_fluid_total_delta = reduced_budget_values[12];
  reduced_budget->new_e_fluid_total = reduced_budget_values[13];
  reduced_budget->e_fluid_total_residual = reduced_budget_values[14];
  reduced_budget->report_line =
      "mpi_aggregated_budget_mass_residual=" + std::to_string(reduced_budget->mass_residual) +
      "; mpi_aggregated_budget_mom_r_residual=" + std::to_string(reduced_budget->mom_r_residual) +
      "; mpi_aggregated_budget_e_residual=" + std::to_string(reduced_budget->e_fluid_total_residual);
}

}  // namespace dec3d::testsupport
