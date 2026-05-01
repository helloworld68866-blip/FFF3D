#pragma once

#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "hydro_mpi_real_case_support.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace dec3d::testsupport {

inline void SeedMacroZoningRadialWave(
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
    const auto& base = inside ? inner_conservative : outer_conservative;
    const double electron_pressure = inside ? 0.4 : 0.04;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      const double theta_mode = 1.0 + 0.01 * std::cos(0.7 * static_cast<double>(theta));
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const double phi_mode = 1.0 + 0.01 * std::sin(0.5 * static_cast<double>(phi));
        const double modifier = theta_mode * phi_mode;

        state.rho(radial, theta, phi) = base.rho * modifier;
        state.mom_r(radial, theta, phi) = base.mom_r;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = base.e_fluid_total * modifier;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure * modifier);

        before.rho(radial, theta, phi) = state.rho(radial, theta, phi);
        before.mom_r(radial, theta, phi) = state.mom_r(radial, theta, phi);
        before.mom_theta(radial, theta, phi) = state.mom_theta(radial, theta, phi);
        before.mom_phi(radial, theta, phi) = state.mom_phi(radial, theta, phi);
        before.e_fluid_total(radial, theta, phi) = state.e_fluid_total(radial, theta, phi);
        before.e_electron(radial, theta, phi) = state.e_electron(radial, theta, phi);
      }
    }
  }
}

inline dec3d::hydro::StaticGridHydroOptions MacroZoningActiveOptions() noexcept {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = false;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_macro_zoning = true;
  options.macro_zoning_coarse_factor = 0.5;
  return options;
}

inline double VelocityMagnitudeAt(
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

inline dec3d::hydro::SphericalRadialWaveMpiParityCheck BuildMacroZoningParitySummary(
    const dec3d::state::CanonicalState& single_rank_after,
    const dec3d::state::CanonicalState& multi_rank_after,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    double parity_tolerance,
    double seam_tolerance) {
  dec3d::hydro::SphericalRadialWaveMpiParityCheck summary;
  summary.rank_count = decomposition.slices.size();
  summary.rank_decomposition_valid = decomposition.is_valid();
  summary.macro_zoning_executed = true;
  summary.reconstruction_mode = "macro_zoning_mpi_real_first_order_angular";
  summary.macro_zoning_mode = "mpi_real_first_order_angular_macro";

  for (std::size_t radial = 0; radial < single_rank_after.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < single_rank_after.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < single_rank_after.layout.phi_cells; ++phi) {
        summary.max_rho_difference = std::max(
            summary.max_rho_difference,
            std::abs(single_rank_after.rho(radial, theta, phi) -
                     multi_rank_after.rho(radial, theta, phi)));
        summary.max_mom_r_difference = std::max(
            summary.max_mom_r_difference,
            std::abs(single_rank_after.mom_r(radial, theta, phi) -
                     multi_rank_after.mom_r(radial, theta, phi)));
        summary.max_e_fluid_total_difference = std::max(
            summary.max_e_fluid_total_difference,
            std::abs(single_rank_after.e_fluid_total(radial, theta, phi) -
                     multi_rank_after.e_fluid_total(radial, theta, phi)));
      }
    }
  }

  for (std::size_t slice_index = 0; slice_index + 1u < decomposition.slices.size(); ++slice_index) {
    const auto& left_slice = decomposition.slices[slice_index];
    const std::size_t left_radial = left_slice.end_index - 1u;
    const std::size_t right_radial = left_slice.end_index;
    for (std::size_t theta = 0; theta < single_rank_after.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < single_rank_after.layout.phi_cells; ++phi) {
        const double single_rho_jump =
            std::abs(single_rank_after.rho(left_radial, theta, phi) -
                     single_rank_after.rho(right_radial, theta, phi));
        const double multi_rho_jump =
            std::abs(multi_rank_after.rho(left_radial, theta, phi) -
                     multi_rank_after.rho(right_radial, theta, phi));
        summary.max_seam_rho_jump = std::max(
            summary.max_seam_rho_jump,
            std::abs(multi_rho_jump - single_rho_jump));

        const double single_velocity_jump =
            std::abs(VelocityMagnitudeAt(single_rank_after, left_radial, theta, phi) -
                     VelocityMagnitudeAt(single_rank_after, right_radial, theta, phi));
        const double multi_velocity_jump =
            std::abs(VelocityMagnitudeAt(multi_rank_after, left_radial, theta, phi) -
                     VelocityMagnitudeAt(multi_rank_after, right_radial, theta, phi));
        summary.max_seam_velocity_jump = std::max(
            summary.max_seam_velocity_jump,
            std::abs(multi_velocity_jump - single_velocity_jump));
      }
    }
  }

  summary.parity_within_tolerance =
      summary.max_rho_difference <= parity_tolerance &&
      summary.max_mom_r_difference <= parity_tolerance &&
      summary.max_e_fluid_total_difference <= parity_tolerance;
  summary.seam_jumps_bounded =
      summary.max_seam_rho_jump <= seam_tolerance &&
      summary.max_seam_velocity_jump <= seam_tolerance;
  summary.success = summary.rank_decomposition_valid &&
                    summary.parity_within_tolerance &&
                    summary.seam_jumps_bounded;

  if (!summary.success) {
    if (!summary.rank_decomposition_valid) {
      summary.failure_reason = "macro-zoning radial ownership decomposition is invalid";
    } else if (!summary.parity_within_tolerance) {
      summary.failure_reason = "macro-zoning MPI parity drift exceeded tolerance";
    } else {
      summary.failure_reason = "macro-zoning seam mismatch exceeded tolerance";
    }
  }

  std::ostringstream report;
  report << "macro_zoning_mpi_parity_success=" << (summary.success ? "true" : "false")
         << "; rank_count=" << summary.rank_count
         << "; rank_decomposition_valid=" << (summary.rank_decomposition_valid ? "true" : "false")
         << "; parity_within_tolerance=" << (summary.parity_within_tolerance ? "true" : "false")
         << "; seam_jumps_bounded=" << (summary.seam_jumps_bounded ? "true" : "false")
         << "; max_rho_difference=" << summary.max_rho_difference
         << "; max_mom_r_difference=" << summary.max_mom_r_difference
         << "; max_e_fluid_total_difference=" << summary.max_e_fluid_total_difference
         << "; max_seam_rho_jump_mismatch=" << summary.max_seam_rho_jump
         << "; max_seam_velocity_jump_mismatch=" << summary.max_seam_velocity_jump
         << "; macro_zoning_executed=true"
         << "; reconstruction_mode=" << summary.reconstruction_mode;
  if (!summary.failure_reason.empty()) {
    report << "; failure_reason=" << summary.failure_reason;
  }
  summary.report_line = report.str();
  return summary;
}

}  // namespace dec3d::testsupport
