#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace {

struct SymmetryMetrics {
  double max_rho_angular_relative_spread{0.0};
  double max_abs_mom_theta{0.0};
  double max_abs_mom_phi{0.0};
  double tangential_to_radial_momentum_ratio{0.0};
};

void SeedNohState(dec3d::state::CanonicalState& state) {
  constexpr double kRho0 = 1.0;
  constexpr double kPressure0 = 1.0e-6;
  constexpr double kRadialVelocity0 = -1.0;
  constexpr double kElectronEnergyFraction = 0.5;
  const double e_internal = kPressure0 / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double e_total =
      e_internal + (0.5 * kRho0 * kRadialVelocity0 * kRadialVelocity0);
  const double e_electron =
      dec3d::state::ElectronEnergyDensityFromPressure(kElectronEnergyFraction * kPressure0);

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = kRho0;
        state.mom_r(radial, theta, phi) = kRho0 * kRadialVelocity0;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = e_total;
        state.e_electron(radial, theta, phi) = e_electron;
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::RadialGhostOverride BuildNohOuterInflowOverride(
    std::size_t theta_cells,
    std::size_t phi_cells,
    std::size_t ghost_layers) {
  dec3d::hydro::RadialGhostOverride radial_override;
  radial_override.has_outer_neighbor = true;
  radial_override.ghost_layers = ghost_layers;
  radial_override.outer_ghost_reuses_boundary_partition = true;
  radial_override.outer_ghost_states =
      dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
          ghost_layers,
          theta_cells,
          phi_cells);

  const dec3d::hydro::HydroPrimitiveState primitive{
      1.0,
      -1.0,
      0.0,
      0.0,
      1.0e-6,
      dec3d::state::ChiEFromElectronPressure(0.5e-6)};
  const auto inflow = dec3d::hydro::MakeConservativeState(primitive);
  for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        radial_override.outer_ghost_states(ghost, theta, phi) = inflow;
      }
    }
  }
  radial_override.report_line = "short_run_noh_outer_inflow=true";
  return radial_override;
}

[[nodiscard]] const dec3d::core::DiagnosticMessage* FindDiagnostic(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return &entry;
    }
  }
  return nullptr;
}

[[nodiscard]] bool HasDiagnosticMessageFragment(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code,
    const char* fragment) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code &&
        entry.message.find(fragment) != std::string::npos) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] SymmetryMetrics MeasureSymmetry(const dec3d::state::CanonicalState& state) {
  SymmetryMetrics metrics;
  double tangential_sum = 0.0;
  double radial_sum = 0.0;

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    double min_rho = state.rho(radial, 0, 0);
    double max_rho = state.rho(radial, 0, 0);
    double rho_sum = 0.0;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const double rho = state.rho(radial, theta, phi);
        min_rho = std::min(min_rho, rho);
        max_rho = std::max(max_rho, rho);
        rho_sum += rho;
        metrics.max_abs_mom_theta =
            std::max(metrics.max_abs_mom_theta, std::abs(state.mom_theta(radial, theta, phi)));
        metrics.max_abs_mom_phi =
            std::max(metrics.max_abs_mom_phi, std::abs(state.mom_phi(radial, theta, phi)));
        tangential_sum +=
            std::abs(state.mom_theta(radial, theta, phi)) +
            std::abs(state.mom_phi(radial, theta, phi));
        radial_sum += std::abs(state.mom_r(radial, theta, phi));
      }
    }
    const double mean_rho =
        rho_sum / static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);
    metrics.max_rho_angular_relative_spread = std::max(
        metrics.max_rho_angular_relative_spread,
        (max_rho - min_rho) / std::abs(mean_rho));
  }

  metrics.tangential_to_radial_momentum_ratio =
      radial_sum > 0.0 ? tangential_sum / radial_sum : 0.0;
  return metrics;
}

[[nodiscard]] SymmetryMetrics RunShortNoh(
    const dec3d::hydro::StaticGridHydroOptions& options,
    std::size_t step_count) {
  constexpr std::size_t kRadialCells = 32u;
  constexpr std::size_t kThetaCells = 8u;
  constexpr std::size_t kPhiCells = 8u;
  constexpr std::size_t kGhostLayers = 3u;
  constexpr double kDt = 1.0e-4;

  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
  SeedNohState(state);
  const auto geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
  DEC3D_CHECK(geometry.is_valid());

  auto radial_override =
      BuildNohOuterInflowOverride(kThetaCells, kPhiCells, kGhostLayers);

  for (std::size_t step = 0; step < step_count; ++step) {
    auto hydro_view = dec3d::state::BuildHydroStateView(state);
    DEC3D_CHECK(hydro_view.is_complete());
    const auto result = dec3d::hydro::AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        kDt,
        radial_override,
        options);
    DEC3D_CHECK(result.is_complete());
  }

  return MeasureSymmetry(state);
}

[[nodiscard]] dec3d::hydro::StaticGridHydroResult RunSingleNohStep(
    const dec3d::hydro::StaticGridHydroOptions& options) {
  constexpr std::size_t kRadialCells = 32u;
  constexpr std::size_t kThetaCells = 8u;
  constexpr std::size_t kPhiCells = 8u;
  constexpr std::size_t kGhostLayers = 3u;
  constexpr double kDt = 1.0e-4;

  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
  SeedNohState(state);
  const auto geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
  DEC3D_CHECK(geometry.is_valid());

  auto radial_override =
      BuildNohOuterInflowOverride(kThetaCells, kPhiCells, kGhostLayers);
  auto hydro_view = dec3d::state::BuildHydroStateView(state);
  DEC3D_CHECK(hydro_view.is_complete());
  return dec3d::hydro::AdvanceStaticGridHydro(
      hydro_view,
      geometry,
      kDt,
      radial_override,
      options);
}

}  // namespace

int main() {
  try {
    dec3d::hydro::StaticGridHydroOptions fine_ppm;
    fine_ppm.use_ppm_reconstruction = true;
    fine_ppm.reconstruction_ghost_layers = 3u;

    auto radial_only = fine_ppm;
    radial_only.apply_geometric_source = false;
    radial_only.apply_theta_sweep = false;
    radial_only.apply_phi_sweep = false;
    const auto radial_metrics = RunShortNoh(radial_only, 20u);

    auto angular_only = fine_ppm;
    angular_only.apply_radial_sweep = false;
    const auto angular_metrics = RunShortNoh(angular_only, 20u);

    const auto full_metrics = RunShortNoh(fine_ppm, 20u);

    auto macro_radial_only = radial_only;
    macro_radial_only.use_macro_zoning = true;
    const auto macro_radial_metrics = RunShortNoh(macro_radial_only, 1u);

    auto macro_angular_only = angular_only;
    macro_angular_only.use_macro_zoning = true;
    const auto macro_angular_metrics = RunShortNoh(macro_angular_only, 1u);

    auto macro_full = fine_ppm;
    macro_full.use_macro_zoning = true;
    macro_full.debug_angular_stage_diagnostics = true;
    const auto macro_full_metrics = RunShortNoh(macro_full, 20u);
    const auto macro_full_single_step = RunSingleNohStep(macro_full);
    DEC3D_CHECK(macro_full_single_step.is_complete());

    const auto* bootstrap_diagnostic = FindDiagnostic(
        macro_full_single_step.diagnostics,
        "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed");
    DEC3D_CHECK(bootstrap_diagnostic != nullptr);
    DEC3D_CHECK(
        bootstrap_diagnostic->message.find("outer_topology_source=boundary_interior_layer") !=
        std::string::npos);
    DEC3D_CHECK(HasDiagnosticMessageFragment(
        macro_full_single_step.diagnostics,
        "p1.hydro.debug.angular.stage",
        "stage=post_macro_geometric_source_coarse_prediction"));

    DEC3D_CHECK(radial_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(radial_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);
    DEC3D_CHECK(angular_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(angular_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);
    DEC3D_CHECK(full_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(full_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);
    DEC3D_CHECK(full_metrics.max_abs_mom_theta < 1.0e-12);
    DEC3D_CHECK(full_metrics.max_abs_mom_phi < 1.0e-12);
    DEC3D_CHECK(macro_radial_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(macro_radial_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);
    DEC3D_CHECK(macro_angular_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(macro_angular_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);
    DEC3D_CHECK(macro_full_metrics.max_rho_angular_relative_spread < 1.0e-12);
    DEC3D_CHECK(macro_full_metrics.tangential_to_radial_momentum_ratio < 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
