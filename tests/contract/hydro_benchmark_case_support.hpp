#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace dec3d::testsupport {

struct SedovBaselineRunResult {
  bool success{false};
  dec3d::state::CanonicalState before;
  dec3d::state::CanonicalState after;
  std::vector<dec3d::hydro::SedovShockRadiusSample> shock_radius_history;
  std::vector<dec3d::hydro::SedovRadialProfileSample> profile_history;
  dec3d::hydro::HydroBudgetResidualSummary final_budget;
  double final_time_s{0.0};
  std::string failure_reason;
};

inline void SeedSedovSphericalBlastCase(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double ambient_density,
    double ambient_pressure,
    double blast_energy,
    double blast_radius,
    double electron_energy_fraction) {
  const double base_electron_pressure = electron_energy_fraction * ambient_pressure;
  const double base_electron_energy =
      dec3d::state::ElectronEnergyDensityFromPressure(base_electron_pressure);

  double deposited_volume = 0.0;
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const double outer_face = geometry.radial_faces[radial + 1u];
    if (outer_face <= blast_radius + 1.0e-14) {
      for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
          const std::size_t linear_index =
              ((radial * state.layout.theta_cells) + theta) * state.layout.phi_cells + phi;
          deposited_volume += geometry.cell_volumes[linear_index];
        }
      }
    }
  }

  const double blast_energy_density =
      deposited_volume > 0.0 ? blast_energy / deposited_volume : 0.0;
  const double base_internal_energy =
      ambient_pressure / (dec3d::state::HydroIdealGasGamma() - 1.0);

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const double outer_face = geometry.radial_faces[radial + 1u];
    const bool in_blast = outer_face <= blast_radius + 1.0e-14;
    const double internal_energy = base_internal_energy + (in_blast ? blast_energy_density : 0.0);
    const double electron_energy =
        base_electron_energy + (in_blast ? electron_energy_fraction * blast_energy_density : 0.0);

    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = ambient_density;
        state.mom_r(radial, theta, phi) = 0.0;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = internal_energy;
        state.e_electron(radial, theta, phi) = electron_energy;

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

inline dec3d::hydro::StaticGridHydroOptions SedovBaselineOptions() noexcept {
  dec3d::hydro::StaticGridHydroOptions options{};
  options.apply_geometric_source = true;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  options.apply_radial_ale_flux_correction = false;
  options.request_radial_ale_proposal = false;
  options.use_macro_zoning = false;
  return options;
}

inline SedovBaselineRunResult RunSedovBaseline(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells,
    double ambient_density,
    double ambient_pressure,
    double blast_energy,
    double blast_radius,
    double electron_energy_fraction,
    double final_time_s,
    std::size_t max_steps,
    const std::filesystem::path* progress_path = nullptr,
    const dec3d::hydro::StaticGridHydroOptions* options_override = nullptr) {
  SedovBaselineRunResult run;
  run.before = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{radial_cells, theta_cells, phi_cells, 0u});
  run.after = dec3d::state::CanonicalState::Create(run.before.layout);

  if (!geometry.is_valid()) {
    run.failure_reason = "Sedov baseline requires valid spherical geometry";
    return run;
  }
  if (!(final_time_s > 0.0) || max_steps == 0u) {
    run.failure_reason = "Sedov baseline requires positive final time and max steps";
    return run;
  }

  SeedSedovSphericalBlastCase(
      run.after,
      run.before,
      geometry,
      ambient_density,
      ambient_pressure,
      blast_energy,
      blast_radius,
      electron_energy_fraction);

  std::ofstream progress_output;
  std::filesystem::path live_output_directory;
  if (progress_path != nullptr) {
    std::error_code error;
    std::filesystem::create_directories(progress_path->parent_path(), error);
    if (error) {
      run.failure_reason = "Sedov baseline could not create progress directory";
      return run;
    }
    live_output_directory = progress_path->parent_path();
    progress_output.open(*progress_path, std::ios::out | std::ios::trunc);
    if (!progress_output.is_open()) {
      run.failure_reason = "Sedov baseline could not open dec3d.out";
      return run;
    }
    progress_output << "# kind=run_progress\n";
    progress_output << "# columns=step time_s status\n";
    progress_output.flush();
  }

  auto write_profile_snapshot = [&](const dec3d::hydro::SedovRadialProfileSample& sample) -> bool {
    run.profile_history.push_back(sample);
    if (live_output_directory.empty()) {
      return true;
    }

    std::ostringstream file_name;
    file_name << "radial_profile_step"
              << std::setw(6) << std::setfill('0') << sample.step_index
              << ".txt";
    std::ostringstream label;
    label << "step_" << sample.step_index;
    return dec3d::hydro::WriteSedovRadialProfileSampleFile(
        live_output_directory / file_name.str(),
        label.str().c_str(),
        sample);
  };

  auto initial_profile = dec3d::hydro::BuildSedovRadialProfileSample(run.after, geometry, 0u, 0.0);
  initial_profile.reference_shock_radius = 0.0;
  initial_profile.numerical_shock_radius = 0.0;
  if (!write_profile_snapshot(initial_profile)) {
    run.failure_reason = "Sedov baseline could not write initial radial profile snapshot";
    return run;
  }

  const auto options = options_override != nullptr ? *options_override : SedovBaselineOptions();
  std::size_t step = 0u;
  double time_s = 0.0;

  while (time_s + 1.0e-16 < final_time_s && step < max_steps) {
    dec3d::hydro::HydroOperator hydro;
    hydro.SetStaticGridOptions(options);

    const dec3d::core::StageContext probe_context{
        time_s,
        1.0e-12,
        static_cast<std::uint64_t>(step + 1u),
        dec3d::core::PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.sedov",
        "ownership:rank0",
        "diagnostics:p1.hydro.sedov"};
    if (!hydro.bind(probe_context, geometry, run.after)) {
      run.failure_reason = "Sedov baseline could not bind hydro operator";
      return run;
    }

    const auto dt_advice = hydro.estimate_dt();
    if (!dt_advice.is_complete() || !(dt_advice.hard_cap_dt > 0.0) ||
        !std::isfinite(dt_advice.hard_cap_dt)) {
      run.failure_reason = "Sedov baseline could not derive a finite dt";
      return run;
    }

    const double step_dt = std::min(final_time_s - time_s, dt_advice.hard_cap_dt);
    const dec3d::core::StageContext step_context{
        time_s,
        step_dt,
        static_cast<std::uint64_t>(step + 1u),
        dec3d::core::PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.sedov",
        "ownership:rank0",
        "diagnostics:p1.hydro.sedov"};
    if (!hydro.bind(step_context, geometry, run.after)) {
      run.failure_reason = "Sedov baseline could not rebind hydro operator for advance";
      return run;
    }

    const auto result = hydro.advance();
    if (!result.is_semantically_complete()) {
      run.failure_reason = result.failure_reason.empty() ? "Sedov baseline hydro advance failed"
                                                         : result.failure_reason;
      return run;
    }

    if (!dec3d::hydro::TryExtractHydroBudgetResidualSummary(
            result.diagnostics,
            &run.final_budget)) {
      run.failure_reason = "Sedov baseline could not extract budget summary";
      return run;
    }

    time_s += step_dt;
    run.final_time_s = time_s;
    const auto check = dec3d::hydro::EvaluateSedovSphericalBlast(
        run.before,
        run.after,
        geometry,
        blast_energy,
        ambient_density,
        time_s,
        1.0,
        1.0,
        1.0);
    run.shock_radius_history.push_back(
        {time_s,
         check.numerical_shock_radius,
         check.reference_shock_radius,
         check.shock_radius_relative_error});
    if (((step + 1u) % 20u) == 0u || time_s + 1.0e-16 >= final_time_s) {
      auto profile = dec3d::hydro::BuildSedovRadialProfileSample(
          run.after,
          geometry,
          step + 1u,
          time_s);
      profile.reference_shock_radius = check.reference_shock_radius;
      profile.numerical_shock_radius = check.numerical_shock_radius;
      if (!write_profile_snapshot(profile)) {
        run.failure_reason = "Sedov baseline could not write radial profile snapshot";
        return run;
      }
    }
    if (progress_output.is_open()) {
      progress_output << (step + 1u) << ' '
                      << time_s << ' '
                      << "pending\n";
      progress_output.flush();
    }
    ++step;
  }

  if (time_s + 1.0e-16 < final_time_s) {
    run.failure_reason = "Sedov baseline hit max_steps before final_time";
    return run;
  }

  run.success = true;
  return run;
}

}  // namespace dec3d::testsupport
