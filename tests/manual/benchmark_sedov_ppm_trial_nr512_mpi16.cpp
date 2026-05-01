#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "hydro_benchmark_case_support.hpp"
#include "hydro_mpi_real_case_support.hpp"

#include <mpi.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifndef DEC3D_SEDOV_MPI16_RADIAL_CELLS
#define DEC3D_SEDOV_MPI16_RADIAL_CELLS 512u
#endif

#ifndef DEC3D_SEDOV_MPI16_USE_MACRO_ZONING
#define DEC3D_SEDOV_MPI16_USE_MACRO_ZONING 0
#endif

#ifndef DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT
#define DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT 16
#endif

#ifndef DEC3D_SEDOV_MPI16_AMBIENT_PRESSURE
#define DEC3D_SEDOV_MPI16_AMBIENT_PRESSURE 1.0e-5
#endif

#ifndef DEC3D_SEDOV_MPI16_BLAST_RADIUS
#define DEC3D_SEDOV_MPI16_BLAST_RADIUS 2.0e-2
#endif

#ifndef DEC3D_SEDOV_MPI16_OUTPUT_DIR_LITERAL
#define DEC3D_SEDOV_MPI16_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_sedov_spherical_ppm_trial_nr512_mpi16"
#endif

#ifndef DEC3D_SEDOV_MPI16_CASE_NAME_LITERAL
#define DEC3D_SEDOV_MPI16_CASE_NAME_LITERAL "case_sedov_spherical_ppm_trial_nr512_mpi16"
#endif

#ifndef DEC3D_SEDOV_MPI16_MESH_HANDLE_LITERAL
#define DEC3D_SEDOV_MPI16_MESH_HANDLE_LITERAL "mesh:p1.sedov.mpi16"
#endif

#ifndef DEC3D_SEDOV_MPI16_DIAGNOSTICS_HANDLE_LITERAL
#define DEC3D_SEDOV_MPI16_DIAGNOSTICS_HANDLE_LITERAL "diagnostics:p1.hydro.sedov.mpi16"
#endif

namespace {

[[nodiscard]] std::vector<double> BuildLocalShellAverage(
    const dec3d::core::Array3D<double>& local_field) {
  std::vector<double> values(local_field.extent_r(), 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_field.extent_theta() * local_field.extent_phi());
  for (std::size_t radial = 0; radial < local_field.extent_r(); ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_field.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < local_field.extent_phi(); ++phi) {
        sum += local_field(radial, theta, phi);
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

[[nodiscard]] std::vector<double> BuildLocalHydroOnlyTeAverage(
    const dec3d::state::CanonicalState& local_state) {
  std::vector<double> values(local_state.layout.radial_cells, 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_state.layout.theta_cells * local_state.layout.phi_cells);
  for (std::size_t radial = 0; radial < local_state.layout.radial_cells; ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < local_state.layout.phi_cells; ++phi) {
        const double rho = local_state.rho(radial, theta, phi);
        const double pe = dec3d::state::ElectronPressureFromElectronEnergyDensity(
            local_state.e_electron(radial, theta, phi));
        sum += (rho > 0.0) ? pe / rho : 0.0;
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

[[nodiscard]] std::vector<double> BuildLocalVelocityMagnitudeAverage(
    const dec3d::state::CanonicalState& local_state) {
  std::vector<double> values(local_state.layout.radial_cells, 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_state.layout.theta_cells * local_state.layout.phi_cells);
  for (std::size_t radial = 0; radial < local_state.layout.radial_cells; ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < local_state.layout.phi_cells; ++phi) {
        const double rho = local_state.rho(radial, theta, phi);
        if (rho > 0.0) {
          const double vr = local_state.mom_r(radial, theta, phi) / rho;
          const double vt = local_state.mom_theta(radial, theta, phi) / rho;
          const double vp = local_state.mom_phi(radial, theta, phi) / rho;
          sum += std::sqrt(vr * vr + vt * vt + vp * vp);
        }
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

void GatherRadialAverageToRoot(
    const std::vector<double>& local_values,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    int rank,
    std::vector<double>* global_values) {
  std::vector<int> recvcounts;
  std::vector<int> displacements;
  if (rank == 0) {
    recvcounts.resize(decomposition.slices.size(), 0);
    displacements.resize(decomposition.slices.size(), 0);
    for (std::size_t index = 0; index < decomposition.slices.size(); ++index) {
      const auto& slice = decomposition.slices[index];
      recvcounts[index] = static_cast<int>(slice.local_cell_count());
      displacements[index] = static_cast<int>(slice.begin_index);
    }
  }

  MPI_Gatherv(
      const_cast<double*>(local_values.data()),
      static_cast<int>(local_values.size()),
      MPI_DOUBLE,
      rank == 0 ? global_values->data() : nullptr,
      rank == 0 ? recvcounts.data() : nullptr,
      rank == 0 ? displacements.data() : nullptr,
      MPI_DOUBLE,
      0,
      MPI_COMM_WORLD);
}

[[nodiscard]] double EstimateShockRadiusFromShellAverage(
    const std::vector<double>& rho_profile,
    const std::vector<double>& radial_faces) {
  if (rho_profile.size() < 2u || radial_faces.size() < rho_profile.size() + 1u) {
    return 0.0;
  }

  std::size_t shock_index = 0u;
  double max_positive_jump = -std::numeric_limits<double>::infinity();
  for (std::size_t radial = 0; radial + 1u < rho_profile.size(); ++radial) {
    const double jump = rho_profile[radial + 1u] - rho_profile[radial];
    if (jump > max_positive_jump) {
      max_positive_jump = jump;
      shock_index = radial;
    }
  }

  return 0.5 * (radial_faces[shock_index + 1u] + radial_faces[shock_index + 2u]);
}

[[nodiscard]] dec3d::hydro::SedovRadialProfileSample BuildSampleFromShellAverages(
    const std::vector<double>& radial_faces,
    const std::vector<double>& rho_avg,
    const std::vector<double>& hydro_only_te_avg,
    const std::vector<double>& velocity_magnitude_avg,
    const std::vector<double>& mom_r_avg,
    const std::vector<double>& e_fluid_total_avg,
    std::size_t step_index,
    double time_s,
    double reference_shock_radius,
    double numerical_shock_radius) {
  dec3d::hydro::SedovRadialProfileSample sample;
  sample.step_index = step_index;
  sample.time_s = time_s;
  sample.reference_shock_radius = reference_shock_radius;
  sample.numerical_shock_radius = numerical_shock_radius;
  sample.radial_centers.reserve(rho_avg.size());
  for (std::size_t radial = 0; radial < rho_avg.size(); ++radial) {
    sample.radial_centers.push_back(0.5 * (radial_faces[radial] + radial_faces[radial + 1u]));
  }
  sample.rho_avg = rho_avg;
  sample.hydro_only_te_avg = hydro_only_te_avg;
  sample.velocity_magnitude_avg = velocity_magnitude_avg;
  sample.mom_r_avg = mom_r_avg;
  sample.e_fluid_total_avg = e_fluid_total_avg;
  return sample;
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

    if (rank_count != DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT) {
      if (rank == 0) {
        std::cerr << "expected " << DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT << " MPI ranks\n";
      }
      MPI_Finalize();
      return 1;
    }

    constexpr std::size_t kRadialCells = DEC3D_SEDOV_MPI16_RADIAL_CELLS;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kAmbientDensity = 1.0;
    constexpr double kAmbientPressure = DEC3D_SEDOV_MPI16_AMBIENT_PRESSURE;
    constexpr double kBlastEnergy = 1.0;
    constexpr double kBlastRadius = DEC3D_SEDOV_MPI16_BLAST_RADIUS;
    constexpr double kElectronEnergyFraction = 0.5;
    constexpr double kFinalTime = 2.0e-3;
    constexpr std::size_t kMaxSteps = 1000000u;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
    if (!geometry.is_valid()) {
      if (rank == 0) {
        std::cerr << "invalid geometry\n";
      }
      MPI_Finalize();
      return 1;
    }

    auto global_before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto global_seed = dec3d::state::CanonicalState::Create(global_before.layout);
    dec3d::testsupport::SeedSedovSphericalBlastCase(
        global_seed,
        global_before,
        geometry,
        kAmbientDensity,
        kAmbientPressure,
        kBlastEnergy,
        kBlastRadius,
        kElectronEnergyFraction);

    const auto decomposition = dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    if (!decomposition.is_valid()) {
      if (rank == 0) {
        std::cerr << decomposition.failure_reason << '\n';
      }
      MPI_Finalize();
      return 1;
    }

    const auto& local_slice = decomposition.slices[static_cast<std::size_t>(rank)];
    const auto local_geometry =
        dec3d::testsupport::BuildLocalGeometrySlice(geometry, local_slice, kThetaCells, kPhiCells);
    auto local_state = dec3d::testsupport::BuildLocalStateSlice(global_seed, local_slice);

    dec3d::hydro::StaticGridHydroOptions ppm_options{};
    ppm_options.apply_geometric_source = true;
    ppm_options.apply_radial_sweep = true;
    ppm_options.apply_theta_sweep = true;
    ppm_options.apply_phi_sweep = true;
    ppm_options.use_ppm_reconstruction = true;
    ppm_options.reconstruction_ghost_layers = kGhostLayers;
    ppm_options.apply_radial_ale_flux_correction = false;
    ppm_options.request_radial_ale_proposal = false;
    ppm_options.use_macro_zoning = DEC3D_SEDOV_MPI16_USE_MACRO_ZONING != 0;
    ppm_options.macro_zoning_coarse_factor = 0.5;

    std::vector<dec3d::hydro::SedovShockRadiusSample> shock_history;
    std::vector<dec3d::hydro::SedovRadialProfileSample> profile_history;
    dec3d::hydro::HydroBudgetResidualSummary reduced_budget{};
    std::ofstream progress_output;
    const std::filesystem::path output_dir =
        DEC3D_SEDOV_MPI16_OUTPUT_DIR_LITERAL;
    if (rank == 0) {
      std::filesystem::remove_all(output_dir);
      std::filesystem::create_directories(output_dir);
      progress_output.open(output_dir / "dec3d.out", std::ios::out | std::ios::trunc);
      progress_output << "# kind=run_progress\n";
      progress_output << "# columns=step time_s status\n";
      progress_output.flush();
    }

    double time_s = 0.0;
    std::size_t step = 0u;

    auto record_profile_snapshot = [&](std::size_t step_index, double sample_time) -> bool {
      auto local_rho = BuildLocalShellAverage(local_state.rho);
      auto local_te = BuildLocalHydroOnlyTeAverage(local_state);
      auto local_velocity = BuildLocalVelocityMagnitudeAverage(local_state);
      auto local_mom_r = BuildLocalShellAverage(local_state.mom_r);
      auto local_e_total = BuildLocalShellAverage(local_state.e_fluid_total);

      std::vector<double> global_rho;
      std::vector<double> global_te;
      std::vector<double> global_velocity;
      std::vector<double> global_mom_r;
      std::vector<double> global_e_total;
      if (rank == 0) {
        global_rho.resize(kRadialCells);
        global_te.resize(kRadialCells);
        global_velocity.resize(kRadialCells);
        global_mom_r.resize(kRadialCells);
        global_e_total.resize(kRadialCells);
      }

      GatherRadialAverageToRoot(local_rho, decomposition, rank, &global_rho);
      GatherRadialAverageToRoot(local_te, decomposition, rank, &global_te);
      GatherRadialAverageToRoot(local_velocity, decomposition, rank, &global_velocity);
      GatherRadialAverageToRoot(local_mom_r, decomposition, rank, &global_mom_r);
      GatherRadialAverageToRoot(local_e_total, decomposition, rank, &global_e_total);

      if (rank == 0) {
        const double reference_shock_radius =
            1.15 * std::pow(kBlastEnergy * sample_time * sample_time / kAmbientDensity, 0.2);
        const double numerical_shock_radius =
            sample_time > 0.0 ? EstimateShockRadiusFromShellAverage(global_rho, geometry.radial_faces) : 0.0;
        profile_history.push_back(BuildSampleFromShellAverages(
            geometry.radial_faces,
            global_rho,
            global_te,
            global_velocity,
            global_mom_r,
            global_e_total,
            step_index,
            sample_time,
            reference_shock_radius,
            numerical_shock_radius));
      }
      return true;
    };

    if (!record_profile_snapshot(0u, 0.0)) {
      if (rank == 0) {
        std::cerr << "failed to record initial profile\n";
      }
      MPI_Finalize();
      return 1;
    }

    while (time_s + 1.0e-16 < kFinalTime && step < kMaxSteps) {
      dec3d::hydro::HydroOperator dt_hydro;
      dt_hydro.SetStaticGridOptions(ppm_options);
      const dec3d::core::StageContext probe_context{
          time_s,
          1.0e-12,
          static_cast<std::uint64_t>(step + 1u),
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          DEC3D_SEDOV_MPI16_MESH_HANDLE_LITERAL,
          "ownership:rank" + std::to_string(rank),
          DEC3D_SEDOV_MPI16_DIAGNOSTICS_HANDLE_LITERAL};
      if (!dt_hydro.bind(probe_context, local_geometry, local_state)) {
        if (rank == 0) {
          std::cerr << "bind failed for dt estimate\n";
        }
        MPI_Finalize();
        return 1;
      }

      const auto local_dt = dt_hydro.estimate_dt();
      double step_dt = local_dt.hard_cap_dt;
      MPI_Allreduce(MPI_IN_PLACE, &step_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      step_dt = std::min(step_dt, kFinalTime - time_s);

      const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.rho, kGhostLayers);
      const auto mom_r_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, kGhostLayers);
      const auto mom_theta_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, kGhostLayers);
      const auto mom_phi_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, kGhostLayers);
      const auto e_total_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, kGhostLayers);
      const auto e_electron_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, kGhostLayers);

      auto local_hydro_view = dec3d::state::BuildHydroWorkView(local_state);
      if (!local_hydro_view.is_complete()) {
        if (rank == 0) {
          std::cerr << local_hydro_view.failure_reason << '\n';
        }
        MPI_Finalize();
        return 1;
      }

      const auto local_override = dec3d::testsupport::BuildRadialGhostOverride(
          rho_exchange,
          mom_r_exchange,
          mom_theta_exchange,
          mom_phi_exchange,
          e_total_exchange,
          e_electron_exchange);

      const auto local_result = dec3d::hydro::AdvanceStaticGridHydro(
          local_hydro_view,
          local_geometry,
          step_dt,
          local_override,
          ppm_options);
      if (!local_result.success) {
        if (rank == 0) {
          std::cerr << local_result.failure_reason << '\n';
        }
        MPI_Finalize();
        return 1;
      }

      const auto writeback = dec3d::state::CommitHydroWriteback(
          local_state,
          local_hydro_view,
          dec3d::state::BuildHydroAuthorizedWriteMask());
      if (!writeback.success) {
        if (rank == 0) {
          std::cerr << writeback.failure_reason << '\n';
        }
        MPI_Finalize();
        return 1;
      }

      dec3d::testsupport::ReduceBudgetSummaryToRoot(local_result.budget, rank, &reduced_budget);

      time_s += step_dt;
      if (rank == 0) {
        shock_history.push_back(
            {time_s, 0.0, 1.15 * std::pow(kBlastEnergy * time_s * time_s / kAmbientDensity, 0.2), 0.0});
        progress_output << (step + 1u) << ' ' << std::setprecision(17) << time_s << " pending\n";
        progress_output.flush();
      }

      if (((step + 1u) % 20u) == 0u || time_s + 1.0e-16 >= kFinalTime) {
        if (!record_profile_snapshot(step + 1u, time_s)) {
          if (rank == 0) {
            std::cerr << "failed to record profile snapshot\n";
          }
          MPI_Finalize();
          return 1;
        }
      }

      ++step;
    }

    if (time_s + 1.0e-16 < kFinalTime) {
      if (rank == 0) {
        std::cerr << "Sedov MPI baseline hit max_steps before final_time\n";
      }
      MPI_Finalize();
      return 1;
    }

    auto gathered_after = dec3d::state::CanonicalState::Create(global_before.layout);
    dec3d::testsupport::GatherFieldToRoot(local_state.rho, decomposition, rank, gathered_after.rho);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_r, decomposition, rank, gathered_after.mom_r);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_theta, decomposition, rank, gathered_after.mom_theta);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_phi, decomposition, rank, gathered_after.mom_phi);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_fluid_total, decomposition, rank, gathered_after.e_fluid_total);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_electron, decomposition, rank, gathered_after.e_electron);

    if (rank == 0) {
      for (auto& sample : shock_history) {
        const auto matching = std::find_if(
            profile_history.begin(),
            profile_history.end(),
            [&](const dec3d::hydro::SedovRadialProfileSample& profile) {
              return std::abs(profile.time_s - sample.time_s) <= 1.0e-15;
            });
        if (matching != profile_history.end()) {
          sample.numerical_shock_radius = matching->numerical_shock_radius;
          sample.relative_error =
              matching->reference_shock_radius > 0.0
                  ? std::abs(matching->numerical_shock_radius - matching->reference_shock_radius) /
                        matching->reference_shock_radius
                  : 0.0;
        }
      }

      const auto summary = dec3d::hydro::EvaluateSedovSphericalBlast(
          global_before,
          gathered_after,
          geometry,
          kBlastEnergy,
          kAmbientDensity,
          time_s,
          0.25,
          1.0e-6,
          1.0e-8);
      if (!summary.is_complete()) {
        std::cerr << "sedov mpi summary incomplete\n";
        MPI_Finalize();
        return 1;
      }

      if (!dec3d::hydro::WriteSedovSphericalBlastOutputs(
              global_before,
              gathered_after,
              geometry,
              output_dir,
              summary,
              shock_history,
              profile_history,
              &reduced_budget,
              DEC3D_SEDOV_MPI16_CASE_NAME_LITERAL)) {
        std::cerr << "failed to write sedov mpi outputs\n";
        MPI_Finalize();
        return 1;
      }

      std::cout << summary.report_line << '\n';
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::cerr << error.what() << '\n';
    return 1;
  }
}
