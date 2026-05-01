#include "benchmarks/p3_radiation_benchmark.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/distributed_multigroup_radiation.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace {

constexpr double kOuterRadiusCm = 1.0e-2;
constexpr double kInterfaceRadiusCm = 5.0e-3;
constexpr double kRhoInner = 50.0;
constexpr double kRhoOuter = 100.0;
constexpr double kTeInnerKeV = 5.0;
constexpr double kTeOuterKeV = 0.5;
constexpr const char* kTableRoot = "F:/dec3d/data/opacities/tops_dt_2026_04_27";

std::string JoinPath(const std::string& root, const std::string& leaf) {
  return (std::filesystem::path(root) / leaf).generic_string();
}

dec3d::radiation::TopsOpacityProviderOptions MakeProviderOptions() {
  dec3d::radiation::TopsOpacityProviderOptions options;
  options.table_root = kTableRoot;
  options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  options.density_clip_policy = dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return options;
}

dec3d::mesh::SphericalGeometryMetadata MakeGeometry(std::size_t radial_cells) {
  return dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{radial_cells, 1u, 2u, 0.0, kOuterRadiusCm});
}

dec3d::transport::GenericDiffusionBoundaryPolicy MakeBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::periodic};
}

void FillLocalState(
    dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  for (std::size_t local_r = 0u; local_r < state.layout.radial_cells; ++local_r) {
    const std::size_t global_r = ownership.global_radial_begin + local_r;
    const double radius =
        0.5 * (geometry.radial_faces[global_r] + geometry.radial_faces[global_r + 1u]);
    const bool inner = radius < kInterfaceRadiusCm;
    const double rho = inner ? kRhoInner : kRhoOuter;
    const double te = dec3d::physics::ErgFromKeV(inner ? kTeInnerKeV : kTeOuterKeV);
    const double ti = te;
    const double ni = rho / mean_ion_mass;
    const double ne = ni;
    const double e_e = ne * te / gamma_minus_one;
    const double e_i = ni * ti / gamma_minus_one;
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        state.rho(local_r, theta, phi) = rho;
        state.mom_r(local_r, theta, phi) = 0.0;
        state.mom_theta(local_r, theta, phi) = 0.0;
        state.mom_phi(local_r, theta, phi) = 0.0;
        state.e_electron(local_r, theta, phi) = e_e;
        state.e_fluid_total(local_r, theta, phi) = e_e + e_i;
      }
    }
  }
}

void InitializeRadiationFromProvider(
    dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationCoefficientProviderResult& provider) {
  for (std::size_t group = 0u; group < state.radiation_groups.size(); ++group) {
    for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
      for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
          state.radiation_groups[group](r, theta, phi) =
              provider.coefficients.B_erg_per_cm3[group](r, theta, phi);
        }
      }
    }
  }
}

bool WriteLocalProfile(
    const std::string& output_dir,
    std::size_t step,
    int rank,
    const dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  std::filesystem::create_directories(output_dir);
  std::ostringstream name;
  name << "profile_step" << std::setw(6) << std::setfill('0') << step
       << "_rank" << std::setw(5) << std::setfill('0') << rank << ".csv";
  std::ofstream out(JoinPath(output_dir, name.str()));
  if (!out) {
    return false;
  }
  out << std::setprecision(17)
      << "step,rank,global_r,radius_cm,rho_g_cm3,Te_keV,sum_Ug_erg_cm3\n";
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  for (std::size_t local_r = 0u; local_r < state.layout.radial_cells; ++local_r) {
    const std::size_t global_r = ownership.global_radial_begin + local_r;
    const double radius =
        0.5 * (geometry.radial_faces[global_r] + geometry.radial_faces[global_r + 1u]);
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double rho = state.rho(local_r, theta, phi);
        const double ne = rho / mean_ion_mass;
        const double te_erg = gamma_minus_one * state.e_electron(local_r, theta, phi) / ne;
        double sum_ug = 0.0;
        for (const auto& group : state.radiation_groups) {
          sum_ug += group(local_r, theta, phi);
        }
        out << step << ',' << rank << ',' << global_r << ',' << radius << ','
            << rho << ',' << dec3d::physics::KeVFromErg(te_erg) << ','
            << sum_ug << '\n';
      }
    }
  }
  return static_cast<bool>(out);
}

double ElectronTemperatureKeV(
    const dec3d::state::CanonicalState& state,
    std::size_t local_r,
    std::size_t theta,
    std::size_t phi) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  const double rho = state.rho(local_r, theta, phi);
  const double ne = rho / mean_ion_mass;
  const double te_erg = gamma_minus_one * state.e_electron(local_r, theta, phi) / ne;
  return dec3d::physics::KeVFromErg(te_erg);
}

double OuterRadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta,
    std::size_t phi) {
  const double radius = geometry.radial_faces.back();
  const double polar_factor =
      std::cos(geometry.theta_faces[theta]) - std::cos(geometry.theta_faces[theta + 1u]);
  const double dphi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return radius * radius * polar_factor * dphi;
}

double MarshakFaceConductanceSpeed(double Dbar_cm2_s, double center_to_face_distance_cm) {
  const double h = 0.5 * dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double d_over_distance = Dbar_cm2_s / center_to_face_distance_cm;
  return (h * d_over_distance) / (h + d_over_distance);
}

dec3d::radiation::RadiationFluxLimiterModel ParseRadiationFluxLimiterModel(
    const std::string& value) {
  using dec3d::radiation::RadiationFluxLimiterModel;
  if (value == "disabled" || value == "none" || value == "no_limiter") {
    return RadiationFluxLimiterModel::disabled;
  }
  if (value == "harmonic" || value == "harmonic_eq_5_209") {
    return RadiationFluxLimiterModel::harmonic_eq_5_209;
  }
  if (value == "larsen" || value == "larsen_n2" || value == "larsen_eq_5_210_n2") {
    return RadiationFluxLimiterModel::larsen_eq_5_210_n2;
  }
  if (value == "minmax" || value == "min_max" || value == "minmax_eq_5_211") {
    return RadiationFluxLimiterModel::minmax_eq_5_211;
  }
  return RadiationFluxLimiterModel::unsupported;
}

bool WriteOuterGroupDiagnostics(
    const std::string& output_dir,
    std::size_t step,
    int rank,
    const dec3d::state::CanonicalState& old_state,
    const dec3d::state::CanonicalState& new_state,
    const dec3d::radiation::RadiationCoefficientProviderResult& old_time_provider,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) {
  if (ownership.global_radial_end != ownership.global_radial_cells ||
      new_state.layout.radial_cells == 0u) {
    return true;
  }
  std::filesystem::create_directories(output_dir);
  std::ostringstream name;
  name << "outer_group_diagnostics_rank" << std::setw(5) << std::setfill('0') << rank
       << ".csv";
  const auto path = JoinPath(output_dir, name.str());
  const bool write_header = !std::filesystem::exists(path);
  std::ofstream out(path, std::ios::app);
  if (!out) {
    return false;
  }
  if (write_header) {
    out << "step,rank,group,theta,phi,global_r,radius_um,Te_old_keV,Te_new_keV,"
           "U_old_erg_cm3,U_new_erg_cm3,B_erg_cm3,kappaP_cm_inv,Dbar_cm2_s,"
           "source_gain_density_erg_cm3,delta_U_density_erg_cm3,"
           "budget_boundary_density_erg_cm3,marshak_G_face_cm_s,"
           "marshak_loss_rate_s_inv,marshak_loss_density_est_erg_cm3,"
           "U_over_B_new\n";
  }

  const std::size_t local_r = new_state.layout.radial_cells - 1u;
  const std::size_t global_r = ownership.global_radial_cells - 1u;
  const double radius_center =
      0.5 * (geometry.radial_faces[global_r] + geometry.radial_faces[global_r + 1u]);
  const double center_to_outer_face = geometry.radial_faces.back() - radius_center;
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;

  out << std::setprecision(17);
  for (std::size_t theta = 0u; theta < new_state.layout.theta_cells; ++theta) {
    for (std::size_t phi = 0u; phi < new_state.layout.phi_cells; ++phi) {
      const std::size_t global_cell_index =
          ((global_r * ownership.global_theta_cells) + theta) * ownership.global_phi_cells +
          phi;
      const double volume = geometry.cell_volumes.at(global_cell_index);
      const double area = OuterRadialFaceArea(geometry, theta, phi);
      for (std::size_t group = 0u; group < new_state.radiation_groups.size(); ++group) {
        const double old_u = old_state.radiation_groups[group](local_r, theta, phi);
        const double new_u = new_state.radiation_groups[group](local_r, theta, phi);
        const double B =
            old_time_provider.coefficients.B_erg_per_cm3[group](local_r, theta, phi);
        const double kappaP =
            old_time_provider.coefficients.kappaP_cm_inv[group](local_r, theta, phi);
        const double Dbar =
            old_time_provider.coefficients.Dbar_cm2_per_s[group](local_r, theta, phi);
        const double source_gain = dt_s * c * kappaP * (B - new_u);
        const double delta_u = new_u - old_u;
        const double budget_boundary = source_gain - delta_u;
        const double G = MarshakFaceConductanceSpeed(Dbar, center_to_outer_face);
        const double loss_rate = (area / volume) * G;
        const double loss_density_est = dt_s * loss_rate * new_u;
        out << step << ',' << rank << ',' << group << ',' << theta << ',' << phi << ','
            << global_r << ',' << radius_center * 1.0e4 << ','
            << ElectronTemperatureKeV(old_state, local_r, theta, phi) << ','
            << ElectronTemperatureKeV(new_state, local_r, theta, phi) << ','
            << old_u << ',' << new_u << ',' << B << ',' << kappaP << ',' << Dbar << ','
            << source_gain << ',' << delta_u << ',' << budget_boundary << ',' << G << ','
            << loss_rate << ',' << loss_density_est << ','
            << (B > 0.0 ? new_u / B : 0.0) << '\n';
      }
    }
  }
  return static_cast<bool>(out);
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const std::size_t radial_cells =
      argc > 1 ? static_cast<std::size_t>(std::strtoull(argv[1], nullptr, 10)) : 256u;
  const std::size_t steps =
      argc > 2 ? static_cast<std::size_t>(std::strtoull(argv[2], nullptr, 10)) : 5u;
  const double dt_s = argc > 3 ? std::strtod(argv[3], nullptr) : 1.0e-10;
  const std::string output_dir =
      argc > 4 ? argv[4]
               : "F:/dec3d/analysis/output/p3_radiation_benchmarks/"
                 "p3_b1_mpi24_nr256";
  const std::string limiter_arg = argc > 5 ? argv[5] : "disabled";
  const auto limiter_model = ParseRadiationFluxLimiterModel(limiter_arg);
  const int max_iterations =
      argc > 6 ? static_cast<int>(std::strtol(argv[6], nullptr, 10)) : 200;
  if (max_iterations <= 0) {
    if (rank == 0) {
      std::cerr << "max_iterations must be positive\n";
    }
    MPI_Finalize();
    return 1;
  }
  if (limiter_model == dec3d::radiation::RadiationFluxLimiterModel::unsupported) {
    if (rank == 0) {
      std::cerr << "unsupported radiation flux limiter model: " << limiter_arg
                << "\nvalid values: disabled, harmonic, larsen, minmax\n";
    }
    MPI_Finalize();
    return 1;
  }

  if (rank == 0) {
    std::cout << "p3_b1_mpi=true"
              << "; mpi_ranks=" << size
              << "; radial_cells=" << radial_cells
              << "; theta_cells=1"
               << "; phi_cells=2"
               << "; steps=" << steps
               << "; dt_s=" << std::setprecision(17) << dt_s
               << "; radiation_flux_limiter="
               << dec3d::radiation::RadiationFluxLimiterModelName(limiter_model)
               << "; gmres_max_iterations=" << max_iterations
               << "; output_dir=" << output_dir << '\n';
  }

  const auto group_layout = dec3d::benchmarks::MakeP3B1TwelveGroupLayout();
  const auto geometry = MakeGeometry(radial_cells);
  const auto ownership =
      dec3d::transport::BuildDistributedDiffusionRowOwnership(
          MPI_COMM_WORLD, radial_cells, 1u, 2u);
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{
          ownership.global_radial_end - ownership.global_radial_begin,
          1u,
          2u,
          group_layout.group_count});
  FillLocalState(state, ownership, geometry);

  const auto provider_options = MakeProviderOptions();
  const auto table_load_begin = std::chrono::steady_clock::now();
  const auto table = dec3d::radiation::LoadTopsOpacityTable(provider_options);
  const auto table_load_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - table_load_begin)
          .count();
  if (!table.success) {
    if (rank == 0) {
      std::cerr << table.failure_reason << '\n' << table.failure_diagnostics << '\n';
    }
    MPI_Finalize();
    return 1;
  }

  auto initial_provider = dec3d::radiation::BuildRadiationCoefficientArrays(
      state, geometry, group_layout, table.table, provider_options);
  if (!initial_provider.success) {
    if (rank == 0) {
      std::cerr << initial_provider.failure_reason << '\n'
                << initial_provider.failure_diagnostics << '\n';
    }
    MPI_Finalize();
    return 1;
  }
  InitializeRadiationFromProvider(state, initial_provider);

  WriteLocalProfile(output_dir, 0u, rank, state, ownership, geometry);
  WriteOuterGroupDiagnostics(
      output_dir, 0u, rank, state, state, initial_provider, ownership, geometry, 0.0);

  if (rank == 0) {
    std::ofstream timing(JoinPath(output_dir, "timing_rank0.csv"));
    timing << "step,max_step_wall_seconds,success,delta_radiation_total,"
              "delta_electron_total,boundary_leak_total,budget_residual\n";
  }

  for (std::size_t step = 1u; step <= steps; ++step) {
    const auto state_before_step = state;
    auto provider_before_step = dec3d::radiation::BuildRadiationCoefficientArrays(
        state_before_step, geometry, group_layout, table.table, provider_options);
    if (!provider_before_step.success) {
      if (rank == 0) {
        std::cerr << provider_before_step.failure_reason << '\n'
                  << provider_before_step.failure_diagnostics << '\n';
      }
      MPI_Finalize();
      return 1;
    }

    dec3d::radiation::DistributedMultigroupRadiationProblem problem;
    problem.ownership = ownership;
    problem.local_state = state;
    problem.global_geometry = geometry;
    problem.boundary_policy = MakeBoundary();
    problem.group_layout = group_layout;
    problem.opacity_table = &table.table;
    problem.dt_s = dt_s;

    dec3d::radiation::DistributedMultigroupRadiationOptions options;
    options.provider_options = provider_options;
    options.boundary_model =
        dec3d::radiation::DistributedRadiationBoundaryModel::thesis_marshak_vacuum;
    options.solve_options.communicator = MPI_COMM_WORLD;
    options.solve_options.relative_tolerance = 1.0e-8;
    options.solve_options.max_iterations = max_iterations;
    options.recovery_options.electron_energy_floor = 0.0;
    options.recovery_options.ion_energy_floor = 0.0;
    options.radiation_energy_floor_erg_per_cm3 = 0.0;
    options.electron_energy_floor_erg_per_cm3 = 0.0;
    options.ion_energy_floor_erg_per_cm3 = 0.0;
    options.radiation_flux_limiter.model = limiter_model;
    options.table_reused_across_steps = true;

    MPI_Barrier(MPI_COMM_WORLD);
    const auto step_begin = std::chrono::steady_clock::now();
    const auto result =
        dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(problem, options);
    const double local_seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - step_begin)
            .count();
    double max_seconds = 0.0;
    MPI_Reduce(&local_seconds, &max_seconds, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    int local_success = result.success ? 1 : 0;
    int global_success = 0;
    MPI_Allreduce(&local_success, &global_success, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    if (global_success == 0) {
      std::cerr << "rank=" << rank << "; step=" << step
                << "; failure_reason=" << result.failure_reason
                << "; diagnostics=" << result.failure_diagnostics << '\n';
      MPI_Finalize();
      return 2;
    }
    state = problem.local_state;
    WriteLocalProfile(output_dir, step, rank, state, ownership, geometry);
    if (!WriteOuterGroupDiagnostics(
            output_dir,
            step,
            rank,
            state_before_step,
            state,
            provider_before_step,
            ownership,
            geometry,
            dt_s)) {
      std::cerr << "rank=" << rank
                << "; failure_reason=failed to write outer group diagnostics\n";
      MPI_Finalize();
      return 3;
    }

    if (rank == 0) {
      std::ofstream timing(JoinPath(output_dir, "timing_rank0.csv"), std::ios::app);
      timing << std::setprecision(17)
             << step << ',' << max_seconds << ",true,"
             << result.delta_radiation_total_all_groups << ','
             << result.delta_electron_total << ','
             << result.boundary_leak_total_all_groups << ','
             << result.global_radiation_plus_electron_residual << '\n';
      std::cout << "step=" << step
                << "; max_step_wall_seconds=" << max_seconds
                << "; delta_radiation_total=" << result.delta_radiation_total_all_groups
                << "; delta_electron_total=" << result.delta_electron_total
                << "; boundary_leak_total=" << result.boundary_leak_total_all_groups
                << "; budget_residual=" << result.global_radiation_plus_electron_residual
                << '\n';
    }
  }

  if (rank == 0) {
    std::ofstream metadata(JoinPath(output_dir, "metadata_mpi.json"));
    metadata << std::setprecision(17)
             << "{\n"
             << "  \"case_id\": \"p3_b1_woo_5_25_slab_no_limiter_mpi\",\n"
             << "  \"mpi_ranks\": " << size << ",\n"
             << "  \"radial_cells\": " << radial_cells << ",\n"
             << "  \"theta_cells\": 1,\n"
             << "  \"phi_cells\": 2,\n"
             << "  \"steps\": " << steps << ",\n"
              << "  \"dt_s\": " << dt_s << ",\n"
              << "  \"gmres_max_iterations\": " << max_iterations << ",\n"
              << "  \"table_load_seconds_rank0\": " << table_load_seconds << ",\n"
              << "  \"radiation_boundary_model\": \"thesis_marshak_vacuum\",\n"
              << "  \"radiation_flux_limiter\": \""
              << dec3d::radiation::RadiationFluxLimiterModelName(limiter_model)
              << "\",\n"
              << "  \"parity_claim_allowed\": false\n"
              << "}\n";
    std::cout << "metadata=" << JoinPath(output_dir, "metadata_mpi.json") << '\n';
  }

  MPI_Finalize();
  return 0;
}
