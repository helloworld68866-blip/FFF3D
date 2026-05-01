#include "benchmarks/p2_thermal_benchmark_runner.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/electron_ion_equilibration.hpp"
#include "state/thermodynamics/electron_ion_tau.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"
#include "transport/thermal/distributed_thermal_conduction.hpp"
#include "transport/thermal/thermal_conductivity.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace dec3d::benchmarks {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

[[nodiscard]] dec3d::transport::ThermalConductionKappaModel ToTransportKappa(
    BenchmarkThermalKappaModel model) noexcept {
  switch (model) {
    case BenchmarkThermalKappaModel::spitzer_no_degeneracy:
      return dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy;
    case BenchmarkThermalKappaModel::lee_more_with_degeneracy:
      return dec3d::transport::ThermalConductionKappaModel::lee_more_with_degeneracy;
    case BenchmarkThermalKappaModel::constant_kappa_exact:
      return dec3d::transport::ThermalConductionKappaModel::constant_user_supplied;
  }
  return dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy;
}

[[nodiscard]] dec3d::transport::ThermalConductivityModel ToProviderKappa(
    BenchmarkThermalKappaModel model) noexcept {
  switch (model) {
    case BenchmarkThermalKappaModel::spitzer_no_degeneracy:
      return dec3d::transport::ThermalConductivityModel::spitzer_no_degeneracy;
    case BenchmarkThermalKappaModel::lee_more_with_degeneracy:
      return dec3d::transport::ThermalConductivityModel::lee_more_with_degeneracy;
    case BenchmarkThermalKappaModel::constant_kappa_exact:
      return dec3d::transport::ThermalConductivityModel::spitzer_no_degeneracy;
  }
  return dec3d::transport::ThermalConductivityModel::spitzer_no_degeneracy;
}

[[nodiscard]] std::string BuildFailure(
    const char* reason,
    const WooThermalBenchmarkDescriptor& descriptor) {
  std::ostringstream out;
  out << "diagnostic_id=p2.thermal_benchmark.failure"
      << "; case_id=" << descriptor.case_id
      << "; failure_reason=" << reason
      << "; canonical_state_mutated=false";
  return out.str();
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildBenchmarkGeometry(
    const WooThermalBenchmarkDescriptor& descriptor,
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  geometry.theta_faces.resize(theta_cells + 1u);
  geometry.phi_faces.resize(phi_cells + 1u);
  geometry.cell_volumes.resize(radial_cells * theta_cells * phi_cells, 0.0);

  for (std::size_t r = 0; r <= radial_cells; ++r) {
    geometry.radial_faces[r] =
        descriptor.r_max_cm * static_cast<double>(r) / static_cast<double>(radial_cells);
  }
  constexpr double theta_lower = 0.25 * kPi;
  constexpr double theta_upper = 0.75 * kPi;
  for (std::size_t theta = 0; theta <= theta_cells; ++theta) {
    geometry.theta_faces[theta] =
        theta_lower +
        (theta_upper - theta_lower) * static_cast<double>(theta) /
            static_cast<double>(theta_cells);
  }
  for (std::size_t phi = 0; phi <= phi_cells; ++phi) {
    geometry.phi_faces[phi] =
        2.0 * kPi * static_cast<double>(phi) / static_cast<double>(phi_cells);
  }

  geometry.global_volume = 0.0;
  std::size_t linear = 0;
  for (std::size_t r = 0; r < radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) / 3.0;
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      const double polar_factor =
          std::cos(geometry.theta_faces[theta]) -
          std::cos(geometry.theta_faces[theta + 1u]);
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const double azimuthal_factor =
            geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
        geometry.cell_volumes[linear] =
            radial_factor * polar_factor * azimuthal_factor;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy Thesis1dBoundary() noexcept {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

void FillBenchmarkSlabState(
    dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const WooThermalBenchmarkDescriptor& descriptor) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();

  for (std::size_t lr = 0; lr < state.layout.radial_cells; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    const double radial_center =
        0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]);
    const bool inner = radial_center < descriptor.r0_cm;
    const double rho = inner ? descriptor.rho_inner_g_cm3 : descriptor.rho_outer_g_cm3;
    const double te = dec3d::physics::ErgFromKeV(
        inner ? descriptor.te_inner_keV : descriptor.te_outer_keV);
    const double ti =
        descriptor.initial_ti_policy == "ion_uniform_0p5keV"
            ? dec3d::physics::ErgFromKeV(0.5)
            : te;
    const double ni = rho / mean_ion_mass;
    const double ne = ni * dec3d::physics::DefaultZbar();
    const double e_electron = ne * te / gamma_minus_one;
    const double e_ion = ni * ti / gamma_minus_one;

    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(lr, theta, phi) = rho;
        state.mom_r(lr, theta, phi) = 0.0;
        state.mom_theta(lr, theta, phi) = 0.0;
        state.mom_phi(lr, theta, phi) = 0.0;
        state.e_electron(lr, theta, phi) = e_electron;
        state.e_fluid_total(lr, theta, phi) = e_electron + e_ion;
      }
    }
  }
}

[[nodiscard]] dec3d::transport::DistributedThermalConductionProblem BuildBenchmarkProblem(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalBenchmarkRunOptions& options) {
  dec3d::transport::DistributedThermalConductionProblem problem;
  problem.ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
      communicator,
      options.global_radial_cells,
      options.theta_cells,
      options.phi_cells);
  const std::size_t local_radial =
      problem.ownership.global_radial_end >= problem.ownership.global_radial_begin
          ? problem.ownership.global_radial_end - problem.ownership.global_radial_begin
          : 0u;
  problem.local_state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{
          local_radial,
          options.theta_cells,
          options.phi_cells,
          0});
  problem.global_geometry = BuildBenchmarkGeometry(
      descriptor,
      options.global_radial_cells,
      options.theta_cells,
      options.phi_cells);
  problem.boundary_policy = Thesis1dBoundary();
  FillBenchmarkSlabState(
      problem.local_state,
      problem.ownership,
      problem.global_geometry,
      descriptor);
  return problem;
}

struct LocalShellProfiles {
  std::vector<double> r_cm;
  std::vector<double> rho_g_cm3;
  std::vector<double> te_keV;
  std::vector<double> ti_keV;
  std::vector<double> kappa_e_cm_inv_s;
  std::vector<double> kappa_i_cm_inv_s;
  std::vector<double> tau_ei_s;
};

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  return geometry.cell_volumes.at(((radial * theta_cells) + theta) * phi_cells + phi);
}

[[nodiscard]] double LocalThermalEnergy(
    const dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const auto recovered = dec3d::state::RecoverThermodynamicState(state);
  if (!recovered.success) {
    return 0.0;
  }

  double thermal_energy = 0.0;
  for (std::size_t lr = 0; lr < state.layout.radial_cells; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const auto& cell = recovered.cells(lr, theta, phi);
        thermal_energy +=
            CellVolume(
                geometry,
                gr,
                theta,
                phi,
                state.layout.theta_cells,
                state.layout.phi_cells) *
            (cell.e_electron_erg_per_cm3 + cell.e_ion_erg_per_cm3);
      }
    }
  }
  return thermal_energy;
}

[[nodiscard]] LocalShellProfiles BuildLocalShellProfiles(
    const dec3d::transport::DistributedThermalConductionProblem& problem,
    const WooThermalBenchmarkDescriptor& descriptor) {
  LocalShellProfiles profiles;
  const auto& state = problem.local_state;
  const std::size_t local_radial = state.layout.radial_cells;
  profiles.r_cm.resize(local_radial, 0.0);
  profiles.rho_g_cm3.resize(local_radial, 0.0);
  profiles.te_keV.resize(local_radial, 0.0);
  profiles.ti_keV.resize(local_radial, 0.0);
  profiles.kappa_e_cm_inv_s.resize(local_radial, 0.0);
  profiles.kappa_i_cm_inv_s.resize(local_radial, 0.0);
  profiles.tau_ei_s.resize(local_radial, 0.0);

  const auto recovered = dec3d::state::RecoverThermodynamicState(state);
  if (!recovered.success) {
    return profiles;
  }

  const auto conductivity_model = ToProviderKappa(descriptor.thermal_kappa_model);
  for (std::size_t lr = 0; lr < local_radial; ++lr) {
    const std::size_t gr = problem.ownership.global_radial_begin + lr;
    profiles.r_cm[lr] =
        0.5 * (problem.global_geometry.radial_faces[gr] +
               problem.global_geometry.radial_faces[gr + 1u]);

    double angular_count = 0.0;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const auto& cell = recovered.cells(lr, theta, phi);
        profiles.rho_g_cm3[lr] += cell.rho_g_per_cm3;
        profiles.te_keV[lr] += cell.t_e_keV;
        profiles.ti_keV[lr] += cell.t_i_keV;

        const auto conductivity = dec3d::transport::ComputeThermalConductivity(
            dec3d::transport::ThermalConductivityInput{
                cell.t_e_erg_per_particle,
                cell.t_i_erg_per_particle,
                cell.n_e_cm3,
                cell.n_i_cm3,
                cell.zbar,
                cell.mean_ion_mass_g,
                conductivity_model});
        if (conductivity.success) {
          profiles.kappa_e_cm_inv_s[lr] += conductivity.kappa_e_cm_inv_s;
          profiles.kappa_i_cm_inv_s[lr] += conductivity.kappa_i_cm_inv_s;
        }

        if (descriptor.tau_model == BenchmarkTauModel::thesis_spitzer_eq_5_241) {
          const auto tau = dec3d::state::ComputeElectronIonTau(
              dec3d::state::ElectronIonTauInput{
                  dec3d::state::ElectronIonTauModel::thesis_spitzer_eq_5_241,
                  cell.t_e_erg_per_particle,
                  cell.t_i_erg_per_particle,
                  cell.n_e_cm3,
                  cell.n_i_cm3,
                  cell.zbar,
                  cell.mean_ion_mass_g,
                  0.0});
          if (tau.success) {
            profiles.tau_ei_s[lr] += tau.tau_ei_s;
          }
        }

        angular_count += 1.0;
      }
    }

    if (angular_count > 0.0) {
      profiles.rho_g_cm3[lr] /= angular_count;
      profiles.te_keV[lr] /= angular_count;
      profiles.ti_keV[lr] /= angular_count;
      profiles.kappa_e_cm_inv_s[lr] /= angular_count;
      profiles.kappa_i_cm_inv_s[lr] /= angular_count;
      profiles.tau_ei_s[lr] /= angular_count;
    }
  }
  return profiles;
}

[[nodiscard]] std::vector<double> GatherShellVector(
    MPI_Comm communicator,
    const std::vector<double>& local_values) {
  int rank = 0;
  int rank_count = 1;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &rank_count);

  const int local_count = static_cast<int>(local_values.size());
  std::vector<int> counts(static_cast<std::size_t>(rank_count), 0);
  MPI_Gather(
      &local_count,
      1,
      MPI_INT,
      counts.data(),
      1,
      MPI_INT,
      0,
      communicator);

  std::vector<int> displacements(static_cast<std::size_t>(rank_count), 0);
  int total_count = 0;
  if (rank == 0) {
    for (int r = 0; r < rank_count; ++r) {
      displacements[static_cast<std::size_t>(r)] = total_count;
      total_count += counts[static_cast<std::size_t>(r)];
    }
  }

  std::vector<double> gathered;
  if (rank == 0) {
    gathered.resize(static_cast<std::size_t>(total_count), 0.0);
  }
  MPI_Gatherv(
      local_values.empty() ? nullptr : local_values.data(),
      local_count,
      MPI_DOUBLE,
      gathered.empty() ? nullptr : gathered.data(),
      counts.data(),
      displacements.data(),
      MPI_DOUBLE,
      0,
      communicator);
  return gathered;
}

void WriteTwoColumnProfile(
    const std::filesystem::path& path,
    const char* value_name,
    const std::vector<double>& r_cm,
    const std::vector<double>& values) {
  std::ofstream out(path);
  out << std::setprecision(17)
      << "# r_cm " << value_name << '\n';
  for (std::size_t i = 0; i < r_cm.size() && i < values.size(); ++i) {
    out << r_cm[i] << ' ' << values[i] << '\n';
  }
}

[[nodiscard]] bool WriteArtifacts(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalReferencePolicy& reference_policy,
    const WooThermalBenchmarkRunOptions& options,
    const dec3d::state::CanonicalState& initial_state,
    const dec3d::transport::DistributedThermalConductionProblem& problem,
    const dec3d::transport::DistributedThermalConductionResult& thermal,
    const std::string& benchmark_report) noexcept {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);

  const double local_initial_energy =
      LocalThermalEnergy(initial_state, problem.ownership, problem.global_geometry);
  const double local_final_energy =
      LocalThermalEnergy(problem.local_state, problem.ownership, problem.global_geometry);
  double global_initial_energy = 0.0;
  double global_final_energy = 0.0;
  MPI_Reduce(
      &local_initial_energy,
      &global_initial_energy,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      communicator);
  MPI_Reduce(
      &local_final_energy,
      &global_final_energy,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      0,
      communicator);

  const auto local_profiles = BuildLocalShellProfiles(problem, descriptor);
  const auto r_cm = GatherShellVector(communicator, local_profiles.r_cm);
  const auto rho = GatherShellVector(communicator, local_profiles.rho_g_cm3);
  const auto te = GatherShellVector(communicator, local_profiles.te_keV);
  const auto ti = GatherShellVector(communicator, local_profiles.ti_keV);
  const auto kappa_e = GatherShellVector(communicator, local_profiles.kappa_e_cm_inv_s);
  const auto kappa_i = GatherShellVector(communicator, local_profiles.kappa_i_cm_inv_s);
  const auto tau = GatherShellVector(communicator, local_profiles.tau_ei_s);

  int artifact_ok = 1;
  if (rank == 0) {
    try {
      const auto manifest = BuildWooThermalArtifactManifest(descriptor, options.output_root);
      const std::filesystem::path output_dir(manifest.output_directory);
      std::filesystem::create_directories(output_dir);

      {
        std::ofstream out(output_dir / "dec3d.out");
        out << std::setprecision(17)
            << "step=1\n"
            << "time_s=" << descriptor.comparison_time_s << '\n'
            << "pending=false\n";
      }
      {
        std::ofstream out(output_dir / "benchmark_summary.md");
        out << "# P2-8 Woo Thermal Benchmark\n\n"
            << "- case_id: `" << descriptor.case_id << "`\n"
            << "- woo_figure: `" << descriptor.woo_figure << "`\n"
            << "- boundary_model: `" << ToString(descriptor.boundary_model) << "`\n"
            << "- lilac_reference_available: `"
            << (reference_policy.lilac_reference_available ? "true" : "false") << "`\n"
            << "- parity_claim_allowed: `"
            << (reference_policy.parity_claim_allowed ? "true" : "false") << "`\n";
      }
      {
        std::ofstream out(output_dir / "benchmark_diagnostics.txt");
        out << benchmark_report << '\n'
            << thermal.report_line << '\n';
      }
      {
        std::ofstream out(output_dir / "energy_budget_vs_time.txt");
        out << std::setprecision(17)
            << "# time_s thermal_energy_erg\n"
            << "0 " << global_initial_energy << '\n'
            << descriptor.comparison_time_s << ' ' << global_final_energy << '\n';
      }
      {
        std::ofstream out(output_dir / "comparison_metrics.txt");
        out << std::setprecision(17)
            << "diagnostic_id=p2.thermal_benchmark.comparison\n"
            << "comparison_mode=" << reference_policy.comparison_mode << '\n'
            << "lilac_reference_available="
            << (reference_policy.lilac_reference_available ? "true" : "false") << '\n'
            << "parity_claim_allowed="
            << (reference_policy.parity_claim_allowed ? "true" : "false") << '\n'
            << "comparison_time_s=" << descriptor.comparison_time_s << '\n';
      }
      {
        std::ofstream out(output_dir / "profile_step000001.txt");
        out << std::setprecision(17)
            << "# r_cm rho_g_cm3 Te_keV Ti_keV kappa_e_cm_inv_s kappa_i_cm_inv_s tau_ei_s\n";
        for (std::size_t i = 0; i < r_cm.size(); ++i) {
          out << r_cm[i] << ' '
              << (i < rho.size() ? rho[i] : 0.0) << ' '
              << (i < te.size() ? te[i] : 0.0) << ' '
              << (i < ti.size() ? ti[i] : 0.0) << ' '
              << (i < kappa_e.size() ? kappa_e[i] : 0.0) << ' '
              << (i < kappa_i.size() ? kappa_i[i] : 0.0) << ' '
              << (i < tau.size() ? tau[i] : 0.0) << '\n';
        }
      }
      WriteTwoColumnProfile(output_dir / "te_profile_step000001.txt", "Te_keV", r_cm, te);
      WriteTwoColumnProfile(output_dir / "ti_profile_step000001.txt", "Ti_keV", r_cm, ti);
      WriteTwoColumnProfile(output_dir / "rho_profile_step000001.txt", "rho_g_cm3", r_cm, rho);
      WriteTwoColumnProfile(
          output_dir / "kappa_e_profile_step000001.txt",
          "kappa_e_cm_inv_s",
          r_cm,
          kappa_e);
      WriteTwoColumnProfile(
          output_dir / "kappa_i_profile_step000001.txt",
          "kappa_i_cm_inv_s",
          r_cm,
          kappa_i);
      WriteTwoColumnProfile(output_dir / "tau_ei_profile_step000001.txt", "tau_ei_s", r_cm, tau);
    } catch (...) {
      artifact_ok = 0;
    }
  }
  MPI_Bcast(&artifact_ok, 1, MPI_INT, 0, communicator);
  return artifact_ok != 0;
}

struct ExactConstants {
  double rho0_g_cm3{1.0};
  double te_bar_erg{dec3d::physics::ErgFromKeV(5.0)};
  double te_amp_erg{dec3d::physics::ErgFromKeV(0.5)};
  double ti_bar_erg{dec3d::physics::ErgFromKeV(3.0)};
  double ti_amp_erg{dec3d::physics::ErgFromKeV(0.3)};
  double kappa_e_cm_inv_s{1.0e8};
  double kappa_i_cm_inv_s{2.0e7};
  double beta{2.081575977818};
  int mode_l{1};
  int mode_m{1};
};

struct ExactErrorSummary {
  double te_l1_keV{0.0};
  double te_linf_keV{0.0};
  double ti_l1_keV{0.0};
  double ti_linf_keV{0.0};
  double e_electron_l1_erg_cm3{0.0};
  double e_fluid_total_l1_erg_cm3{0.0};
};

struct ExactProfileSnapshot {
  double time_s{0.0};
  std::vector<double> r_cm;
  std::vector<double> te_numeric_keV;
  std::vector<double> te_exact_keV;
  std::vector<double> ti_numeric_keV;
  std::vector<double> ti_exact_keV;
};

struct ExactSliceSnapshot {
  double time_s{0.0};
  std::vector<double> x_um;
  std::vector<double> z_um;
  std::vector<double> te_numeric_keV;
  std::vector<double> te_exact_keV;
  std::vector<double> ti_numeric_keV;
  std::vector<double> ti_exact_keV;
};

struct ExactRunData {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
  std::string report_line;
  ExactErrorSummary errors;
  bool origin_remap_used{false};
  bool pole_remap_used{false};
  std::size_t global_radial_seam_coupling_count{0};
  std::vector<ExactProfileSnapshot> profiles;
  std::vector<ExactSliceSnapshot> slices;
};

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy FullSphereExactBoundary()
    noexcept {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::periodic};
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildFullSphereGeometry(
    double r_max_cm,
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  geometry.theta_faces.resize(theta_cells + 1u);
  geometry.phi_faces.resize(phi_cells + 1u);
  geometry.cell_volumes.resize(radial_cells * theta_cells * phi_cells, 0.0);

  for (std::size_t r = 0; r <= radial_cells; ++r) {
    geometry.radial_faces[r] =
        r_max_cm * static_cast<double>(r) / static_cast<double>(radial_cells);
  }
  for (std::size_t theta = 0; theta <= theta_cells; ++theta) {
    geometry.theta_faces[theta] =
        kPi * static_cast<double>(theta) / static_cast<double>(theta_cells);
  }
  for (std::size_t phi = 0; phi <= phi_cells; ++phi) {
    geometry.phi_faces[phi] =
        2.0 * kPi * static_cast<double>(phi) / static_cast<double>(phi_cells);
  }

  geometry.global_volume = 0.0;
  std::size_t linear = 0;
  for (std::size_t r = 0; r < radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) /
        3.0;
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      const double polar_factor =
          std::cos(geometry.theta_faces[theta]) -
          std::cos(geometry.theta_faces[theta + 1u]);
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const double azimuthal_factor =
            geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
        geometry.cell_volumes[linear] =
            radial_factor * polar_factor * azimuthal_factor;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

[[nodiscard]] double SphericalBesselJ1(double x) noexcept {
  if (std::abs(x) < 1.0e-6) {
    return x / 3.0 - x * x * x / 30.0;
  }
  return std::sin(x) / (x * x) - std::cos(x) / x;
}

[[nodiscard]] double RadialJ1Primitive(double alpha, double r) noexcept {
  if (r == 0.0) {
    return -2.0 / (alpha * alpha * alpha);
  }
  return -(r * std::sin(alpha * r)) / (alpha * alpha) -
         (2.0 * std::cos(alpha * r)) / (alpha * alpha * alpha);
}

[[nodiscard]] double RadialJ0Primitive(double alpha, double r) noexcept {
  if (r == 0.0) {
    return 0.0;
  }
  return -(r * std::cos(alpha * r)) / (alpha * alpha) +
         std::sin(alpha * r) / (alpha * alpha * alpha);
}

[[nodiscard]] double ExactCellAverageModeShape(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    double r_max_cm,
    std::size_t global_radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double alpha = constants.beta / r_max_cm;
  const double r0 = geometry.radial_faces[global_radial];
  const double r1 = geometry.radial_faces[global_radial + 1u];
  const double radial_factor = (std::pow(r1, 3) - std::pow(r0, 3)) / 3.0;
  if (constants.mode_l == 0) {
    const double radial_integral =
        RadialJ0Primitive(alpha, r1) - RadialJ0Primitive(alpha, r0);
    return radial_factor > 0.0 ? radial_integral / radial_factor : 0.0;
  }

  const double radial_integral =
      RadialJ1Primitive(alpha, r1) - RadialJ1Primitive(alpha, r0);

  const double theta0 = geometry.theta_faces[theta];
  const double theta1 = geometry.theta_faces[theta + 1u];
  const double phi0 = geometry.phi_faces[phi];
  const double phi1 = geometry.phi_faces[phi + 1u];
  const double theta_integral =
      0.5 * (theta1 - theta0) -
      0.25 * (std::sin(2.0 * theta1) - std::sin(2.0 * theta0));
  const double phi_integral = std::sin(phi1) - std::sin(phi0);
  const double solid_angle = (std::cos(theta0) - std::cos(theta1)) * (phi1 - phi0);
  if (radial_factor <= 0.0 || solid_angle == 0.0) {
    return 0.0;
  }
  return (radial_integral * theta_integral * phi_integral) /
         (radial_factor * solid_angle);
}

[[nodiscard]] ExactConstants MakeExactConstants(
    const WooThermalBenchmarkDescriptor& descriptor) noexcept {
  ExactConstants constants;
  if (descriptor.case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
    constants.beta = 4.4934094579090642;
    constants.mode_l = 0;
    constants.mode_m = 0;
  }
  return constants;
}

[[nodiscard]] const char* ExactSolutionName(const ExactConstants& constants) noexcept {
  return constants.mode_l == 0
             ? "l0_spherical_bessel_neumann"
             : "l1m1_spherical_bessel_neumann";
}

[[nodiscard]] double ExactTemperatureErg(
    bool electron,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    double r_max_cm,
    double number_density,
    std::size_t global_radial,
    std::size_t theta,
    std::size_t phi,
    double time_s) noexcept {
  const double kappa = electron ? constants.kappa_e_cm_inv_s : constants.kappa_i_cm_inv_s;
  const double chi =
      (dec3d::state::HydroIdealGasGamma() - 1.0) * kappa / number_density;
  const double decay =
      std::exp(-chi * (constants.beta / r_max_cm) * (constants.beta / r_max_cm) * time_s);
  const double base = electron ? constants.te_bar_erg : constants.ti_bar_erg;
  const double amplitude = electron ? constants.te_amp_erg : constants.ti_amp_erg;
  return base +
         amplitude *
             ExactCellAverageModeShape(geometry, constants, r_max_cm, global_radial, theta, phi) *
             decay;
}

[[nodiscard]] std::pair<std::size_t, std::size_t> SelectPositiveModeProfileLine(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  double best_value = -std::numeric_limits<double>::infinity();
  std::pair<std::size_t, std::size_t> best{0u, 0u};
  for (std::size_t theta = 0; theta < theta_cells; ++theta) {
    const double theta_center =
        0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
    for (std::size_t phi = 0; phi < phi_cells; ++phi) {
      const double phi_center =
          0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
      const double value = std::sin(theta_center) * std::cos(phi_center);
      if (value > best_value) {
        best_value = value;
        best = {theta, phi};
      }
    }
  }
  return best;
}

[[nodiscard]] dec3d::transport::DistributedGenericDiffusionProblem BuildExactScalarProblem(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::transport::GenericDiffusionBoundaryPolicy& boundary,
    const dec3d::core::Array3D<double>& scalar_old,
    double coefficient_a,
    double coefficient_d,
    double dt_s,
    dec3d::transport::DiffusionPoleMetricMode pole_metric_mode) {
  dec3d::transport::DistributedGenericDiffusionProblem problem;
  problem.ownership = ownership;
  problem.global_geometry = geometry;
  problem.dt_s = dt_s;
  problem.boundary_policy = boundary;
  problem.pole_metric_mode = pole_metric_mode;
  const std::size_t local_radial =
      ownership.global_radial_end - ownership.global_radial_begin;
  problem.local_coefficient_A = dec3d::core::Array3D<double>(
      local_radial,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      coefficient_a);
  problem.local_coefficient_D = dec3d::core::Array3D<double>(
      local_radial,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      coefficient_d);
  problem.local_coefficient_C = dec3d::core::Array3D<double>(
      local_radial,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      0.0);
  problem.local_coefficient_B = dec3d::core::Array3D<double>(
      local_radial,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      0.0);
  problem.local_scalar_old = scalar_old;
  return problem;
}

void FillExactInitialScalars(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    double ne,
    double ni,
    dec3d::core::Array3D<double>& te,
    dec3d::core::Array3D<double>& ti) {
  for (std::size_t lr = 0; lr < te.extent_r(); ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t theta = 0; theta < ownership.global_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < ownership.global_phi_cells; ++phi) {
        te(lr, theta, phi) = ExactTemperatureErg(
            true,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ne,
            gr,
            theta,
            phi,
            0.0);
        ti(lr, theta, phi) = ExactTemperatureErg(
            false,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ni,
            gr,
            theta,
            phi,
            0.0);
      }
    }
  }
}

void CopySolveToScalar(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const std::vector<double>& scalar_new,
    dec3d::core::Array3D<double>& out) {
  for (std::size_t lr = 0; lr < out.extent_r(); ++lr) {
    for (std::size_t theta = 0; theta < ownership.global_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < ownership.global_phi_cells; ++phi) {
        const std::size_t local_index =
            ((lr * ownership.global_theta_cells) + theta) * ownership.global_phi_cells + phi;
        out(lr, theta, phi) = scalar_new.at(local_index);
      }
    }
  }
}

[[nodiscard]] ExactErrorSummary ComputeExactErrors(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    const dec3d::core::Array3D<double>& te,
    const dec3d::core::Array3D<double>& ti,
    double ne,
    double ni,
    double time_s) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  double local_volume = 0.0;
  double local_te_l1 = 0.0;
  double local_ti_l1 = 0.0;
  double local_ee_l1 = 0.0;
  double local_et_l1 = 0.0;
  double local_te_linf = 0.0;
  double local_ti_linf = 0.0;

  for (std::size_t lr = 0; lr < te.extent_r(); ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t theta = 0; theta < ownership.global_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < ownership.global_phi_cells; ++phi) {
        const double volume = CellVolume(
            geometry,
            gr,
            theta,
            phi,
            ownership.global_theta_cells,
            ownership.global_phi_cells);
        const double exact_te = ExactTemperatureErg(
            true,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ne,
            gr,
            theta,
            phi,
            time_s);
        const double exact_ti = ExactTemperatureErg(
            false,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ni,
            gr,
            theta,
            phi,
            time_s);
        const double te_abs_keV = dec3d::physics::KeVFromErg(std::abs(te(lr, theta, phi) - exact_te));
        const double ti_abs_keV = dec3d::physics::KeVFromErg(std::abs(ti(lr, theta, phi) - exact_ti));
        const double ee_error =
            std::abs(ne * (te(lr, theta, phi) - exact_te) / gamma_minus_one);
        const double ion_error =
            std::abs(ni * (ti(lr, theta, phi) - exact_ti) / gamma_minus_one);
        local_volume += volume;
        local_te_l1 += volume * te_abs_keV;
        local_ti_l1 += volume * ti_abs_keV;
        local_ee_l1 += volume * ee_error;
        local_et_l1 += volume * (ee_error + ion_error);
        local_te_linf = std::max(local_te_linf, te_abs_keV);
        local_ti_linf = std::max(local_ti_linf, ti_abs_keV);
      }
    }
  }

  double global_volume = 0.0;
  double global_te_l1 = 0.0;
  double global_ti_l1 = 0.0;
  double global_ee_l1 = 0.0;
  double global_et_l1 = 0.0;
  double global_te_linf = 0.0;
  double global_ti_linf = 0.0;
  MPI_Allreduce(&local_volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM, communicator);
  MPI_Allreduce(&local_te_l1, &global_te_l1, 1, MPI_DOUBLE, MPI_SUM, communicator);
  MPI_Allreduce(&local_ti_l1, &global_ti_l1, 1, MPI_DOUBLE, MPI_SUM, communicator);
  MPI_Allreduce(&local_ee_l1, &global_ee_l1, 1, MPI_DOUBLE, MPI_SUM, communicator);
  MPI_Allreduce(&local_et_l1, &global_et_l1, 1, MPI_DOUBLE, MPI_SUM, communicator);
  MPI_Allreduce(&local_te_linf, &global_te_linf, 1, MPI_DOUBLE, MPI_MAX, communicator);
  MPI_Allreduce(&local_ti_linf, &global_ti_linf, 1, MPI_DOUBLE, MPI_MAX, communicator);

  ExactErrorSummary errors;
  const double inv_volume = global_volume > 0.0 ? 1.0 / global_volume : 0.0;
  errors.te_l1_keV = global_te_l1 * inv_volume;
  errors.ti_l1_keV = global_ti_l1 * inv_volume;
  errors.e_electron_l1_erg_cm3 = global_ee_l1 * inv_volume;
  errors.e_fluid_total_l1_erg_cm3 = global_et_l1 * inv_volume;
  errors.te_linf_keV = global_te_linf;
  errors.ti_linf_keV = global_ti_linf;
  return errors;
}

[[nodiscard]] ExactProfileSnapshot BuildExactProfileSnapshot(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    const dec3d::core::Array3D<double>& te,
    const dec3d::core::Array3D<double>& ti,
    double ne,
    double ni,
    double time_s) {
  const auto [theta_line, phi_line] = SelectPositiveModeProfileLine(
      geometry,
      ownership.global_theta_cells,
      ownership.global_phi_cells);

  std::vector<double> local_r;
  std::vector<double> local_te;
  std::vector<double> local_te_exact;
  std::vector<double> local_ti;
  std::vector<double> local_ti_exact;
  local_r.reserve(te.extent_r());
  local_te.reserve(te.extent_r());
  local_te_exact.reserve(te.extent_r());
  local_ti.reserve(te.extent_r());
  local_ti_exact.reserve(te.extent_r());
  for (std::size_t lr = 0; lr < te.extent_r(); ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    local_r.push_back(0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]));
    local_te.push_back(dec3d::physics::KeVFromErg(te(lr, theta_line, phi_line)));
    local_ti.push_back(dec3d::physics::KeVFromErg(ti(lr, theta_line, phi_line)));
    local_te_exact.push_back(dec3d::physics::KeVFromErg(ExactTemperatureErg(
        true,
        geometry,
        constants,
        geometry.radial_faces.back(),
        ne,
        gr,
        theta_line,
        phi_line,
        time_s)));
    local_ti_exact.push_back(dec3d::physics::KeVFromErg(ExactTemperatureErg(
        false,
        geometry,
        constants,
        geometry.radial_faces.back(),
        ni,
        gr,
        theta_line,
        phi_line,
        time_s)));
  }

  ExactProfileSnapshot snapshot;
  snapshot.time_s = time_s;
  snapshot.r_cm = GatherShellVector(communicator, local_r);
  snapshot.te_numeric_keV = GatherShellVector(communicator, local_te);
  snapshot.te_exact_keV = GatherShellVector(communicator, local_te_exact);
  snapshot.ti_numeric_keV = GatherShellVector(communicator, local_ti);
  snapshot.ti_exact_keV = GatherShellVector(communicator, local_ti_exact);
  return snapshot;
}

[[nodiscard]] ExactSliceSnapshot BuildExactSliceSnapshot(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ExactConstants& constants,
    const dec3d::core::Array3D<double>& te,
    const dec3d::core::Array3D<double>& ti,
    double ne,
    double ni,
    double time_s) {
  const std::size_t phi_zero = 0u;
  const std::size_t phi_pi = ownership.global_phi_cells / 2u;
  const std::size_t phi_count = te.extent_phi();
  std::vector<double> local_x;
  std::vector<double> local_z;
  std::vector<double> local_te;
  std::vector<double> local_te_exact;
  std::vector<double> local_ti;
  std::vector<double> local_ti_exact;
  local_x.reserve(te.extent_r() * te.extent_theta() * 2u);
  local_z.reserve(te.extent_r() * te.extent_theta() * 2u);
  local_te.reserve(te.extent_r() * te.extent_theta() * 2u);
  local_te_exact.reserve(te.extent_r() * te.extent_theta() * 2u);
  local_ti.reserve(te.extent_r() * te.extent_theta() * 2u);
  local_ti_exact.reserve(te.extent_r() * te.extent_theta() * 2u);

  const std::size_t phi_indices[2] = {phi_zero, phi_pi};
  for (std::size_t lr = 0; lr < te.extent_r(); ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    const double radius =
        0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]);
    for (std::size_t theta = 0; theta < te.extent_theta(); ++theta) {
      const double theta_center =
          0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
      for (const std::size_t phi_candidate : phi_indices) {
        const std::size_t phi = phi_candidate < phi_count ? phi_candidate : phi_count - 1u;
        const double phi_center =
            0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
        local_x.push_back(radius * std::sin(theta_center) * std::cos(phi_center) * 1.0e4);
        local_z.push_back(radius * std::cos(theta_center) * 1.0e4);
        local_te.push_back(dec3d::physics::KeVFromErg(te(lr, theta, phi)));
        local_ti.push_back(dec3d::physics::KeVFromErg(ti(lr, theta, phi)));
        local_te_exact.push_back(dec3d::physics::KeVFromErg(ExactTemperatureErg(
            true,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ne,
            gr,
            theta,
            phi,
            time_s)));
        local_ti_exact.push_back(dec3d::physics::KeVFromErg(ExactTemperatureErg(
            false,
            geometry,
            constants,
            geometry.radial_faces.back(),
            ni,
            gr,
            theta,
            phi,
            time_s)));
      }
    }
  }

  ExactSliceSnapshot snapshot;
  snapshot.time_s = time_s;
  snapshot.x_um = GatherShellVector(communicator, local_x);
  snapshot.z_um = GatherShellVector(communicator, local_z);
  snapshot.te_numeric_keV = GatherShellVector(communicator, local_te);
  snapshot.te_exact_keV = GatherShellVector(communicator, local_te_exact);
  snapshot.ti_numeric_keV = GatherShellVector(communicator, local_ti);
  snapshot.ti_exact_keV = GatherShellVector(communicator, local_ti_exact);
  return snapshot;
}

[[nodiscard]] bool WriteExactArtifacts(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalBenchmarkRunOptions& options,
    const ExactRunData& data) noexcept {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  int artifact_ok = 1;
  if (rank == 0) {
    try {
      const auto manifest = BuildWooThermalArtifactManifest(descriptor, options.output_root);
      const std::filesystem::path output_dir(manifest.output_directory);
      std::filesystem::create_directories(output_dir);
      {
        std::ofstream out(output_dir / "dec3d.out");
        out << std::setprecision(17)
            << "step=" << data.profiles.size() << '\n'
            << "time_s=" << descriptor.comparison_time_s << '\n'
            << "pending=false\n";
      }
      {
        std::ofstream out(output_dir / "exact_benchmark_summary.md");
        out << "# P2-8e Exact Spherical Harmonic Thermal Benchmark\n\n"
            << "- exact_solution: `l1m1_spherical_bessel_neumann`\n"
            << "- boundary_model: `full_sphere_scalar_remap`\n"
            << "- te_l1_keV: `" << data.errors.te_l1_keV << "`\n"
            << "- ti_l1_keV: `" << data.errors.ti_l1_keV << "`\n";
      }
      {
        std::ofstream out(output_dir / "benchmark_diagnostics.txt");
        out << data.report_line << '\n';
      }
      for (std::size_t step = 0; step < data.profiles.size(); ++step) {
        const auto& profile = data.profiles[step];
        std::ostringstream name;
        name << "profile_exact_vs_numeric_step"
             << std::setw(6) << std::setfill('0') << (step + 1u) << ".txt";
        std::ofstream out(output_dir / name.str());
        out << std::setprecision(17)
            << "# time_s " << profile.time_s << '\n'
            << "# r_cm Te_numeric_keV Te_exact_keV Ti_numeric_keV Ti_exact_keV\n";
        for (std::size_t i = 0; i < profile.r_cm.size(); ++i) {
          out << profile.r_cm[i] << ' '
              << (i < profile.te_numeric_keV.size() ? profile.te_numeric_keV[i] : 0.0) << ' '
              << (i < profile.te_exact_keV.size() ? profile.te_exact_keV[i] : 0.0) << ' '
              << (i < profile.ti_numeric_keV.size() ? profile.ti_numeric_keV[i] : 0.0) << ' '
              << (i < profile.ti_exact_keV.size() ? profile.ti_exact_keV[i] : 0.0) << '\n';
        }
      }
      for (std::size_t step = 0; step < data.slices.size(); ++step) {
        const auto& slice = data.slices[step];
        std::ostringstream name;
        name << "slice_xz_exact_vs_numeric_step"
             << std::setw(6) << std::setfill('0') << (step + 1u) << ".txt";
        std::ofstream out(output_dir / name.str());
        out << std::setprecision(17)
            << "# time_s " << slice.time_s << '\n'
            << "# x_um z_um Te_numeric_keV Te_exact_keV Ti_numeric_keV Ti_exact_keV\n";
        for (std::size_t i = 0; i < slice.x_um.size(); ++i) {
          out << slice.x_um[i] << ' '
              << (i < slice.z_um.size() ? slice.z_um[i] : 0.0) << ' '
              << (i < slice.te_numeric_keV.size() ? slice.te_numeric_keV[i] : 0.0) << ' '
              << (i < slice.te_exact_keV.size() ? slice.te_exact_keV[i] : 0.0) << ' '
              << (i < slice.ti_numeric_keV.size() ? slice.ti_numeric_keV[i] : 0.0) << ' '
              << (i < slice.ti_exact_keV.size() ? slice.ti_exact_keV[i] : 0.0) << '\n';
        }
      }
    } catch (...) {
      artifact_ok = 0;
    }
  }
  MPI_Bcast(&artifact_ok, 1, MPI_INT, 0, communicator);
  return artifact_ok != 0;
}

[[nodiscard]] ExactRunData RunExactL1M1Benchmark(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalBenchmarkRunOptions& options) noexcept {
  ExactRunData data;
  const ExactConstants constants = MakeExactConstants(descriptor);
  const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
      communicator,
      options.global_radial_cells,
      options.theta_cells,
      options.phi_cells);
  if (!ownership.success) {
    data.failure_reason = "exact benchmark ownership failed";
    data.failure_diagnostics = ownership.failure_diagnostics;
    return data;
  }

  const auto geometry = BuildFullSphereGeometry(
      descriptor.r_max_cm,
      options.global_radial_cells,
      options.theta_cells,
      options.phi_cells);
  const auto boundary = FullSphereExactBoundary();
  const std::size_t local_radial =
      ownership.global_radial_end - ownership.global_radial_begin;
  dec3d::core::Array3D<double> te(
      local_radial,
      options.theta_cells,
      options.phi_cells,
      0.0);
  dec3d::core::Array3D<double> ti(
      local_radial,
      options.theta_cells,
      options.phi_cells,
      0.0);
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  const double ni = constants.rho0_g_cm3 / mean_ion_mass;
  const double ne = dec3d::physics::DefaultZbar() * ni;
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  FillExactInitialScalars(ownership, geometry, constants, ne, ni, te, ti);

  double previous_time = 0.0;
  dec3d::transport::GenericDiffusionHypreSolveOptions solve_options;
  solve_options.communicator = communicator;
  solve_options.relative_tolerance = 1.0e-11;
  solve_options.max_iterations = 300;

  for (const double output_time : descriptor.output_times_s) {
    const double dt_s = output_time - previous_time;
    if (!std::isfinite(dt_s) || dt_s <= 0.0) {
      data.failure_reason = "exact benchmark output_times_s must be strictly increasing";
      data.failure_diagnostics = BuildFailure(data.failure_reason.c_str(), descriptor);
      return data;
    }

    auto electron_problem = BuildExactScalarProblem(
        ownership,
        geometry,
        boundary,
        te,
        ne / gamma_minus_one,
        constants.kappa_e_cm_inv_s,
        dt_s,
        options.pole_metric_mode);
    auto ion_problem = BuildExactScalarProblem(
        ownership,
        geometry,
        boundary,
        ti,
        ni / gamma_minus_one,
        constants.kappa_i_cm_inv_s,
        dt_s,
        options.pole_metric_mode);
    const auto electron_assembly =
        dec3d::transport::AssembleDistributedGenericDiffusionSystem(electron_problem);
    if (!electron_assembly.success) {
      data.failure_reason = "electron exact diffusion assembly failed";
      data.failure_diagnostics = electron_assembly.failure_diagnostics;
      return data;
    }
    const auto ion_assembly =
        dec3d::transport::AssembleDistributedGenericDiffusionSystem(ion_problem);
    if (!ion_assembly.success) {
      data.failure_reason = "ion exact diffusion assembly failed";
      data.failure_diagnostics = ion_assembly.failure_diagnostics;
      return data;
    }

    const auto electron_solve =
        dec3d::transport::SolveDistributedGenericDiffusionHypre(electron_assembly, solve_options);
    if (!electron_solve.success) {
      data.failure_reason = "electron exact diffusion solve failed";
      data.failure_diagnostics = electron_solve.failure_diagnostics;
      return data;
    }
    const auto ion_solve =
        dec3d::transport::SolveDistributedGenericDiffusionHypre(ion_assembly, solve_options);
    if (!ion_solve.success) {
      data.failure_reason = "ion exact diffusion solve failed";
      data.failure_diagnostics = ion_solve.failure_diagnostics;
      return data;
    }

    CopySolveToScalar(ownership, electron_solve.local_scalar_new, te);
    CopySolveToScalar(ownership, ion_solve.local_scalar_new, ti);
    data.origin_remap_used = data.origin_remap_used ||
                             electron_assembly.global_origin_remap_used ||
                             ion_assembly.global_origin_remap_used;
    data.pole_remap_used = data.pole_remap_used ||
                           electron_assembly.global_pole_remap_used ||
                           ion_assembly.global_pole_remap_used;
    data.global_radial_seam_coupling_count +=
        electron_assembly.global_radial_seam_coupling_count +
        ion_assembly.global_radial_seam_coupling_count;
    data.profiles.push_back(BuildExactProfileSnapshot(
        communicator,
        ownership,
        geometry,
        constants,
        te,
        ti,
        ne,
        ni,
        output_time));
    data.slices.push_back(BuildExactSliceSnapshot(
        communicator,
        ownership,
        geometry,
        constants,
        te,
        ti,
        ne,
        ni,
        output_time));
    previous_time = output_time;
  }

  data.errors = ComputeExactErrors(
      communicator,
      ownership,
      geometry,
      constants,
      te,
      ti,
      ne,
      ni,
      descriptor.comparison_time_s);

  std::ostringstream out;
  out << std::setprecision(17)
      << "; exact_solution=" << ExactSolutionName(constants)
      << "; exact_projection=cell_volume_average"
      << "; exact_beta=" << constants.beta
      << "; exact_kappa_e_cm_inv_s=" << constants.kappa_e_cm_inv_s
      << "; exact_kappa_i_cm_inv_s=" << constants.kappa_i_cm_inv_s
      << "; exact_rho0_g_cm3=" << constants.rho0_g_cm3
      << "; exact_te_bar_keV=" << dec3d::physics::KeVFromErg(constants.te_bar_erg)
      << "; exact_te_amp_keV=" << dec3d::physics::KeVFromErg(constants.te_amp_erg)
      << "; exact_ti_bar_keV=" << dec3d::physics::KeVFromErg(constants.ti_bar_erg)
      << "; exact_ti_amp_keV=" << dec3d::physics::KeVFromErg(constants.ti_amp_erg)
      << "; exact_mode_l=" << constants.mode_l
      << "; exact_mode_m=" << constants.mode_m
      << "; coefficient_time_level=constant"
      << "; pole_metric_mode="
      << (options.pole_metric_mode ==
                  dec3d::transport::DiffusionPoleMetricMode::axis_regular_polar_phi
              ? "axis_regular_polar_phi"
              : "cell_centered_spherical")
      << "; origin_remap_used=" << (data.origin_remap_used ? "true" : "false")
      << "; pole_remap_used=" << (data.pole_remap_used ? "true" : "false")
      << "; radial_seam_coupling_count=" << data.global_radial_seam_coupling_count
      << "; te_l1_keV=" << data.errors.te_l1_keV
      << "; te_linf_keV=" << data.errors.te_linf_keV
      << "; ti_l1_keV=" << data.errors.ti_l1_keV
      << "; ti_linf_keV=" << data.errors.ti_linf_keV
      << "; e_electron_l1_erg_cm3=" << data.errors.e_electron_l1_erg_cm3
      << "; e_fluid_total_l1_erg_cm3=" << data.errors.e_fluid_total_l1_erg_cm3;
  data.report_line = out.str();
  data.success = true;
  return data;
}

}  // namespace

WooThermalBenchmarkRunResult RunWooThermalBenchmarkCase(
    MPI_Comm communicator,
    const WooThermalBenchmarkDescriptor& descriptor,
    const WooThermalReferencePolicy& reference_policy,
    const WooThermalBenchmarkRunOptions& options) noexcept {
  WooThermalBenchmarkRunResult result;
  const auto descriptor_validation = ValidateWooThermalBenchmarkDescriptor(descriptor);
  if (!descriptor_validation.success) {
    result.failure_reason = descriptor_validation.failure_reason;
    result.failure_diagnostics = descriptor_validation.failure_diagnostics;
    return result;
  }

  if (descriptor.case_kind == WooThermalBenchmarkCase::exact_l1m1_constant_kappa ||
      descriptor.case_kind == WooThermalBenchmarkCase::exact_l0_radial_constant_kappa) {
    const auto exact = RunExactL1M1Benchmark(communicator, descriptor, options);
    result.thermal_stage_executed = exact.success;
    result.equilibration_stage_executed = false;
    result.momentum_update_executed = false;
    result.output_directory =
        BuildWooThermalArtifactManifest(descriptor, options.output_root).output_directory;
    const auto base = BuildWooThermalBenchmarkDiagnosticsLine(descriptor, reference_policy);
    std::ostringstream out;
    out << std::setprecision(17)
        << base
        << "; thermal_stage_executed=" << (exact.success ? "true" : "false")
        << "; equilibration_stage_executed=false"
        << "; momentum_update_executed=false"
        << "; thermal_stage_report_present=" << (exact.success ? "true" : "false");
    if (exact.success) {
      out << exact.report_line;
    }
    if (!exact.success) {
      result.failure_reason = exact.failure_reason;
      result.failure_diagnostics = exact.failure_diagnostics;
      result.report_line = out.str();
      return result;
    }
    if (options.write_artifacts) {
      const bool artifacts_written =
          WriteExactArtifacts(communicator, descriptor, options, exact);
      result.profile_artifacts_written = artifacts_written;
      result.benchmark_summary_written = artifacts_written;
      result.benchmark_diagnostics_written = artifacts_written;
      result.comparison_metrics_written = artifacts_written;
      result.native_profiles_written = artifacts_written;
      result.energy_budget_written = false;
      out << "; artifact_gather_used=true"
          << "; artifact_gather_is_not_production_writeback=true"
          << "; profile_artifacts_written=" << (artifacts_written ? "true" : "false")
          << "; energy_budget_written=false"
          << "; benchmark_summary_written=" << (artifacts_written ? "true" : "false")
          << "; benchmark_diagnostics_written=" << (artifacts_written ? "true" : "false")
          << "; comparison_metrics_written=" << (artifacts_written ? "true" : "false")
          << "; native_profiles_written=" << (artifacts_written ? "true" : "false");
      if (!artifacts_written) {
        result.failure_reason = "exact artifact writing failed";
        result.failure_diagnostics = BuildFailure(result.failure_reason.c_str(), descriptor);
        result.report_line = out.str();
        return result;
      }
    } else {
      result.profile_artifacts_written = false;
      result.energy_budget_written = false;
      result.benchmark_summary_written = false;
      result.benchmark_diagnostics_written = false;
      result.comparison_metrics_written = false;
      result.native_profiles_written = false;
      out << "; artifact_gather_used=false"
          << "; artifact_gather_is_not_production_writeback=true"
          << "; profile_artifacts_written=false"
          << "; energy_budget_written=false";
    }
    result.report_line = out.str();
    result.success = ValidateWooThermalBenchmarkDiagnostics(result.report_line);
    if (!result.success) {
      result.failure_reason = "exact benchmark diagnostics incomplete";
      result.failure_diagnostics = BuildFailure(result.failure_reason.c_str(), descriptor);
    }
    return result;
  }

  auto problem = BuildBenchmarkProblem(communicator, descriptor, options);
  const auto initial_state = problem.local_state;
  if (descriptor.case_kind == WooThermalBenchmarkCase::ei_0d_thesis_spitzer) {
    auto equilibration_options = dec3d::state::ElectronIonEquilibrationOptions{
        descriptor.comparison_time_s,
        dec3d::state::ElectronIonTauModel::thesis_spitzer_eq_5_241,
        0.0,
        0.0,
        0.0};
    const auto equilibration = dec3d::state::ApplyLocalElectronIonEquilibration(
        problem.local_state,
        equilibration_options);
    if (!equilibration.success) {
      result.failure_reason = "0D electron-ion equilibration failed";
      result.failure_diagnostics = equilibration.failure_diagnostics;
      return result;
    }
    result.thermal_stage_executed = false;
    result.equilibration_stage_executed = true;
    result.momentum_update_executed = false;
    result.output_directory =
        BuildWooThermalArtifactManifest(descriptor, options.output_root).output_directory;
    const auto base = BuildWooThermalBenchmarkDiagnosticsLine(descriptor, reference_policy);
    std::ostringstream out;
    out << std::setprecision(17)
        << base
        << "; stage_order=E"
        << "; equilibration_stage_report_present=true"
        << "; tau_model_requested=thesis_spitzer_eq_5_241"
        << "; tau_model_executed=thesis_spitzer_eq_5_241"
        << "; max_deltaT_before_erg=" << equilibration.max_deltaT_before_erg
        << "; max_deltaT_after_erg=" << equilibration.max_deltaT_after_erg
        << "; max_dt_over_tau=" << equilibration.max_dt_over_tau
        << "; max_expected_decay_residual_erg="
        << equilibration.max_expected_decay_residual_erg
        << "; max_total_thermal_energy_residual="
        << equilibration.max_total_thermal_energy_residual
        << "; artifact_gather_used=false"
        << "; artifact_gather_is_not_production_writeback=true"
        << "; profile_artifacts_written=false"
        << "; energy_budget_written=false";
    result.report_line = out.str();
    result.profile_artifacts_written = false;
    result.energy_budget_written = false;
    result.benchmark_summary_written = false;
    result.benchmark_diagnostics_written = false;
    result.comparison_metrics_written = false;
    result.native_profiles_written = false;
    result.success = ValidateWooThermalBenchmarkDiagnostics(result.report_line);
    if (!result.success) {
      result.failure_reason = "0D e-i benchmark diagnostics incomplete";
      result.failure_diagnostics = BuildFailure(result.failure_reason.c_str(), descriptor);
    }
    return result;
  }

  auto thermal_options = dec3d::transport::DistributedThermalConductionOptions{
      descriptor.comparison_time_s,
      ToTransportKappa(descriptor.thermal_kappa_model),
      dec3d::transport::ElectronFluxLimiterModel::disabled,
      0.0};
  if (descriptor.case_kind == WooThermalBenchmarkCase::flux_limiter_stress) {
    thermal_options.electron_flux_limiter_model =
        dec3d::transport::ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa;
    thermal_options.electron_flux_limiter_alpha_e = 1.0e-3;
  }
  const auto thermal = dec3d::transport::ApplyDistributedVariableKappaThermalConduction(
      problem,
      thermal_options);
  if (!thermal.success) {
    result.failure_reason = "thermal stage failed";
    result.failure_diagnostics = thermal.failure_diagnostics;
    return result;
  }

  result.thermal_stage_executed = true;
  result.equilibration_stage_executed = descriptor.equilibration_stage_enabled;
  result.momentum_update_executed = false;
  result.output_directory =
      BuildWooThermalArtifactManifest(descriptor, options.output_root).output_directory;

  dec3d::state::ElectronIonEquilibrationResult equilibration;
  if (descriptor.equilibration_stage_enabled) {
    equilibration = dec3d::state::ApplyLocalElectronIonEquilibration(
        problem.local_state,
        dec3d::state::ElectronIonEquilibrationOptions{
            descriptor.comparison_time_s,
            dec3d::state::ElectronIonTauModel::thesis_spitzer_eq_5_241,
            0.0,
            0.0,
            0.0});
    if (!equilibration.success) {
      result.failure_reason = "electron-ion equilibration stage failed";
      result.failure_diagnostics = equilibration.failure_diagnostics;
      return result;
    }
  }

  const auto base = BuildWooThermalBenchmarkDiagnosticsLine(descriptor, reference_policy);
  std::ostringstream out;
  out << std::setprecision(17)
      << base
      << "; thermal_stage_executed=true"
      << "; equilibration_stage_executed="
      << (result.equilibration_stage_executed ? "true" : "false")
      << "; momentum_update_executed=false"
      << "; thermal_stage_report_present=true";
  if (descriptor.equilibration_stage_enabled) {
    out << "; stage_order=H,T,E"
        << "; tau_model_requested=thesis_spitzer_eq_5_241"
        << "; tau_model_executed=thesis_spitzer_eq_5_241"
        << "; equilibration_stage_report_present=true"
        << "; max_dt_over_tau=" << equilibration.max_dt_over_tau
        << "; max_expected_decay_residual_erg="
        << equilibration.max_expected_decay_residual_erg;
  }
  if (descriptor.case_kind == WooThermalBenchmarkCase::flux_limiter_stress) {
    out << "; electron_flux_limiter_required=true"
        << "; flux_limiter_report_present=true"
        << "; flux_limiter_report={" << thermal.flux_limiter_report << "}";
  }
  if (options.write_artifacts) {
    const bool artifacts_written = WriteArtifacts(
        communicator,
        descriptor,
        reference_policy,
        options,
        initial_state,
        problem,
        thermal,
        out.str());
    result.profile_artifacts_written = artifacts_written;
    result.energy_budget_written = artifacts_written;
    result.benchmark_summary_written = artifacts_written;
    result.benchmark_diagnostics_written = artifacts_written;
    result.comparison_metrics_written = artifacts_written;
    result.native_profiles_written = artifacts_written;
    out << "; artifact_gather_used=true"
        << "; artifact_gather_is_not_production_writeback=true"
        << "; profile_artifacts_written=" << (artifacts_written ? "true" : "false")
        << "; energy_budget_written=" << (artifacts_written ? "true" : "false")
        << "; benchmark_summary_written=" << (artifacts_written ? "true" : "false")
        << "; benchmark_diagnostics_written=" << (artifacts_written ? "true" : "false")
        << "; comparison_metrics_written=" << (artifacts_written ? "true" : "false")
        << "; native_profiles_written=" << (artifacts_written ? "true" : "false");
    if (!artifacts_written) {
      result.failure_reason = "artifact writing failed";
      result.failure_diagnostics = BuildFailure(result.failure_reason.c_str(), descriptor);
      result.report_line = out.str();
      return result;
    }
  } else {
    result.profile_artifacts_written = false;
    result.energy_budget_written = false;
    result.benchmark_summary_written = false;
    result.benchmark_diagnostics_written = false;
    result.comparison_metrics_written = false;
    result.native_profiles_written = false;
    out << "; artifact_gather_used=false"
        << "; artifact_gather_is_not_production_writeback=true"
        << "; profile_artifacts_written=false"
        << "; energy_budget_written=false";
  }
  result.report_line = out.str();
  result.success = ValidateWooThermalBenchmarkDiagnostics(result.report_line);
  if (!result.success) {
    result.failure_reason = "benchmark diagnostics incomplete";
    result.failure_diagnostics = BuildFailure(result.failure_reason.c_str(), descriptor);
  }
  return result;
}

}  // namespace dec3d::benchmarks
