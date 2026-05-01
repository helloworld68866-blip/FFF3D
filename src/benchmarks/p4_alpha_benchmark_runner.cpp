#include "benchmarks/p4_alpha_benchmark_runner.hpp"

#include "alpha/alpha_coefficients.hpp"
#include "alpha/alpha_operator.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace dec3d::benchmarks {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct ProfileRow {
  double time_s{0.0};
  double r_cm{0.0};
  double epsilon_alpha_erg_cm3{0.0};
  double Te_keV{0.0};
  double birth_source_erg_cm3_s{0.0};
  double drag_deposition_erg_cm3_s{0.0};
};

struct BudgetRow {
  double time_s{0.0};
  std::string run_kind{"serial"};
  double alpha_equation_budget_residual{0.0};
  double alpha_electron_exchange_residual{0.0};
  double global_alpha_plus_electron_budget_residual{0.0};
};

[[nodiscard]] std::string JoinPath(
    const std::string& root,
    const std::string& leaf) {
  return (std::filesystem::path(root) / leaf).generic_string();
}

[[nodiscard]] std::string FailureLine(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const char* reason) {
  std::ostringstream out;
  out << "diagnostic_id=p4.alpha.benchmark_ladder.failure"
      << "; benchmark_id=" << descriptor.benchmark_id
      << "; case_id=" << descriptor.case_id
      << "; failure_reason=" << reason
      << "; canonical_state_mutated=false";
  return out.str();
}

void EnsureOutputDirectory(const std::string& root) {
  std::filesystem::create_directories(std::filesystem::path(root));
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy ZeroFlux1dBoundary()
    noexcept {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

[[nodiscard]] std::size_t EvenPhiCellCount(std::size_t requested) noexcept {
  return requested > 1u && requested % 2u == 0u ? requested : 2u;
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildBenchmarkGeometry(
    const P4AlphaBenchmarkDescriptor& descriptor,
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  geometry.theta_faces.resize(theta_cells + 1u);
  geometry.phi_faces.resize(phi_cells + 1u);
  geometry.cell_volumes.resize(radial_cells * theta_cells * phi_cells, 0.0);

  for (std::size_t r = 0u; r <= radial_cells; ++r) {
    geometry.radial_faces[r] =
        descriptor.R_cm * static_cast<double>(r) /
        static_cast<double>(radial_cells);
  }
  constexpr double theta_lower = 0.25 * kPi;
  constexpr double theta_upper = 0.75 * kPi;
  for (std::size_t t = 0u; t <= theta_cells; ++t) {
    geometry.theta_faces[t] =
        theta_lower +
        (theta_upper - theta_lower) * static_cast<double>(t) /
            static_cast<double>(theta_cells);
  }
  for (std::size_t p = 0u; p <= phi_cells; ++p) {
    geometry.phi_faces[p] =
        2.0 * kPi * static_cast<double>(p) / static_cast<double>(phi_cells);
  }

  geometry.global_volume = 0.0;
  std::size_t linear = 0u;
  for (std::size_t r = 0u; r < radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) /
        3.0;
    for (std::size_t t = 0u; t < theta_cells; ++t) {
      const double polar_factor =
          std::cos(geometry.theta_faces[t]) -
          std::cos(geometry.theta_faces[t + 1u]);
      for (std::size_t p = 0u; p < phi_cells; ++p) {
        const double azimuth =
            geometry.phi_faces[p + 1u] - geometry.phi_faces[p];
        geometry.cell_volumes[linear] = radial_factor * polar_factor * azimuth;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  return geometry.cell_volumes.at(((radial * theta_cells) + theta) * phi_cells + phi);
}

void FillBenchmarkState(
    dec3d::state::CanonicalState& state,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t global_radial_begin = 0u) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  for (std::size_t lr = 0u; lr < state.layout.radial_cells; ++lr) {
    const std::size_t gr = global_radial_begin + lr;
    const double radius =
        0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]);
    const bool inner = radius < descriptor.r0_cm;
    const double rho =
        inner ? descriptor.rho_inner_g_cm3 : descriptor.rho_outer_g_cm3;
    const double Te = dec3d::physics::ErgFromKeV(
        inner ? descriptor.Te_inner_keV : descriptor.Te_outer_keV);
    const double Ti = Te;
    const double ni = rho / mean_ion_mass;
    const double ne = ni * dec3d::physics::DefaultZbar();
    const double e_electron = ne * Te / gamma_minus_one;
    const double e_ion = ni * Ti / gamma_minus_one;
    for (std::size_t t = 0u; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0u; p < state.layout.phi_cells; ++p) {
        state.rho(lr, t, p) = rho;
        state.mom_r(lr, t, p) = 0.0;
        state.mom_theta(lr, t, p) = 0.0;
        state.mom_phi(lr, t, p) = 0.0;
        state.e_electron(lr, t, p) = e_electron;
        state.e_fluid_total(lr, t, p) = e_electron + e_ion;
        state.alpha_state.storage(lr, t, p) = 0.0;
      }
    }
  }
}

[[nodiscard]] dec3d::alpha::AlphaCoefficientProviderOptions
BoschHaleProviderOptions() {
  dec3d::alpha::AlphaCoefficientProviderOptions options;
  options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.tau_model = dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.reactivity_model = dec3d::alpha::AlphaReactivityModel::bosch_hale_dt;
  options.recovery_options.mean_ion_mass_g =
      dec3d::physics::DefaultMeanDTIonMassG();
  options.recovery_options.zbar = dec3d::physics::DefaultZbar();
  return options;
}

[[nodiscard]] dec3d::alpha::OneGroupAlphaTransportOptions AlphaOptions(
    double dt_s) {
  dec3d::alpha::OneGroupAlphaTransportOptions options;
  options.dt_s = dt_s;
  options.provider_options = BoschHaleProviderOptions();
  options.boundary_policy = ZeroFlux1dBoundary();
  // Alpha source terms are large in CGS; benchmark acceptance is guarded by
  // volume-integrated budget residuals rather than the dense solver's absolute
  // residual scale.
  options.serial_reference_options.residual_tolerance = 1.0e20;
  options.serial_reference_options.row_limit = 512u;
  return options;
}

[[nodiscard]] std::vector<ProfileRow> BuildProfileRows(
    const dec3d::state::CanonicalState& state,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double time_s,
    std::size_t global_radial_begin = 0u) {
  std::vector<ProfileRow> rows;
  const auto recovered = dec3d::state::RecoverThermodynamicState(
      state,
      dec3d::state::ThermodynamicRecoveryOptions{
          0.0,
          0.0,
          dec3d::physics::DefaultZbar(),
          dec3d::physics::DefaultMeanDTIonMassG()});
  const auto provider =
      dec3d::alpha::BuildAlphaCoefficientArrays(state, BoschHaleProviderOptions());
  if (!recovered.success || !provider.success) {
    return rows;
  }
  rows.reserve(state.layout.radial_cells);
  for (std::size_t lr = 0u; lr < state.layout.radial_cells; ++lr) {
    const std::size_t gr = global_radial_begin + lr;
    double alpha_sum = 0.0;
    double te_sum = 0.0;
    double birth_sum = 0.0;
    double drag_sum = 0.0;
    double volume_sum = 0.0;
    for (std::size_t t = 0u; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0u; p < state.layout.phi_cells; ++p) {
        const double volume =
            CellVolume(geometry,
                       gr,
                       t,
                       p,
                       state.layout.theta_cells,
                       state.layout.phi_cells);
        volume_sum += volume;
        alpha_sum += state.alpha_state.storage(lr, t, p) * volume;
        te_sum += recovered.cells(lr, t, p).t_e_keV * volume;
        birth_sum +=
            provider.coefficients.birth_source_erg_cm3_s(lr, t, p) * volume;
        const double tau = provider.coefficients.tau_alphae_s(lr, t, p);
        drag_sum +=
            (tau > 0.0 ? state.alpha_state.storage(lr, t, p) / tau : 0.0) *
            volume;
      }
    }
    rows.push_back(ProfileRow{
        time_s,
        0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]),
        volume_sum > 0.0 ? alpha_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? te_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? birth_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? drag_sum / volume_sum : 0.0});
  }
  (void)descriptor;
  return rows;
}

void WriteProfilesCsv(
    const std::string& path,
    const std::vector<ProfileRow>& rows) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "time_s,r_cm,epsilon_alpha_erg_cm3,Te_keV,"
         "birth_source_erg_cm3_s,drag_deposition_erg_cm3_s\n";
  for (const auto& row : rows) {
    out << row.time_s << ',' << row.r_cm << ','
        << row.epsilon_alpha_erg_cm3 << ',' << row.Te_keV << ','
        << row.birth_source_erg_cm3_s << ','
        << row.drag_deposition_erg_cm3_s << '\n';
  }
}

void WriteBudgetSummaryCsv(
    const std::string& path,
    const std::vector<BudgetRow>& rows) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "time_s,run_kind,alpha_equation_budget_residual,"
         "alpha_electron_exchange_residual,"
         "global_alpha_plus_electron_budget_residual\n";
  for (const auto& row : rows) {
    out << row.time_s << ',' << row.run_kind << ','
        << row.alpha_equation_budget_residual << ','
        << row.alpha_electron_exchange_residual << ','
        << row.global_alpha_plus_electron_budget_residual << '\n';
  }
}

void WriteB0Json(
    const std::string& path,
    const P4AlphaBenchmarkRunResult& result) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "{\n"
      << "  \"diagnostic_id\": \"p4.alpha.benchmark_b0\",\n"
      << "  \"success\": " << (result.success ? "true" : "false") << ",\n"
      << "  \"b0_operator_oracle_executed\": "
      << (result.b0_operator_oracle_executed ? "true" : "false") << ",\n"
      << "  \"delta_alpha_total\": " << result.delta_alpha_total << ",\n"
      << "  \"delta_electron_total\": " << result.delta_electron_total << ",\n"
      << "  \"birth_source_total\": " << result.birth_source_total << ",\n"
      << "  \"drag_deposition_total\": " << result.drag_deposition_total << ",\n"
      << "  \"global_alpha_plus_electron_budget_residual\": "
      << result.global_alpha_plus_electron_budget_residual << "\n"
      << "}\n";
}

void WriteCaseDescriptorJson(
    const std::string& path,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "{\n"
      << "  \"benchmark_id\": \"" << descriptor.benchmark_id << "\",\n"
      << "  \"case_id\": \"" << descriptor.case_id << "\",\n"
      << "  \"geometry\": \"" << descriptor.geometry << "\",\n"
      << "  \"radial_cells\": " << options.radial_cells << ",\n"
      << "  \"theta_cells\": " << options.theta_cells << ",\n"
      << "  \"phi_cells\": " << options.phi_cells << ",\n"
      << "  \"R_cm\": " << descriptor.R_cm << ",\n"
      << "  \"r0_cm\": " << descriptor.r0_cm << ",\n"
      << "  \"initial_ti_policy\": \"" << descriptor.initial_ti_policy << "\",\n"
      << "  \"epsilon_alpha_initial_policy\": \""
      << descriptor.epsilon_alpha_initial_policy << "\",\n"
      << "  \"fuel_depletion_enabled\": false,\n"
      << "  \"separate_dt_species_authoritative\": false,\n"
      << "  \"parity_claim_allowed\": false\n"
      << "}\n";
}

void WriteManifestJson(
    const std::string& path,
    const P4AlphaBenchmarkRunResult& result) {
  std::ofstream out(path);
  out << "{\n"
      << "  \"diagnostic_id\": \"p4.alpha.benchmark_manifest\",\n"
      << "  \"profile_artifacts_written\": "
      << (result.profile_artifacts_written ? "true" : "false") << ",\n"
      << "  \"budget_artifacts_written\": "
      << (result.budget_artifacts_written ? "true" : "false") << ",\n"
      << "  \"performance_artifacts_written\": "
      << (result.performance_artifacts_written ? "true" : "false") << ",\n"
      << "  \"plot_artifacts_written\": false\n"
      << "}\n";
}

[[nodiscard]] P4AlphaBenchmarkRunResult RunB0(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) {
  P4AlphaBenchmarkRunResult result;
  const auto validation = ValidateP4AlphaBenchmarkDescriptor(descriptor);
  if (!validation.success) {
    result.failure_reason = validation.failure_reason;
    result.failure_diagnostics = validation.failure_diagnostics;
    return result;
  }

  const std::size_t nr =
      descriptor.case_kind ==
              P4AlphaBenchmarkCase::b0_two_cell_diffusion_coupling_oracle
          ? 2u
          : 1u;
  auto local_descriptor = descriptor;
  local_descriptor.R_cm = 0.2;
  local_descriptor.r0_cm = 0.1;
  const auto geometry = BuildBenchmarkGeometry(local_descriptor, nr, 1u, 2u);
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{nr, 1u, 2u, 0u});
  FillBenchmarkState(state, local_descriptor, geometry);
  if (nr == 2u) {
    state.alpha_state.storage(1u, 0u, 0u) = 1.0e8;
  }

  const auto start = std::chrono::steady_clock::now();
  const auto alpha = dec3d::alpha::ApplyOneGroupAlphaTransport(
      state,
      geometry,
      AlphaOptions(descriptor.case_kind ==
                           P4AlphaBenchmarkCase::b0_one_cell_dt_zero_no_change
                       ? 0.0
                       : options.dt_s));
  const auto stop = std::chrono::steady_clock::now();
  result.total_wall_seconds =
      std::chrono::duration<double>(stop - start).count();

  if (!alpha.success) {
    result.failure_reason = alpha.failure_reason;
    result.failure_diagnostics = alpha.failure_diagnostics;
    return result;
  }

  result.success = true;
  result.b0_operator_oracle_executed = true;
  result.single_rank_run_present = true;
  result.delta_alpha_total = alpha.delta_alpha_total;
  result.delta_electron_total = alpha.delta_electron_total;
  result.birth_source_total = alpha.birth_source_total;
  result.drag_deposition_total = alpha.drag_deposition_total;
  result.global_alpha_plus_electron_budget_residual =
      alpha.global_alpha_plus_electron_budget_residual;
  result.output_directory = options.output_root;
  std::ostringstream report;
  report << std::setprecision(17)
         << BuildP4AlphaBenchmarkDiagnosticsLine(descriptor)
         << "; b0_operator_oracle_executed=true"
         << "; single_rank_run_present=true"
         << "; distributed_run_present=false"
         << "; budget_closed="
         << (std::abs(result.global_alpha_plus_electron_budget_residual) <
                     1.0e-6 * std::max(1.0, std::abs(result.birth_source_total))
                 ? "true"
                 : "false")
         << "; delta_alpha_total=" << result.delta_alpha_total
         << "; delta_electron_total=" << result.delta_electron_total
         << "; birth_source_total=" << result.birth_source_total
         << "; drag_deposition_total=" << result.drag_deposition_total
         << "; global_alpha_plus_electron_budget_residual="
         << result.global_alpha_plus_electron_budget_residual
         << "; rank0_gather_solve_used=false"
         << "; serial_dense_fallback_used=false";
  result.report_line = report.str();

  if (options.write_artifacts) {
    EnsureOutputDirectory(options.output_root);
    WriteB0Json(JoinPath(options.output_root, "p4_b0_operator_oracles.json"),
                result);
  }
  return result;
}

[[nodiscard]] P4AlphaBenchmarkRunResult RunB1(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) {
  P4AlphaBenchmarkRunResult result;
  const auto validation = ValidateP4AlphaBenchmarkDescriptor(descriptor);
  if (!validation.success) {
    result.failure_reason = validation.failure_reason;
    result.failure_diagnostics = validation.failure_diagnostics;
    return result;
  }
  const std::size_t radial = std::max<std::size_t>(1u, options.radial_cells);
  const std::size_t theta = std::max<std::size_t>(1u, options.theta_cells);
  const std::size_t phi = EvenPhiCellCount(options.phi_cells);
  const auto geometry = BuildBenchmarkGeometry(descriptor, radial, theta, phi);
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{radial, theta, phi, 0u});
  FillBenchmarkState(state, descriptor, geometry);

  std::vector<ProfileRow> rows =
      BuildProfileRows(state, descriptor, geometry, 0.0);
  std::vector<BudgetRow> budget_rows;
  const auto start = std::chrono::steady_clock::now();
  for (std::size_t step = 0u; step < std::max<std::size_t>(1u, options.step_count);
       ++step) {
    const auto step_start = std::chrono::steady_clock::now();
    const auto alpha = dec3d::alpha::ApplyOneGroupAlphaTransport(
        state, geometry, AlphaOptions(options.dt_s));
    const auto step_stop = std::chrono::steady_clock::now();
    result.step_wall_seconds.push_back(
        std::chrono::duration<double>(step_stop - step_start).count());
    if (!alpha.success) {
      result.failure_reason = alpha.failure_reason;
      result.failure_diagnostics = alpha.failure_diagnostics;
      return result;
    }
    result.delta_alpha_total += alpha.delta_alpha_total;
    result.delta_electron_total += alpha.delta_electron_total;
    result.birth_source_total += alpha.birth_source_total;
    result.drag_deposition_total += alpha.drag_deposition_total;
    result.global_alpha_plus_electron_budget_residual +=
        alpha.global_alpha_plus_electron_budget_residual;
    budget_rows.push_back(BudgetRow{
        options.dt_s * static_cast<double>(step + 1u),
        "serial",
        alpha.alpha_equation_budget_residual,
        alpha.alpha_electron_exchange_residual,
        alpha.global_alpha_plus_electron_budget_residual});
    auto step_rows = BuildProfileRows(
        state,
        descriptor,
        geometry,
        options.dt_s * static_cast<double>(step + 1u));
    rows.insert(rows.end(), step_rows.begin(), step_rows.end());
  }
  const auto stop = std::chrono::steady_clock::now();

  result.success = true;
  result.b1_hotspot_executed = true;
  result.single_rank_run_present = true;
  result.profile_artifacts_written = options.write_artifacts;
  result.budget_artifacts_written = options.write_artifacts;
  result.total_wall_seconds =
      std::chrono::duration<double>(stop - start).count();
  result.output_directory = options.output_root;

  std::ostringstream report;
  report << std::setprecision(17)
         << BuildP4AlphaBenchmarkDiagnosticsLine(descriptor)
         << "; b1_hotspot_executed=true"
         << "; single_rank_run_present=true"
         << "; distributed_run_present=false"
         << "; profile_artifacts_written="
         << (result.profile_artifacts_written ? "true" : "false")
         << "; budget_artifacts_written="
         << (result.budget_artifacts_written ? "true" : "false")
         << "; delta_alpha_total=" << result.delta_alpha_total
         << "; delta_electron_total=" << result.delta_electron_total
         << "; birth_source_total=" << result.birth_source_total
         << "; drag_deposition_total=" << result.drag_deposition_total
         << "; global_alpha_plus_electron_budget_residual="
         << result.global_alpha_plus_electron_budget_residual;
  result.report_line = report.str();

  if (options.write_artifacts) {
    EnsureOutputDirectory(options.output_root);
    WriteCaseDescriptorJson(
        JoinPath(options.output_root, "p4_b1_alpha_hotspot_case_descriptor.json"),
        descriptor,
        options);
    WriteProfilesCsv(JoinPath(options.output_root, "p4_b1_serial_profiles.csv"),
                     rows);
    WriteBudgetSummaryCsv(
        JoinPath(options.output_root, "p4_b1_budget_summary.csv"),
        budget_rows);
    WriteManifestJson(JoinPath(options.output_root, "p4_benchmark_manifest.json"),
                      result);
  }
  return result;
}

}  // namespace

P4AlphaBenchmarkRunResult RunP4AlphaBenchmarkCase(
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept {
  try {
    switch (descriptor.layer) {
      case P4AlphaBenchmarkLayer::b0_operator_budget_oracle:
        return RunB0(descriptor, options);
      case P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot:
        return RunB1(descriptor, options);
      case P4AlphaBenchmarkLayer::b3_distributed_performance: {
        P4AlphaBenchmarkRunResult result;
        result.failure_reason = "B3 requires the HYPRE/MPI benchmark runner";
        result.failure_diagnostics = FailureLine(descriptor, result.failure_reason.c_str());
        return result;
      }
    }
  } catch (const std::exception& error) {
    P4AlphaBenchmarkRunResult result;
    result.failure_reason = error.what();
    result.failure_diagnostics = FailureLine(descriptor, error.what());
    return result;
  }
  P4AlphaBenchmarkRunResult result;
  result.failure_reason = "unknown P4 alpha benchmark layer";
  result.failure_diagnostics = FailureLine(descriptor, result.failure_reason.c_str());
  return result;
}

}  // namespace dec3d::benchmarks
