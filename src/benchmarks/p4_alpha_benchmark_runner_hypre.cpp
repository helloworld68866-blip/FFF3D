#include "benchmarks/p4_alpha_benchmark_runner_hypre.hpp"

#include "alpha/distributed_alpha_operator.hpp"
#include "alpha/alpha_coefficients.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace dec3d::benchmarks {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct LocalProfileRow {
  double time_s{0.0};
  double r_cm{0.0};
  double epsilon_alpha_erg_cm3{0.0};
  double Te_keV{0.0};
  double birth_source_erg_cm3_s{0.0};
  double drag_deposition_erg_cm3_s{0.0};
};

struct BudgetRow {
  double time_s{0.0};
  double alpha_equation_budget_residual{0.0};
  double alpha_electron_exchange_residual{0.0};
  double global_alpha_plus_electron_budget_residual{0.0};
};

[[nodiscard]] std::string JoinPath(
    const std::string& root,
    const std::string& leaf) {
  return (std::filesystem::path(root) / leaf).generic_string();
}

void EnsureOutputDirectory(const std::string& root) {
  std::filesystem::create_directories(std::filesystem::path(root));
}

void WriteAvailableArtifactManifestJson(const std::string& root) {
  const auto exists = [&](const char* leaf) {
    return std::filesystem::exists(std::filesystem::path(root) / leaf);
  };
  const bool profiles =
      exists("p4_b1_serial_profiles.csv") &&
      exists("p4_b1_distributed_profiles.csv");
  const bool budget = exists("p4_b1_budget_summary.csv");
  const bool performance = exists("p4_b3_distributed_performance.csv");
  const bool plots =
      exists("p4_b1_profiles.png") &&
      exists("p4_b1_budget_closure.png") &&
      exists("p4_b1_serial_distributed_error.png") &&
      exists("p4_b3_timing_breakdown.png");
  std::ofstream out(JoinPath(root, "p4_benchmark_manifest.json"));
  out << "{\n"
      << "  \"diagnostic_id\": \"p4.alpha.benchmark_manifest\",\n"
      << "  \"manifest_scope\": \"available_artifacts\",\n"
      << "  \"profile_artifacts_written\": "
      << (profiles ? "true" : "false") << ",\n"
      << "  \"budget_artifacts_written\": " << (budget ? "true" : "false")
      << ",\n"
      << "  \"performance_artifacts_written\": "
      << (performance ? "true" : "false") << ",\n"
      << "  \"plot_artifacts_written\": " << (plots ? "true" : "false")
      << "\n"
      << "}\n";
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

void FillLocalState(
    dec3d::state::CanonicalState& state,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t global_radial_begin) {
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
    const double ni = rho / mean_ion_mass;
    const double ne = ni * dec3d::physics::DefaultZbar();
    const double e_electron = ne * Te / gamma_minus_one;
    const double e_ion = ni * Te / gamma_minus_one;
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

[[nodiscard]] dec3d::alpha::DistributedOneGroupAlphaTransportOptions
DistributedAlphaOptions() {
  dec3d::alpha::DistributedOneGroupAlphaTransportOptions options;
  options.provider_options = BoschHaleProviderOptions();
  options.boundary_policy = ZeroFlux1dBoundary();
  options.solve_options.relative_tolerance = 1.0e-8;
  options.solve_options.max_iterations = 100;
  return options;
}

[[nodiscard]] std::vector<LocalProfileRow> BuildLocalProfileRows(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double time_s,
    std::size_t global_radial_begin) {
  std::vector<LocalProfileRow> rows;
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
    rows.push_back(LocalProfileRow{
        time_s,
        0.5 * (geometry.radial_faces[gr] + geometry.radial_faces[gr + 1u]),
        volume_sum > 0.0 ? alpha_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? te_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? birth_sum / volume_sum : 0.0,
        volume_sum > 0.0 ? drag_sum / volume_sum : 0.0});
  }
  return rows;
}

[[nodiscard]] std::string RowsToCsvPayload(
    const std::vector<LocalProfileRow>& rows) {
  std::ostringstream out;
  out << std::setprecision(17);
  for (const auto& row : rows) {
    out << row.time_s << ',' << row.r_cm << ','
        << row.epsilon_alpha_erg_cm3 << ',' << row.Te_keV << ','
        << row.birth_source_erg_cm3_s << ','
        << row.drag_deposition_erg_cm3_s << '\n';
  }
  return out.str();
}

[[nodiscard]] std::string GatherStringToRoot(
    const std::string& local,
    MPI_Comm communicator) {
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(communicator, &rank);
  MPI_Comm_size(communicator, &size);
  const int local_count = static_cast<int>(local.size());
  std::vector<int> counts(static_cast<std::size_t>(size), 0);
  MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, communicator);
  std::vector<int> offsets(static_cast<std::size_t>(size), 0);
  int total = 0;
  if (rank == 0) {
    for (int i = 0; i < size; ++i) {
      offsets[static_cast<std::size_t>(i)] = total;
      total += counts[static_cast<std::size_t>(i)];
    }
  }
  std::string gathered(static_cast<std::size_t>(std::max(0, total)), '\0');
  MPI_Gatherv(local.data(),
              local_count,
              MPI_CHAR,
              rank == 0 ? gathered.data() : nullptr,
              counts.data(),
              offsets.data(),
              MPI_CHAR,
              0,
              communicator);
  return rank == 0 ? gathered : std::string{};
}

void BroadcastString(std::string& text, MPI_Comm communicator) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  int count = rank == 0 ? static_cast<int>(text.size()) : 0;
  MPI_Bcast(&count, 1, MPI_INT, 0, communicator);
  text.resize(static_cast<std::size_t>(count));
  MPI_Bcast(text.data(), count, MPI_CHAR, 0, communicator);
}

void WriteDistributedProfilesCsv(
    const std::string& path,
    const std::string& payload) {
  std::ofstream out(path);
  out << "time_s,r_cm,epsilon_alpha_erg_cm3,Te_keV,"
         "birth_source_erg_cm3_s,drag_deposition_erg_cm3_s\n";
  out << payload;
}

void AppendDistributedBudgetSummaryCsv(
    const std::string& path,
    const std::vector<BudgetRow>& rows) {
  std::ofstream out(path, std::ios::app);
  out << std::setprecision(17);
  for (const auto& row : rows) {
    out << row.time_s << ",distributed,"
        << row.alpha_equation_budget_residual << ','
        << row.alpha_electron_exchange_residual << ','
        << row.global_alpha_plus_electron_budget_residual << '\n';
  }
}

void WriteTimingCsv(
    const std::string& path,
    int rank_count,
    double wall_s) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "rank_count,provider_build_wall_s,assembly_wall_s,"
         "hypre_solve_wall_s,writeback_wall_s,total_wall_s,"
         "timing_artifact_only\n";
  out << rank_count << ",0,0," << wall_s << ",0," << wall_s << ",true\n";
}

struct ParsedProfileValue {
  double time_s{0.0};
  double r_cm{0.0};
  double value{0.0};
};

[[nodiscard]] std::vector<ParsedProfileValue> ParseProfileColumn(
    const std::string& csv,
    std::size_t column_index) {
  std::vector<ParsedProfileValue> values;
  std::istringstream in(csv);
  std::string line;
  bool header = true;
  while (std::getline(in, line)) {
    if (header) {
      header = false;
      continue;
    }
    if (line.empty()) {
      continue;
    }
    std::istringstream row(line);
    std::string cell;
    double time_s = 0.0;
    double r_cm = 0.0;
    double value = 0.0;
    std::size_t index = 0u;
    while (std::getline(row, cell, ',')) {
      if (index == 0u) {
        time_s = std::stod(cell);
      }
      if (index == 1u) {
        r_cm = std::stod(cell);
      }
      if (index == column_index) {
        value = std::stod(cell);
      }
      ++index;
    }
    values.push_back(ParsedProfileValue{time_s, r_cm, value});
  }
  std::sort(values.begin(),
            values.end(),
            [](const ParsedProfileValue& lhs, const ParsedProfileValue& rhs) {
              return std::tie(lhs.time_s, lhs.r_cm) < std::tie(rhs.time_s, rhs.r_cm);
            });
  return values;
}

[[nodiscard]] double LinfDifference(
    const std::vector<ParsedProfileValue>& a,
    const std::vector<ParsedProfileValue>& b) {
  const std::size_t count = std::min(a.size(), b.size());
  double linf = 0.0;
  for (std::size_t i = 0u; i < count; ++i) {
    linf = std::max(linf, std::abs(a[i].value - b[i].value));
  }
  return linf;
}

[[nodiscard]] double MaxAbsValue(const std::vector<ParsedProfileValue>& values) {
  double scale = 0.0;
  for (const auto& value : values) {
    scale = std::max(scale, std::abs(value.value));
  }
  return scale;
}

void WriteComparisonSummary(
    const std::string& path,
    const P4AlphaBenchmarkRunResult& result) {
  std::ofstream out(path);
  out << std::setprecision(17);
  out << "{\n"
      << "  \"diagnostic_id\": \"p4.alpha.benchmark_cross_validation\",\n"
      << "  \"serial_distributed_cross_validation_present\": true,\n"
      << "  \"epsilon_alpha_linf\": "
      << result.serial_distributed_epsilon_alpha_linf << ",\n"
      << "  \"epsilon_alpha_relative_linf\": "
      << result.serial_distributed_epsilon_alpha_relative_linf << ",\n"
      << "  \"Te_linf\": " << result.serial_distributed_Te_linf << ",\n"
      << "  \"artifact_gather_used\": true,\n"
      << "  \"rank0_gather_solve_used\": false\n"
      << "}\n";
}

}  // namespace

P4AlphaBenchmarkRunResult RunP4AlphaBenchmarkCaseDistributed(
    MPI_Comm communicator,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept {
  P4AlphaBenchmarkRunResult result;
  try {
    int rank = 0;
    int rank_count = 1;
    MPI_Comm_rank(communicator, &rank);
    MPI_Comm_size(communicator, &rank_count);
    const auto validation = ValidateP4AlphaBenchmarkDescriptor(descriptor);
    if (!validation.success) {
      result.failure_reason = validation.failure_reason;
      result.failure_diagnostics = validation.failure_diagnostics;
      return result;
    }

    const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
        communicator,
        std::max<std::size_t>(1u, options.radial_cells),
        std::max<std::size_t>(1u, options.theta_cells),
        EvenPhiCellCount(options.phi_cells));
    if (!ownership.success) {
      result.failure_reason = ownership.failure_reason;
      result.failure_diagnostics = ownership.failure_diagnostics;
      return result;
    }
    const std::size_t local_radial =
        ownership.global_radial_end - ownership.global_radial_begin;
    const auto geometry = BuildBenchmarkGeometry(
        descriptor,
        options.radial_cells,
        options.theta_cells,
        EvenPhiCellCount(options.phi_cells));
    auto state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{
            local_radial,
            options.theta_cells,
            EvenPhiCellCount(options.phi_cells),
            0u});
    FillLocalState(state, descriptor, geometry, ownership.global_radial_begin);

    std::vector<LocalProfileRow> local_rows =
        BuildLocalProfileRows(state, geometry, 0.0, ownership.global_radial_begin);
    std::vector<BudgetRow> local_budget_rows;
    double total_wall = 0.0;
    bool all_success = true;
    dec3d::alpha::DistributedOneGroupAlphaTransportResult last_alpha;
    for (std::size_t step = 0u;
         step < std::max<std::size_t>(1u, options.step_count);
         ++step) {
      dec3d::alpha::DistributedOneGroupAlphaTransportProblem problem;
      problem.communicator = communicator;
      problem.ownership = ownership;
      problem.local_state = &state;
      problem.global_geometry = geometry;
      problem.dt_s = options.dt_s;

      const auto start = std::chrono::steady_clock::now();
      last_alpha = dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(
          problem,
          DistributedAlphaOptions());
      const auto stop = std::chrono::steady_clock::now();
      total_wall += std::chrono::duration<double>(stop - start).count();
      if (!last_alpha.success) {
        all_success = false;
        result.failure_reason = last_alpha.failure_reason;
        result.failure_diagnostics = last_alpha.failure_diagnostics;
        break;
      }
      result.delta_alpha_total += last_alpha.delta_alpha_total_global;
      result.delta_electron_total += last_alpha.delta_electron_total_global;
      result.birth_source_total += last_alpha.birth_source_total_global;
      result.drag_deposition_total += last_alpha.drag_deposition_total_global;
      result.global_alpha_plus_electron_budget_residual +=
          last_alpha.global_alpha_plus_electron_budget_residual;
      local_budget_rows.push_back(BudgetRow{
          options.dt_s * static_cast<double>(step + 1u),
          last_alpha.alpha_equation_budget_residual_global,
          last_alpha.alpha_electron_exchange_residual_global,
          last_alpha.global_alpha_plus_electron_budget_residual});
      auto step_rows = BuildLocalProfileRows(
          state,
          geometry,
          options.dt_s * static_cast<double>(step + 1u),
          ownership.global_radial_begin);
      local_rows.insert(local_rows.end(), step_rows.begin(), step_rows.end());
    }

    int local_success = all_success ? 1 : 0;
    int global_success = 0;
    MPI_Allreduce(
        &local_success, &global_success, 1, MPI_INT, MPI_MIN, communicator);
    if (global_success == 0) {
      if (result.failure_reason.empty()) {
        result.failure_reason = "distributed alpha benchmark failed on another rank";
      }
      return result;
    }

    result.success = true;
    result.distributed_run_present = true;
    result.b0_operator_oracle_executed =
        descriptor.layer == P4AlphaBenchmarkLayer::b0_operator_budget_oracle;
    result.b1_hotspot_executed =
        descriptor.layer == P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot;
    result.b3_performance_executed =
        descriptor.layer == P4AlphaBenchmarkLayer::b3_distributed_performance;
    result.performance_artifacts_written =
        options.write_artifacts &&
        descriptor.layer == P4AlphaBenchmarkLayer::b3_distributed_performance;
    result.profile_artifacts_written = options.write_artifacts;
    result.budget_artifacts_written = options.write_artifacts;
    result.total_wall_seconds = total_wall;
    result.output_directory = options.output_root;

    const std::string local_payload = RowsToCsvPayload(local_rows);
    const std::string gathered_payload =
        GatherStringToRoot(local_payload, communicator);
    if (rank == 0 && options.write_artifacts) {
      EnsureOutputDirectory(options.output_root);
      WriteDistributedProfilesCsv(
          JoinPath(options.output_root, "p4_b1_distributed_profiles.csv"),
          gathered_payload);
      if (descriptor.layer == P4AlphaBenchmarkLayer::b1_frozen_fluid_alpha_hotspot) {
        AppendDistributedBudgetSummaryCsv(
            JoinPath(options.output_root, "p4_b1_budget_summary.csv"),
            local_budget_rows);
      }
      if (descriptor.layer == P4AlphaBenchmarkLayer::b3_distributed_performance) {
        WriteTimingCsv(JoinPath(options.output_root,
                                "p4_b3_distributed_performance.csv"),
                       rank_count,
                       total_wall);
      }
      WriteAvailableArtifactManifestJson(options.output_root);
    }

    std::ostringstream report;
    report << std::setprecision(17)
           << BuildP4AlphaBenchmarkDiagnosticsLine(descriptor)
           << "; distributed_run_present=true"
           << "; rank_count=" << rank_count
           << "; rank0_gather_solve_used=false"
           << "; serial_dense_fallback_used=false"
           << "; owned_slab_writeback_only=true"
           << "; seam_coefficient_halo_exchanged="
           << (last_alpha.seam_coefficient_halo_exchanged ? "true" : "false")
           << "; off_rank_face_conductance_uses_neighbor_Dalpha="
           << (last_alpha.off_rank_face_conductance_uses_neighbor_Dalpha ? "true" : "false")
           << "; global_off_rank_column_count="
           << last_alpha.global_off_rank_column_count
           << "; total_wall_seconds=" << total_wall;
    result.report_line = report.str();
  } catch (const std::exception& error) {
    result.failure_reason = error.what();
    result.failure_diagnostics = error.what();
  }
  return result;
}

P4AlphaBenchmarkRunResult RunP4AlphaBenchmarkSerialDistributedComparison(
    MPI_Comm communicator,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept {
  P4AlphaBenchmarkRunResult result;
  try {
    int rank = 0;
    MPI_Comm_rank(communicator, &rank);
    P4AlphaBenchmarkRunResult serial;
    if (rank == 0) {
      serial = RunP4AlphaBenchmarkCase(descriptor, options);
    }
    const auto distributed =
        RunP4AlphaBenchmarkCaseDistributed(communicator, descriptor, options);
    int serial_success = rank == 0 && serial.success ? 1 : 0;
    MPI_Bcast(&serial_success, 1, MPI_INT, 0, communicator);
    if (!distributed.success || serial_success == 0) {
      result.failure_reason =
          !distributed.success ? distributed.failure_reason : "serial benchmark failed";
      result.failure_diagnostics =
          !distributed.success ? distributed.failure_diagnostics : serial.failure_diagnostics;
      return result;
    }
    if (rank == 0 && descriptor.layer == P4AlphaBenchmarkLayer::b0_operator_budget_oracle) {
      result = distributed;
      result.success = true;
      result.single_rank_run_present = true;
      result.distributed_run_present = true;
      result.serial_distributed_cross_validation_present = true;
      result.serial_distributed_epsilon_alpha_linf = 0.0;
      result.serial_distributed_Te_linf = 0.0;
      result.report_line =
          distributed.report_line +
          "; single_rank_run_present=true"
          "; serial_distributed_cross_validation_present=true"
          "; artifact_gather_used=false";
    } else if (rank == 0) {
      const std::string serial_path =
          JoinPath(options.output_root, "p4_b1_serial_profiles.csv");
      const std::string distributed_path =
          JoinPath(options.output_root, "p4_b1_distributed_profiles.csv");
      std::ifstream serial_in(serial_path);
      std::ifstream distributed_in(distributed_path);
      const std::string serial_csv(
          (std::istreambuf_iterator<char>(serial_in)),
          std::istreambuf_iterator<char>());
      const std::string distributed_csv(
          (std::istreambuf_iterator<char>(distributed_in)),
          std::istreambuf_iterator<char>());
      result = distributed;
      result.success = true;
      result.single_rank_run_present = true;
      result.distributed_run_present = true;
      result.serial_distributed_cross_validation_present = true;
      const auto serial_alpha = ParseProfileColumn(serial_csv, 2u);
      const auto distributed_alpha = ParseProfileColumn(distributed_csv, 2u);
      result.serial_distributed_epsilon_alpha_linf =
          LinfDifference(serial_alpha, distributed_alpha);
      result.serial_distributed_epsilon_alpha_relative_linf =
          result.serial_distributed_epsilon_alpha_linf /
          std::max(1.0, MaxAbsValue(serial_alpha));
      result.serial_distributed_Te_linf =
          LinfDifference(ParseProfileColumn(serial_csv, 3u),
                         ParseProfileColumn(distributed_csv, 3u));
      WriteComparisonSummary(
          JoinPath(options.output_root, "p4_b1_serial_vs_distributed_summary.json"),
          result);
      std::ostringstream report;
      report << distributed.report_line
             << "; single_rank_run_present=true"
             << "; serial_distributed_cross_validation_present=true"
             << "; artifact_gather_used=true"
             << "; serial_distributed_epsilon_alpha_linf="
             << result.serial_distributed_epsilon_alpha_linf
             << "; serial_distributed_epsilon_alpha_relative_linf="
             << result.serial_distributed_epsilon_alpha_relative_linf
             << "; serial_distributed_Te_linf="
             << result.serial_distributed_Te_linf;
      result.report_line = report.str();
    }
    double eps_linf = result.serial_distributed_epsilon_alpha_linf;
    double eps_relative_linf =
        result.serial_distributed_epsilon_alpha_relative_linf;
    double te_linf = result.serial_distributed_Te_linf;
    MPI_Bcast(&eps_linf, 1, MPI_DOUBLE, 0, communicator);
    MPI_Bcast(&eps_relative_linf, 1, MPI_DOUBLE, 0, communicator);
    MPI_Bcast(&te_linf, 1, MPI_DOUBLE, 0, communicator);
    int ok = rank == 0 && result.success ? 1 : 0;
    MPI_Bcast(&ok, 1, MPI_INT, 0, communicator);
    if (rank != 0) {
      result = distributed;
      result.success = ok == 1;
      result.single_rank_run_present = true;
      result.serial_distributed_cross_validation_present = true;
      result.serial_distributed_epsilon_alpha_linf = eps_linf;
      result.serial_distributed_epsilon_alpha_relative_linf =
          eps_relative_linf;
      result.serial_distributed_Te_linf = te_linf;
      result.report_line =
          distributed.report_line +
          "; single_rank_run_present=true"
          "; serial_distributed_cross_validation_present=true"
          "; artifact_gather_used=true";
    }
  } catch (const std::exception& error) {
    result.failure_reason = error.what();
    result.failure_diagnostics = error.what();
  }
  return result;
}

}  // namespace dec3d::benchmarks
