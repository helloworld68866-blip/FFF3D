#include "alpha/distributed_alpha_operator.hpp"

#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::alpha {
namespace {

constexpr double kBudgetTolerance = 1.0e-6;

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] int Rank(MPI_Comm communicator) noexcept {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  return rank;
}

[[nodiscard]] int AllreduceMinInt(MPI_Comm communicator, int value) noexcept {
  int global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global;
}

[[nodiscard]] double AllreduceSumDouble(MPI_Comm communicator, double value) noexcept {
  double global = 0.0;
  MPI_Allreduce(&value, &global, 1, MPI_DOUBLE, MPI_SUM, communicator);
  return global;
}

[[nodiscard]] double AllreduceMinDouble(MPI_Comm communicator, double value) noexcept {
  double global = 0.0;
  MPI_Allreduce(&value, &global, 1, MPI_DOUBLE, MPI_MIN, communicator);
  return global;
}

[[nodiscard]] double AllreduceMaxDouble(MPI_Comm communicator, double value) noexcept {
  double global = 0.0;
  MPI_Allreduce(&value, &global, 1, MPI_DOUBLE, MPI_MAX, communicator);
  return global;
}

[[nodiscard]] std::size_t LocalCellCount(const dec3d::state::CanonicalState& state) noexcept {
  return state.rho.extent_r() * state.rho.extent_theta() * state.rho.extent_phi();
}

[[nodiscard]] double BudgetScale(double value) noexcept {
  return std::max(1.0, std::abs(value));
}

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] std::string UpdatedFieldsLabel(std::uint32_t mask) {
  if (mask == kDistributedAlphaWritesNoFields) {
    return "none";
  }
  std::string label;
  if ((mask & kDistributedAlphaWritesAlphaState) != 0u) {
    label += "alpha_state";
  }
  if ((mask & kDistributedAlphaWritesElectronEnergy) != 0u) {
    if (!label.empty()) {
      label += ",";
    }
    label += "e_electron";
  }
  if ((mask & kDistributedAlphaWritesFluidTotalEnergy) != 0u) {
    if (!label.empty()) {
      label += ",";
    }
    label += "e_fluid_total";
  }
  return label;
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask DistributedAlphaWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

[[nodiscard]] dec3d::core::Array3D<double> FilledLike(
    const dec3d::state::CanonicalState& state,
    double value) {
  return dec3d::core::Array3D<double>(
      state.rho.extent_r(),
      state.rho.extent_theta(),
      state.rho.extent_phi(),
      value);
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    std::size_t local_r,
    std::size_t theta,
    std::size_t phi) noexcept {
  const std::size_t global_r = ownership.global_radial_begin + local_r;
  const std::size_t index =
      ((global_r * ownership.global_theta_cells) + theta) * ownership.global_phi_cells + phi;
  return geometry.cell_volumes[index];
}

[[nodiscard]] dec3d::transport::DistributedGenericDiffusionProblem BuildAlphaDiffusionProblem(
    const DistributedOneGroupAlphaTransportProblem& problem,
    const AlphaCoefficientProviderResult& provider,
    const DistributedOneGroupAlphaTransportOptions& options) {
  dec3d::transport::DistributedGenericDiffusionProblem diffusion;
  diffusion.ownership = problem.ownership;
  diffusion.global_geometry = problem.global_geometry;
  diffusion.dt_s = problem.dt_s;
  diffusion.local_coefficient_A = FilledLike(*problem.local_state, 1.0);
  diffusion.local_coefficient_D = provider.coefficients.D_alpha_cm2_s;
  diffusion.local_coefficient_C = FilledLike(*problem.local_state, 0.0);
  diffusion.local_coefficient_B = provider.coefficients.birth_source_erg_cm3_s;
  diffusion.local_scalar_old = problem.local_state->alpha_state.storage;
  diffusion.boundary_policy = options.boundary_policy;
  for (std::size_t r = 0; r < diffusion.local_coefficient_C.extent_r(); ++r) {
    for (std::size_t t = 0; t < diffusion.local_coefficient_C.extent_theta(); ++t) {
      for (std::size_t p = 0; p < diffusion.local_coefficient_C.extent_phi(); ++p) {
        diffusion.local_coefficient_C(r, t, p) =
            -1.0 / provider.coefficients.tau_alphae_s(r, t, p);
      }
    }
  }
  return diffusion;
}

void BuildReport(
    DistributedOneGroupAlphaTransportResult& result,
    const char* diagnostic_id) {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=" << diagnostic_id
         << "; phase_id=P4"
         << "; stage_id=A"
         << "; alpha_transport_model=atzeni_one_group"
         << "; distributed_alpha_solve=" << BoolToken(result.distributed_alpha_solve_executed)
         << "; backend_requested=" << result.backend_requested
         << "; backend_executed=" << result.backend_executed
         << "; coefficient_time_level=old_time_lagged"
         << "; coefficient_provider_wall_s=" << result.coefficient_provider_wall_s
         << "; assembly_wall_s=" << result.assembly_wall_s
         << "; hypre_setup_wall_s=" << result.hypre_setup_wall_s
         << "; hypre_solve_wall_s=" << result.hypre_solve_wall_s
         << "; writeback_wall_s=" << result.writeback_wall_s
         << "; solver_iterations=" << result.solver_iterations
         << "; provider_report_present=" << BoolToken(!result.provider_report.empty())
         << "; assembly_report_present=" << BoolToken(!result.assembly_report.empty())
         << "; hypre_solve_report_present=" << BoolToken(!result.hypre_solve_report.empty())
         << "; thermodynamic_recovery_report_present="
         << BoolToken(!result.thermodynamic_recovery_report.empty())
         << "; composition_model=equimolar_dt_from_p2_recovery"
         << "; separate_dt_species_authoritative=false"
         << "; fuel_depletion_enabled=false"
         << "; tau_alphae_model=thesis_spitzer_eq_5_262"
         << "; D_alpha_model=atzeni_drag_one_group"
         << "; reactivity_model_requested=bosch_hale_dt"
         << "; reactivity_model_executed=bosch_hale_dt"
         << "; bosch_hale_coefficients_source=project_owner_2026_04_28"
         << "; bosch_hale_temperature_source=Ti_old"
         << "; dt_reactivity_internal_unit=cm3_s"
         << "; owned_slab_coefficient_build_only=true"
         << "; coefficient_scope=owned_slab"
         << "; seam_coefficient_halo_exchanged=" << BoolToken(result.seam_coefficient_halo_exchanged)
         << "; off_rank_face_conductance_uses_neighbor_Dalpha="
         << BoolToken(result.off_rank_face_conductance_uses_neighbor_Dalpha)
         << "; global_off_rank_column_count=" << result.global_off_rank_column_count
         << "; owned_slab_writeback_only=" << BoolToken(result.owned_slab_writeback_only)
         << "; local_owned_cell_count=" << result.local_owned_cell_count
         << "; global_owned_cell_count=" << result.global_owned_cell_count
         << "; local_stage_ok=" << BoolToken(result.local_stage_ok)
         << "; global_stage_ok=" << BoolToken(result.global_stage_ok)
         << "; local_publishable=" << BoolToken(result.local_publishable)
         << "; global_publish_ok=" << BoolToken(result.global_publish_ok)
         << "; dense_fallback_used=false"
         << "; serial_dense_fallback_used=false"
         << "; rank0_gather_solve_used=false"
         << "; fallback_used=false"
         << "; solve_count=" << result.solve_count
         << "; hypre_solve_executed=" << BoolToken(result.hypre_solve_executed)
         << "; updated_fields=" << result.updated_fields_label
         << "; canonical_state_mutated=" << BoolToken(result.canonical_state_mutated)
         << "; metadata_written=" << BoolToken(result.metadata_written)
         << "; delta_alpha_total_global=" << result.delta_alpha_total_global
         << "; delta_electron_total_global=" << result.delta_electron_total_global
         << "; birth_source_total_global=" << result.birth_source_total_global
         << "; drag_deposition_total_global=" << result.drag_deposition_total_global
         << "; alpha_equation_budget_residual_global="
         << result.alpha_equation_budget_residual_global
         << "; alpha_electron_exchange_residual_global="
         << result.alpha_electron_exchange_residual_global
         << "; global_alpha_plus_electron_budget_residual="
         << result.global_alpha_plus_electron_budget_residual
         << "; global_budget_semantics=mpi_allreduced_volume_integral"
         << "; min_epsilon_alpha_before_local=" << result.min_epsilon_alpha_before_local
         << "; max_epsilon_alpha_before_local=" << result.max_epsilon_alpha_before_local
         << "; min_epsilon_alpha_after_local=" << result.min_epsilon_alpha_after_local
         << "; max_epsilon_alpha_after_local=" << result.max_epsilon_alpha_after_local
         << "; min_epsilon_alpha_after_global=" << result.min_epsilon_alpha_after_global
         << "; max_epsilon_alpha_after_global=" << result.max_epsilon_alpha_after_global
         << "; min_e_electron_after_global=" << result.min_e_electron_after_global
         << "; min_e_ion_after_global=" << result.min_e_ion_after_global
         << "; min_D_alpha_cm2_s=" << result.min_D_alpha_cm2_s
         << "; max_D_alpha_cm2_s=" << result.max_D_alpha_cm2_s
         << "; min_tau_alphae_s=" << result.min_tau_alphae_s
         << "; max_tau_alphae_s=" << result.max_tau_alphae_s
         << "; min_birth_source_erg_cm3_s=" << result.min_birth_source_erg_cm3_s
         << "; max_birth_source_erg_cm3_s=" << result.max_birth_source_erg_cm3_s
         << "; min_dt_reactivity_cm3_s=" << result.min_dt_reactivity_cm3_s
         << "; max_dt_reactivity_cm3_s=" << result.max_dt_reactivity_cm3_s;
  result.report_line = report.str();
}

DistributedOneGroupAlphaTransportResult Fail(
    DistributedOneGroupAlphaTransportResult result,
    const std::string& reason,
    MPI_Comm communicator,
    bool local_failed) {
  result.success = false;
  result.canonical_state_mutated = false;
  result.metadata_written = false;
  result.updated_fields = kDistributedAlphaWritesNoFields;
  result.updated_fields_label = "none";
  result.failure_reason = reason;
  const int local_rank = local_failed ? Rank(communicator) : std::numeric_limits<int>::max();
  result.first_failing_rank = AllreduceMinInt(communicator, local_rank);
  if (result.first_failing_rank == std::numeric_limits<int>::max()) {
    result.first_failing_rank = -1;
  }
  BuildReport(result, "p4.alpha.distributed_one_group_operator.failure");
  result.failure_diagnostics =
      result.report_line +
      "; first_failing_rank=" + std::to_string(result.first_failing_rank) +
      "; failure_reason=" + reason;
  return result;
}

struct LocalBudget {
  double delta_alpha{0.0};
  double delta_electron{0.0};
  double birth_source_dt{0.0};
  double drag_deposition{0.0};
  double min_alpha_before{std::numeric_limits<double>::infinity()};
  double max_alpha_before{-std::numeric_limits<double>::infinity()};
  double min_alpha_after{std::numeric_limits<double>::infinity()};
  double max_alpha_after{-std::numeric_limits<double>::infinity()};
};

}  // namespace

bool DistributedOneGroupAlphaTransportResult::is_complete() const noexcept {
  return success &&
         !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p4.alpha.distributed_one_group_operator") &&
         Contains(report_line, "phase_id=P4") &&
         Contains(report_line, "stage_id=A") &&
         Contains(report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(report_line, "updated_fields=") &&
         Contains(report_line, "global_budget_semantics=mpi_allreduced_volume_integral");
}

DistributedOneGroupAlphaTransportResult ApplyDistributedOneGroupAlphaTransport(
    DistributedOneGroupAlphaTransportProblem& problem,
    const DistributedOneGroupAlphaTransportOptions& options) noexcept {
  DistributedOneGroupAlphaTransportResult result;
  result.dt_s = problem.dt_s;
  result.backend_requested = "hypre_parcsr_gmres_boomeramg";

  const int local_has_state = problem.local_state == nullptr ? 0 : 1;
  const int global_has_state = AllreduceMinInt(problem.communicator, local_has_state);
  if (problem.local_state != nullptr) {
    result.local_owned_cell_count = LocalCellCount(*problem.local_state);
  }
  result.global_owned_cell_count =
      static_cast<std::size_t>(AllreduceSumDouble(
          problem.communicator,
          static_cast<double>(result.local_owned_cell_count)));
  if (global_has_state == 0) {
    return Fail(result, "distributed alpha local_state pointer is null", problem.communicator,
                local_has_state == 0);
  }

  if (problem.dt_s == 0.0) {
    result.success = true;
    result.local_stage_ok = true;
    result.global_stage_ok = AllreduceMinInt(problem.communicator, 1) != 0;
    result.local_publishable = true;
    result.global_publish_ok = AllreduceMinInt(problem.communicator, 1) != 0;
    result.canonical_state_mutated = false;
    result.metadata_written = false;
    result.updated_fields = kDistributedAlphaWritesNoFields;
    result.updated_fields_label = "none";
    result.backend_executed = "none";
    result.solve_count = 0u;
    result.distributed_alpha_solve_executed = false;
    result.hypre_solve_executed = false;
    result.provider_report =
        "diagnostic_id=p4.alpha.coefficients.provider_skipped; "
        "provider_skipped_for_no_change=true";
    result.thermodynamic_recovery_report =
        "diagnostic_id=p4.alpha.distributed_no_change_recovery; "
        "thermodynamic_recovery_skipped_for_no_change=true";
    BuildReport(result, "p4.alpha.distributed_one_group_operator");
    result.report_line += "; provider_skipped_for_no_change=true";
    return result;
  }

  auto provider_options = options.provider_options;
  provider_options.alpha_energy_floor_erg_cm3 =
      options.alpha_energy_floor_erg_per_cm3;
  auto recovery_options_for_provider = provider_options.recovery_options;
  recovery_options_for_provider.electron_energy_floor =
      options.electron_energy_floor_erg_per_cm3;
  recovery_options_for_provider.ion_energy_floor =
      options.ion_energy_floor_erg_per_cm3;
  provider_options.recovery_options = recovery_options_for_provider;

  const auto provider_timer = std::chrono::steady_clock::now();
  const auto provider = BuildAlphaCoefficientArrays(*problem.local_state, provider_options);
  result.coefficient_provider_wall_s = ElapsedSecondsSince(provider_timer);
  result.provider_report = provider.success ? provider.report_line : provider.failure_diagnostics;
  result.min_D_alpha_cm2_s = provider.min_D_alpha_cm2_s;
  result.max_D_alpha_cm2_s = provider.max_D_alpha_cm2_s;
  result.min_tau_alphae_s = provider.min_tau_alphae_s;
  result.max_tau_alphae_s = provider.max_tau_alphae_s;
  result.min_birth_source_erg_cm3_s = provider.min_birth_source_erg_cm3_s;
  result.max_birth_source_erg_cm3_s = provider.max_birth_source_erg_cm3_s;
  result.min_dt_reactivity_cm3_s = provider.min_dt_reactivity_cm3_s;
  result.max_dt_reactivity_cm3_s = provider.max_dt_reactivity_cm3_s;

  const bool local_provider_ok =
      provider.success && ValidateAlphaCoefficientDiagnostics(provider);
  const bool global_provider_ok =
      AllreduceMinInt(problem.communicator, local_provider_ok ? 1 : 0) != 0;
  if (!global_provider_ok) {
    return Fail(result, "distributed alpha provider failed", problem.communicator,
                !local_provider_ok);
  }

  auto diffusion = BuildAlphaDiffusionProblem(problem, provider, options);
  const auto assembly_timer = std::chrono::steady_clock::now();
  const auto assembly = dec3d::transport::AssembleDistributedGenericDiffusionSystem(diffusion);
  result.assembly_wall_s = ElapsedSecondsSince(assembly_timer);
  result.assembly_report = assembly.success ? assembly.report_line : assembly.failure_diagnostics;
  result.seam_coefficient_halo_exchanged =
      assembly.coefficient_D_halo_lower_received ||
      assembly.coefficient_D_halo_upper_received ||
      problem.ownership.rank_count == 1;
  result.off_rank_face_conductance_uses_neighbor_Dalpha =
      assembly.global_radial_seam_coupling_count > 0u ||
      problem.ownership.rank_count == 1;
  result.global_off_rank_column_count = assembly.global_off_rank_column_count;
  const bool local_assembly_ok =
      assembly.success &&
      dec3d::transport::ValidateDistributedGenericDiffusionAssemblyDiagnostics(assembly);
  const bool global_assembly_ok =
      AllreduceMinInt(problem.communicator, local_assembly_ok ? 1 : 0) != 0;
  if (!global_assembly_ok) {
    return Fail(result, "distributed alpha assembly failed", problem.communicator,
                !local_assembly_ok);
  }

  const auto solve =
      dec3d::transport::SolveDistributedGenericDiffusionHypre(assembly, options.solve_options);
  result.hypre_solve_report = solve.success ? solve.report_line : solve.failure_diagnostics;
  result.backend_executed = solve.backend;
  result.solve_count = 1u;
  result.distributed_alpha_solve_executed = true;
  result.hypre_solve_executed = true;
  result.hypre_setup_wall_s = solve.hypre_setup_wall_s;
  result.hypre_solve_wall_s = solve.hypre_solve_wall_s;
  result.solver_iterations = solve.gmres_iterations;
  const bool local_solve_ok =
      solve.success &&
      dec3d::transport::ValidateDistributedGenericDiffusionHypreSolveDiagnostics(solve);
  const bool global_solve_ok =
      AllreduceMinInt(problem.communicator, local_solve_ok ? 1 : 0) != 0;
  if (!global_solve_ok) {
    return Fail(result, "distributed alpha HYPRE solve failed", problem.communicator,
                !local_solve_ok);
  }

  const auto writeback_timer = std::chrono::steady_clock::now();
  dec3d::core::Array3D<double> staged_alpha = problem.local_state->alpha_state.storage;
  dec3d::core::Array3D<double> staged_electron = problem.local_state->e_electron;
  dec3d::core::Array3D<double> staged_total = problem.local_state->e_fluid_total;
  LocalBudget local_budget;
  bool local_staged_finite = true;
  std::size_t flat = 0u;
  for (std::size_t r = 0; r < staged_alpha.extent_r(); ++r) {
    for (std::size_t t = 0; t < staged_alpha.extent_theta(); ++t) {
      for (std::size_t p = 0; p < staged_alpha.extent_phi(); ++p) {
        const double old_alpha = problem.local_state->alpha_state.storage(r, t, p);
        const double new_alpha = solve.local_scalar_new[flat++];
        const double tau = provider.coefficients.tau_alphae_s(r, t, p);
        const double birth = provider.coefficients.birth_source_erg_cm3_s(r, t, p);
        const double drag = problem.dt_s * new_alpha / tau;
        const double volume = CellVolume(problem.global_geometry, problem.ownership, r, t, p);

        staged_alpha(r, t, p) = new_alpha;
        staged_electron(r, t, p) = problem.local_state->e_electron(r, t, p) + drag;
        staged_total(r, t, p) = problem.local_state->e_fluid_total(r, t, p) + drag;

        local_budget.delta_alpha += (new_alpha - old_alpha) * volume;
        local_budget.delta_electron += drag * volume;
        local_budget.birth_source_dt += problem.dt_s * birth * volume;
        local_budget.drag_deposition += drag * volume;
        local_budget.min_alpha_before = std::min(local_budget.min_alpha_before, old_alpha);
        local_budget.max_alpha_before = std::max(local_budget.max_alpha_before, old_alpha);
        local_budget.min_alpha_after = std::min(local_budget.min_alpha_after, new_alpha);
        local_budget.max_alpha_after = std::max(local_budget.max_alpha_after, new_alpha);

        local_staged_finite =
            local_staged_finite &&
            std::isfinite(new_alpha) &&
            new_alpha >= options.alpha_energy_floor_erg_per_cm3 &&
            std::isfinite(drag) &&
            std::isfinite(volume) &&
            volume > 0.0 &&
            std::isfinite(staged_electron(r, t, p)) &&
            staged_electron(r, t, p) >= options.electron_energy_floor_erg_per_cm3 &&
            std::isfinite(staged_total(r, t, p));
      }
    }
  }

  result.delta_alpha_total_global =
      AllreduceSumDouble(problem.communicator, local_budget.delta_alpha);
  result.delta_electron_total_global =
      AllreduceSumDouble(problem.communicator, local_budget.delta_electron);
  result.birth_source_total_global =
      AllreduceSumDouble(problem.communicator, local_budget.birth_source_dt);
  result.drag_deposition_total_global =
      AllreduceSumDouble(problem.communicator, local_budget.drag_deposition);
  result.alpha_equation_budget_residual_global =
      std::abs(result.delta_alpha_total_global -
               result.birth_source_total_global +
               result.drag_deposition_total_global);
  result.alpha_electron_exchange_residual_global =
      std::abs(result.delta_electron_total_global -
               result.drag_deposition_total_global);
  result.global_alpha_plus_electron_budget_residual =
      std::abs(result.delta_alpha_total_global +
               result.delta_electron_total_global -
               result.birth_source_total_global);
  result.min_epsilon_alpha_before_local = local_budget.min_alpha_before;
  result.max_epsilon_alpha_before_local = local_budget.max_alpha_before;
  result.min_epsilon_alpha_after_local = local_budget.min_alpha_after;
  result.max_epsilon_alpha_after_local = local_budget.max_alpha_after;
  result.min_epsilon_alpha_after_global =
      AllreduceMinDouble(problem.communicator, local_budget.min_alpha_after);
  result.max_epsilon_alpha_after_global =
      AllreduceMaxDouble(problem.communicator, local_budget.max_alpha_after);

  auto staged_state = *problem.local_state;
  staged_state.alpha_state.storage = staged_alpha;
  staged_state.e_electron = staged_electron;
  staged_state.e_fluid_total = staged_total;
  auto recovery_options = options.recovery_options;
  recovery_options.electron_energy_floor = options.electron_energy_floor_erg_per_cm3;
  recovery_options.ion_energy_floor = options.ion_energy_floor_erg_per_cm3;
  const auto recovered =
      dec3d::state::RecoverThermodynamicState(staged_state, recovery_options);
  result.thermodynamic_recovery_report =
      recovered.success ? recovered.recovery_diagnostics : recovered.failure_diagnostics;
  const bool local_recovery_ok =
      recovered.success &&
      dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovered);
  double local_min_electron = std::numeric_limits<double>::infinity();
  double local_min_ion = std::numeric_limits<double>::infinity();
  if (recovered.success) {
    for (std::size_t r = 0; r < recovered.cells.extent_r(); ++r) {
      for (std::size_t t = 0; t < recovered.cells.extent_theta(); ++t) {
        for (std::size_t p = 0; p < recovered.cells.extent_phi(); ++p) {
          local_min_electron =
              std::min(local_min_electron, recovered.cells(r, t, p).e_electron_erg_per_cm3);
          local_min_ion =
              std::min(local_min_ion, recovered.cells(r, t, p).e_ion_erg_per_cm3);
        }
      }
    }
  }
  result.min_e_electron_after_global =
      AllreduceMinDouble(problem.communicator, local_min_electron);
  result.min_e_ion_after_global =
      AllreduceMinDouble(problem.communicator, local_min_ion);

  const bool local_budget_ok =
      result.alpha_equation_budget_residual_global <=
          kBudgetTolerance * BudgetScale(result.birth_source_total_global) &&
      result.alpha_electron_exchange_residual_global <=
          kBudgetTolerance * BudgetScale(result.drag_deposition_total_global) &&
      result.global_alpha_plus_electron_budget_residual <=
          kBudgetTolerance * BudgetScale(result.birth_source_total_global);
  result.local_stage_ok = local_staged_finite && local_recovery_ok && local_budget_ok;
  result.global_stage_ok =
      AllreduceMinInt(problem.communicator, result.local_stage_ok ? 1 : 0) != 0;
  if (!result.global_stage_ok) {
    return Fail(result,
                local_recovery_ok ? "distributed alpha staged state is not feasible"
                                  : "distributed alpha staged thermodynamic recovery failed",
                problem.communicator,
                !result.local_stage_ok);
  }

  result.local_publishable =
      !options.force_publish_preflight_failure_for_test &&
      staged_alpha.extent_r() == problem.local_state->alpha_state.storage.extent_r() &&
      staged_electron.extent_r() == problem.local_state->e_electron.extent_r() &&
      staged_total.extent_r() == problem.local_state->e_fluid_total.extent_r();
  result.global_publish_ok =
      AllreduceMinInt(problem.communicator, result.local_publishable ? 1 : 0) != 0;
  if (!result.global_publish_ok) {
    return Fail(result,
                "distributed alpha publish preflight failed",
                problem.communicator,
                !result.local_publishable);
  }

  problem.local_state->alpha_state.storage = std::move(staged_alpha);
  problem.local_state->e_electron = std::move(staged_electron);
  problem.local_state->e_fluid_total = std::move(staged_total);
  problem.local_state->ApplyAuthoritativeWrite(DistributedAlphaWriteMask());
  result.writeback_wall_s = ElapsedSecondsSince(writeback_timer);

  result.success = true;
  result.canonical_state_mutated = true;
  result.metadata_written = true;
  result.updated_fields =
      kDistributedAlphaWritesAlphaState |
      kDistributedAlphaWritesElectronEnergy |
      kDistributedAlphaWritesFluidTotalEnergy;
  result.updated_fields_label = UpdatedFieldsLabel(result.updated_fields);
  BuildReport(result, "p4.alpha.distributed_one_group_operator");
  return result;
}

bool ValidateDistributedOneGroupAlphaDiagnostics(
    const DistributedOneGroupAlphaTransportResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p4.alpha.distributed_one_group_operator") &&
         Contains(line, "phase_id=P4") &&
         Contains(line, "stage_id=A") &&
         Contains(line, "alpha_transport_model=atzeni_one_group") &&
         Contains(line, "distributed_alpha_solve=") &&
         Contains(line, "backend_requested=hypre_parcsr_gmres_boomeramg") &&
         Contains(line, "coefficient_provider_wall_s=") &&
         Contains(line, "assembly_wall_s=") &&
         Contains(line, "hypre_setup_wall_s=") &&
         Contains(line, "hypre_solve_wall_s=") &&
         Contains(line, "writeback_wall_s=") &&
         Contains(line, "solver_iterations=") &&
         Contains(line, "composition_model=equimolar_dt_from_p2_recovery") &&
         Contains(line, "separate_dt_species_authoritative=false") &&
         Contains(line, "fuel_depletion_enabled=false") &&
         Contains(line, "tau_alphae_model=thesis_spitzer_eq_5_262") &&
         Contains(line, "D_alpha_model=atzeni_drag_one_group") &&
         Contains(line, "reactivity_model_requested=bosch_hale_dt") &&
         Contains(line, "reactivity_model_executed=bosch_hale_dt") &&
         Contains(line, "owned_slab_coefficient_build_only=true") &&
         Contains(line, "coefficient_scope=owned_slab") &&
         Contains(line, "seam_coefficient_halo_exchanged=") &&
         Contains(line, "off_rank_face_conductance_uses_neighbor_Dalpha=") &&
         Contains(line, "global_off_rank_column_count=") &&
         Contains(line, "local_owned_cell_count=") &&
         Contains(line, "global_owned_cell_count=") &&
         Contains(line, "local_stage_ok=") &&
         Contains(line, "global_stage_ok=") &&
         Contains(line, "local_publishable=") &&
         Contains(line, "global_publish_ok=") &&
         Contains(line, "updated_fields=") &&
         Contains(line, "canonical_state_mutated=") &&
         Contains(line, "metadata_written=") &&
         Contains(line, "solve_count=") &&
         Contains(line, "hypre_solve_executed=") &&
         Contains(line, "dense_fallback_used=false") &&
         Contains(line, "serial_dense_fallback_used=false") &&
         Contains(line, "rank0_gather_solve_used=false") &&
         Contains(line, "fallback_used=false") &&
         Contains(line, "global_budget_semantics=mpi_allreduced_volume_integral") &&
         line.find("per_group_solve_count") == std::string::npos;
}

}  // namespace dec3d::alpha
