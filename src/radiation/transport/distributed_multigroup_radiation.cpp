#include "radiation/transport/distributed_multigroup_radiation.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace dec3d::radiation {
namespace {

constexpr const char* kBackend = "hypre_parcsr_gmres_boomeramg";
constexpr std::size_t kNoFailingGroup = std::numeric_limits<std::size_t>::max();

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] const char* BoolToken(bool value) noexcept {
  return value ? "true" : "false";
}

[[nodiscard]] std::size_t FlatIndex(
    const dec3d::state::CanonicalStateLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return (radial * layout.theta_cells + theta) * layout.phi_cells + phi;
}

[[nodiscard]] int AllreduceMinBool(MPI_Comm communicator, bool value) noexcept {
  int local = value ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global;
}

[[nodiscard]] unsigned long long AllreduceSumUll(
    MPI_Comm communicator,
    unsigned long long value) noexcept {
  unsigned long long global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
  return global;
}

[[nodiscard]] double AllreduceSumDouble(MPI_Comm communicator, double value) noexcept {
  double global = 0.0;
  MPI_Allreduce(&value, &global, 1, MPI_DOUBLE, MPI_SUM, communicator);
  return global;
}

[[nodiscard]] bool AllreduceOrBool(MPI_Comm communicator, bool value) noexcept {
  int local = value ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, communicator);
  return global != 0;
}

[[nodiscard]] std::string UpdatedFieldsLabel(std::uint32_t mask) {
  if (mask == kDistributedRadiationWritesNoFields) {
    return "none";
  }
  const std::uint32_t expected =
      kDistributedRadiationWritesRadiationGroups |
      kDistributedRadiationWritesElectronEnergy |
      kDistributedRadiationWritesFluidTotalEnergy;
  return mask == expected ? "radiation_groups,e_electron,e_fluid_total" : "partial";
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask AuthoritativeMaskForWriteback() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy EffectiveBoundaryPolicy(
    const DistributedMultigroupRadiationProblem& problem,
    const DistributedMultigroupRadiationOptions& options) noexcept {
  auto policy = problem.boundary_policy;
  if (options.boundary_model == DistributedRadiationBoundaryModel::thesis_marshak_vacuum) {
    policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::radiation_marshak_vacuum;
  }
  return policy;
}

void BuildReport(DistributedMultigroupRadiationResult& result, const char* diagnostic_id) {
  const bool provider_reports_present =
      result.per_group_provider_reports.size() == result.group_count;
  const bool matrix_reports_present =
      result.per_group_matrix_reports.size() == result.group_count ||
      result.per_group_solve_count == 0u;
  const bool solve_reports_present =
      result.per_group_solve_reports.size() == result.group_count ||
      result.per_group_solve_count == 0u;
  const bool limiter_reports_present =
      !result.radiation_flux_limiter_enabled ||
      (!result.per_group_limiter_reports.empty() &&
       result.per_group_limiter_reports.size() == result.group_count);
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=" << diagnostic_id
      << "; stage_id=R"
      << "; radiation_model=multigroup_gray"
      << "; frozen_hydro=true"
      << "; provider_fed=true"
      << "; opacity_provider=tops_dt_tabulated"
      << "; opacity_data_source=TOPS"
      << "; opacity_source_thesis_exact_match=false"
      << "; group_layout_mode=explicit_frequency_groups"
      << "; frequency_edges_unit=Hz"
      << "; group_layout_report_present=true"
      << "; group_count=" << result.group_count
      << "; backend_requested=" << kBackend
      << "; backend_executed=" << (result.per_group_solve_count == 0u ? "none" : kBackend)
      << "; serial_for_groups=true"
      << "; distributed_per_group=true"
      << "; group_parallelism=disabled"
      << "; coefficient_time_level=old_time_lagged"
      << "; coefficient_provider_wall_s=" << result.coefficient_provider_wall_s
      << "; flux_limiter_wall_s=" << result.flux_limiter_wall_s
      << "; assembly_wall_s=" << result.assembly_wall_s
      << "; hypre_setup_wall_s=" << result.hypre_setup_wall_s
      << "; hypre_solve_wall_s=" << result.hypre_solve_wall_s
      << "; writeback_wall_s=" << result.writeback_wall_s
      << "; solver_iterations=" << result.solver_iterations
      << "; radiation_flux_limiter_enabled="
      << BoolToken(result.radiation_flux_limiter_enabled)
      << "; radiation_flux_limiter_model=" << result.radiation_flux_limiter_model
      << "; face_limiter_time_level=old_time_lagged"
      << "; U_floor_is_numerical_regularization=true"
      << "; limited_face_count_local=" << result.limited_face_count_local
      << "; limited_face_count_global=" << result.limited_face_count_global
      << "; min_limiter_scale=" << result.min_limiter_scale
      << "; max_flux_ratio_before_limit=" << result.max_flux_ratio_before_limit
      << "; per_group_limiter_report_present=" << BoolToken(limiter_reports_present)
      << "; outer_marshak_uses_face_effective_D="
      << BoolToken(result.outer_marshak_uses_face_effective_D)
      << "; table_cache_scope=explicit_handle"
      << "; table_handle_supplied=" << BoolToken(result.table_handle_supplied)
      << "; table_reused_across_steps=" << BoolToken(result.table_reused_across_steps)
      << "; per_step_csv_io=" << BoolToken(result.per_step_csv_io)
      << "; per_lookup_full_table_scan=" << BoolToken(result.per_lookup_full_table_scan)
      << "; coefficient_scope=owned_slab"
      << "; seam_coefficient_halo_exchanged=" << BoolToken(result.seam_coefficient_halo_exchanged)
      << "; off_rank_face_conductance_uses_neighbor_Dbar="
      << BoolToken(result.off_rank_face_conductance_uses_neighbor_Dbar)
      << "; global_off_rank_column_count=" << result.global_off_rank_column_count
      << "; per_group_provider_report_present=" << BoolToken(provider_reports_present)
      << "; per_group_hypre_matrix_report_present=" << BoolToken(matrix_reports_present)
      << "; per_group_hypre_solve_report_present=" << BoolToken(solve_reports_present)
      << "; per_group_budget_count=" << result.delta_radiation_total_by_group.size()
      << "; thermodynamic_recovery_report_present="
      << BoolToken(!result.thermodynamic_recovery_report.empty())
      << "; local_stage_ok=" << BoolToken(result.local_stage_ok)
      << "; global_stage_ok=" << BoolToken(result.global_stage_ok)
      << "; local_publishable=" << BoolToken(result.local_publishable)
      << "; global_publish_ok=" << BoolToken(result.global_publish_ok)
      << "; owned_slab_writeback_only=" << BoolToken(result.owned_slab_writeback_only)
      << "; updated_fields=" << result.updated_fields_label
      << "; canonical_state_mutated=" << BoolToken(result.canonical_state_mutated)
      << "; marshak_enabled=" << BoolToken(result.marshak_enabled)
      << "; marshak_only_on_outer_global_rank=true"
      << "; parity_claim_allowed=false"
      << "; multigroup_benchmark_claim_allowed=false"
      << "; per_group_solve_count=" << result.per_group_solve_count
      << "; local_owned_cell_count=" << result.local_owned_cell_count
      << "; global_owned_cell_count=" << result.global_owned_cell_count
      << "; first_failing_group_index="
      << (result.first_failing_group_index == kNoFailingGroup
              ? std::string{"none"}
              : std::to_string(result.first_failing_group_index))
      << "; first_failing_rank="
      << (result.first_failing_rank < 0 ? std::string{"none"}
                                        : std::to_string(result.first_failing_rank))
      << "; delta_radiation_total_all_groups="
      << result.delta_radiation_total_all_groups
      << "; delta_electron_total=" << result.delta_electron_total
      << "; source_gain_radiation_total_all_groups="
      << result.source_gain_radiation_total_all_groups
      << "; boundary_leak_total_all_groups=" << result.boundary_leak_total_all_groups
      << "; radiation_electron_exchange_residual="
      << result.radiation_electron_exchange_residual
      << "; global_radiation_plus_electron_residual="
      << result.global_radiation_plus_electron_residual
      << "; rank0_gather_solve_used=false"
      << "; fallback_used=false";
  result.report_line = out.str();
}

[[nodiscard]] DistributedMultigroupRadiationResult Fail(
    DistributedMultigroupRadiationResult result,
    const std::string& reason,
    MPI_Comm communicator,
    bool local_failed) {
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  const int local_failure_rank = local_failed ? rank : std::numeric_limits<int>::max();
  int first_rank = std::numeric_limits<int>::max();
  MPI_Allreduce(&local_failure_rank, &first_rank, 1, MPI_INT, MPI_MIN, communicator);
  result.first_failing_rank =
      first_rank == std::numeric_limits<int>::max() ? -1 : first_rank;
  result.success = false;
  result.canonical_state_mutated = false;
  result.updated_fields = kDistributedRadiationWritesNoFields;
  result.updated_fields_label = "none";
  result.failure_reason = reason;
  BuildReport(result, "p3.radiation.distributed_provider_fed_multigroup.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  if (!result.per_group_matrix_reports.empty()) {
    result.failure_diagnostics +=
        "; nested_matrix_diagnostics=" + result.per_group_matrix_reports.back();
  }
  if (!result.per_group_solve_reports.empty()) {
    result.failure_diagnostics +=
        "; nested_solve_diagnostics=" + result.per_group_solve_reports.back();
  }
  if (!result.per_group_limiter_reports.empty()) {
    result.failure_diagnostics +=
        "; nested_limiter_diagnostics=" + result.per_group_limiter_reports.back();
  }
  return result;
}

[[nodiscard]] dec3d::transport::DistributedGenericDiffusionProblem BuildGroupProblem(
    const DistributedMultigroupRadiationProblem& problem,
    const DistributedMultigroupRadiationOptions& options,
    const MultigroupGrayRadiationCoefficients& coefficients,
    std::size_t group,
    const dec3d::transport::GenericDiffusionFaceEffectiveCoefficients* face_effective) {
  dec3d::transport::DistributedGenericDiffusionProblem diffusion;
  diffusion.ownership = problem.ownership;
  diffusion.global_geometry = problem.global_geometry;
  diffusion.dt_s = problem.dt_s;
  diffusion.boundary_policy = EffectiveBoundaryPolicy(problem, options);
  diffusion.local_scalar_old = problem.local_state.radiation_groups[group];
  const auto& layout = problem.local_state.layout;
  diffusion.local_coefficient_A =
      dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 1.0);
  diffusion.local_coefficient_D = coefficients.Dbar_cm2_per_s[group];
  diffusion.local_coefficient_C =
      dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  diffusion.local_coefficient_B =
      dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const double kappa = coefficients.kappaP_cm_inv[group](r, t, p);
        const double b = coefficients.B_erg_per_cm3[group](r, t, p);
        diffusion.local_coefficient_C(r, t, p) = -c * kappa;
        diffusion.local_coefficient_B(r, t, p) = c * kappa * b;
      }
    }
  }
  if (face_effective != nullptr && face_effective->enabled) {
    diffusion.face_effective_coefficients = *face_effective;
  }
  return diffusion;
}

[[nodiscard]] bool LocalStateShapeMatches(
    const DistributedMultigroupRadiationProblem& problem) noexcept {
  const std::size_t expected_radial =
      problem.ownership.global_radial_end - problem.ownership.global_radial_begin;
  return problem.local_state.layout.radial_cells == expected_radial &&
         problem.local_state.layout.theta_cells == problem.ownership.global_theta_cells &&
         problem.local_state.layout.phi_cells == problem.ownership.global_phi_cells &&
         problem.local_state.radiation_groups.size() == problem.group_layout.group_count;
}

}  // namespace

bool DistributedMultigroupRadiationResult::is_complete() const noexcept {
  return success && !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p3.radiation.distributed_provider_fed_multigroup") &&
         Contains(report_line, "stage_id=R") &&
         Contains(report_line, "provider_fed=true") &&
         Contains(report_line, "owned_slab_writeback_only=true") &&
         Contains(report_line, "rank0_gather_solve_used=false");
}

DistributedMultigroupRadiationResult ApplyDistributedProviderFedMultigroupRadiation(
    DistributedMultigroupRadiationProblem& problem,
    const DistributedMultigroupRadiationOptions& options) noexcept {
  DistributedMultigroupRadiationResult result;
  const auto& ownership = problem.ownership;
  const MPI_Comm communicator = ownership.communicator;
  result.group_count = problem.group_layout.group_count;
  result.table_handle_supplied = problem.opacity_table != nullptr;
  result.table_reused_across_steps = options.table_reused_across_steps;
  result.per_step_csv_io = false;
  result.per_lookup_full_table_scan = false;
  result.marshak_enabled =
      options.boundary_model == DistributedRadiationBoundaryModel::thesis_marshak_vacuum;
  result.radiation_flux_limiter_enabled =
      options.radiation_flux_limiter.model != RadiationFluxLimiterModel::disabled;
  result.radiation_flux_limiter_model =
      RadiationFluxLimiterModelName(options.radiation_flux_limiter.model);
  result.outer_marshak_uses_face_effective_D =
      result.radiation_flux_limiter_enabled && result.marshak_enabled;
  result.local_owned_cell_count = ownership.local_row_count;
  result.global_owned_cell_count = static_cast<std::size_t>(
      AllreduceSumUll(communicator, static_cast<unsigned long long>(ownership.local_row_count)));
  result.first_failing_group_index = kNoFailingGroup;

  bool local_preflight =
      ownership.success &&
      problem.opacity_table != nullptr &&
      LocalStateShapeMatches(problem) &&
      problem.global_geometry.is_valid() &&
      problem.group_layout.mode == RadiationGroupMode::explicit_frequency_groups &&
      std::isfinite(problem.dt_s) &&
      problem.dt_s >= 0.0;
  result.local_stage_ok = local_preflight;
  result.global_stage_ok = AllreduceMinBool(communicator, local_preflight) != 0;
  if (!result.global_stage_ok) {
    return Fail(result,
                "distributed radiation preflight failed",
                communicator,
                !local_preflight);
  }

  if (problem.dt_s == 0.0) {
    result.success = true;
    result.local_stage_ok = true;
    result.global_stage_ok = true;
    result.local_publishable = true;
    result.global_publish_ok = true;
    result.owned_slab_writeback_only = true;
    result.updated_fields = kDistributedRadiationWritesNoFields;
    result.updated_fields_label = "none";
    result.canonical_state_mutated = false;
    result.per_group_solve_count = 0u;
    result.per_group_provider_reports.assign(
        result.group_count, "provider_skipped=true; reason=no_change_contract_path");
    BuildReport(result, "p3.radiation.distributed_provider_fed_multigroup");
    return result;
  }

  const auto provider_timer = std::chrono::steady_clock::now();
  auto provider = BuildRadiationCoefficientArrays(
      problem.local_state,
      problem.global_geometry,
      problem.group_layout,
      *problem.opacity_table,
      options.provider_options);
  result.coefficient_provider_wall_s = ElapsedSecondsSince(provider_timer);
  result.coefficient_report = provider.success ? provider.report_line
                                               : provider.failure_diagnostics;
  result.per_group_provider_reports.assign(result.group_count, result.coefficient_report);
  bool local_provider_ok = provider.success &&
                           provider.coefficients.Dbar_cm2_per_s.size() == result.group_count;
  if (!local_provider_ok && provider.first_bad_group != std::numeric_limits<std::size_t>::max()) {
    result.first_failing_group_index = provider.first_bad_group;
  }
  result.local_stage_ok = local_provider_ok;
  result.global_stage_ok = AllreduceMinBool(communicator, local_provider_ok) != 0;
  if (!result.global_stage_ok) {
    return Fail(result,
                provider.failure_reason.empty() ? "distributed radiation coefficient stage failed"
                                                : provider.failure_reason,
                communicator,
                !local_provider_ok);
  }

  std::vector<dec3d::core::Array3D<double>> staged_groups =
      problem.local_state.radiation_groups;
  result.per_group_matrix_reports.clear();
  result.per_group_solve_reports.clear();

  bool any_halo = false;
  std::size_t max_global_off_rank_columns = 0u;
  bool local_group_ok = true;

  for (std::size_t group = 0; group < result.group_count; ++group) {
    dec3d::transport::GenericDiffusionFaceEffectiveCoefficients face_effective;
    const dec3d::transport::GenericDiffusionFaceEffectiveCoefficients* face_ptr = nullptr;
    if (result.radiation_flux_limiter_enabled) {
      const auto limiter_timer = std::chrono::steady_clock::now();
      const auto limited_faces = BuildRadiationFluxLimitedFaceCoefficients(
          ownership,
          problem.global_geometry,
          problem.local_state.radiation_groups[group],
          provider.coefficients.Dbar_cm2_per_s[group],
          options.radiation_flux_limiter);
      result.flux_limiter_wall_s += ElapsedSecondsSince(limiter_timer);
      const bool local_limiter_ok =
          limited_faces.success &&
          ValidateRadiationFluxLimitedFaceDiagnostics(limited_faces);
      const bool global_limiter_ok =
          AllreduceMinBool(communicator, local_limiter_ok) != 0;
      if (!global_limiter_ok) {
        result.first_failing_group_index = group;
        result.global_stage_ok = false;
        result.per_group_limiter_reports.push_back(
            limited_faces.success ? limited_faces.report_line
                                  : limited_faces.failure_diagnostics);
        return Fail(result,
                    limited_faces.failure_reason.empty()
                        ? "distributed radiation flux limiter failed"
                        : limited_faces.failure_reason,
                    communicator,
                    !local_limiter_ok);
      }
      result.per_group_limiter_reports.push_back(limited_faces.report_line);
      result.limited_face_count_local += limited_faces.limited_face_count_local;
      result.min_limiter_scale =
          std::min(result.min_limiter_scale, limited_faces.min_limiter_scale);
      result.max_flux_ratio_before_limit =
          std::max(result.max_flux_ratio_before_limit,
                   limited_faces.max_flux_ratio_before_limit);
      face_effective = limited_faces.face_coefficients;
      face_ptr = &face_effective;
    }

    auto diffusion =
        BuildGroupProblem(problem, options, provider.coefficients, group, face_ptr);
    const auto assembly_timer = std::chrono::steady_clock::now();
    const auto assembly = dec3d::transport::AssembleDistributedGenericDiffusionSystem(diffusion);
    result.assembly_wall_s += ElapsedSecondsSince(assembly_timer);
    local_group_ok = assembly.success &&
                     dec3d::transport::ValidateDistributedGenericDiffusionAssemblyDiagnostics(assembly);
    const bool global_assembly_ok = AllreduceMinBool(communicator, local_group_ok) != 0;
    if (!global_assembly_ok) {
      result.first_failing_group_index = group;
      result.global_stage_ok = false;
      result.per_group_matrix_reports.push_back(
          assembly.success ? assembly.report_line : assembly.failure_diagnostics);
      return Fail(result,
                  "distributed radiation group assembly failed",
                  communicator,
                  !local_group_ok);
    }
    result.per_group_matrix_reports.push_back(assembly.report_line);
    any_halo = any_halo ||
               assembly.coefficient_D_halo_lower_received ||
               assembly.coefficient_D_halo_upper_received;
    max_global_off_rank_columns =
        std::max(max_global_off_rank_columns, assembly.global_off_rank_column_count);

    auto solve_options = options.solve_options;
    solve_options.communicator = communicator;
    const auto solve =
        dec3d::transport::SolveDistributedGenericDiffusionHypre(assembly, solve_options);
    result.hypre_setup_wall_s += solve.hypre_setup_wall_s;
    result.hypre_solve_wall_s += solve.hypre_solve_wall_s;
    result.solver_iterations += solve.gmres_iterations;
    local_group_ok = solve.success &&
                     dec3d::transport::ValidateDistributedGenericDiffusionHypreSolveDiagnostics(solve) &&
                     solve.local_scalar_new.size() == ownership.local_row_count;
    const bool global_solve_ok = AllreduceMinBool(communicator, local_group_ok) != 0;
    if (!global_solve_ok) {
      result.first_failing_group_index = group;
      result.global_stage_ok = false;
      result.per_group_solve_reports.push_back(
          solve.success ? solve.report_line : solve.failure_diagnostics);
      return Fail(result,
                  "distributed radiation group HYPRE solve failed",
                  communicator,
                  !local_group_ok);
    }
    result.per_group_solve_reports.push_back(solve.report_line);
    ++result.per_group_solve_count;

    std::size_t offset = 0u;
    for (std::size_t r = 0; r < problem.local_state.layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < problem.local_state.layout.theta_cells; ++t) {
        for (std::size_t p = 0; p < problem.local_state.layout.phi_cells; ++p) {
          const double value = solve.local_scalar_new[offset++];
          if (!std::isfinite(value) || value < options.radiation_energy_floor_erg_per_cm3) {
            result.first_failing_group_index = group;
            result.global_stage_ok = false;
            return Fail(result,
                        "distributed radiation group solution violated floor",
                        communicator,
                        true);
          }
          staged_groups[group](r, t, p) = value;
        }
      }
    }
  }

  result.seam_coefficient_halo_exchanged = AllreduceOrBool(communicator, any_halo);
  result.global_off_rank_column_count = max_global_off_rank_columns;
  result.off_rank_face_conductance_uses_neighbor_Dbar =
      result.seam_coefficient_halo_exchanged &&
      result.global_off_rank_column_count > 0u;
  result.limited_face_count_global = static_cast<std::size_t>(
      AllreduceSumUll(communicator,
                      static_cast<unsigned long long>(
                          result.limited_face_count_local)));
  double global_min_limiter_scale = result.min_limiter_scale;
  double global_max_flux_ratio = result.max_flux_ratio_before_limit;
  MPI_Allreduce(&result.min_limiter_scale,
                &global_min_limiter_scale,
                1,
                MPI_DOUBLE,
                MPI_MIN,
                communicator);
  MPI_Allreduce(&result.max_flux_ratio_before_limit,
                &global_max_flux_ratio,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                communicator);
  result.min_limiter_scale = global_min_limiter_scale;
  result.max_flux_ratio_before_limit = global_max_flux_ratio;

  const auto writeback_timer = std::chrono::steady_clock::now();
  auto staged_e = problem.local_state.e_electron;
  auto staged_total = problem.local_state.e_fluid_total;
  bool local_staged_finite = true;
  std::vector<double> local_delta_radiation_by_group(result.group_count, 0.0);
  std::vector<double> local_source_gain_by_group(result.group_count, 0.0);
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  for (std::size_t r = 0; r < problem.local_state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < problem.local_state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < problem.local_state.layout.phi_cells; ++p) {
        const std::size_t local = FlatIndex(problem.local_state.layout, r, t, p);
        const std::size_t global = ownership.local_row_begin + local;
        const double volume = problem.global_geometry.cell_volumes[global];
        double delta_electron = 0.0;
        for (std::size_t group = 0; group < result.group_count; ++group) {
          const double old_u = problem.local_state.radiation_groups[group](r, t, p);
          const double new_u = staged_groups[group](r, t, p);
          const double source = problem.dt_s * c *
                                provider.coefficients.kappaP_cm_inv[group](r, t, p) *
                                (provider.coefficients.B_erg_per_cm3[group](r, t, p) - new_u);
          local_delta_radiation_by_group[group] += (new_u - old_u) * volume;
          local_source_gain_by_group[group] += source * volume;
          delta_electron -= source;
        }
        staged_e(r, t, p) = problem.local_state.e_electron(r, t, p) + delta_electron;
        staged_total(r, t, p) = problem.local_state.e_fluid_total(r, t, p) + delta_electron;
        local_staged_finite =
            local_staged_finite &&
            std::isfinite(staged_e(r, t, p)) &&
            std::isfinite(staged_total(r, t, p)) &&
            staged_e(r, t, p) >= options.electron_energy_floor_erg_per_cm3;
      }
    }
  }

  result.delta_radiation_total_by_group.assign(result.group_count, 0.0);
  result.source_gain_radiation_total_by_group.assign(result.group_count, 0.0);
  result.boundary_leak_total_by_group.assign(result.group_count, 0.0);
  for (std::size_t group = 0; group < result.group_count; ++group) {
    result.delta_radiation_total_by_group[group] =
        AllreduceSumDouble(communicator, local_delta_radiation_by_group[group]);
    result.source_gain_radiation_total_by_group[group] =
        AllreduceSumDouble(communicator, local_source_gain_by_group[group]);
    result.boundary_leak_total_by_group[group] =
        result.source_gain_radiation_total_by_group[group] -
        result.delta_radiation_total_by_group[group];
    result.delta_radiation_total_all_groups +=
        result.delta_radiation_total_by_group[group];
    result.source_gain_radiation_total_all_groups +=
        result.source_gain_radiation_total_by_group[group];
    result.boundary_leak_total_all_groups +=
        result.boundary_leak_total_by_group[group];
  }
  result.delta_electron_total = -result.source_gain_radiation_total_all_groups;
  result.radiation_electron_exchange_residual =
      std::abs(result.delta_electron_total +
               result.source_gain_radiation_total_all_groups);
  result.global_radiation_plus_electron_residual =
      std::abs(result.delta_radiation_total_all_groups +
               result.delta_electron_total +
               result.boundary_leak_total_all_groups);

  auto staged_state = problem.local_state;
  staged_state.radiation_groups = staged_groups;
  staged_state.e_electron = staged_e;
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
  const bool local_stage_ok = local_staged_finite && local_recovery_ok;
  result.local_stage_ok = local_stage_ok;
  result.global_stage_ok = AllreduceMinBool(communicator, local_stage_ok) != 0;
  if (!result.global_stage_ok) {
    return Fail(result,
                local_recovery_ok ? "distributed radiation staged state is not feasible"
                                  : "distributed radiation staged thermodynamic recovery failed",
                communicator,
                !local_stage_ok);
  }

  result.local_publishable =
      staged_groups.size() == result.group_count &&
      staged_e.extent_r() == problem.local_state.e_electron.extent_r() &&
      staged_total.extent_r() == problem.local_state.e_fluid_total.extent_r();
  result.global_publish_ok = AllreduceMinBool(communicator, result.local_publishable) != 0;
  if (!result.global_publish_ok) {
    return Fail(result,
                "distributed radiation publish preflight failed",
                communicator,
                !result.local_publishable);
  }

  problem.local_state.radiation_groups = std::move(staged_groups);
  problem.local_state.e_electron = std::move(staged_e);
  problem.local_state.e_fluid_total = std::move(staged_total);
  problem.local_state.ApplyAuthoritativeWrite(AuthoritativeMaskForWriteback());
  result.writeback_wall_s = ElapsedSecondsSince(writeback_timer);

  result.success = true;
  result.canonical_state_mutated = true;
  result.updated_fields =
      kDistributedRadiationWritesRadiationGroups |
      kDistributedRadiationWritesElectronEnergy |
      kDistributedRadiationWritesFluidTotalEnergy;
  result.updated_fields_label = UpdatedFieldsLabel(result.updated_fields);
  result.owned_slab_writeback_only = true;
  BuildReport(result, "p3.radiation.distributed_provider_fed_multigroup");
  return result;
}

bool ValidateDistributedMultigroupRadiationDiagnostics(
    const DistributedMultigroupRadiationResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.distributed_provider_fed_multigroup") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "radiation_model=multigroup_gray") &&
         Contains(line, "frozen_hydro=true") &&
         Contains(line, "provider_fed=true") &&
         Contains(line, "opacity_provider=tops_dt_tabulated") &&
         Contains(line, "opacity_source_thesis_exact_match=false") &&
         Contains(line, "group_layout_mode=explicit_frequency_groups") &&
         Contains(line, "frequency_edges_unit=Hz") &&
         Contains(line, "group_layout_report_present=true") &&
         Contains(line, "backend_requested=hypre_parcsr_gmres_boomeramg") &&
         Contains(line, "coefficient_provider_wall_s=") &&
         Contains(line, "flux_limiter_wall_s=") &&
         Contains(line, "assembly_wall_s=") &&
         Contains(line, "hypre_setup_wall_s=") &&
         Contains(line, "hypre_solve_wall_s=") &&
         Contains(line, "writeback_wall_s=") &&
         Contains(line, "solver_iterations=") &&
         Contains(line, "radiation_flux_limiter_enabled=") &&
         Contains(line, "radiation_flux_limiter_model=") &&
         Contains(line, "face_limiter_time_level=old_time_lagged") &&
         Contains(line, "U_floor_is_numerical_regularization=true") &&
         Contains(line, "limited_face_count_global=") &&
         Contains(line, "min_limiter_scale=") &&
         Contains(line, "max_flux_ratio_before_limit=") &&
         Contains(line, "per_group_limiter_report_present=") &&
         Contains(line, "outer_marshak_uses_face_effective_D=") &&
         Contains(line, "serial_for_groups=true") &&
         Contains(line, "distributed_per_group=true") &&
         Contains(line, "table_cache_scope=explicit_handle") &&
         Contains(line, "table_handle_supplied=true") &&
         Contains(line, "table_reused_across_steps=") &&
         Contains(line, "per_step_csv_io=false") &&
         Contains(line, "per_lookup_full_table_scan=false") &&
         Contains(line, "coefficient_scope=owned_slab") &&
         Contains(line, "seam_coefficient_halo_exchanged=") &&
         Contains(line, "off_rank_face_conductance_uses_neighbor_Dbar=") &&
         Contains(line, "per_group_provider_report_present=true") &&
         Contains(line, "per_group_hypre_matrix_report_present=true") &&
         Contains(line, "per_group_hypre_solve_report_present=true") &&
         Contains(line, "per_group_budget_count=") &&
         Contains(line, "thermodynamic_recovery_report_present=true") &&
         Contains(line, "local_stage_ok=true") &&
         Contains(line, "global_stage_ok=true") &&
         Contains(line, "local_publishable=true") &&
         Contains(line, "global_publish_ok=true") &&
         Contains(line, "owned_slab_writeback_only=true") &&
         Contains(line, "updated_fields=") &&
         Contains(line, "canonical_state_mutated=") &&
         Contains(line, "marshak_only_on_outer_global_rank=true") &&
         Contains(line, "parity_claim_allowed=false") &&
         Contains(line, "multigroup_benchmark_claim_allowed=false") &&
         Contains(line, "rank0_gather_solve_used=false") &&
         Contains(line, "fallback_used=false");
}

}  // namespace dec3d::radiation
