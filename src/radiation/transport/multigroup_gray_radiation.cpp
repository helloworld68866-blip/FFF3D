#include "radiation/transport/multigroup_gray_radiation.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace dec3d::radiation {
namespace {

constexpr const char* kImplementationId =
    "p3.radiation.multigroup_gray_matter_coupling_v1";
constexpr std::size_t kNoFailingGroup = std::numeric_limits<std::size_t>::max();

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] std::string UpdatedFieldsName(
    dec3d::core::AuthoritativeFieldMask mask) {
  return mask == MultigroupRadiationMatterWriteMask()
             ? "radiation_groups,e_electron,e_fluid_total"
             : "none";
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy EffectiveBoundaryPolicy(
    const MultigroupGrayRadiationMatterCouplingOptions& options) noexcept {
  auto policy = options.boundary_policy;
  if (options.radiation_boundary_model ==
      OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum) {
    policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::radiation_marshak_vacuum;
  }
  return policy;
}

[[nodiscard]] bool SameShape(
    const dec3d::core::Array3D<double>& values,
    const dec3d::state::CanonicalStateLayout& layout) noexcept {
  return values.extent_r() == layout.radial_cells &&
         values.extent_theta() == layout.theta_cells &&
         values.extent_phi() == layout.phi_cells;
}

[[nodiscard]] bool CoefficientsHaveGroupCount(
    const MultigroupGrayRadiationCoefficients& coefficients,
    std::size_t group_count) noexcept {
  return coefficients.Dbar_cm2_per_s.size() == group_count &&
         coefficients.kappaP_cm_inv.size() == group_count &&
         coefficients.B_erg_per_cm3.size() == group_count;
}

[[nodiscard]] bool CoefficientGroupIsValid(
    const MultigroupGrayRadiationCoefficients& coefficients,
    std::size_t group,
    const dec3d::state::CanonicalStateLayout& layout) noexcept {
  if (!SameShape(coefficients.Dbar_cm2_per_s[group], layout) ||
      !SameShape(coefficients.kappaP_cm_inv[group], layout) ||
      !SameShape(coefficients.B_erg_per_cm3[group], layout)) {
    return false;
  }
  for (std::size_t r = 0u; r < layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < layout.phi_cells; ++phi) {
        const double d = coefficients.Dbar_cm2_per_s[group](r, theta, phi);
        const double k = coefficients.kappaP_cm_inv[group](r, theta, phi);
        const double b = coefficients.B_erg_per_cm3[group](r, theta, phi);
        if (!std::isfinite(d) || !std::isfinite(k) || !std::isfinite(b) ||
            d < 0.0 || k < 0.0) {
          return false;
        }
      }
    }
  }
  return true;
}

[[nodiscard]] double BudgetScale(double value) noexcept {
  return std::max(1.0, std::abs(value));
}

void BuildReport(
    MultigroupGrayRadiationMatterCouplingResult& result,
    bool canonical_state_mutated,
    const char* diagnostic_id =
        "p3.radiation.multigroup_gray_matter_coupling") {
  const bool per_group_reports_present =
      result.per_group_radiation_reports.size() == result.group_count;
  std::ostringstream out;
  out << "diagnostic_id=" << diagnostic_id
      << "; stage_id=R"
      << "; implementation_id=" << kImplementationId
      << "; radiation_electron_coupling_enabled=true"
      << "; radiation_group_mode=explicit_frequency_groups"
      << "; group_count=" << result.group_count
      << "; frequency_edges_used=true"
      << "; frequency_edges_unit=Hz"
      << "; group_layout_report_present=" << BoolToken(!result.group_layout_report.empty())
      << "; coefficient_time_level=old_time_lagged"
      << "; per_group_coefficients_source=" << result.per_group_coefficients_source
      << "; per_group_solve_count=" << result.per_group_solve_count
      << "; all_groups_updated=" << BoolToken(result.all_groups_updated)
      << "; per_group_matrix_report_present=" << BoolToken(per_group_reports_present)
      << "; per_group_solve_report_present=" << BoolToken(per_group_reports_present)
      << "; per_group_backend_executed=serial_dense_reference"
      << "; backend_requested=" << result.backend_requested
      << "; backend_executed=" << result.backend_executed
      << "; boundary_model="
      << (result.marshak_enabled ? "thesis_marshak_vacuum"
                                 : "contract_zero_flux_or_scalar_remap")
      << "; marshak_enabled=" << BoolToken(result.marshak_enabled)
      << "; electron_writeback=true"
      << "; updated_fields=" << UpdatedFieldsName(result.updated_fields)
      << "; opacity_provider=" << result.opacity_provider
      << "; Bg_source=" << result.Bg_source
      << "; Bg_feedback=disabled"
      << "; radiation_flux_limiter=disabled"
      << "; advection_enabled=false"
      << "; radiation_pressure_work_enabled=false"
      << "; distributed_multigroup_enabled=false"
      << "; group_parallelism_enabled=false"
      << "; parity_claim_allowed=false"
      << "; first_failing_group_index="
      << (result.first_failing_group_index == kNoFailingGroup
              ? std::string{"none"}
              : std::to_string(result.first_failing_group_index))
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
      << "; canonical_state_mutated=" << BoolToken(canonical_state_mutated);
  result.report_line = out.str();
}

[[nodiscard]] MultigroupGrayRadiationMatterCouplingResult Fail(
    MultigroupGrayRadiationMatterCouplingResult result,
    const std::string& reason) {
  const std::string nested = result.failure_diagnostics;
  result.success = false;
  result.updated_fields = 0u;
  result.all_groups_updated = false;
  result.failure_reason = reason;
  BuildReport(result, false, "p3.radiation.multigroup_gray_matter_coupling.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  if (!nested.empty()) {
    result.failure_diagnostics += "; nested_failure_diagnostics={" + nested + "}";
  }
  return result;
}

[[nodiscard]] dec3d::transport::GenericDiffusionProblem BuildGroupProblem(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const MultigroupGrayRadiationMatterCouplingOptions& options,
    std::size_t group) {
  dec3d::transport::GenericDiffusionProblem problem;
  problem.layout = dec3d::transport::DiffusionGridLayout{
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells};
  problem.geometry = geometry;
  problem.dt_s = options.dt_s;
  problem.boundary_policy = EffectiveBoundaryPolicy(options);
  problem.scalar_old = state.radiation_groups[group];
  problem.coefficient_A = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      1.0);
  problem.coefficient_D = options.coefficients.Dbar_cm2_per_s[group];
  problem.coefficient_C = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      0.0);
  problem.coefficient_B = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      0.0);

  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double kappa = options.coefficients.kappaP_cm_inv[group](r, theta, phi);
        const double b = options.coefficients.B_erg_per_cm3[group](r, theta, phi);
        problem.coefficient_C(r, theta, phi) = -c * kappa;
        problem.coefficient_B(r, theta, phi) = c * kappa * b;
      }
    }
  }
  return problem;
}

}  // namespace

bool MultigroupGrayRadiationMatterCouplingResult::is_complete() const noexcept {
  return success && !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p3.radiation.multigroup_gray_matter_coupling") &&
         Contains(report_line, "stage_id=R") &&
         Contains(report_line, "radiation_group_mode=explicit_frequency_groups") &&
         Contains(report_line, "updated_fields=");
}

dec3d::core::AuthoritativeFieldMask MultigroupRadiationMatterWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

MultigroupGrayRadiationMatterCouplingResult
ApplyMultigroupGrayRadiationMatterCoupling(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const MultigroupGrayRadiationMatterCouplingOptions& options) noexcept {
  MultigroupGrayRadiationMatterCouplingResult result;
  result.group_count = options.group_layout.group_count;
  result.backend_requested = "serial_dense_reference";
  result.backend_executed = "serial_dense_reference";
  result.per_group_coefficients_source = options.coefficients.source.empty()
                                             ? std::string{"missing"}
                                             : options.coefficients.source;
  if (result.per_group_coefficients_source == "tops_dt_tabulated") {
    result.opacity_provider = "tops_dt_tabulated";
    result.Bg_source = "blackbody_group_integral";
  }
  result.marshak_enabled =
      options.radiation_boundary_model ==
      OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;

  const auto layout_result = ValidateRadiationGroupStateStorage(options.group_layout, state);
  result.group_layout_report = layout_result.report_line;
  if (!layout_result.success) {
    result.failure_diagnostics = layout_result.report_line;
    return Fail(result, layout_result.failure_reason);
  }
  if (options.group_layout.mode != RadiationGroupMode::explicit_frequency_groups) {
    return Fail(result, "P3-3 requires explicit frequency radiation groups");
  }
  if (!geometry.is_valid()) {
    return Fail(result, "multigroup radiation coupling requires valid spherical geometry");
  }
  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    return Fail(result, "multigroup radiation coupling dt_s must be finite and nonnegative");
  }
  if (!std::isfinite(options.radiation_energy_floor_erg_per_cm3) ||
      !std::isfinite(options.electron_energy_floor_erg_per_cm3) ||
      !std::isfinite(options.ion_energy_floor_erg_per_cm3)) {
    return Fail(result, "multigroup radiation floors must be finite");
  }
  if (!CoefficientsHaveGroupCount(options.coefficients, result.group_count)) {
    return Fail(result, "multigroup radiation coefficient group count mismatch");
  }
  if (options.coefficients.source != "fixed_user_supplied" &&
      options.coefficients.source != "tops_dt_tabulated") {
    return Fail(result, "multigroup radiation coefficient source is unsupported");
  }
  for (std::size_t g = 0u; g < result.group_count; ++g) {
    if (!CoefficientGroupIsValid(options.coefficients, g, state.layout)) {
      result.first_failing_group_index = g;
      return Fail(result, "multigroup radiation coefficient group is invalid");
    }
  }

  if (options.dt_s == 0.0) {
    result.success = true;
    result.updated_fields = 0u;
    result.first_failing_group_index = kNoFailingGroup;
    result.per_group_radiation_reports.assign(
        result.group_count,
        "assembly_skipped=true; solve_skipped=true; reason=no_change_contract_path");
    BuildReport(result, false);
    return result;
  }

  auto staged_state = state;
  result.delta_radiation_total_by_group.assign(result.group_count, 0.0);
  result.source_gain_radiation_total_by_group.assign(result.group_count, 0.0);
  result.boundary_leak_total_by_group.assign(result.group_count, 0.0);
  result.per_group_radiation_reports.clear();
  result.per_group_radiation_reports.reserve(result.group_count);

  std::vector<dec3d::core::Array3D<double>> staged_groups = state.radiation_groups;

  for (std::size_t g = 0u; g < result.group_count; ++g) {
    const auto problem = BuildGroupProblem(state, geometry, options, g);
    const auto assembly = dec3d::transport::AssembleGenericImplicitDiffusionSystem(problem);
    if (!assembly.success ||
        !dec3d::transport::ValidateGenericDiffusionAssemblyDiagnostics(assembly)) {
      result.first_failing_group_index = g;
      result.failure_diagnostics = assembly.failure_diagnostics.empty()
                                       ? assembly.report_line
                                       : assembly.failure_diagnostics;
      return Fail(result, "multigroup radiation group assembly failed");
    }

    const auto solve =
        dec3d::transport::SolveGenericDiffusionReference(assembly, options.serial_reference_options);
    if (!solve.success ||
        !dec3d::transport::ValidateGenericDiffusionSolveDiagnostics(solve)) {
      result.first_failing_group_index = g;
      result.failure_diagnostics = solve.failure_diagnostics.empty()
                                       ? solve.report_line
                                       : solve.failure_diagnostics;
      return Fail(result, "multigroup radiation group solve failed");
    }
    if (solve.scalar_new.size() != state.radiation_groups[g].size()) {
      result.first_failing_group_index = g;
      return Fail(result, "multigroup radiation group solve returned wrong scalar size");
    }

    result.per_group_radiation_reports.push_back(
        assembly.report_line + "; nested_solve_report={" + solve.report_line + "}");
    ++result.per_group_solve_count;

    std::size_t offset = 0u;
    for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
      for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
          const double value = solve.scalar_new[offset++];
          if (!std::isfinite(value) ||
              value < options.radiation_energy_floor_erg_per_cm3) {
            result.first_failing_group_index = g;
            return Fail(result, "multigroup radiation group output violated floor");
          }
          staged_groups[g](r, theta, phi) = value;
        }
      }
    }
  }

  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  bool valid = true;
  std::size_t linear = 0u;

  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double volume = geometry.cell_volumes[linear++];
        double delta_electron_cell = 0.0;
        valid = valid && std::isfinite(volume) && volume > 0.0;

        for (std::size_t g = 0u; g < result.group_count; ++g) {
          const double u_old = state.radiation_groups[g](r, theta, phi);
          const double u_new = staged_groups[g](r, theta, phi);
          const double kappa = options.coefficients.kappaP_cm_inv[g](r, theta, phi);
          const double b = options.coefficients.B_erg_per_cm3[g](r, theta, phi);
          const double source_gain = options.dt_s * c * kappa * (b - u_new);
          const double delta_u = u_new - u_old;

          result.delta_radiation_total_by_group[g] += delta_u * volume;
          result.source_gain_radiation_total_by_group[g] += source_gain * volume;
          delta_electron_cell -= source_gain;
          valid = valid && std::isfinite(u_old) && std::isfinite(u_new) &&
                  std::isfinite(source_gain);
        }

        staged_state.e_electron(r, theta, phi) =
            state.e_electron(r, theta, phi) + delta_electron_cell;
        staged_state.e_fluid_total(r, theta, phi) =
            state.e_fluid_total(r, theta, phi) + delta_electron_cell;
        valid = valid &&
                std::isfinite(staged_state.e_electron(r, theta, phi)) &&
                std::isfinite(staged_state.e_fluid_total(r, theta, phi)) &&
                staged_state.e_electron(r, theta, phi) >=
                    options.electron_energy_floor_erg_per_cm3;
      }
    }
  }

  for (std::size_t g = 0u; g < result.group_count; ++g) {
    result.boundary_leak_total_by_group[g] =
        result.source_gain_radiation_total_by_group[g] -
        result.delta_radiation_total_by_group[g];
    result.delta_radiation_total_all_groups += result.delta_radiation_total_by_group[g];
    result.source_gain_radiation_total_all_groups +=
        result.source_gain_radiation_total_by_group[g];
    result.boundary_leak_total_all_groups += result.boundary_leak_total_by_group[g];
  }
  result.delta_electron_total = -result.source_gain_radiation_total_all_groups;
  result.radiation_electron_exchange_residual =
      std::abs(result.delta_electron_total +
               result.source_gain_radiation_total_all_groups);
  result.global_radiation_plus_electron_residual =
      std::abs(result.delta_radiation_total_all_groups +
               result.delta_electron_total +
               result.boundary_leak_total_all_groups);

  if (!valid) {
    return Fail(result, "multigroup radiation matter coupling produced invalid staged state");
  }

  staged_state.radiation_groups = staged_groups;
  const auto recovery = dec3d::state::RecoverThermodynamicState(
      staged_state,
      dec3d::state::ThermodynamicRecoveryOptions{
          options.electron_energy_floor_erg_per_cm3,
          options.ion_energy_floor_erg_per_cm3});
  result.thermodynamic_recovery_report =
      recovery.success ? recovery.recovery_diagnostics : recovery.failure_diagnostics;
  if (!recovery.success ||
      !dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovery)) {
    result.failure_diagnostics = recovery.failure_diagnostics;
    return Fail(result, "staged multigroup radiation matter thermodynamic recovery failed");
  }

  const double exchange_tolerance =
      1.0e-12 * BudgetScale(result.source_gain_radiation_total_all_groups);
  const double total_tolerance =
      1.0e-12 * BudgetScale(result.delta_radiation_total_all_groups +
                            result.delta_electron_total);
  if (result.radiation_electron_exchange_residual > exchange_tolerance) {
    return Fail(result, "multigroup radiation electron exchange residual exceeds tolerance");
  }
  if (result.global_radiation_plus_electron_residual > total_tolerance) {
    return Fail(result, "multigroup radiation plus electron residual exceeds tolerance");
  }
  if (!result.marshak_enabled &&
      std::abs(result.boundary_leak_total_all_groups) > total_tolerance) {
    return Fail(result, "zero-flux multigroup radiation reported boundary leak");
  }

  state.radiation_groups = staged_state.radiation_groups;
  state.e_electron = staged_state.e_electron;
  state.e_fluid_total = staged_state.e_fluid_total;
  result.updated_fields = MultigroupRadiationMatterWriteMask();
  result.all_groups_updated = true;
  result.first_failing_group_index = kNoFailingGroup;
  state.ApplyAuthoritativeWrite(result.updated_fields);
  result.success = true;
  BuildReport(result, true);
  return result;
}

bool ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(
    const MultigroupGrayRadiationMatterCouplingResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }
  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.multigroup_gray_matter_coupling") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "implementation_id=p3.radiation.multigroup_gray_matter_coupling_v1") &&
         Contains(line, "radiation_electron_coupling_enabled=true") &&
         Contains(line, "radiation_group_mode=explicit_frequency_groups") &&
         Contains(line, "group_count=") &&
         Contains(line, "frequency_edges_used=true") &&
         Contains(line, "frequency_edges_unit=Hz") &&
         Contains(line, "group_layout_report_present=true") &&
         Contains(line, "coefficient_time_level=old_time_lagged") &&
         Contains(line, "per_group_coefficients_source=") &&
         Contains(line, "per_group_solve_count=") &&
         Contains(line, "all_groups_updated=") &&
         Contains(line, "per_group_matrix_report_present=true") &&
         Contains(line, "per_group_solve_report_present=true") &&
         Contains(line, "per_group_backend_executed=serial_dense_reference") &&
         Contains(line, "backend_requested=serial_dense_reference") &&
         Contains(line, "backend_executed=serial_dense_reference") &&
         Contains(line, "boundary_model=") &&
         Contains(line, "marshak_enabled=") &&
         Contains(line, "electron_writeback=true") &&
         Contains(line, "updated_fields=") &&
         Contains(line, "opacity_provider=") &&
         Contains(line, "Bg_source=") &&
         Contains(line, "Bg_feedback=disabled") &&
         Contains(line, "radiation_flux_limiter=disabled") &&
         Contains(line, "advection_enabled=false") &&
         Contains(line, "radiation_pressure_work_enabled=false") &&
         Contains(line, "distributed_multigroup_enabled=false") &&
         Contains(line, "group_parallelism_enabled=false") &&
         Contains(line, "parity_claim_allowed=false") &&
         Contains(line, "first_failing_group_index=") &&
         Contains(line, "delta_radiation_total_all_groups=") &&
         Contains(line, "delta_electron_total=") &&
         Contains(line, "source_gain_radiation_total_all_groups=") &&
         Contains(line, "boundary_leak_total_all_groups=") &&
         Contains(line, "radiation_electron_exchange_residual=") &&
         Contains(line, "global_radiation_plus_electron_residual=") &&
         Contains(line, "canonical_state_mutated=");
}

}  // namespace dec3d::radiation
