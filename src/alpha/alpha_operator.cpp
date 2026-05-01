#include "alpha/alpha_operator.hpp"

#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace dec3d::alpha {
namespace {

constexpr const char* kImplementationId = "p4.alpha.one_group_operator_v1";

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] std::string UpdatedFieldsName(dec3d::core::AuthoritativeFieldMask mask) {
  return mask == AlphaOperatorWriteMask()
             ? "alpha_state,e_electron,e_fluid_total"
             : "none";
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t linear) noexcept {
  return geometry.cell_volumes[linear];
}

[[nodiscard]] double BudgetScale(double value) noexcept {
  return std::max(1.0, std::abs(value));
}

void PopulateProviderDiagnostics(
    OneGroupAlphaTransportResult& result,
    const AlphaCoefficientProviderResult& provider) {
  result.provider_report = provider.report_line;
  result.reactivity_model_requested = provider.reactivity_model_requested;
  result.reactivity_model_executed = provider.reactivity_model_executed;
  result.bosch_hale_coefficients_source = provider.bosch_hale_coefficients_source;
  result.bosch_hale_temperature_source = provider.bosch_hale_temperature_source;
  result.min_dt_reactivity_cm3_s = provider.min_dt_reactivity_cm3_s;
  result.max_dt_reactivity_cm3_s = provider.max_dt_reactivity_cm3_s;
  result.thesis_reactivity_claim_allowed = provider.thesis_reactivity_claim_allowed;
}

void BuildReport(
    OneGroupAlphaTransportResult& result,
    bool canonical_state_mutated,
    const char* diagnostic_id = "p4.alpha.one_group_operator") {
  const bool metadata_written = result.updated_fields != 0u;
  std::ostringstream out;
  out << "diagnostic_id=" << diagnostic_id
      << "; phase_id=P4"
      << "; stage_id=A"
      << "; implementation_id=" << kImplementationId
      << "; alpha_transport_model=atzeni_one_group"
      << "; unknown=epsilon_alpha"
      << "; epsilon_alpha_unit=erg_per_cm3"
      << "; coefficient_time_level=old_time_lagged"
      << "; composition_model=equimolar_dt_from_p2_recovery"
      << "; equimolar_dt_assumption=true"
      << "; tau_alphae_model=thesis_spitzer_eq_5_262"
      << "; D_alpha_model=atzeni_drag_one_group"
      << "; birth_source_model=nD_nT_dt_reactivity_Ealpha0"
      << "; drag_deposition_time_level=implicit_epsilon_alpha_new"
      << "; reactivity_model_requested=" << result.reactivity_model_requested
      << "; reactivity_model_executed=" << result.reactivity_model_executed
      << "; bosch_hale_coefficients_source=" << result.bosch_hale_coefficients_source
      << "; bosch_hale_temperature_source=" << result.bosch_hale_temperature_source
      << "; dt_reactivity_internal_unit=cm3_s"
      << "; min_dt_reactivity_cm3_s=" << result.min_dt_reactivity_cm3_s
      << "; max_dt_reactivity_cm3_s=" << result.max_dt_reactivity_cm3_s
      << "; thesis_reactivity_claim_allowed="
      << BoolToken(result.thesis_reactivity_claim_allowed)
      << "; backend_requested=" << result.backend_requested
      << "; backend_executed=" << result.backend_executed
      << "; provider_report_present=" << BoolToken(!result.provider_report.empty())
      << "; assembly_report_present=" << BoolToken(!result.assembly_report.empty())
      << "; matrix_report_present=" << BoolToken(!result.assembly_report.empty())
      << "; solve_report_present=" << BoolToken(!result.solve_report.empty())
      << "; thermodynamic_recovery_report_present="
      << BoolToken(!result.thermodynamic_recovery_report.empty())
      << "; cell_count=" << result.cell_count
      << "; dt_s=" << result.dt_s
      << "; delta_alpha_total=" << result.delta_alpha_total
      << "; delta_electron_total=" << result.delta_electron_total
      << "; birth_source_total=" << result.birth_source_total
      << "; drag_deposition_total=" << result.drag_deposition_total
      << "; alpha_electron_exchange_residual="
      << result.alpha_electron_exchange_residual
      << "; alpha_equation_budget_residual="
      << result.alpha_equation_budget_residual
      << "; global_alpha_plus_electron_budget_residual="
      << result.global_alpha_plus_electron_budget_residual
      << "; global_budget_semantics=single_rank_whole_domain_volume_integral"
      << "; min_epsilon_alpha_before=" << result.min_epsilon_alpha_before
      << "; max_epsilon_alpha_before=" << result.max_epsilon_alpha_before
      << "; min_epsilon_alpha_after=" << result.min_alpha_after
      << "; max_epsilon_alpha_after=" << result.max_epsilon_alpha_after
      << "; min_alpha_after=" << result.min_alpha_after
      << "; min_e_electron_after=" << result.min_e_electron_after
      << "; min_e_ion_after=" << result.min_e_ion_after
      << "; alpha_hydro_advection_enabled=false"
      << "; distributed_backend_enabled=false"
      << "; production_A_registered=false"
      << "; updated_fields=" << UpdatedFieldsName(result.updated_fields)
      << "; metadata_written=" << BoolToken(metadata_written)
      << "; canonical_state_mutated=" << BoolToken(canonical_state_mutated)
      << "; fallback_used=false";
  result.report_line = out.str();
}

[[nodiscard]] OneGroupAlphaTransportResult Fail(
    OneGroupAlphaTransportResult result,
    const std::string& reason) {
  const std::string nested = result.failure_diagnostics;
  result.success = false;
  result.updated_fields = 0u;
  result.failure_reason = reason;
  BuildReport(result, false, "p4.alpha.one_group_operator.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  if (!nested.empty()) {
    result.failure_diagnostics += "; nested_failure_diagnostics={" + nested + "}";
  }
  return result;
}

[[nodiscard]] bool FloorsAreValid(const OneGroupAlphaTransportOptions& options) noexcept {
  return std::isfinite(options.electron_energy_floor_erg_per_cm3) &&
         options.electron_energy_floor_erg_per_cm3 >= 0.0 &&
         std::isfinite(options.ion_energy_floor_erg_per_cm3) &&
         options.ion_energy_floor_erg_per_cm3 >= 0.0 &&
         std::isfinite(options.alpha_energy_floor_erg_per_cm3) &&
         options.alpha_energy_floor_erg_per_cm3 >= 0.0;
}

}  // namespace

bool OneGroupAlphaTransportResult::is_complete() const noexcept {
  return success &&
         !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p4.alpha.one_group_operator") &&
         Contains(report_line, "phase_id=P4") &&
         Contains(report_line, "stage_id=A") &&
         Contains(report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(report_line, "updated_fields=") &&
         Contains(report_line, "alpha_equation_budget_residual=") &&
         Contains(report_line, "global_budget_semantics=single_rank_whole_domain_volume_integral");
}

dec3d::core::AuthoritativeFieldMask AlphaOperatorWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

OneGroupAlphaTransportResult ApplyOneGroupAlphaTransport(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupAlphaTransportOptions& options) noexcept {
  OneGroupAlphaTransportResult result;
  result.dt_s = options.dt_s;

  if (!geometry.is_valid()) {
    return Fail(result, "alpha transport requires valid spherical geometry");
  }
  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    return Fail(result, "alpha transport dt_s must be finite and nonnegative");
  }
  if (!FloorsAreValid(options)) {
    return Fail(result, "alpha transport floors must be finite and nonnegative");
  }
  if (!state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state)) {
    return Fail(result, "alpha transport requires authoritative alpha_state storage");
  }

  auto provider_options = options.provider_options;
  provider_options.alpha_energy_floor_erg_cm3 = options.alpha_energy_floor_erg_per_cm3;
  const auto provider = BuildAlphaCoefficientArrays(state, provider_options);
  PopulateProviderDiagnostics(result, provider);
  result.cell_count = provider.cell_count;
  if (!provider.success || !ValidateAlphaCoefficientDiagnostics(provider)) {
    result.failure_diagnostics = provider.failure_diagnostics.empty()
                                     ? provider.report_line
                                     : provider.failure_diagnostics;
    return Fail(result, "alpha coefficient provider failed before transport");
  }

  result.min_epsilon_alpha_before = std::numeric_limits<double>::infinity();
  result.max_epsilon_alpha_before = -std::numeric_limits<double>::infinity();
  result.min_alpha_after = std::numeric_limits<double>::infinity();
  result.max_epsilon_alpha_after = -std::numeric_limits<double>::infinity();
  result.min_e_electron_after = std::numeric_limits<double>::infinity();
  result.min_e_ion_after = std::numeric_limits<double>::infinity();

  if (options.dt_s == 0.0) {
    result.success = true;
    result.updated_fields = 0u;
    result.assembly_report = "assembly_skipped=true; reason=no_change_contract_path";
    result.solve_report = "solve_skipped=true; reason=no_change_contract_path";
    result.thermodynamic_recovery_report =
        "recovery_skipped=true; reason=no_change_contract_path";
    for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
      for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
          result.min_epsilon_alpha_before =
              std::min(result.min_epsilon_alpha_before,
                       state.alpha_state.storage(r, theta, phi));
          result.max_epsilon_alpha_before =
              std::max(result.max_epsilon_alpha_before,
                       state.alpha_state.storage(r, theta, phi));
          result.min_alpha_after =
              std::min(result.min_alpha_after, state.alpha_state.storage(r, theta, phi));
          result.max_epsilon_alpha_after =
              std::max(result.max_epsilon_alpha_after,
                       state.alpha_state.storage(r, theta, phi));
          result.min_e_electron_after =
              std::min(result.min_e_electron_after, state.e_electron(r, theta, phi));
          result.min_e_ion_after =
              std::min(result.min_e_ion_after,
                       provider.recovered.cells(r, theta, phi).e_ion_erg_per_cm3);
        }
      }
    }
    BuildReport(result, false);
    return result;
  }

  dec3d::transport::GenericDiffusionProblem problem;
  problem.layout = dec3d::transport::DiffusionGridLayout{
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells};
  problem.geometry = geometry;
  problem.dt_s = options.dt_s;
  problem.boundary_policy = options.boundary_policy;
  problem.scalar_old = state.alpha_state.storage;
  problem.coefficient_A = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      1.0);
  problem.coefficient_D = provider.coefficients.D_alpha_cm2_s;
  problem.coefficient_C = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      0.0);
  problem.coefficient_B = provider.coefficients.birth_source_erg_cm3_s;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        problem.coefficient_C(r, theta, phi) =
            -1.0 / provider.coefficients.tau_alphae_s(r, theta, phi);
      }
    }
  }

  const auto assembly = dec3d::transport::AssembleGenericImplicitDiffusionSystem(problem);
  result.assembly_report = assembly.report_line;
  if (!assembly.success ||
      !dec3d::transport::ValidateGenericDiffusionAssemblyDiagnostics(assembly)) {
    result.failure_diagnostics = assembly.failure_diagnostics;
    const std::string reason = assembly.failure_reason.empty()
                                   ? "alpha diffusion assembly diagnostics missing"
                                   : assembly.failure_reason;
    return Fail(result, reason);
  }

  const auto solve =
      dec3d::transport::SolveGenericDiffusionReference(assembly, options.serial_reference_options);
  result.solve_report = solve.report_line;
  result.backend_executed = solve.backend;
  if (!solve.success || !dec3d::transport::ValidateGenericDiffusionSolveDiagnostics(solve)) {
    result.failure_diagnostics = solve.failure_diagnostics;
    const std::string reason = solve.failure_reason.empty()
                                   ? "alpha diffusion solve diagnostics missing"
                                   : solve.failure_reason;
    return Fail(result, reason);
  }
  if (solve.scalar_new.size() != state.alpha_state.storage.size()) {
    return Fail(result, "alpha solve returned wrong scalar size");
  }

  auto staged_state = state;
  std::size_t offset = 0u;
  std::size_t linear = 0u;
  bool staged_valid = true;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double volume = CellVolume(geometry, linear++);
        const double alpha_old = state.alpha_state.storage(r, theta, phi);
        const double alpha_new = solve.scalar_new[offset++];
        const double tau = provider.coefficients.tau_alphae_s(r, theta, phi);
        const double birth = provider.coefficients.birth_source_erg_cm3_s(r, theta, phi);
        const double deposition = options.dt_s * alpha_new / tau;

        result.min_epsilon_alpha_before =
            std::min(result.min_epsilon_alpha_before, alpha_old);
        result.max_epsilon_alpha_before =
            std::max(result.max_epsilon_alpha_before, alpha_old);

        staged_state.alpha_state.storage(r, theta, phi) = alpha_new;
        staged_state.e_electron(r, theta, phi) =
            state.e_electron(r, theta, phi) + deposition;
        staged_state.e_fluid_total(r, theta, phi) =
            state.e_fluid_total(r, theta, phi) + deposition;

        result.delta_alpha_total += (alpha_new - alpha_old) * volume;
        result.delta_electron_total += deposition * volume;
        result.birth_source_total += options.dt_s * birth * volume;
        result.drag_deposition_total += deposition * volume;
        result.min_alpha_after = std::min(result.min_alpha_after, alpha_new);
        result.max_epsilon_alpha_after =
            std::max(result.max_epsilon_alpha_after, alpha_new);
        result.min_e_electron_after =
            std::min(result.min_e_electron_after, staged_state.e_electron(r, theta, phi));

        staged_valid = staged_valid && std::isfinite(volume) && volume > 0.0 &&
                       std::isfinite(alpha_new) &&
                       alpha_new >= options.alpha_energy_floor_erg_per_cm3 &&
                       std::isfinite(deposition) &&
                       std::isfinite(staged_state.e_electron(r, theta, phi)) &&
                       std::isfinite(staged_state.e_fluid_total(r, theta, phi)) &&
                       staged_state.e_electron(r, theta, phi) >=
                           options.electron_energy_floor_erg_per_cm3;
      }
    }
  }
  if (!staged_valid) {
    return Fail(result, "alpha transport produced nonphysical staged state");
  }

  auto recovery_options = options.provider_options.recovery_options;
  recovery_options.electron_energy_floor = options.electron_energy_floor_erg_per_cm3;
  recovery_options.ion_energy_floor = options.ion_energy_floor_erg_per_cm3;
  const auto recovery =
      dec3d::state::RecoverThermodynamicState(staged_state, recovery_options);
  result.thermodynamic_recovery_report =
      recovery.success ? recovery.recovery_diagnostics : recovery.failure_diagnostics;
  if (!recovery.success ||
      !dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovery)) {
    result.failure_diagnostics = recovery.failure_diagnostics;
    return Fail(result, "staged alpha thermodynamic recovery failed");
  }
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        result.min_e_ion_after =
            std::min(result.min_e_ion_after,
                     recovery.cells(r, theta, phi).e_ion_erg_per_cm3);
      }
    }
  }

  result.alpha_electron_exchange_residual =
      std::abs(result.delta_electron_total - result.drag_deposition_total);
  result.alpha_equation_budget_residual =
      std::abs(result.delta_alpha_total - result.birth_source_total +
               result.drag_deposition_total);
  result.global_alpha_plus_electron_budget_residual =
      std::abs(result.delta_alpha_total + result.delta_electron_total -
               result.birth_source_total);
  const double exchange_tolerance =
      1.0e-12 * BudgetScale(result.drag_deposition_total);
  const double budget_tolerance =
      1.0e-10 * BudgetScale(result.birth_source_total);
  if (result.alpha_electron_exchange_residual > exchange_tolerance) {
    return Fail(result, "alpha electron exchange residual exceeds tolerance");
  }
  if (result.alpha_equation_budget_residual > budget_tolerance) {
    return Fail(result, "alpha equation budget residual exceeds tolerance");
  }
  if (result.global_alpha_plus_electron_budget_residual > budget_tolerance) {
    return Fail(result, "alpha plus electron budget residual exceeds tolerance");
  }

  state.alpha_state = staged_state.alpha_state;
  state.e_electron = staged_state.e_electron;
  state.e_fluid_total = staged_state.e_fluid_total;
  result.updated_fields = AlphaOperatorWriteMask();
  state.ApplyAuthoritativeWrite(result.updated_fields);
  result.success = true;
  BuildReport(result, true);
  return result;
}

bool ValidateAlphaOneGroupOperatorDiagnostics(
    const OneGroupAlphaTransportResult& result) noexcept {
  return result.is_complete() &&
         Contains(result.report_line, "composition_model=equimolar_dt_from_p2_recovery") &&
         Contains(result.report_line, "tau_alphae_model=thesis_spitzer_eq_5_262") &&
         Contains(result.report_line, "D_alpha_model=atzeni_drag_one_group") &&
         Contains(result.report_line, "reactivity_model_requested=") &&
         Contains(result.report_line, "reactivity_model_executed=") &&
         Contains(result.report_line, "bosch_hale_coefficients_source=") &&
         Contains(result.report_line, "bosch_hale_temperature_source=") &&
         Contains(result.report_line, "dt_reactivity_internal_unit=cm3_s") &&
         Contains(result.report_line, "min_dt_reactivity_cm3_s=") &&
         Contains(result.report_line, "max_dt_reactivity_cm3_s=") &&
         Contains(result.report_line, "thesis_reactivity_claim_allowed=") &&
         Contains(result.report_line, "provider_report_present=true") &&
         Contains(result.report_line, "assembly_report_present=true") &&
         Contains(result.report_line, "matrix_report_present=true") &&
         Contains(result.report_line, "solve_report_present=true") &&
         Contains(result.report_line, "thermodynamic_recovery_report_present=true") &&
         Contains(result.report_line, "min_epsilon_alpha_before=") &&
         Contains(result.report_line, "max_epsilon_alpha_before=") &&
         Contains(result.report_line, "min_epsilon_alpha_after=") &&
         Contains(result.report_line, "max_epsilon_alpha_after=") &&
         Contains(result.report_line, "alpha_equation_budget_residual=") &&
         Contains(result.report_line, "drag_deposition_time_level=implicit_epsilon_alpha_new") &&
         Contains(result.report_line, "global_budget_semantics=single_rank_whole_domain_volume_integral") &&
         Contains(result.report_line, "fallback_used=false");
}

}  // namespace dec3d::alpha
