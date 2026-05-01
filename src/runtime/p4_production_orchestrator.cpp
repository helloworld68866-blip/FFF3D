#include "runtime/p4_production_orchestrator.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

namespace dec3d::runtime {
namespace {

constexpr const char* kPostHGeometry = "post_H_committed_ALE_geometry";
constexpr const char* kPostEState = "post_E_committed";
constexpr const char* kPostRState = "post_R_committed";
constexpr const char* kPostHRadiationGroups = "post_H_committed";
constexpr const char* kPostHAlphaState = "post_H_committed";

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] bool StageComplete(
    const P4ProductionStageResult& stage,
    const char* expected_id,
    double dt_s) noexcept {
  return !stage.stage_id.empty() &&
         stage.stage_id == expected_id &&
         std::isfinite(stage.dt_s) &&
         stage.dt_s == dt_s &&
         (stage.success ? !stage.report_line.empty() : !stage.failure_diagnostics.empty());
}

[[nodiscard]] const char* RuntimeStageName(dec3d::core::StageId stage_id) noexcept {
  switch (stage_id) {
    case dec3d::core::StageId::hydro:
      return "H";
    case dec3d::core::StageId::thermal:
      return "T";
    case dec3d::core::StageId::equilibration:
      return "E";
    case dec3d::core::StageId::radiation:
      return "R";
    case dec3d::core::StageId::alpha:
      return "A";
  }
  return "?";
}

[[nodiscard]] dec3d::core::StageId StageIdFromChar(const char* stage_id) noexcept {
  if (stage_id != nullptr && stage_id[0] == 'T') {
    return dec3d::core::StageId::thermal;
  }
  if (stage_id != nullptr && stage_id[0] == 'E') {
    return dec3d::core::StageId::equilibration;
  }
  if (stage_id != nullptr && stage_id[0] == 'R') {
    return dec3d::core::StageId::radiation;
  }
  if (stage_id != nullptr && stage_id[0] == 'A') {
    return dec3d::core::StageId::alpha;
  }
  return dec3d::core::StageId::hydro;
}

[[nodiscard]] const char* CommittedPriorStages(const char* failure_stage) noexcept {
  if (failure_stage == nullptr) {
    return "none";
  }
  if (failure_stage[0] == 'T') {
    return "H";
  }
  if (failure_stage[0] == 'E') {
    return "H,T";
  }
  if (failure_stage[0] == 'R') {
    return "H,T,E";
  }
  if (failure_stage[0] == 'A') {
    return "H,T,E,R";
  }
  return "none";
}

void Fail(
    P4ProductionStepResult& result,
    const char* reason,
    const char* failure_stage,
    bool h_executed,
    bool t_executed,
    bool e_executed,
    bool r_executed,
    bool a_executed,
    const std::string& nested_diagnostics = {}) {
  result.success = false;
  result.advance_to_next_timestep = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p4.production.step.failure"
      << "; phase_id=p4"
      << "; failure_reason=" << reason
      << "; failure_stage=" << failure_stage
      << "; stage_order=H,T,E,R,A"
      << "; h_stage_executed=" << (h_executed ? "true" : "false")
      << "; t_stage_executed=" << (t_executed ? "true" : "false")
      << "; e_stage_executed=" << (e_executed ? "true" : "false")
      << "; r_stage_executed=" << (r_executed ? "true" : "false")
      << "; a_stage_executed=" << (a_executed ? "true" : "false")
      << "; r_stage_registered=true"
      << "; a_stage_registered=true"
      << "; committed_prior_stages=" << CommittedPriorStages(failure_stage)
      << "; failed_stage_published=false"
      << "; later_stages_executed=false"
      << "; stage_local_atomic=true"
      << "; step_atomic_across_H_T_E_R_A=false"
      << "; advance_to_next_timestep=false";
  if (failure_stage != nullptr && failure_stage[0] == 'A') {
    out << "; a_stage_published=false";
  }
  if (failure_stage != nullptr && failure_stage[0] == 'R') {
    out << "; r_stage_published=false";
  }
  if (!nested_diagnostics.empty()) {
    out << "; nested_failure_diagnostics={" << nested_diagnostics << "}";
  }
  result.failure_diagnostics = out.str();
}

[[nodiscard]] bool HydroDiagnosticsCompleteForP4(
    const P4ProductionStageResult& hydro) noexcept {
  return Contains(hydro.report_line, "alpha_hydro_terms_report_present=true") &&
         Contains(hydro.report_line, "alpha_hydro_terms_enabled=true") &&
         Contains(
             hydro.report_line,
             "alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure");
}

[[nodiscard]] bool RadiationDiagnosticsComplete(
    const P4ProductionStageResult& radiation) noexcept {
  return Contains(radiation.report_line, "radiation_stage_report_present=true") &&
         Contains(radiation.report_line, "radiation_hydro_terms_report_present=true") &&
         Contains(radiation.report_line, "radiation_groups_epoch=post_H_committed") &&
         Contains(radiation.report_line, "electron_thermal_state_epoch=post_E_committed") &&
         Contains(
             radiation.report_line,
             "radiation_geometry_epoch=post_H_committed_ALE_geometry") &&
         Contains(
             radiation.report_line,
             "radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure") &&
         Contains(radiation.report_line, "radiation_updated_fields=") &&
         Contains(radiation.report_line, "radiation_boundary_model=") &&
         Contains(radiation.report_line, "radiation_flux_limiter_model=") &&
         Contains(radiation.report_line, "opacity_source_thesis_exact_match=false") &&
         Contains(radiation.report_line, "parity_claim_allowed=false");
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask AlphaRequiredWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask AlphaForbiddenWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::rho) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_r) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_theta) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_phi) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
}

[[nodiscard]] bool AlphaWriteSetValid(
    const P4ProductionStageResult& alpha,
    double dt_s) noexcept {
  if ((alpha.updated_fields & AlphaForbiddenWriteMask()) != 0u) {
    return false;
  }
  if (dt_s == 0.0 && alpha.updated_fields == 0u &&
      Contains(alpha.report_line, "alpha_updated_fields=none")) {
    return true;
  }
  return alpha.updated_fields == AlphaRequiredWriteMask() &&
         Contains(alpha.report_line, "alpha_updated_fields=alpha_state,e_electron,e_fluid_total");
}

[[nodiscard]] bool AlphaDiagnosticsComplete(
    const P4ProductionStageResult& alpha,
    double dt_s) noexcept {
  if (!Contains(alpha.report_line, "alpha_stage_report_present=true") ||
      !Contains(alpha.report_line, "alpha_geometry_epoch=post_H_committed_ALE_geometry") ||
      !Contains(alpha.report_line, "alpha_state_epoch=post_H_committed") ||
      !Contains(alpha.report_line, "alpha_thermal_state_epoch=post_R_committed") ||
      !Contains(alpha.report_line, "alpha_transport_model=atzeni_one_group") ||
      !Contains(alpha.report_line, "reactivity_model_executed=bosch_hale_dt") ||
      !Contains(alpha.report_line, "coefficient_time_level=old_time_lagged_to_A_stage_entry") ||
      !Contains(alpha.report_line, "fuel_depletion_enabled=false") ||
      !Contains(alpha.report_line, "separate_dt_species_authoritative=false") ||
      !Contains(alpha.report_line, "owned_slab_writeback_only=true") ||
      !Contains(alpha.report_line, "rank0_gather_solve_used=false") ||
      !Contains(alpha.report_line, "serial_dense_fallback_used=false") ||
      !Contains(alpha.report_line, "parity_claim_allowed=false")) {
    return false;
  }
  if (dt_s == 0.0 && Contains(alpha.report_line, "alpha_updated_fields=none")) {
    return Contains(alpha.report_line, "canonical_state_mutated=false") &&
           Contains(alpha.report_line, "metadata_written=false");
  }
  return Contains(alpha.report_line, "alpha_backend_executed=hypre_parcsr_gmres_boomeramg") &&
         Contains(alpha.report_line, "alpha_updated_fields=alpha_state,e_electron,e_fluid_total");
}

[[nodiscard]] dec3d::core::StageResult ToRuntimeStageResult(
    const P4ProductionStageResult& stage,
    const char* implementation_id) {
  dec3d::core::DiagnosticsPayload diagnostics;
  diagnostics.entries.push_back({
      std::string{"p4.production.stage."} + stage.stage_id,
      stage.report_line});
  dec3d::core::ExecutionEvidence evidence;
  evidence.entered_stage = true;
  evidence.implementation_id = implementation_id;
  evidence.touched_cell_count = 1u;
  return dec3d::core::StageResult::Successful(
      stage.updated_fields,
      std::move(diagnostics),
      std::move(evidence));
}

[[nodiscard]] bool RegisterP4Stages(RuntimeScaffold& runtime) noexcept {
  const bool hydro_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::hydro) ||
      RegisterStage(runtime, dec3d::core::StageId::hydro);
  const bool thermal_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::thermal) ||
      RegisterStage(runtime, dec3d::core::StageId::thermal);
  const bool equilibration_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::equilibration) ||
      RegisterStage(runtime, dec3d::core::StageId::equilibration);
  const bool radiation_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::radiation) ||
      RegisterStage(runtime, dec3d::core::StageId::radiation);
  const bool alpha_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::alpha) ||
      RegisterStage(runtime, dec3d::core::StageId::alpha);
  return hydro_registered && thermal_registered && equilibration_registered &&
         radiation_registered && alpha_registered;
}

}  // namespace

P4ProductionStepResult ExecuteP4ProductionStep(
    const P4ProductionOptions& options,
    const P4ProductionOperatorHooks& hooks) noexcept {
  P4ProductionStepResult result;
  result.dt_s = options.dt_s;

  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    Fail(result, "dt_s must be finite and non-negative", "preflight", false, false, false, false, false);
    return result;
  }
  if (!hooks.hydro || !hooks.thermal || !hooks.equilibration || !hooks.radiation ||
      !hooks.alpha) {
    Fail(result, "missing A hook or required H/T/E/R/A hook", "preflight", false, false, false, false, false);
    return result;
  }

  result.hydro = hooks.hydro(options.dt_s);
  if (!StageComplete(result.hydro, "H", options.dt_s) || !result.hydro.success) {
    Fail(result, "hydro stage failed", "H", true, false, false, false, false);
    return result;
  }
  if (!HydroDiagnosticsCompleteForP4(result.hydro)) {
    Fail(result, "hydro stage did not report alpha hydro terms", "H", true, false, false, false, false);
    return result;
  }

  result.thermal = hooks.thermal(options.dt_s);
  if (!StageComplete(result.thermal, "T", options.dt_s) || !result.thermal.success) {
    Fail(result, "thermal stage failed", "T", true, true, false, false, false);
    return result;
  }
  if (result.hydro.ale_geometry_committed) {
    if (result.hydro.geometry_epoch.empty() ||
        result.thermal.geometry_epoch != result.hydro.geometry_epoch) {
      Fail(result, "thermal stage did not consume committed ALE geometry", "T", true, true, false, false, false);
      return result;
    }
  }

  result.equilibration = hooks.equilibration(options.dt_s);
  if (!StageComplete(result.equilibration, "E", options.dt_s) ||
      !result.equilibration.success) {
    Fail(result, "equilibration stage failed", "E", true, true, true, false, false);
    return result;
  }

  result.radiation = hooks.radiation(options.dt_s);
  if (!StageComplete(result.radiation, "R", options.dt_s) ||
      !result.radiation.success) {
    Fail(result, "radiation stage failed", "R", true, true, true, true, false);
    return result;
  }
  if (result.hydro.ale_geometry_committed &&
      result.radiation.geometry_epoch != result.hydro.geometry_epoch) {
    Fail(result, "radiation stage did not consume committed ALE geometry", "R", true, true, true, true, false);
    return result;
  }
  if (result.radiation.geometry_epoch != kPostHGeometry) {
    Fail(result, "radiation stage did not report post-H geometry epoch", "R", true, true, true, true, false);
    return result;
  }
  if (result.radiation.radiation_groups_epoch != kPostHRadiationGroups) {
    Fail(result, "radiation stage did not report post-H radiation groups epoch", "R", true, true, true, true, false);
    return result;
  }
  if (result.radiation.electron_thermal_state_epoch != kPostEState) {
    Fail(result, "radiation stage did not report post-E electron thermal state epoch", "R", true, true, true, true, false);
    return result;
  }
  if (!RadiationDiagnosticsComplete(result.radiation)) {
    Fail(result, "radiation stage diagnostics are incomplete", "R", true, true, true, true, false);
    return result;
  }

  result.alpha = hooks.alpha(options.dt_s);
  if (!StageComplete(result.alpha, "A", options.dt_s) || !result.alpha.success) {
    Fail(
        result,
        "alpha stage failed",
        "A",
        true,
        true,
        true,
        true,
        true,
        result.alpha.failure_diagnostics);
    return result;
  }
  if (result.hydro.ale_geometry_committed &&
      result.alpha.alpha_geometry_epoch != result.hydro.geometry_epoch) {
    Fail(result, "alpha stage did not consume committed ALE geometry", "A", true, true, true, true, true);
    return result;
  }
  if (result.alpha.alpha_geometry_epoch != kPostHGeometry) {
    Fail(result, "alpha stage did not report post-H geometry epoch", "A", true, true, true, true, true);
    return result;
  }
  if (result.alpha.alpha_state_epoch != kPostHAlphaState) {
    Fail(result, "alpha stage did not report post-H alpha state epoch", "A", true, true, true, true, true);
    return result;
  }
  if (result.alpha.alpha_thermal_state_epoch != kPostRState) {
    Fail(result, "alpha stage did not report post-R thermal state epoch", "A", true, true, true, true, true);
    return result;
  }
  if (!AlphaWriteSetValid(result.alpha, options.dt_s)) {
    Fail(result, "alpha stage write set is invalid", "A", true, true, true, true, true);
    return result;
  }
  if (!AlphaDiagnosticsComplete(result.alpha, options.dt_s)) {
    Fail(result, "alpha stage diagnostics are incomplete", "A", true, true, true, true, true);
    return result;
  }

  result.success = true;
  result.advance_to_next_timestep = true;
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p4.production.step"
      << "; phase_id=p4"
      << "; stage_order=H,T,E,R,A"
      << "; same_dt_for_H_T_E_R_A=true"
      << "; stage_local_atomic=true"
      << "; step_atomic_across_H_T_E_R_A=false"
      << "; h_stage_executed=true"
      << "; t_stage_executed=true"
      << "; e_stage_executed=true"
      << "; r_stage_executed=true"
      << "; a_stage_executed=true"
      << "; r_stage_registered=true"
      << "; a_stage_registered=true"
      << "; radiation_stage_report_present=true"
      << "; radiation_hydro_terms_report_present=true"
      << "; radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"
      << "; alpha_stage_report_present=true"
      << "; alpha_hydro_terms_report_present=true"
      << "; alpha_hydro_terms_enabled=true"
      << "; alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure"
      << "; radiation_groups_epoch=post_H_committed"
      << "; electron_thermal_state_epoch=post_E_committed"
      << "; radiation_geometry_epoch=post_H_committed_ALE_geometry"
      << "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
      << "; alpha_state_epoch=post_H_committed"
      << "; alpha_thermal_state_epoch=post_R_committed"
      << "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
      << "; alpha_transport_model=atzeni_one_group"
      << "; reactivity_model_executed=bosch_hale_dt"
      << "; fuel_depletion_enabled=false"
      << "; separate_dt_species_authoritative=false"
      << (options.dt_s == 0.0 ? "; alpha_updated_fields=none"
                              : "; alpha_updated_fields=alpha_state,e_electron,e_fluid_total")
      << "; opacity_source_thesis_exact_match=false"
      << "; alpha_parity_claim_allowed=false"
      << "; radiation_parity_claim_allowed=false"
      << "; fallback_used=false"
      << "; advance_to_next_timestep=true"
      << "; p4_production_orchestrator_backend_dependency=hook_only"
      << "; hypre_link_required_by_runtime_orchestrator=false"
      << "; ale_geometry_committed=" << (result.hydro.ale_geometry_committed ? "true" : "false")
      << "; hydro_geometry_epoch=" << result.hydro.geometry_epoch
      << "; thermal_geometry_epoch=" << result.thermal.geometry_epoch
      << "; dt_s=" << options.dt_s;
  if (options.dt_s == 0.0) {
    out << "; canonical_state_mutated=false; metadata_written=false";
  }
  result.report_line = out.str();
  return result;
}

P4ProductionStepResult ExecuteP4ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P4ProductionOptions& options,
    const P4ProductionOperatorHooks& hooks) noexcept {
  P4ProductionStepResult result;
  if (!runtime.initialized || runtime.phase_id != dec3d::core::PhaseId::p4) {
    Fail(result, "runtime scaffold must be initialized for P4", "preflight", false, false, false, false, false);
    return result;
  }
  if (!RegisterP4Stages(runtime)) {
    Fail(result, "failed to register P4 H/T/E/R/A stages", "runtime_registration", false, false, false, false, false);
    return result;
  }

  result = ExecuteP4ProductionStep(options, hooks);
  if (!result.success) {
    return result;
  }

  const auto hydro_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::hydro,
      ToRuntimeStageResult(result.hydro, "p4.production.stage.hydro"));
  const auto thermal_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::thermal,
      ToRuntimeStageResult(result.thermal, "p4.production.stage.thermal"));
  const auto equilibration_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::equilibration,
      ToRuntimeStageResult(result.equilibration, "p4.production.stage.equilibration"));
  const auto radiation_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::radiation,
      ToRuntimeStageResult(result.radiation, "p4.production.stage.radiation"));
  const auto alpha_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::alpha,
      ToRuntimeStageResult(result.alpha, "p4.production.stage.alpha"));

  if (!hydro_ingested || !thermal_ingested || !equilibration_ingested ||
      !radiation_ingested || !alpha_ingested) {
    Fail(result, "runtime stage-result ingestion failed", "runtime", true, true, true, true, true);
    return result;
  }

  result.report_line += "; runtime_scaffold_registered=true; runtime_scaffold_ingested=true";
  return result;
}

bool ValidateP4ProductionStepDiagnostics(
    const P4ProductionStepResult& result) noexcept {
  if (!result.success) {
    return !result.failure_diagnostics.empty() &&
           Contains(result.failure_diagnostics, "diagnostic_id=p4.production.step.failure") &&
           Contains(result.failure_diagnostics, "failure_stage=") &&
           Contains(result.failure_diagnostics, "stage_order=H,T,E,R,A") &&
           Contains(result.failure_diagnostics, "step_atomic_across_H_T_E_R_A=false") &&
           Contains(result.failure_diagnostics, "advance_to_next_timestep=false");
  }
  const bool alpha_fields_ok =
      Contains(result.report_line, "alpha_updated_fields=alpha_state,e_electron,e_fluid_total") ||
      Contains(result.report_line, "alpha_updated_fields=none");
  return Contains(result.report_line, "diagnostic_id=p4.production.step") &&
         Contains(result.report_line, "phase_id=p4") &&
         Contains(result.report_line, "stage_order=H,T,E,R,A") &&
         Contains(result.report_line, "same_dt_for_H_T_E_R_A=true") &&
         Contains(result.report_line, "stage_local_atomic=true") &&
         Contains(result.report_line, "step_atomic_across_H_T_E_R_A=false") &&
         Contains(result.report_line, "h_stage_executed=true") &&
         Contains(result.report_line, "t_stage_executed=true") &&
         Contains(result.report_line, "e_stage_executed=true") &&
         Contains(result.report_line, "r_stage_executed=true") &&
         Contains(result.report_line, "a_stage_executed=true") &&
         Contains(result.report_line, "r_stage_registered=true") &&
         Contains(result.report_line, "a_stage_registered=true") &&
         Contains(result.report_line, "radiation_stage_report_present=true") &&
         Contains(result.report_line, "radiation_hydro_terms_report_present=true") &&
         Contains(result.report_line, "radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure") &&
         Contains(result.report_line, "alpha_stage_report_present=true") &&
         Contains(result.report_line, "alpha_hydro_terms_report_present=true") &&
         Contains(result.report_line, "alpha_hydro_terms_enabled=true") &&
         Contains(result.report_line, "alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure") &&
         Contains(result.report_line, "radiation_groups_epoch=post_H_committed") &&
         Contains(result.report_line, "electron_thermal_state_epoch=post_E_committed") &&
         Contains(result.report_line, "radiation_geometry_epoch=post_H_committed_ALE_geometry") &&
         Contains(result.report_line, "alpha_geometry_epoch=post_H_committed_ALE_geometry") &&
         Contains(result.report_line, "alpha_state_epoch=post_H_committed") &&
         Contains(result.report_line, "alpha_thermal_state_epoch=post_R_committed") &&
         Contains(result.report_line, "coefficient_time_level=old_time_lagged_to_A_stage_entry") &&
         Contains(result.report_line, "alpha_transport_model=atzeni_one_group") &&
         Contains(result.report_line, "reactivity_model_executed=bosch_hale_dt") &&
         Contains(result.report_line, "fuel_depletion_enabled=false") &&
         Contains(result.report_line, "separate_dt_species_authoritative=false") &&
         alpha_fields_ok &&
         Contains(result.report_line, "opacity_source_thesis_exact_match=false") &&
         Contains(result.report_line, "alpha_parity_claim_allowed=false") &&
         Contains(result.report_line, "radiation_parity_claim_allowed=false") &&
         Contains(result.report_line, "fallback_used=false") &&
         Contains(result.report_line, "advance_to_next_timestep=true") &&
         Contains(result.report_line, "p4_production_orchestrator_backend_dependency=hook_only") &&
         Contains(result.report_line, "hypre_link_required_by_runtime_orchestrator=false");
}

P4ProductionStageResult MakeP4ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success) {
  P4ProductionStageResult result;
  result.success = success;
  result.stage_id = stage_id;
  result.dt_s = dt_s;
  const auto runtime_stage_id = StageIdFromChar(stage_id);
  switch (runtime_stage_id) {
    case dec3d::core::StageId::hydro:
      result.updated_fields =
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::rho) |
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
      result.geometry_epoch = kPostHGeometry;
      result.alpha_state_epoch = kPostHAlphaState;
      break;
    case dec3d::core::StageId::thermal:
      result.updated_fields =
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
      result.geometry_epoch = kPostHGeometry;
      break;
    case dec3d::core::StageId::equilibration:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron);
      result.state_epoch = kPostEState;
      break;
    case dec3d::core::StageId::radiation:
      result.updated_fields =
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups) |
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
      result.geometry_epoch = kPostHGeometry;
      result.state_epoch = kPostRState;
      result.radiation_groups_epoch = kPostHRadiationGroups;
      result.electron_thermal_state_epoch = kPostEState;
      result.radiation_boundary_model = "thesis_marshak_vacuum";
      result.radiation_flux_limiter_model = "harmonic";
      break;
    case dec3d::core::StageId::alpha:
      result.updated_fields = AlphaRequiredWriteMask();
      result.geometry_epoch = kPostHGeometry;
      result.state_epoch = kPostRState;
      result.alpha_geometry_epoch = kPostHGeometry;
      result.alpha_state_epoch = kPostHAlphaState;
      result.alpha_thermal_state_epoch = kPostRState;
      break;
  }

  if (success) {
    std::ostringstream out;
    out << "diagnostic_id=p4.production.stage"
        << "; stage_id=" << stage_id
        << "; runtime_stage_id=" << RuntimeStageName(runtime_stage_id)
        << "; success=true";
    result.report_line = out.str();
  } else {
    std::ostringstream out;
    out << "diagnostic_id=p4.production.stage.failure"
        << "; stage_id=" << stage_id
        << "; success=false";
    result.failure_diagnostics = out.str();
  }
  return result;
}

}  // namespace dec3d::runtime
