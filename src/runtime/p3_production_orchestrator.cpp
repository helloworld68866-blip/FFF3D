#include "runtime/p3_production_orchestrator.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

namespace dec3d::runtime {
namespace {

constexpr const char* kPostHGeometry = "post_H_committed_ALE_geometry";
constexpr const char* kPostEState = "post_E_committed";
constexpr const char* kPostHRadiationGroups = "post_H_committed";

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] bool StageComplete(
    const P3ProductionStageResult& stage,
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
  return "none";
}

void Fail(
    P3ProductionStepResult& result,
    const char* reason,
    const char* failure_stage,
    bool h_executed,
    bool t_executed,
    bool e_executed,
    bool r_executed) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.production.step.failure"
      << "; failure_reason=" << reason
      << "; failure_stage=" << failure_stage
      << "; stage_order=H,T,E,R"
      << "; h_stage_executed=" << (h_executed ? "true" : "false")
      << "; t_stage_executed=" << (t_executed ? "true" : "false")
      << "; e_stage_executed=" << (e_executed ? "true" : "false")
      << "; r_stage_executed=" << (r_executed ? "true" : "false")
      << "; r_stage_registered=true"
      << "; a_stage_registered=false"
      << "; committed_prior_stages=" << CommittedPriorStages(failure_stage)
      << "; failed_stage_published=false"
      << "; later_stages_executed=false"
      << "; stage_local_atomic=true"
      << "; step_atomic_across_H_T_E_R=false";
  if (failure_stage != nullptr && failure_stage[0] == 'R') {
    out << "; r_stage_published=false";
  }
  result.failure_diagnostics = out.str();
}

[[nodiscard]] bool RadiationDiagnosticsComplete(
    const P3ProductionStageResult& radiation) noexcept {
  return Contains(radiation.report_line, "radiation_stage_report_present=true") &&
         Contains(radiation.report_line, "radiation_groups_epoch=post_H_committed") &&
         Contains(radiation.report_line, "electron_thermal_state_epoch=post_E_committed") &&
         Contains(
             radiation.report_line,
             "radiation_geometry_epoch=post_H_committed_ALE_geometry") &&
         Contains(
             radiation.report_line,
             "radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry") &&
         Contains(radiation.report_line, "radiation_updated_fields=") &&
         Contains(radiation.report_line, "radiation_boundary_model=") &&
         Contains(radiation.report_line, "radiation_flux_limiter_model=") &&
         Contains(radiation.report_line, "opacity_source_thesis_exact_match=false") &&
         Contains(radiation.report_line, "parity_claim_allowed=false");
}

[[nodiscard]] dec3d::core::StageResult ToRuntimeStageResult(
    const P3ProductionStageResult& stage,
    const char* implementation_id) {
  dec3d::core::DiagnosticsPayload diagnostics;
  diagnostics.entries.push_back({
      std::string{"p3.production.stage."} + stage.stage_id,
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

[[nodiscard]] bool RegisterP3Stages(RuntimeScaffold& runtime) noexcept {
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
  return hydro_registered && thermal_registered && equilibration_registered &&
         radiation_registered &&
         !IsStageRegistered(runtime, dec3d::core::StageId::alpha);
}

}  // namespace

P3ProductionStepResult ExecuteP3ProductionStep(
    const P3ProductionOptions& options,
    const P3ProductionOperatorHooks& hooks) noexcept {
  P3ProductionStepResult result;
  result.dt_s = options.dt_s;

  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    Fail(result, "dt_s must be finite and non-negative", "preflight", false, false, false, false);
    return result;
  }
  if (options.register_a_stage) {
    Fail(result, "forbidden alpha registration", "preflight", false, false, false, false);
    return result;
  }
  if (!hooks.hydro || !hooks.thermal || !hooks.equilibration || !hooks.radiation) {
    Fail(result, "H/T/E/R hooks are required", "preflight", false, false, false, false);
    return result;
  }

  result.hydro = hooks.hydro(options.dt_s);
  if (!StageComplete(result.hydro, "H", options.dt_s) || !result.hydro.success) {
    Fail(result, "hydro stage failed", "H", true, false, false, false);
    return result;
  }

  result.thermal = hooks.thermal(options.dt_s);
  if (!StageComplete(result.thermal, "T", options.dt_s) || !result.thermal.success) {
    Fail(result, "thermal stage failed", "T", true, true, false, false);
    return result;
  }
  if (result.hydro.ale_geometry_committed) {
    if (result.hydro.geometry_epoch.empty() ||
        result.thermal.geometry_epoch != result.hydro.geometry_epoch) {
      Fail(
          result,
          "thermal stage did not consume committed ALE geometry",
          "T",
          true,
          true,
          false,
          false);
      return result;
    }
  }

  result.equilibration = hooks.equilibration(options.dt_s);
  if (!StageComplete(result.equilibration, "E", options.dt_s) ||
      !result.equilibration.success) {
    Fail(result, "equilibration stage failed", "E", true, true, true, false);
    return result;
  }

  result.radiation = hooks.radiation(options.dt_s);
  if (!StageComplete(result.radiation, "R", options.dt_s) ||
      !result.radiation.success) {
    Fail(result, "radiation stage failed", "R", true, true, true, true);
    return result;
  }
  if (result.hydro.ale_geometry_committed &&
      result.radiation.geometry_epoch != result.hydro.geometry_epoch) {
    Fail(
        result,
        "radiation stage did not consume committed ALE geometry",
        "R",
        true,
        true,
        true,
        true);
    return result;
  }
  if (result.radiation.geometry_epoch != kPostHGeometry) {
    Fail(
        result,
        "radiation stage did not report post-H geometry epoch",
        "R",
        true,
        true,
        true,
        true);
    return result;
  }
  if (result.radiation.radiation_groups_epoch != kPostHRadiationGroups) {
    Fail(
        result,
        "radiation stage did not report post-H radiation groups epoch",
        "R",
        true,
        true,
        true,
        true);
    return result;
  }
  if (result.radiation.electron_thermal_state_epoch != kPostEState) {
    Fail(
        result,
        "radiation stage did not report post-E electron thermal state epoch",
        "R",
        true,
        true,
        true,
        true);
    return result;
  }
  if (!RadiationDiagnosticsComplete(result.radiation)) {
    Fail(result, "radiation stage diagnostics are incomplete", "R", true, true, true, true);
    return result;
  }

  result.success = true;
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.production.step"
      << "; phase_id=p3"
      << "; stage_order=H,T,E,R"
      << "; same_dt_for_H_T_E_R=true"
      << "; stage_local_atomic=true"
      << "; step_atomic_across_H_T_E_R=false"
      << "; h_stage_executed=true"
      << "; t_stage_executed=true"
      << "; e_stage_executed=true"
      << "; r_stage_executed=true"
      << "; r_stage_registered=true"
      << "; a_stage_registered=false"
      << "; radiation_stage_report_present=true"
      << "; radiation_groups_epoch=post_H_committed"
      << "; electron_thermal_state_epoch=post_E_committed"
      << "; radiation_geometry_epoch=post_H_committed_ALE_geometry"
      << "; radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry"
      << "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
      << "; radiation_boundary_model=" << result.radiation.radiation_boundary_model
      << "; radiation_flux_limiter_model=" << result.radiation.radiation_flux_limiter_model
      << "; opacity_source_thesis_exact_match=false"
      << "; parity_claim_allowed=false"
      << "; ale_geometry_committed=" << (result.hydro.ale_geometry_committed ? "true" : "false")
      << "; hydro_geometry_epoch=" << result.hydro.geometry_epoch
      << "; thermal_geometry_epoch=" << result.thermal.geometry_epoch
      << "; dt_s=" << options.dt_s
      << "; fallback_used=false";
  result.report_line = out.str();
  return result;
}

P3ProductionStepResult ExecuteP3ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P3ProductionOptions& options,
    const P3ProductionOperatorHooks& hooks) noexcept {
  P3ProductionStepResult result;
  if (!runtime.initialized || runtime.phase_id != dec3d::core::PhaseId::p3) {
    Fail(result, "runtime scaffold must be initialized for P3", "preflight", false, false, false, false);
    return result;
  }
  if (options.register_a_stage) {
    Fail(result, "forbidden alpha registration", "preflight", false, false, false, false);
    return result;
  }
  if (!RegisterP3Stages(runtime)) {
    Fail(
        result,
        "failed to register P3 H/T/E/R stages",
        "runtime_registration",
        false,
        false,
        false,
        false);
    return result;
  }

  result = ExecuteP3ProductionStep(options, hooks);
  if (!result.success) {
    return result;
  }

  const auto hydro_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::hydro,
      ToRuntimeStageResult(result.hydro, "p3.production.stage.hydro"));
  const auto thermal_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::thermal,
      ToRuntimeStageResult(result.thermal, "p3.production.stage.thermal"));
  const auto equilibration_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::equilibration,
      ToRuntimeStageResult(result.equilibration, "p3.production.stage.equilibration"));
  const auto radiation_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::radiation,
      ToRuntimeStageResult(result.radiation, "p3.production.stage.radiation"));

  if (!hydro_ingested || !thermal_ingested || !equilibration_ingested ||
      !radiation_ingested) {
    Fail(result, "runtime stage-result ingestion failed", "runtime", true, true, true, true);
    return result;
  }

  result.report_line += "; runtime_scaffold_registered=true; runtime_scaffold_ingested=true";
  return result;
}

bool ValidateP3ProductionStepDiagnostics(
    const P3ProductionStepResult& result) noexcept {
  if (!result.success) {
    return !result.failure_diagnostics.empty() &&
           Contains(result.failure_diagnostics, "diagnostic_id=p3.production.step.failure") &&
           Contains(result.failure_diagnostics, "failure_stage=") &&
           Contains(result.failure_diagnostics, "step_atomic_across_H_T_E_R=false");
  }
  return Contains(result.report_line, "diagnostic_id=p3.production.step") &&
         Contains(result.report_line, "phase_id=p3") &&
         Contains(result.report_line, "stage_order=H,T,E,R") &&
         Contains(result.report_line, "same_dt_for_H_T_E_R=true") &&
         Contains(result.report_line, "stage_local_atomic=true") &&
         Contains(result.report_line, "step_atomic_across_H_T_E_R=false") &&
         Contains(result.report_line, "r_stage_registered=true") &&
         Contains(result.report_line, "a_stage_registered=false") &&
         Contains(result.report_line, "radiation_stage_report_present=true") &&
         Contains(result.report_line, "radiation_groups_epoch=post_H_committed") &&
         Contains(result.report_line, "electron_thermal_state_epoch=post_E_committed") &&
         Contains(result.report_line, "radiation_geometry_epoch=post_H_committed_ALE_geometry") &&
         Contains(
             result.report_line,
             "radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry") &&
         Contains(result.report_line, "radiation_updated_fields=") &&
         Contains(result.report_line, "radiation_boundary_model=") &&
         Contains(result.report_line, "radiation_flux_limiter_model=") &&
         Contains(result.report_line, "opacity_source_thesis_exact_match=false") &&
         Contains(result.report_line, "parity_claim_allowed=false") &&
         Contains(result.report_line, "fallback_used=false");
}

P3ProductionStageResult MakeP3ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success) {
  P3ProductionStageResult result;
  result.success = success;
  result.stage_id = stage_id;
  result.dt_s = dt_s;
  const auto runtime_stage_id = StageIdFromChar(stage_id);
  switch (runtime_stage_id) {
    case dec3d::core::StageId::hydro:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::rho);
      result.geometry_epoch = kPostHGeometry;
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
      result.state_epoch = kPostEState;
      result.radiation_groups_epoch = kPostHRadiationGroups;
      result.electron_thermal_state_epoch = kPostEState;
      result.radiation_boundary_model = "thesis_marshak_vacuum";
      result.radiation_flux_limiter_model = "harmonic";
      break;
    case dec3d::core::StageId::alpha:
      result.updated_fields = 0u;
      break;
  }

  if (success) {
    std::ostringstream out;
    out << "diagnostic_id=p3.production.stage"
        << "; stage_id=" << stage_id
        << "; runtime_stage_id=" << RuntimeStageName(runtime_stage_id)
        << "; success=true";
    result.report_line = out.str();
  } else {
    std::ostringstream out;
    out << "diagnostic_id=p3.production.stage.failure"
        << "; stage_id=" << stage_id
        << "; success=false";
    result.failure_diagnostics = out.str();
  }
  return result;
}

}  // namespace dec3d::runtime
