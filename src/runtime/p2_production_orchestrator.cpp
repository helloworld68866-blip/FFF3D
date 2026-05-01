#include "runtime/p2_production_orchestrator.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <utility>

namespace dec3d::runtime {
namespace {

[[nodiscard]] bool StageComplete(
    const P2ProductionStageResult& stage,
    const char* expected_id,
    double dt_s) noexcept {
  return !stage.stage_id.empty() &&
         stage.stage_id == expected_id &&
         std::isfinite(stage.dt_s) &&
         stage.dt_s == dt_s &&
         (stage.success ? !stage.report_line.empty() : !stage.failure_diagnostics.empty());
}

void Fail(
    P2ProductionStepResult& result,
    const char* reason,
    const char* stage_id,
    bool e_stage_executed) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p2.production.step.failure"
      << "; failure_reason=" << reason
      << "; stage_id=" << stage_id
      << "; stage_order=H,T,E"
      << "; e_stage_executed=" << (e_stage_executed ? "true" : "false")
      << "; r_stage_registered=false"
      << "; a_stage_registered=false"
      << "; canonical_state_mutated=false";
  result.failure_diagnostics = out.str();
}

[[nodiscard]] dec3d::core::StageId StageIdFromChar(const char* stage_id) noexcept {
  if (stage_id != nullptr && stage_id[0] == 'T') {
    return dec3d::core::StageId::thermal;
  }
  if (stage_id != nullptr && stage_id[0] == 'E') {
    return dec3d::core::StageId::equilibration;
  }
  return dec3d::core::StageId::hydro;
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

[[nodiscard]] dec3d::core::StageResult ToRuntimeStageResult(
    const P2ProductionStageResult& stage,
    const char* implementation_id) {
  dec3d::core::DiagnosticsPayload diagnostics;
  diagnostics.entries.push_back({
      std::string{"p2.production.stage."} + stage.stage_id,
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

[[nodiscard]] bool RegisterP2Stages(RuntimeScaffold& runtime) noexcept {
  const bool hydro_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::hydro) ||
      RegisterStage(runtime, dec3d::core::StageId::hydro);
  const bool thermal_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::thermal) ||
      RegisterStage(runtime, dec3d::core::StageId::thermal);
  const bool equilibration_registered =
      IsStageRegistered(runtime, dec3d::core::StageId::equilibration) ||
      RegisterStage(runtime, dec3d::core::StageId::equilibration);
  return hydro_registered && thermal_registered && equilibration_registered &&
         !IsStageRegistered(runtime, dec3d::core::StageId::radiation) &&
         !IsStageRegistered(runtime, dec3d::core::StageId::alpha);
}

}  // namespace

P2ProductionStepResult ExecuteP2ProductionStep(
    const P2ProductionOptions& options,
    const P2ProductionOperatorHooks& hooks) noexcept {
  P2ProductionStepResult result;
  result.dt_s = options.dt_s;

  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    Fail(result, "dt_s must be finite and non-negative", "preflight", false);
    return result;
  }
  if (options.register_r_stage || options.register_a_stage) {
    Fail(result, "R/A stages are not registered in P2", "preflight", false);
    return result;
  }
  if (!hooks.hydro || !hooks.thermal || !hooks.equilibration) {
    Fail(result, "H/T/E hooks are required", "preflight", false);
    return result;
  }

  result.hydro = hooks.hydro(options.dt_s);
  if (!StageComplete(result.hydro, "H", options.dt_s) || !result.hydro.success) {
    Fail(result, "hydro stage failed", "H", false);
    return result;
  }

  result.thermal = hooks.thermal(options.dt_s);
  if (!StageComplete(result.thermal, "T", options.dt_s) || !result.thermal.success) {
    Fail(result, "thermal stage failed", "T", false);
    return result;
  }
  if (result.hydro.ale_geometry_committed) {
    if (result.hydro.geometry_epoch.empty() ||
        result.thermal.geometry_epoch != result.hydro.geometry_epoch) {
      Fail(
          result,
          "thermal stage did not consume committed ALE geometry",
          "T",
          false);
      return result;
    }
  }

  result.equilibration = hooks.equilibration(options.dt_s);
  if (!StageComplete(result.equilibration, "E", options.dt_s) ||
      !result.equilibration.success) {
    Fail(result, "equilibration stage failed", "E", true);
    return result;
  }

  result.success = true;
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.production.step"
      << "; stage_order=H,T,E"
      << "; h_stage_executed=true"
      << "; t_stage_executed=true"
      << "; e_stage_executed=true"
      << "; r_stage_registered=false"
      << "; a_stage_registered=false"
      << "; ale_geometry_committed="
      << (result.hydro.ale_geometry_committed ? "true" : "false")
      << "; hydro_geometry_epoch=" << result.hydro.geometry_epoch
      << "; thermal_geometry_epoch=" << result.thermal.geometry_epoch
      << "; dt_s=" << options.dt_s
      << "; fallback_used=false";
  result.report_line = out.str();
  return result;
}

P2ProductionStepResult ExecuteP2ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P2ProductionOptions& options,
    const P2ProductionOperatorHooks& hooks) noexcept {
  P2ProductionStepResult result;
  if (!runtime.initialized || runtime.phase_id != dec3d::core::PhaseId::p2) {
    Fail(result, "runtime scaffold must be initialized for P2", "preflight", false);
    return result;
  }
  if (!RegisterP2Stages(runtime)) {
    Fail(result, "failed to register P2 H/T/E stages", "preflight", false);
    return result;
  }

  result = ExecuteP2ProductionStep(options, hooks);
  if (!result.success) {
    return result;
  }

  const auto hydro_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::hydro,
      ToRuntimeStageResult(result.hydro, "p2.production.stage.hydro"));
  const auto thermal_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::thermal,
      ToRuntimeStageResult(result.thermal, "p2.production.stage.thermal"));
  const auto equilibration_ingested = IngestStageResult(
      runtime,
      dec3d::core::StageId::equilibration,
      ToRuntimeStageResult(result.equilibration, "p2.production.stage.equilibration"));

  if (!hydro_ingested || !thermal_ingested || !equilibration_ingested) {
    Fail(result, "runtime stage-result ingestion failed", "runtime", true);
    return result;
  }

  result.report_line += "; runtime_scaffold_registered=true; runtime_scaffold_ingested=true";
  return result;
}

bool ValidateP2ProductionStepDiagnostics(
    const P2ProductionStepResult& result) noexcept {
  if (!result.success) {
    return !result.failure_diagnostics.empty() &&
           result.failure_diagnostics.find("diagnostic_id=p2.production.step.failure") !=
               std::string::npos &&
           result.failure_diagnostics.find("stage_id=") != std::string::npos;
  }
  return result.report_line.find("diagnostic_id=p2.production.step") != std::string::npos &&
         result.report_line.find("stage_order=H,T,E") != std::string::npos &&
         result.report_line.find("r_stage_registered=false") != std::string::npos &&
         result.report_line.find("a_stage_registered=false") != std::string::npos &&
         result.report_line.find("fallback_used=false") != std::string::npos;
}

P2ProductionStageResult MakeP2ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success) {
  P2ProductionStageResult result;
  result.success = success;
  result.stage_id = stage_id;
  result.dt_s = dt_s;
  const auto runtime_stage_id = StageIdFromChar(stage_id);
  switch (runtime_stage_id) {
    case dec3d::core::StageId::hydro:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::rho);
      break;
    case dec3d::core::StageId::thermal:
      result.updated_fields =
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
          dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
      break;
    case dec3d::core::StageId::equilibration:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron);
      break;
    case dec3d::core::StageId::radiation:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
      break;
    case dec3d::core::StageId::alpha:
      result.updated_fields = dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
      break;
  }
  if (success) {
    std::ostringstream out;
    out << "diagnostic_id=p2.production.stage"
        << "; stage_id=" << stage_id
        << "; runtime_stage_id=" << RuntimeStageName(runtime_stage_id)
        << "; success=true";
    result.report_line = out.str();
  } else {
    std::ostringstream out;
    out << "diagnostic_id=p2.production.stage.failure"
        << "; stage_id=" << stage_id
        << "; success=false";
    result.failure_diagnostics = out.str();
  }
  return result;
}

}  // namespace dec3d::runtime
