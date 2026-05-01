#include "runtime/runtime_scaffold.hpp"

#include "mesh/ale/radial_ale.hpp"

#include <sstream>

namespace dec3d::runtime {

[[maybe_unused]] static constexpr const char* kRuntimeRole = "runtime.substrate.scaffold";

namespace {

[[nodiscard]] std::uint32_t StageBit(dec3d::core::StageId stage_id) noexcept {
  return 1u << static_cast<std::uint32_t>(stage_id);
}

[[nodiscard]] const char* StageName(dec3d::core::StageId stage_id) noexcept {
  using dec3d::core::StageId;

  switch (stage_id) {
    case StageId::hydro:
      return "hydro";
    case StageId::thermal:
      return "thermal";
    case StageId::equilibration:
      return "equilibration";
    case StageId::radiation:
      return "radiation";
    case StageId::alpha:
      return "alpha";
  }

  return "unknown";
}

[[nodiscard]] bool IsKnownPhase(dec3d::core::PhaseId phase_id) noexcept {
  using dec3d::core::PhaseId;

  switch (phase_id) {
    case PhaseId::p0:
    case PhaseId::p1:
    case PhaseId::p2:
    case PhaseId::p3:
    case PhaseId::p4:
    case PhaseId::p5:
      return true;
  }

  return false;
}

[[nodiscard]] std::uint32_t RequiredStageMaskForPhase(dec3d::core::PhaseId phase_id) noexcept {
  using dec3d::core::PhaseId;
  using dec3d::core::StageId;

  switch (phase_id) {
    case PhaseId::p0:
      return 0u;
    case PhaseId::p1:
      return StageBit(StageId::hydro);
    case PhaseId::p2:
      return StageBit(StageId::hydro) |
             StageBit(StageId::thermal) |
             StageBit(StageId::equilibration);
    case PhaseId::p3:
      return StageBit(StageId::hydro) |
             StageBit(StageId::thermal) |
             StageBit(StageId::equilibration) |
             StageBit(StageId::radiation);
    case PhaseId::p4:
    case PhaseId::p5:
      return StageBit(StageId::hydro) |
             StageBit(StageId::thermal) |
             StageBit(StageId::equilibration) |
             StageBit(StageId::radiation) |
             StageBit(StageId::alpha);
  }

  return 0u;
}

[[nodiscard]] const char* PhaseName(dec3d::core::PhaseId phase_id) noexcept {
  using dec3d::core::PhaseId;

  switch (phase_id) {
    case PhaseId::p0:
      return "p0";
    case PhaseId::p1:
      return "p1";
    case PhaseId::p2:
      return "p2";
    case PhaseId::p3:
      return "p3";
    case PhaseId::p4:
      return "p4";
    case PhaseId::p5:
      return "p5";
  }

  return "unknown";
}

[[nodiscard]] std::size_t CountBits(std::uint32_t mask) noexcept {
  std::size_t count = 0;
  while (mask != 0u) {
    count += mask & 1u;
    mask >>= 1u;
  }
  return count;
}

[[nodiscard]] StageRegistryEntry MakeStageRegistryEntry(
    dec3d::core::StageId stage_id) noexcept {
  StageRegistryEntry entry;
  entry.stage_id = stage_id;
  entry.stage_name = StageName(stage_id);
  entry.binding_id = std::string{"p0.runtime.stage_binding."} + entry.stage_name;
  return entry;
}

[[nodiscard]] std::uint32_t BuildIngestedStageMask(
    const RuntimeDiagnosticsSink& diagnostics_sink) noexcept {
  std::uint32_t mask = 0u;
  for (const auto& ingested_result : diagnostics_sink.ingested_stage_results) {
    mask |= StageBit(ingested_result.stage_id);
  }

  return mask;
}

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }

  return false;
}

}  // namespace

bool StageRegistryEntry::is_complete() const noexcept {
  return !stage_name.empty() &&
         !binding_id.empty();
}

bool IngestedStageResult::is_complete() const noexcept {
  return stage_result.is_semantically_complete();
}

bool RuntimeDiagnosticsSink::is_complete() const noexcept {
  return initialized &&
         !sink_id.empty();
}

bool RuntimeMeshCommitState::is_complete() const noexcept {
  return !report_line.empty();
}

bool RuntimeExecutionEvidence::is_present() const noexcept {
  return stage_registry_initialized &&
         required_stage_validation_performed &&
         stage_machine_initialization_attempted &&
         diagnostics_sink_initialized &&
         !implementation_id.empty();
}

bool StageResultIngestionResult::is_complete() const noexcept {
  return attempt_performed &&
         !report_line.empty();
}

bool StageMachineInitializationResult::is_complete() const noexcept {
  return attempt_performed &&
         execution_evidence.is_present() &&
         !report_line.empty();
}

bool StageReportAggregation::is_complete() const noexcept {
  return aggregation_performed &&
         diagnostics_complete &&
         execution_evidence_summary.is_present() &&
         !report_line.empty();
}

bool RuntimeSubstrateReport::is_complete() const noexcept {
  return initialized &&
         phase_id_valid &&
         validation_performed &&
         stage_machine_initialization_attempted &&
         stage_report_aggregation_performed &&
         diagnostics_complete &&
         execution_evidence.is_present() &&
         !report_line.empty();
}

bool RuntimeAcceptanceSnapshot::is_complete() const noexcept {
  return built &&
         !report_line.empty();
}

bool MeshCommitResult::is_complete() const noexcept {
  return attempted &&
         !report_line.empty();
}

bool MeshSubstrateRuntimeReport::is_complete() const noexcept {
  return mesh_geometry_built &&
         ownership_built &&
         ghost_topology_initialized &&
         !implementation_id.empty();
}

MeshSubstrateRuntimeReport BuildMeshSubstrateRuntimeReport(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& ownership,
    const dec3d::mesh::GhostTopology& topology,
    const dec3d::mesh::MeshBuildEvidence& evidence) noexcept {
  MeshSubstrateRuntimeReport report;
  report.mesh_geometry_built = geometry.is_valid() && evidence.is_present();
  report.ownership_built = ownership.is_valid();
  report.ghost_topology_initialized = topology.is_initialized();

  if (evidence.is_present()) {
    report.implementation_id = evidence.implementation_id;
  }

  return report;
}

RuntimeScaffold CreateRuntimeScaffold(
    dec3d::core::PhaseId phase_id,
    std::string contract_version) noexcept {
  RuntimeScaffold runtime;
  runtime.initialized = true;
  runtime.phase_id = phase_id;
  runtime.contract_version = std::move(contract_version);
  runtime.diagnostics_sink.initialized = true;
  runtime.diagnostics_sink.sink_id = "p0.runtime.diagnostics_sink";
  runtime.mesh_commit_state.report_line =
      "mesh_proposal_present=false; mesh_commit_required=false; mesh_commit_performed=false";
  return runtime;
}

bool RegisterStage(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept {
  if (!runtime.initialized) {
    return false;
  }

  if (IsStageRegistered(runtime, stage_id)) {
    return false;
  }

  runtime.registered_stage_mask |= StageBit(stage_id);
  runtime.stage_registry_entries.push_back(MakeStageRegistryEntry(stage_id));
  return true;
}

bool IsStageRegistered(
    const RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept {
  return (runtime.registered_stage_mask & StageBit(stage_id)) != 0u;
}

bool IngestStageResult(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept {
  return TryIngestStageResult(runtime, stage_id, stage_result).accepted;
}

StageResultIngestionResult TryIngestStageResult(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept {
  StageResultIngestionResult result;
  result.attempt_performed = true;
  result.registered_stage_mask = runtime.registered_stage_mask;
  result.ingested_stage_mask = BuildIngestedStageMask(runtime.diagnostics_sink);
  result.diagnostics_sink_initialized = runtime.diagnostics_sink.is_complete();

  if (!runtime.initialized) {
    result.failure_reason = "runtime scaffold is not initialized";
  } else if (!result.diagnostics_sink_initialized) {
    result.failure_reason = "runtime diagnostics sink is not initialized";
  } else if (!IsKnownPhase(runtime.phase_id)) {
    result.failure_reason = "invalid phase id";
  } else {
    result.required_stage_mask = RequiredStageMaskForPhase(runtime.phase_id);
    result.stage_registered = IsStageRegistered(runtime, stage_id);
    result.stage_required_for_phase =
        (result.required_stage_mask & StageBit(stage_id)) != 0u;
    result.stage_result_complete = stage_result.is_semantically_complete();
    result.duplicate_stage_result = HasIngestedStageResult(runtime, stage_id);

    if (!result.stage_registered) {
      result.failure_reason = "stage is not registered";
    } else if (!result.stage_required_for_phase) {
      result.failure_reason = "stage is forbidden for the active phase";
    } else if (!result.stage_result_complete) {
      result.failure_reason = "stage result is semantically incomplete";
    } else if (result.duplicate_stage_result) {
      result.failure_reason = "stage result has already been ingested";
    } else {
      runtime.diagnostics_sink.ingested_stage_results.push_back({stage_id, stage_result});
      result.accepted = true;
      result.ingested_stage_mask = BuildIngestedStageMask(runtime.diagnostics_sink);
      if (stage_result.mesh_update_proposal.has_value() &&
          stage_result.mesh_update_proposal->requested) {
        runtime.mesh_commit_state.proposal_present = true;
        runtime.mesh_commit_state.commit_required = true;
        runtime.mesh_commit_state.commit_performed = false;
        runtime.mesh_commit_state.proposal_summary = stage_result.mesh_update_proposal->summary;
        runtime.mesh_commit_state.failure_reason = "mesh update proposal requires runtime commit";
        runtime.mesh_commit_state.report_line =
            "mesh_proposal_present=true; mesh_commit_required=true; mesh_commit_performed=false";
      }
    }
  }

  std::ostringstream report;
  report << "registered_stage_mask=" << result.registered_stage_mask
         << "; required_stage_mask=" << result.required_stage_mask
         << "; ingested_stage_mask=" << result.ingested_stage_mask
         << "; stage_registered=" << (result.stage_registered ? "true" : "false")
         << "; stage_required_for_phase=" << (result.stage_required_for_phase ? "true" : "false")
         << "; stage_result_complete=" << (result.stage_result_complete ? "true" : "false")
         << "; duplicate_stage_result=" << (result.duplicate_stage_result ? "true" : "false")
         << "; accepted=" << (result.accepted ? "true" : "false");
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

bool HasIngestedStageResult(
    const RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept {
  for (const auto& ingested_result : runtime.diagnostics_sink.ingested_stage_results) {
    if (ingested_result.stage_id == stage_id) {
      return true;
    }
  }

  return false;
}

StageRegistryValidationResult ValidateStageRegistry(
    const RuntimeScaffold& runtime) noexcept {
  StageRegistryValidationResult result;
  if (!runtime.initialized) {
    result.registered_stage_mask = runtime.registered_stage_mask;
    result.failure_reason = "runtime scaffold is not initialized";
    result.report_line = "runtime_uninitialized";
    return result;
  }

  result.registered_stage_mask = runtime.registered_stage_mask;
  if (!IsKnownPhase(runtime.phase_id)) {
    result.forbidden_registered_stage_mask = runtime.registered_stage_mask;
    result.failure_reason = "invalid phase id";

    std::ostringstream report;
    report << "phase=" << PhaseName(runtime.phase_id)
           << "; registered_mask=" << runtime.registered_stage_mask
           << "; missing_required_mask=" << result.missing_required_stage_mask
           << "; forbidden_registered_mask=" << result.forbidden_registered_stage_mask
           << "; invalid_phase=true";
    result.report_line = report.str();
    return result;
  }

  const std::uint32_t required_mask = RequiredStageMaskForPhase(runtime.phase_id);
  result.required_stage_mask = required_mask;
  result.missing_required_stage_mask = required_mask & ~runtime.registered_stage_mask;
  result.forbidden_registered_stage_mask = runtime.registered_stage_mask & ~required_mask;
  result.success =
      result.missing_required_stage_mask == 0u &&
      result.forbidden_registered_stage_mask == 0u;

  std::ostringstream report;
  report << "phase=" << PhaseName(runtime.phase_id)
         << "; registered_mask=" << runtime.registered_stage_mask
         << "; missing_required_mask=" << result.missing_required_stage_mask
         << "; forbidden_registered_mask=" << result.forbidden_registered_stage_mask;
  result.report_line = report.str();

  if (!result.success) {
    if (result.missing_required_stage_mask != 0u) {
      result.failure_reason = "required stages are missing";
    } else {
      result.failure_reason = "forbidden stages are registered";
    }
  }

  return result;
}

StageMachineInitializationResult InitializeStageMachine(
    const RuntimeScaffold& runtime) noexcept {
  StageMachineInitializationResult result;
  if (!runtime.initialized) {
    result.failure_reason = "runtime scaffold is not initialized";
    result.report_line = "stage_machine_uninitialized";
    return result;
  }

  result.attempt_performed = true;
  result.required_stage_mask = 0u;
  result.registered_stage_mask = runtime.registered_stage_mask;
  result.registry_entries = runtime.stage_registry_entries;
  result.registry_entry_count = runtime.stage_registry_entries.size();
  result.zero_physics_mode = runtime.phase_id == dec3d::core::PhaseId::p0 &&
                             runtime.registered_stage_mask == 0u;
  result.execution_evidence.stage_registry_initialized = runtime.initialized;
  result.execution_evidence.required_stage_validation_performed = true;
  result.execution_evidence.stage_machine_initialization_attempted = true;
  result.execution_evidence.diagnostics_sink_initialized = runtime.diagnostics_sink.is_complete();
  result.execution_evidence.implementation_id = "p0.runtime.substrate";

  const auto validation = ValidateStageRegistry(runtime);
  result.required_stage_mask = validation.required_stage_mask;
  result.registered_stage_mask = validation.registered_stage_mask;
  result.missing_required_stage_mask = validation.missing_required_stage_mask;
  result.forbidden_registered_stage_mask = validation.forbidden_registered_stage_mask;

  if (!validation.success) {
    result.failure_reason = validation.failure_reason;
  } else if (result.registry_entry_count != CountBits(result.registered_stage_mask)) {
    result.failure_reason = "stage registry entry count mismatch";
  } else {
    bool entries_are_complete = true;
    for (const auto& entry : result.registry_entries) {
      if (!entry.is_complete()) {
        entries_are_complete = false;
        break;
      }
    }

    if (!entries_are_complete) {
      result.failure_reason = "stage registry entry is incomplete";
    } else {
      result.success = true;
      result.initialized = true;
    }
  }

  result.execution_evidence.stage_machine_initialized = result.initialized;

  std::ostringstream report;
  report << "phase=" << PhaseName(runtime.phase_id)
         << "; zero_physics_mode=" << (result.zero_physics_mode ? "true" : "false")
         << "; registry_entry_count=" << result.registry_entry_count
         << "; required_stage_mask=" << result.required_stage_mask
         << "; registered_stage_mask=" << result.registered_stage_mask
         << "; missing_required_mask=" << result.missing_required_stage_mask
         << "; forbidden_registered_mask=" << result.forbidden_registered_stage_mask
         << "; stage_machine_initialized=" << (result.initialized ? "true" : "false");
  result.report_line = report.str();

  return result;
}

StageReportAggregation BuildStageReportAggregation(
    const RuntimeScaffold& runtime) noexcept {
  StageReportAggregation aggregation;
  if (!runtime.initialized) {
    aggregation.failure_reason = "runtime scaffold is not initialized";
    aggregation.report_line = "stage_report_aggregation_uninitialized";
    return aggregation;
  }

  aggregation.aggregation_performed = true;
  aggregation.registered_stage_count = CountBits(runtime.registered_stage_mask);
  aggregation.registered_stage_digest = runtime.registered_stage_mask;
  aggregation.diagnostics_sink_initialized = runtime.diagnostics_sink.is_complete();
  aggregation.mesh_proposal_present = runtime.mesh_commit_state.proposal_present;
  aggregation.mesh_commit_required = runtime.mesh_commit_state.commit_required;
  aggregation.mesh_commit_performed = runtime.mesh_commit_state.commit_performed;
  aggregation.ingested_stage_result_count = runtime.diagnostics_sink.ingested_stage_results.size();
  aggregation.ingested_stage_mask = BuildIngestedStageMask(runtime.diagnostics_sink);

  const auto validation = ValidateStageRegistry(runtime);
  const auto stage_machine = InitializeStageMachine(runtime);

  aggregation.validation_succeeded = validation.success;
  aggregation.validated_required_stage_digest = validation.required_stage_mask;
  aggregation.stage_machine_initialized = stage_machine.initialized;
  aggregation.execution_evidence_summary = stage_machine.execution_evidence;
  aggregation.validation_failure_reason = validation.failure_reason;
  aggregation.stage_machine_failure_reason = stage_machine.failure_reason;
  aggregation.execution_evidence_summary.diagnostics_sink_initialized =
      runtime.diagnostics_sink.is_complete();
  aggregation.execution_evidence_summary.ingested_stage_result_count =
      aggregation.ingested_stage_result_count;

  for (const auto& ingested_result : runtime.diagnostics_sink.ingested_stage_results) {
    aggregation.authoritative_write_set_digest |= ingested_result.stage_result.updated_fields;
    if (ingested_result.stage_result.execution_evidence.has_value() &&
        ingested_result.stage_result.execution_evidence->is_present()) {
      ++aggregation.ingested_execution_evidence_count;
    }
  }

  aggregation.execution_evidence_summary.ingested_execution_evidence_count =
      aggregation.ingested_execution_evidence_count;
  aggregation.stage_report_count = aggregation.ingested_stage_result_count;
  aggregation.registered_but_not_ingested_mask =
      runtime.registered_stage_mask & ~aggregation.ingested_stage_mask;
  aggregation.missing_stage_report_mask =
      runtime.registered_stage_mask & ~aggregation.ingested_stage_mask;

  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.validation",
      validation.success ? "required-stage validation succeeded"
                         : "required-stage validation failed: " + validation.failure_reason);
  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.stage_machine",
      stage_machine.initialized ? "stage machine initialized"
                                : "stage machine initialization failed: " + stage_machine.failure_reason);
  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.execution_evidence",
      aggregation.execution_evidence_summary.is_present() ? "runtime execution evidence present"
                                                          : "runtime execution evidence missing");
  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.diagnostics_sink",
      aggregation.diagnostics_sink_initialized ? "runtime diagnostics sink initialized"
                                               : "runtime diagnostics sink missing");
  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.ingestion",
      "ingested_stage_results=" + std::to_string(aggregation.ingested_stage_result_count) +
          "; ingested_execution_evidence=" + std::to_string(aggregation.ingested_execution_evidence_count) +
          "; registered_stage_digest=" + std::to_string(aggregation.registered_stage_digest) +
          "; validated_required_stage_digest=" + std::to_string(aggregation.validated_required_stage_digest) +
          "; ingested_stage_mask=" + std::to_string(aggregation.ingested_stage_mask) +
          "; registered_but_not_ingested_mask=" + std::to_string(aggregation.registered_but_not_ingested_mask) +
          "; missing_stage_report_mask=" + std::to_string(aggregation.missing_stage_report_mask) +
          "; authoritative_write_set_digest=" + std::to_string(aggregation.authoritative_write_set_digest));
  AppendDiagnostic(
      aggregation.diagnostics,
      "runtime.stage_report.mesh_commit",
      runtime.mesh_commit_state.commit_required
          ? (runtime.mesh_commit_state.commit_performed
                 ? "runtime committed mesh update proposal"
                 : "runtime mesh commit missing: " + runtime.mesh_commit_state.failure_reason)
          : "no mesh commit required");

  aggregation.diagnostics_complete =
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.validation") &&
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.stage_machine") &&
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.execution_evidence") &&
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.diagnostics_sink") &&
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.ingestion") &&
      HasDiagnosticCode(aggregation.diagnostics, "runtime.stage_report.mesh_commit");

  aggregation.stage_result_ingestion_complete =
      aggregation.ingested_stage_result_count == aggregation.registered_stage_count &&
      aggregation.ingested_stage_mask == runtime.registered_stage_mask &&
      aggregation.registered_but_not_ingested_mask == 0u &&
      aggregation.missing_stage_report_mask == 0u;
  aggregation.execution_evidence_complete =
      aggregation.ingested_execution_evidence_count == aggregation.ingested_stage_result_count;
  aggregation.mesh_commit_complete =
      !aggregation.mesh_commit_required || aggregation.mesh_commit_performed;

  aggregation.success =
      aggregation.validation_succeeded &&
      aggregation.stage_machine_initialized &&
      aggregation.diagnostics_sink_initialized &&
      aggregation.stage_result_ingestion_complete &&
      aggregation.execution_evidence_complete &&
      aggregation.mesh_commit_complete &&
      aggregation.diagnostics_complete &&
      aggregation.execution_evidence_summary.is_present();

  if (!aggregation.validation_succeeded) {
    aggregation.failure_reason = aggregation.validation_failure_reason;
  } else if (!aggregation.stage_machine_initialized) {
    aggregation.failure_reason = aggregation.stage_machine_failure_reason;
  } else if (!aggregation.diagnostics_sink_initialized) {
    aggregation.failure_reason = "runtime diagnostics sink is incomplete";
  } else if (!aggregation.stage_result_ingestion_complete) {
    aggregation.ingestion_failure_reason = "stage-result ingestion is incomplete";
    aggregation.failure_reason = aggregation.ingestion_failure_reason;
  } else if (!aggregation.execution_evidence_complete) {
    aggregation.ingestion_failure_reason = "stage execution evidence aggregation is incomplete";
    aggregation.failure_reason = aggregation.ingestion_failure_reason;
  } else if (!aggregation.mesh_commit_complete) {
    aggregation.mesh_commit_failure_reason = runtime.mesh_commit_state.failure_reason.empty()
                                                ? "mesh update proposal requires runtime commit"
                                                : runtime.mesh_commit_state.failure_reason;
    aggregation.failure_reason = aggregation.mesh_commit_failure_reason;
  } else if (!aggregation.diagnostics_complete) {
    aggregation.failure_reason = "runtime diagnostics are incomplete";
  } else if (!aggregation.execution_evidence_summary.is_present()) {
    aggregation.failure_reason = "runtime execution evidence is incomplete";
  }

  std::ostringstream report;
  report << "registered_stage_count=" << aggregation.registered_stage_count
         << "; ingested_stage_result_count=" << aggregation.ingested_stage_result_count
         << "; ingested_execution_evidence_count=" << aggregation.ingested_execution_evidence_count
         << "; stage_report_count=" << aggregation.stage_report_count
         << "; registered_stage_digest=" << aggregation.registered_stage_digest
         << "; validated_required_stage_digest=" << aggregation.validated_required_stage_digest
         << "; registered_but_not_ingested_mask=" << aggregation.registered_but_not_ingested_mask
         << "; ingested_stage_mask=" << aggregation.ingested_stage_mask
         << "; missing_stage_report_mask=" << aggregation.missing_stage_report_mask
         << "; authoritative_write_set_digest=" << aggregation.authoritative_write_set_digest
         << "; validation_succeeded=" << (aggregation.validation_succeeded ? "true" : "false")
         << "; stage_machine_initialized=" << (aggregation.stage_machine_initialized ? "true" : "false")
         << "; diagnostics_sink_initialized=" << (aggregation.diagnostics_sink_initialized ? "true" : "false")
         << "; stage_result_ingestion_complete=" << (aggregation.stage_result_ingestion_complete ? "true" : "false")
         << "; execution_evidence_complete=" << (aggregation.execution_evidence_complete ? "true" : "false")
         << "; mesh_proposal_present=" << (aggregation.mesh_proposal_present ? "true" : "false")
         << "; mesh_commit_required=" << (aggregation.mesh_commit_required ? "true" : "false")
         << "; mesh_commit_performed=" << (aggregation.mesh_commit_performed ? "true" : "false")
         << "; diagnostics_complete=" << (aggregation.diagnostics_complete ? "true" : "false");
  if (!aggregation.failure_reason.empty()) {
    report << "; failure_reason=" << aggregation.failure_reason;
  }
  aggregation.report_line = report.str();

  return aggregation;
}

RuntimeSubstrateReport BuildRuntimeSubstrateReport(
    const RuntimeScaffold& runtime) noexcept {
  RuntimeSubstrateReport report;
  report.initialized = runtime.initialized;
  report.phase_id = runtime.phase_id;
  report.phase_id_valid = IsKnownPhase(runtime.phase_id);
  report.registered_stage_count = CountBits(runtime.registered_stage_mask);
  report.zero_physics_mode = runtime.phase_id == dec3d::core::PhaseId::p0 &&
                             runtime.registered_stage_mask == 0u;

  const auto validation = ValidateStageRegistry(runtime);
  const auto stage_machine = InitializeStageMachine(runtime);
  const auto aggregation = BuildStageReportAggregation(runtime);
  report.validation_performed = runtime.initialized;
  report.required_stage_contract_satisfied = validation.success;
  report.required_stage_mask = validation.required_stage_mask;
  report.registered_stage_mask = validation.registered_stage_mask;
  report.missing_required_stage_mask = validation.missing_required_stage_mask;
  report.forbidden_registered_stage_mask = validation.forbidden_registered_stage_mask;
  report.registered_stage_digest = aggregation.registered_stage_digest;
  report.validated_required_stage_digest = aggregation.validated_required_stage_digest;
  report.stage_machine_initialization_attempted = stage_machine.attempt_performed;
  report.stage_machine_initialized = stage_machine.initialized;
  report.stage_registry_entry_count = stage_machine.registry_entry_count;
  report.stage_report_aggregation_performed = aggregation.aggregation_performed;
  report.stage_report_aggregation_succeeded = aggregation.success;
  report.diagnostics_sink_initialized = aggregation.diagnostics_sink_initialized;
  report.stage_result_ingestion_complete = aggregation.stage_result_ingestion_complete;
  report.execution_evidence_complete = aggregation.execution_evidence_complete;
  report.stage_report_count = aggregation.stage_report_count;
  report.ingested_stage_result_count = aggregation.ingested_stage_result_count;
  report.ingested_execution_evidence_count = aggregation.ingested_execution_evidence_count;
  report.diagnostics_complete = aggregation.diagnostics_complete;
  report.diagnostics_entry_count = aggregation.diagnostics.entries.size();
  report.registered_but_not_ingested_mask = aggregation.registered_but_not_ingested_mask;
  report.ingested_stage_mask = aggregation.ingested_stage_mask;
  report.missing_stage_report_mask = aggregation.missing_stage_report_mask;
  report.authoritative_write_set_digest = aggregation.authoritative_write_set_digest;
  report.mesh_proposal_present = aggregation.mesh_proposal_present;
  report.mesh_commit_required = aggregation.mesh_commit_required;
  report.mesh_commit_performed = aggregation.mesh_commit_performed;
  report.zero_physics_mode = stage_machine.zero_physics_mode;
  report.execution_evidence = aggregation.execution_evidence_summary;
  report.execution_evidence.stage_registry_initialized = runtime.initialized;
  report.execution_evidence.required_stage_validation_performed = report.validation_performed;
  report.execution_evidence.stage_machine_initialization_attempted =
      stage_machine.attempt_performed;
  report.execution_evidence.stage_machine_initialized = stage_machine.initialized;
  report.execution_evidence.diagnostics_sink_initialized = aggregation.diagnostics_sink_initialized;
  report.execution_evidence.ingested_stage_result_count = aggregation.ingested_stage_result_count;
  report.execution_evidence.ingested_execution_evidence_count = aggregation.ingested_execution_evidence_count;
  report.execution_evidence.implementation_id = "p0.runtime.substrate";
  report.validation_failure_reason = validation.failure_reason;
  report.stage_machine_failure_reason = stage_machine.failure_reason;
  report.ingestion_failure_reason = aggregation.ingestion_failure_reason;
  report.mesh_commit_failure_reason = aggregation.mesh_commit_failure_reason;
  report.stage_report_failure_reason = aggregation.failure_reason;

  std::ostringstream line;
  line << "phase=" << PhaseName(runtime.phase_id)
       << "; phase_id_valid=" << (report.phase_id_valid ? "true" : "false")
       << "; registered_stage_count=" << report.registered_stage_count
       << "; stage_registry_entry_count=" << report.stage_registry_entry_count
       << "; ingested_stage_result_count=" << report.ingested_stage_result_count
       << "; ingested_execution_evidence_count=" << report.ingested_execution_evidence_count
       << "; stage_report_count=" << report.stage_report_count
       << "; registered_stage_digest=" << report.registered_stage_digest
       << "; validated_required_stage_digest=" << report.validated_required_stage_digest
       << "; registered_but_not_ingested_mask=" << report.registered_but_not_ingested_mask
       << "; ingested_stage_mask=" << report.ingested_stage_mask
       << "; missing_stage_report_mask=" << report.missing_stage_report_mask
       << "; diagnostics_entry_count=" << report.diagnostics_entry_count
       << "; authoritative_write_set_digest=" << report.authoritative_write_set_digest
       << "; required_stage_mask=" << report.required_stage_mask
       << "; registered_stage_mask=" << report.registered_stage_mask
       << "; missing_required_mask=" << report.missing_required_stage_mask
       << "; forbidden_registered_mask=" << report.forbidden_registered_stage_mask
       << "; zero_physics_mode=" << (report.zero_physics_mode ? "true" : "false")
       << "; validation_performed=" << (report.validation_performed ? "true" : "false")
        << "; stage_machine_initialization_attempted="
       << (report.stage_machine_initialization_attempted ? "true" : "false")
       << "; stage_machine_initialized=" << (report.stage_machine_initialized ? "true" : "false")
       << "; stage_report_aggregation_performed="
       << (report.stage_report_aggregation_performed ? "true" : "false")
       << "; stage_report_aggregation_succeeded="
       << (report.stage_report_aggregation_succeeded ? "true" : "false")
       << "; diagnostics_sink_initialized=" << (report.diagnostics_sink_initialized ? "true" : "false")
       << "; stage_result_ingestion_complete="
       << (report.stage_result_ingestion_complete ? "true" : "false")
       << "; execution_evidence_complete="
       << (report.execution_evidence_complete ? "true" : "false")
       << "; mesh_proposal_present=" << (report.mesh_proposal_present ? "true" : "false")
       << "; mesh_commit_required=" << (report.mesh_commit_required ? "true" : "false")
       << "; mesh_commit_performed=" << (report.mesh_commit_performed ? "true" : "false")
       << "; diagnostics_complete=" << (report.diagnostics_complete ? "true" : "false")
       << "; required_stage_contract_satisfied="
       << (report.required_stage_contract_satisfied ? "true" : "false");
  if (!report.validation_failure_reason.empty()) {
    line << "; validation_failure_reason=" << report.validation_failure_reason;
  }
  if (!report.stage_machine_failure_reason.empty()) {
    line << "; stage_machine_failure_reason=" << report.stage_machine_failure_reason;
  }
  if (!report.ingestion_failure_reason.empty()) {
    line << "; ingestion_failure_reason=" << report.ingestion_failure_reason;
  }
  if (!report.mesh_commit_failure_reason.empty()) {
    line << "; mesh_commit_failure_reason=" << report.mesh_commit_failure_reason;
  }
  if (!report.stage_report_failure_reason.empty()) {
    line << "; stage_report_failure_reason=" << report.stage_report_failure_reason;
  }
  report.report_line = line.str();

  return report;
}

MeshCommitResult CommitMeshUpdateProposal(
    RuntimeScaffold& runtime,
    dec3d::mesh::SphericalGeometryMetadata& geometry,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept {
  MeshCommitResult result;
  result.attempted = true;
  result.proposal_present =
      stage_result.mesh_update_proposal.has_value() &&
      stage_result.mesh_update_proposal->requested;

  if (!runtime.initialized) {
    result.failure_reason = "runtime scaffold is not initialized";
  } else if (!HasIngestedStageResult(runtime, stage_id)) {
    result.failure_reason = "stage result must be ingested before mesh commit";
  } else if (!result.proposal_present) {
    result.failure_reason = "stage result does not contain a mesh update proposal";
  } else {
    const auto commit = dec3d::mesh::ApplyRadialAleMeshUpdateProposal(
        *stage_result.mesh_update_proposal,
        geometry,
        0u);
    result.commit_performed = commit.success;
    result.success = commit.success;
    if (!commit.success) {
      result.failure_reason = commit.failure_reason;
    }
  }

  runtime.mesh_commit_state.proposal_present = result.proposal_present;
  runtime.mesh_commit_state.commit_required = result.proposal_present;
  runtime.mesh_commit_state.commit_performed = result.commit_performed;
  runtime.mesh_commit_state.proposal_summary =
      result.proposal_present && stage_result.mesh_update_proposal.has_value()
          ? stage_result.mesh_update_proposal->summary
          : std::string{};
  runtime.mesh_commit_state.failure_reason = result.failure_reason;

  std::ostringstream report;
  report << "mesh_proposal_present=" << (result.proposal_present ? "true" : "false")
         << "; mesh_commit_performed=" << (result.commit_performed ? "true" : "false")
         << "; success=" << (result.success ? "true" : "false");
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();
  runtime.mesh_commit_state.report_line = result.report_line;

  return result;
}

RuntimeAcceptanceSnapshot BuildRuntimeAcceptanceSnapshot(
    const RuntimeScaffold& runtime) noexcept {
  RuntimeAcceptanceSnapshot snapshot;
  snapshot.built = true;
  snapshot.phase_id = runtime.phase_id;

  const auto report = BuildRuntimeSubstrateReport(runtime);
  snapshot.runtime_report_complete = report.is_complete();
  snapshot.execution_evidence_present = report.execution_evidence.is_present();
  snapshot.diagnostics_complete = report.diagnostics_complete;
  snapshot.registered_stage_digest = report.registered_stage_digest;
  snapshot.validated_required_stage_digest = report.validated_required_stage_digest;
  snapshot.ingested_stage_digest = report.ingested_stage_mask;
  snapshot.registered_but_not_ingested_mask = report.registered_but_not_ingested_mask;
  snapshot.stage_report_count = report.stage_report_count;
  snapshot.diagnostics_entry_count = report.diagnostics_entry_count;

  if (!report.validation_failure_reason.empty()) {
    snapshot.failure_reason = report.validation_failure_reason;
  } else if (!report.stage_machine_failure_reason.empty()) {
    snapshot.failure_reason = report.stage_machine_failure_reason;
  } else if (!report.ingestion_failure_reason.empty()) {
    snapshot.failure_reason = report.ingestion_failure_reason;
  } else if (!report.stage_report_failure_reason.empty()) {
    snapshot.failure_reason = report.stage_report_failure_reason;
  }

  std::ostringstream artifact;
  artifact << "phase=" << PhaseName(runtime.phase_id)
           << "; runtime_report_complete=" << (snapshot.runtime_report_complete ? "true" : "false")
           << "; execution_evidence_present=" << (snapshot.execution_evidence_present ? "true" : "false")
           << "; diagnostics_complete=" << (snapshot.diagnostics_complete ? "true" : "false")
           << "; registered_stage_digest=" << snapshot.registered_stage_digest
           << "; validated_required_stage_digest=" << snapshot.validated_required_stage_digest
           << "; ingested_stage_digest=" << snapshot.ingested_stage_digest
           << "; registered_but_not_ingested_mask=" << snapshot.registered_but_not_ingested_mask
           << "; stage_report_count=" << snapshot.stage_report_count
           << "; diagnostics_entry_count=" << snapshot.diagnostics_entry_count;
  if (!snapshot.failure_reason.empty()) {
    artifact << "; failure_reason=" << snapshot.failure_reason;
  }
  snapshot.report_line = artifact.str();

  return snapshot;
}

}  // namespace dec3d::runtime
