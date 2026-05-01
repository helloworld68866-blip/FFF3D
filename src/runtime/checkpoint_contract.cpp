#include "runtime/checkpoint_contract.hpp"

#include <sstream>

namespace dec3d::runtime {

namespace {

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

[[nodiscard]] bool RuntimeSnapshotIsUsable(
    const RuntimeAcceptanceSnapshot& snapshot) noexcept {
  return snapshot.built &&
         snapshot.runtime_report_complete &&
         snapshot.execution_evidence_present &&
         snapshot.diagnostics_complete;
}

}  // namespace

bool CheckpointPayload::is_complete() const noexcept {
  const auto matches_layout = [this](const dec3d::core::Array3D<double>& field) {
    return field.extent_r() == layout.radial_cells &&
           field.extent_theta() == layout.theta_cells &&
           field.extent_phi() == layout.phi_cells;
  };

  if (!(serialized &&
        layout.is_valid() &&
        matches_layout(rho) &&
        matches_layout(mom_r) &&
        matches_layout(mom_theta) &&
        matches_layout(mom_phi) &&
        matches_layout(e_fluid_total) &&
        matches_layout(e_electron) &&
        matches_layout(alpha_state.storage) &&
        radiation_groups.size() == layout.radiation_group_count &&
        RuntimeSnapshotIsUsable(runtime_snapshot) &&
        diagnostics.has_entries() &&
        !report_line.empty())) {
    return false;
  }

  for (const auto& radiation_group : radiation_groups) {
    if (!matches_layout(radiation_group)) {
      return false;
    }
  }

  return true;
}

bool RestartValidationResult::is_complete() const noexcept {
  return validation_performed &&
         !report_line.empty();
}

bool CheckpointContractReport::is_complete() const noexcept {
  return built &&
         diagnostics_complete &&
         !report_line.empty();
}

CheckpointPayload BuildCheckpointPayload(
    const dec3d::state::CanonicalState& state,
    const RuntimeScaffold& runtime) {
  CheckpointPayload payload;
  payload.serialized = state.layout.is_valid();
  payload.phase_id = runtime.phase_id;
  payload.contract_version = runtime.contract_version;
  payload.layout = state.layout;
  payload.rho = state.rho;
  payload.mom_r = state.mom_r;
  payload.mom_theta = state.mom_theta;
  payload.mom_phi = state.mom_phi;
  payload.e_fluid_total = state.e_fluid_total;
  payload.e_electron = state.e_electron;
  payload.radiation_groups = state.radiation_groups;
  payload.alpha_state = state.alpha_state;
  payload.runtime_snapshot = BuildRuntimeAcceptanceSnapshot(runtime);

  const auto summary = dec3d::state::BuildCanonicalStateSummary(state);
  payload.authoritative_field_mask = summary.restart_authoritative_mask;
  payload.serialized_cached_field_mask = 0u;

  AppendDiagnostic(
      payload.diagnostics,
      "checkpoint.payload.authoritative_state",
      payload.serialized ? "authoritative state serialized"
                         : "authoritative state serialization unavailable");
  AppendDiagnostic(
      payload.diagnostics,
      "checkpoint.payload.runtime_snapshot",
      payload.runtime_snapshot.runtime_report_complete
          ? "runtime acceptance snapshot captured"
          : "runtime acceptance snapshot incomplete");
  AppendDiagnostic(
      payload.diagnostics,
      "checkpoint.payload.cached_policy",
      payload.serialized_cached_field_mask == 0u
          ? "cached fields omitted from restart payload"
          : "cached fields serialized in restart payload");

  std::ostringstream report;
  report << "phase=" << static_cast<int>(payload.phase_id)
         << "; authoritative_field_mask=" << payload.authoritative_field_mask
         << "; serialized_cached_field_mask=" << payload.serialized_cached_field_mask
         << "; runtime_registered_stage_digest=" << payload.runtime_snapshot.registered_stage_digest
         << "; runtime_validated_required_stage_digest="
         << payload.runtime_snapshot.validated_required_stage_digest
         << "; runtime_stage_report_count=" << payload.runtime_snapshot.stage_report_count;
  payload.report_line = report.str();

  return payload;
}

dec3d::state::CanonicalState RestoreCanonicalStateFromCheckpoint(
    const CheckpointPayload& payload) {
  auto state = dec3d::state::CanonicalState::Create(payload.layout);
  if (!payload.is_complete()) {
    return state;
  }

  state.rho = payload.rho;
  state.mom_r = payload.mom_r;
  state.mom_theta = payload.mom_theta;
  state.mom_phi = payload.mom_phi;
  state.e_fluid_total = payload.e_fluid_total;
  state.e_electron = payload.e_electron;
  state.radiation_groups = payload.radiation_groups;
  state.alpha_state = payload.alpha_state;
  dec3d::state::InvalidateRecoveredAndCachedFieldsForRestart(state);
  return state;
}

RestartValidationResult ValidateRestartFromCheckpoint(
    const CheckpointPayload& payload) {
  RestartValidationResult result;
  result.validation_performed = true;
  result.runtime_snapshot_complete = RuntimeSnapshotIsUsable(payload.runtime_snapshot);
  result.payload_complete = payload.is_complete();
  result.cached_restart_truth_rejected = payload.serialized_cached_field_mask == 0u;

  if (result.payload_complete && result.cached_restart_truth_rejected) {
    const auto restored_state = RestoreCanonicalStateFromCheckpoint(payload);
    const auto summary = dec3d::state::BuildCanonicalStateSummary(restored_state);
    result.authoritative_state_restored = summary.allocated;
    result.authoritative_mask_matches_payload =
        summary.restart_authoritative_mask == payload.authoritative_field_mask;
    result.cached_fields_invalid_after_restart = summary.cached_valid_mask == 0u;
    result.restart_invalidated_cached_mask =
        dec3d::state::RestartInvalidatedCachedMask(restored_state);
  } else {
    const auto restart_state = dec3d::state::CanonicalState::Create(payload.layout);
    result.restart_invalidated_cached_mask =
        dec3d::state::RestartInvalidatedCachedMask(restart_state);
  }

  result.success =
      result.payload_complete &&
      result.runtime_snapshot_complete &&
      result.authoritative_state_restored &&
      result.authoritative_mask_matches_payload &&
      result.cached_restart_truth_rejected &&
      result.cached_fields_invalid_after_restart;

  if (!result.runtime_snapshot_complete) {
    result.failure_reason = "runtime snapshot is incomplete";
  } else if (!result.payload_complete) {
    result.failure_reason = "checkpoint payload is incomplete";
  } else if (!result.cached_restart_truth_rejected) {
    result.failure_reason = "cached fields were serialized as restart truth";
  } else if (!result.authoritative_state_restored) {
    result.failure_reason = "authoritative state restoration failed";
  } else if (!result.authoritative_mask_matches_payload) {
    result.failure_reason = "authoritative restart mask does not match payload";
  } else if (!result.cached_fields_invalid_after_restart) {
    result.failure_reason = "cached fields remained valid after restart";
  }

  std::ostringstream report;
  report << "payload_complete=" << (result.payload_complete ? "true" : "false")
         << "; runtime_snapshot_complete=" << (result.runtime_snapshot_complete ? "true" : "false")
         << "; authoritative_state_restored=" << (result.authoritative_state_restored ? "true" : "false")
         << "; authoritative_mask_matches_payload="
         << (result.authoritative_mask_matches_payload ? "true" : "false")
         << "; cached_restart_truth_rejected="
         << (result.cached_restart_truth_rejected ? "true" : "false")
         << "; cached_fields_invalid_after_restart="
         << (result.cached_fields_invalid_after_restart ? "true" : "false")
         << "; restart_invalidated_cached_mask="
         << result.restart_invalidated_cached_mask;
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

CheckpointContractReport BuildCheckpointContractReport(
    const CheckpointPayload& payload) {
  CheckpointContractReport report;
  report.built = true;
  report.payload_complete = payload.is_complete();
  report.authoritative_field_mask = payload.authoritative_field_mask;
  report.serialized_cached_field_mask = payload.serialized_cached_field_mask;
  report.runtime_snapshot_complete = RuntimeSnapshotIsUsable(payload.runtime_snapshot);
  report.runtime_registered_stage_digest = payload.runtime_snapshot.registered_stage_digest;
  report.runtime_validated_required_stage_digest =
      payload.runtime_snapshot.validated_required_stage_digest;
  report.runtime_stage_report_count = payload.runtime_snapshot.stage_report_count;

  const auto restart_validation = ValidateRestartFromCheckpoint(payload);
  report.restart_validation_succeeded = restart_validation.success;
  report.restart_invalidated_cached_mask = restart_validation.restart_invalidated_cached_mask;

  AppendDiagnostic(
      report.diagnostics,
      "checkpoint.contract.payload",
      report.payload_complete ? "checkpoint payload complete"
                              : "checkpoint payload incomplete");
  AppendDiagnostic(
      report.diagnostics,
      "checkpoint.contract.restart_validation",
      restart_validation.success ? "restart validation succeeded"
                                 : "restart validation failed: " + restart_validation.failure_reason);
  AppendDiagnostic(
      report.diagnostics,
      "checkpoint.contract.runtime_snapshot",
      report.runtime_snapshot_complete ? "runtime snapshot complete"
                                       : "runtime snapshot incomplete");

  report.diagnostics_complete =
      HasDiagnosticCode(report.diagnostics, "checkpoint.contract.payload") &&
      HasDiagnosticCode(report.diagnostics, "checkpoint.contract.restart_validation") &&
      HasDiagnosticCode(report.diagnostics, "checkpoint.contract.runtime_snapshot");

  report.success =
      report.payload_complete &&
      report.restart_validation_succeeded &&
      report.runtime_snapshot_complete &&
      report.diagnostics_complete;

  if (!report.payload_complete) {
    report.failure_reason = "checkpoint payload is incomplete";
  } else if (!report.runtime_snapshot_complete) {
    report.failure_reason = "runtime snapshot is incomplete";
  } else if (!report.restart_validation_succeeded) {
    report.failure_reason = restart_validation.failure_reason;
  } else if (!report.diagnostics_complete) {
    report.failure_reason = "checkpoint diagnostics are incomplete";
  }

  std::ostringstream line;
  line << "payload_complete=" << (report.payload_complete ? "true" : "false")
       << "; restart_validation_succeeded="
       << (report.restart_validation_succeeded ? "true" : "false")
       << "; runtime_snapshot_complete=" << (report.runtime_snapshot_complete ? "true" : "false")
       << "; diagnostics_complete=" << (report.diagnostics_complete ? "true" : "false")
       << "; authoritative_field_mask=" << report.authoritative_field_mask
       << "; serialized_cached_field_mask=" << report.serialized_cached_field_mask
       << "; restart_invalidated_cached_mask=" << report.restart_invalidated_cached_mask
       << "; runtime_registered_stage_digest=" << report.runtime_registered_stage_digest
       << "; runtime_validated_required_stage_digest="
       << report.runtime_validated_required_stage_digest
       << "; runtime_stage_report_count=" << report.runtime_stage_report_count;
  if (!report.failure_reason.empty()) {
    line << "; failure_reason=" << report.failure_reason;
  }
  report.report_line = line.str();

  return report;
}

}  // namespace dec3d::runtime
