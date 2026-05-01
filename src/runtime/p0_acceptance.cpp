#include "runtime/p0_acceptance.hpp"

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

}  // namespace

bool BuildVerificationArtifact::is_complete() const noexcept {
  return recorded &&
         clean_rebuild_verified &&
         !report_line.empty();
}

bool FailingFirstEvidenceArtifact::is_complete() const noexcept {
  return recorded &&
         failing_test_count > 0 &&
         !report_line.empty();
}

bool AssumptionLedgerDeltaArtifact::is_complete() const noexcept {
  const bool count_matches =
      (empty && entry_count == 0) || (!empty && entry_count > 0);
  return recorded &&
         count_matches &&
         !report_line.empty();
}

bool P0AcceptanceSummary::is_complete() const noexcept {
  return built &&
         diagnostics_complete &&
         !report_line.empty();
}

AssumptionLedgerDeltaArtifact BuildEmptyP0AssumptionLedgerDelta() noexcept {
  AssumptionLedgerDeltaArtifact artifact;
  artifact.recorded = true;
  artifact.empty = true;
  artifact.entry_count = 0;
  artifact.report_line = "phase=p0; recorded=true; empty=true; entry_count=0";
  return artifact;
}

P0AcceptanceSummary BuildP0AcceptanceSummary(
    const BuildVerificationArtifact& build_artifact,
    const FailingFirstEvidenceArtifact& failing_first_evidence,
    const RuntimeSubstrateReport& runtime_report,
    const CheckpointContractReport& checkpoint_report,
    const AssumptionLedgerDeltaArtifact& assumption_delta) noexcept {
  P0AcceptanceSummary summary;
  summary.built = true;
  summary.build_artifact_complete = build_artifact.is_complete();
  summary.failing_first_evidence_complete = failing_first_evidence.is_complete();
  summary.runtime_report_complete = runtime_report.is_complete();
  summary.checkpoint_report_complete = checkpoint_report.is_complete();
  summary.assumption_delta_complete = assumption_delta.is_complete();
  summary.assumption_delta_empty = assumption_delta.empty;
  summary.zero_physics_runtime = runtime_report.zero_physics_mode;
  summary.no_physics_stage_registered = runtime_report.registered_stage_digest == 0u;
  summary.runtime_registered_stage_digest = runtime_report.registered_stage_digest;
  summary.runtime_validated_required_stage_digest =
      runtime_report.validated_required_stage_digest;
  summary.runtime_ingested_stage_mask = runtime_report.ingested_stage_mask;
  summary.runtime_missing_stage_report_mask = runtime_report.missing_stage_report_mask;
  summary.runtime_stage_report_count = runtime_report.stage_report_count;
  summary.runtime_authoritative_write_set_digest =
      runtime_report.authoritative_write_set_digest;

  summary.runtime_contract_satisfied =
      summary.runtime_report_complete &&
      runtime_report.stage_report_aggregation_succeeded &&
      runtime_report.required_stage_contract_satisfied &&
      runtime_report.zero_physics_mode &&
      runtime_report.registered_stage_digest == 0u &&
      runtime_report.ingested_stage_mask == 0u &&
      runtime_report.missing_stage_report_mask == 0u &&
      runtime_report.authoritative_write_set_digest == 0u &&
      runtime_report.diagnostics_complete;

  summary.checkpoint_contract_satisfied =
      summary.checkpoint_report_complete &&
      checkpoint_report.success &&
      checkpoint_report.serialized_cached_field_mask == 0u &&
      checkpoint_report.diagnostics_complete;

  AppendDiagnostic(
      summary.diagnostics,
      "p0.acceptance.build",
      summary.build_artifact_complete ? "clean rebuild verified"
                                      : "clean rebuild verification missing");
  AppendDiagnostic(
      summary.diagnostics,
      "p0.acceptance.failing_first",
      summary.failing_first_evidence_complete ? "failing-first evidence present"
                                              : "failing-first evidence missing");
  AppendDiagnostic(
      summary.diagnostics,
      "p0.acceptance.runtime",
      summary.runtime_contract_satisfied ? "runtime substrate gate satisfied"
                                         : "runtime substrate gate failed");
  AppendDiagnostic(
      summary.diagnostics,
      "p0.acceptance.checkpoint",
      summary.checkpoint_contract_satisfied ? "checkpoint contract gate satisfied"
                                            : "checkpoint contract gate failed");
  AppendDiagnostic(
      summary.diagnostics,
      "p0.acceptance.assumption_delta",
      summary.assumption_delta_complete ? "assumption-ledger delta recorded"
                                        : "assumption-ledger delta missing");

  summary.diagnostics_complete =
      HasDiagnosticCode(summary.diagnostics, "p0.acceptance.build") &&
      HasDiagnosticCode(summary.diagnostics, "p0.acceptance.failing_first") &&
      HasDiagnosticCode(summary.diagnostics, "p0.acceptance.runtime") &&
      HasDiagnosticCode(summary.diagnostics, "p0.acceptance.checkpoint") &&
      HasDiagnosticCode(summary.diagnostics, "p0.acceptance.assumption_delta");

  summary.success =
      summary.build_artifact_complete &&
      summary.failing_first_evidence_complete &&
      summary.runtime_contract_satisfied &&
      summary.checkpoint_contract_satisfied &&
      summary.assumption_delta_complete &&
      summary.diagnostics_complete;

  if (!summary.build_artifact_complete) {
    summary.failure_reason = "clean rebuild verification is incomplete";
  } else if (!summary.failing_first_evidence_complete) {
    summary.failure_reason = "failing-first evidence is incomplete";
  } else if (!summary.runtime_contract_satisfied) {
    summary.failure_reason = "runtime substrate gate failed";
  } else if (!summary.checkpoint_contract_satisfied) {
    summary.failure_reason = "checkpoint contract gate failed";
  } else if (!summary.assumption_delta_complete) {
    summary.failure_reason = "assumption-ledger delta is incomplete";
  } else if (!summary.diagnostics_complete) {
    summary.failure_reason = "acceptance diagnostics are incomplete";
  }

  std::ostringstream report;
  report << "build_artifact_complete=" << (summary.build_artifact_complete ? "true" : "false")
         << "; failing_first_evidence_complete="
         << (summary.failing_first_evidence_complete ? "true" : "false")
         << "; runtime_report_complete=" << (summary.runtime_report_complete ? "true" : "false")
         << "; runtime_contract_satisfied="
         << (summary.runtime_contract_satisfied ? "true" : "false")
         << "; checkpoint_report_complete="
         << (summary.checkpoint_report_complete ? "true" : "false")
         << "; checkpoint_contract_satisfied="
         << (summary.checkpoint_contract_satisfied ? "true" : "false")
         << "; assumption_delta_complete="
         << (summary.assumption_delta_complete ? "true" : "false")
         << "; assumption_delta_empty=" << (summary.assumption_delta_empty ? "true" : "false")
         << "; zero_physics_runtime=" << (summary.zero_physics_runtime ? "true" : "false")
         << "; no_physics_stage_registered="
         << (summary.no_physics_stage_registered ? "true" : "false")
         << "; diagnostics_complete=" << (summary.diagnostics_complete ? "true" : "false")
         << "; runtime_registered_stage_digest=" << summary.runtime_registered_stage_digest
         << "; runtime_validated_required_stage_digest="
         << summary.runtime_validated_required_stage_digest
         << "; runtime_ingested_stage_mask=" << summary.runtime_ingested_stage_mask
         << "; runtime_missing_stage_report_mask=" << summary.runtime_missing_stage_report_mask
         << "; runtime_stage_report_count=" << summary.runtime_stage_report_count;
  report << "; runtime_authoritative_write_set_digest="
         << summary.runtime_authoritative_write_set_digest;
  if (!summary.failure_reason.empty()) {
    report << "; failure_reason=" << summary.failure_reason;
  }
  summary.report_line = report.str();

  return summary;
}

}  // namespace dec3d::runtime
