#include "runtime/p1_acceptance.hpp"

#include <fstream>
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

[[nodiscard]] bool EnsureDirectory(const std::filesystem::path& path) noexcept {
  std::error_code error;
  std::filesystem::create_directories(path, error);
  return !error;
}

[[nodiscard]] bool WriteMarkdownArtifact(
    const std::filesystem::path& path,
    const std::string& title,
    const std::string& report_line) noexcept {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }
  output << "# " << title << "\n\n";
  output << "```text\n" << report_line << "\n```\n";
  return true;
}

}  // namespace

bool HydroStageArtifact::is_complete() const noexcept {
  return recorded && runtime_executed && diagnostics_complete && authorized_write_set && !report_line.empty();
}

bool AleRuntimeCommitArtifact::is_complete() const noexcept {
  return recorded && contract_satisfied && diagnostics_complete && !report_line.empty();
}

bool MacroZoningArtifact::is_complete() const noexcept {
  return recorded &&
         substrate_complete &&
         single_rank_coupling_complete &&
         mpi_case2_complete &&
         mpi_case3_complete &&
         mpi_case4_complete &&
         diagnostics_complete &&
         !report_line.empty();
}

bool RadialMpiParityArtifact::is_complete() const noexcept {
  return recorded &&
         case2_complete &&
         case3_complete &&
         case4_complete &&
         diagnostics_complete &&
         !report_line.empty();
}

bool HydroNumericalCheckArtifact::is_complete() const noexcept {
  return recorded &&
         case2_complete &&
         case3_complete &&
         case4_complete &&
         budget_complete &&
         diagnostics_complete &&
         !report_line.empty();
}

bool FieldUpdateSummaryArtifact::is_complete() const noexcept {
  return recorded &&
         authoritative_write_mask_nonzero &&
         diagnostics_complete &&
         !report_line.empty();
}

bool P1AcceptanceSummary::is_complete() const noexcept {
  return built && diagnostics_complete && !report_line.empty();
}

P1AcceptanceSummary BuildP1AcceptanceSummary(
    const BuildVerificationArtifact& build_artifact,
    const FailingFirstEvidenceArtifact& failing_first_evidence,
    const HydroStageArtifact& hydro_stage_artifact,
    const AleRuntimeCommitArtifact& ale_artifact,
    const MacroZoningArtifact& macro_zoning_artifact,
    const RadialMpiParityArtifact& radial_mpi_parity_artifact,
    const HydroNumericalCheckArtifact& hydro_numerical_check_artifact,
    const FieldUpdateSummaryArtifact& field_update_artifact,
    const AssumptionLedgerDeltaArtifact& assumption_delta) noexcept {
  P1AcceptanceSummary summary;
  summary.built = true;
  summary.build_artifact_complete = build_artifact.is_complete();
  summary.failing_first_evidence_complete = failing_first_evidence.is_complete();
  summary.hydro_stage_complete = hydro_stage_artifact.is_complete();
  summary.ale_artifact_complete = ale_artifact.is_complete();
  summary.macro_zoning_artifact_complete = macro_zoning_artifact.is_complete();
  summary.radial_mpi_parity_complete = radial_mpi_parity_artifact.is_complete();
  summary.hydro_numerical_check_complete = hydro_numerical_check_artifact.is_complete();
  summary.field_update_summary_complete = field_update_artifact.is_complete();
  summary.assumption_delta_complete = assumption_delta.is_complete();

  summary.workstream7_complete =
      summary.build_artifact_complete &&
      summary.failing_first_evidence_complete &&
      summary.hydro_stage_complete &&
      summary.ale_artifact_complete &&
      summary.macro_zoning_artifact_complete &&
      summary.radial_mpi_parity_complete &&
      summary.hydro_numerical_check_complete &&
      summary.field_update_summary_complete &&
      summary.assumption_delta_complete;

  summary.full_p1_closeout_ready = false;

  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.build",
      summary.build_artifact_complete ? "clean rebuild verified" : "clean rebuild verification missing");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.failing_first",
      summary.failing_first_evidence_complete ? "failing-first evidence present" : "failing-first evidence missing");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.hydro_stage",
      summary.hydro_stage_complete ? "hydro stage evidence present" : "hydro stage evidence missing");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.ale",
      summary.ale_artifact_complete ? "ALE artifact complete" : "ALE artifact incomplete");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.macro_zoning",
      summary.macro_zoning_artifact_complete ? "macro-zoning artifact complete" : "macro-zoning artifact incomplete");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.radial_mpi",
      summary.radial_mpi_parity_complete ? "radial MPI parity artifact complete" : "radial MPI parity artifact incomplete");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.numerics",
      summary.hydro_numerical_check_complete ? "hydro numerical checks complete" : "hydro numerical checks incomplete");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.field_update",
      summary.field_update_summary_complete ? "field update summary complete" : "field update summary incomplete");
  AppendDiagnostic(
      summary.diagnostics,
      "p1.acceptance.assumption_delta",
      summary.assumption_delta_complete ? "assumption-ledger delta recorded" : "assumption-ledger delta missing");

  summary.diagnostics_complete =
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.build") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.failing_first") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.hydro_stage") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.ale") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.macro_zoning") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.radial_mpi") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.numerics") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.field_update") &&
      HasDiagnosticCode(summary.diagnostics, "p1.acceptance.assumption_delta");

  summary.success = summary.workstream7_complete && summary.diagnostics_complete;

  if (!summary.build_artifact_complete) {
    summary.failure_reason = "clean rebuild verification is incomplete";
  } else if (!summary.failing_first_evidence_complete) {
    summary.failure_reason = "failing-first evidence is incomplete";
  } else if (!summary.hydro_stage_complete) {
    summary.failure_reason = "hydro stage artifact is incomplete";
  } else if (!summary.ale_artifact_complete) {
    summary.failure_reason = "ALE artifact is incomplete";
  } else if (!summary.macro_zoning_artifact_complete) {
    summary.failure_reason = "macro-zoning artifact is incomplete";
  } else if (!summary.radial_mpi_parity_complete) {
    summary.failure_reason = "radial MPI parity artifact is incomplete";
  } else if (!summary.hydro_numerical_check_complete) {
    summary.failure_reason = "hydro numerical check artifact is incomplete";
  } else if (!summary.field_update_summary_complete) {
    summary.failure_reason = "field update summary artifact is incomplete";
  } else if (!summary.assumption_delta_complete) {
    summary.failure_reason = "assumption-ledger delta is incomplete";
  } else if (!summary.diagnostics_complete) {
    summary.failure_reason = "acceptance diagnostics are incomplete";
  }

  std::ostringstream report;
  report << "build_artifact_complete=" << (summary.build_artifact_complete ? "true" : "false")
         << "; failing_first_evidence_complete=" << (summary.failing_first_evidence_complete ? "true" : "false")
         << "; hydro_stage_complete=" << (summary.hydro_stage_complete ? "true" : "false")
         << "; ale_artifact_complete=" << (summary.ale_artifact_complete ? "true" : "false")
         << "; macro_zoning_artifact_complete=" << (summary.macro_zoning_artifact_complete ? "true" : "false")
         << "; radial_mpi_parity_complete=" << (summary.radial_mpi_parity_complete ? "true" : "false")
         << "; hydro_numerical_check_complete=" << (summary.hydro_numerical_check_complete ? "true" : "false")
         << "; field_update_summary_complete=" << (summary.field_update_summary_complete ? "true" : "false")
         << "; assumption_delta_complete=" << (summary.assumption_delta_complete ? "true" : "false")
         << "; workstream7_complete=" << (summary.workstream7_complete ? "true" : "false")
         << "; full_p1_closeout_ready=" << (summary.full_p1_closeout_ready ? "true" : "false")
         << "; diagnostics_complete=" << (summary.diagnostics_complete ? "true" : "false");
  if (!summary.failure_reason.empty()) {
    report << "; failure_reason=" << summary.failure_reason;
  }
  summary.report_line = report.str();
  return summary;
}

bool WriteP1AcceptanceArtifacts(
    const std::filesystem::path& artifact_directory,
    const BuildVerificationArtifact& build_artifact,
    const FailingFirstEvidenceArtifact& failing_first_evidence,
    const HydroStageArtifact& hydro_stage_artifact,
    const AleRuntimeCommitArtifact& ale_artifact,
    const MacroZoningArtifact& macro_zoning_artifact,
    const RadialMpiParityArtifact& radial_mpi_parity_artifact,
    const HydroNumericalCheckArtifact& hydro_numerical_check_artifact,
    const FieldUpdateSummaryArtifact& field_update_artifact,
    const AssumptionLedgerDeltaArtifact& assumption_delta,
    const P1AcceptanceSummary& summary) noexcept {
  if (!EnsureDirectory(artifact_directory)) {
    return false;
  }

  return
      WriteMarkdownArtifact(artifact_directory / "build-report.md", "Build Report", build_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "failing-first-evidence.md", "Failing-First Evidence", failing_first_evidence.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "hydro-stage-report.md", "Hydro Stage Report", hydro_stage_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "ale-runtime-commit-report.md", "ALE Runtime Commit Report", ale_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "macro-zoning-report.md", "Macro-Zoning Report", macro_zoning_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "radial-mpi-parity-report.md", "Radial MPI Parity Report", radial_mpi_parity_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "hydro-numerical-check-summary.md", "Hydro Numerical Check Summary", hydro_numerical_check_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "field-update-summary.md", "Field Update Summary", field_update_artifact.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "assumption-ledger-delta.md", "Assumption Ledger Delta", assumption_delta.report_line) &&
      WriteMarkdownArtifact(artifact_directory / "p1-closeout-summary.md", "P1 Closeout Summary", summary.report_line);
}

}  // namespace dec3d::runtime
