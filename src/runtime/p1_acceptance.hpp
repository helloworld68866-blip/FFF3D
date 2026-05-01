#pragma once

#include "runtime/p0_acceptance.hpp"

#include <filesystem>

namespace dec3d::runtime {

struct HydroStageArtifact {
  bool recorded{false};
  bool runtime_executed{false};
  bool diagnostics_complete{false};
  bool authorized_write_set{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct AleRuntimeCommitArtifact {
  bool recorded{false};
  bool contract_satisfied{false};
  bool diagnostics_complete{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MacroZoningArtifact {
  bool recorded{false};
  bool substrate_complete{false};
  bool single_rank_coupling_complete{false};
  bool mpi_case2_complete{false};
  bool mpi_case3_complete{false};
  bool mpi_case4_complete{false};
  bool diagnostics_complete{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RadialMpiParityArtifact {
  bool recorded{false};
  bool case2_complete{false};
  bool case3_complete{false};
  bool case4_complete{false};
  bool diagnostics_complete{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct HydroNumericalCheckArtifact {
  bool recorded{false};
  bool case2_complete{false};
  bool case3_complete{false};
  bool case4_complete{false};
  bool budget_complete{false};
  bool diagnostics_complete{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct FieldUpdateSummaryArtifact {
  bool recorded{false};
  bool authoritative_write_mask_nonzero{false};
  bool diagnostics_complete{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct P1AcceptanceSummary {
  bool built{false};
  bool success{false};
  bool build_artifact_complete{false};
  bool failing_first_evidence_complete{false};
  bool hydro_stage_complete{false};
  bool ale_artifact_complete{false};
  bool macro_zoning_artifact_complete{false};
  bool radial_mpi_parity_complete{false};
  bool hydro_numerical_check_complete{false};
  bool field_update_summary_complete{false};
  bool assumption_delta_complete{false};
  bool workstream7_complete{false};
  bool full_p1_closeout_ready{false};
  bool diagnostics_complete{false};
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] P1AcceptanceSummary BuildP1AcceptanceSummary(
    const BuildVerificationArtifact& build_artifact,
    const FailingFirstEvidenceArtifact& failing_first_evidence,
    const HydroStageArtifact& hydro_stage_artifact,
    const AleRuntimeCommitArtifact& ale_artifact,
    const MacroZoningArtifact& macro_zoning_artifact,
    const RadialMpiParityArtifact& radial_mpi_parity_artifact,
    const HydroNumericalCheckArtifact& hydro_numerical_check_artifact,
    const FieldUpdateSummaryArtifact& field_update_artifact,
    const AssumptionLedgerDeltaArtifact& assumption_delta) noexcept;

[[nodiscard]] bool WriteP1AcceptanceArtifacts(
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
    const P1AcceptanceSummary& summary) noexcept;

}  // namespace dec3d::runtime
