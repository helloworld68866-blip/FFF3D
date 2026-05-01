#pragma once

#include "runtime/checkpoint_contract.hpp"

namespace dec3d::runtime {

struct BuildVerificationArtifact {
  bool recorded{false};
  bool clean_rebuild_verified{false};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct FailingFirstEvidenceArtifact {
  bool recorded{false};
  std::size_t failing_test_count{0};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct AssumptionLedgerDeltaArtifact {
  bool recorded{false};
  bool empty{false};
  std::size_t entry_count{0};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct P0AcceptanceSummary {
  bool built{false};
  bool success{false};
  bool build_artifact_complete{false};
  bool failing_first_evidence_complete{false};
  bool runtime_report_complete{false};
  bool runtime_contract_satisfied{false};
  bool checkpoint_report_complete{false};
  bool checkpoint_contract_satisfied{false};
  bool diagnostics_complete{false};
  bool assumption_delta_complete{false};
  bool assumption_delta_empty{false};
  bool zero_physics_runtime{false};
  bool no_physics_stage_registered{false};
  std::uint32_t runtime_registered_stage_digest{0};
  std::uint32_t runtime_validated_required_stage_digest{0};
  std::uint32_t runtime_ingested_stage_mask{0};
  std::uint32_t runtime_missing_stage_report_mask{0};
  std::size_t runtime_stage_report_count{0};
  dec3d::core::AuthoritativeFieldMask runtime_authoritative_write_set_digest{0};
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] AssumptionLedgerDeltaArtifact BuildEmptyP0AssumptionLedgerDelta() noexcept;

[[nodiscard]] P0AcceptanceSummary BuildP0AcceptanceSummary(
    const BuildVerificationArtifact& build_artifact,
    const FailingFirstEvidenceArtifact& failing_first_evidence,
    const RuntimeSubstrateReport& runtime_report,
    const CheckpointContractReport& checkpoint_report,
    const AssumptionLedgerDeltaArtifact& assumption_delta) noexcept;

}  // namespace dec3d::runtime
