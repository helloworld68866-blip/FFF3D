#include "runtime/checkpoint_contract.hpp"
#include "runtime/p0_acceptance.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

dec3d::runtime::BuildVerificationArtifact MakeBuildArtifact() {
  dec3d::runtime::BuildVerificationArtifact artifact;
  artifact.recorded = true;
  artifact.clean_rebuild_verified = true;
  artifact.report_line = "clean_rebuild=true; generator=vs2022";
  return artifact;
}

dec3d::runtime::FailingFirstEvidenceArtifact MakeFailingFirstArtifact() {
  dec3d::runtime::FailingFirstEvidenceArtifact artifact;
  artifact.recorded = true;
  artifact.failing_test_count = 1;
  artifact.report_line = "failing_test=dec3d_acceptance_phase_p0; rerun=green";
  return artifact;
}

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::BuildCheckpointContractReport;
    using dec3d::runtime::BuildCheckpointPayload;
    using dec3d::runtime::BuildEmptyP0AssumptionLedgerDelta;
    using dec3d::runtime::BuildP0AcceptanceSummary;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::RegisterStage;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    const auto runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    const auto checkpoint_report =
        BuildCheckpointContractReport(BuildCheckpointPayload(state, runtime));
    const auto runtime_report = BuildRuntimeSubstrateReport(runtime);

    const auto success_summary = BuildP0AcceptanceSummary(
        MakeBuildArtifact(),
        MakeFailingFirstArtifact(),
        runtime_report,
        checkpoint_report,
        BuildEmptyP0AssumptionLedgerDelta());

    DEC3D_CHECK(success_summary.is_complete());
    DEC3D_CHECK(success_summary.success);
    DEC3D_CHECK(success_summary.build_artifact_complete);
    DEC3D_CHECK(success_summary.failing_first_evidence_complete);
    DEC3D_CHECK(success_summary.runtime_report_complete);
    DEC3D_CHECK(success_summary.runtime_contract_satisfied);
    DEC3D_CHECK(success_summary.checkpoint_report_complete);
    DEC3D_CHECK(success_summary.checkpoint_contract_satisfied);
    DEC3D_CHECK(success_summary.diagnostics_complete);
    DEC3D_CHECK(success_summary.assumption_delta_complete);
    DEC3D_CHECK(success_summary.assumption_delta_empty);
    DEC3D_CHECK(success_summary.zero_physics_runtime);
    DEC3D_CHECK(success_summary.no_physics_stage_registered);
    DEC3D_CHECK_EQ(success_summary.runtime_ingested_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(success_summary.runtime_missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(
        success_summary.runtime_authoritative_write_set_digest,
        static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(success_summary.failure_reason.empty());

    auto forbidden_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    DEC3D_CHECK(RegisterStage(forbidden_runtime, StageId::hydro));
    const auto forbidden_runtime_report = BuildRuntimeSubstrateReport(forbidden_runtime);
    const auto forbidden_checkpoint_report =
        BuildCheckpointContractReport(BuildCheckpointPayload(state, forbidden_runtime));

    const auto forbidden_summary = BuildP0AcceptanceSummary(
        MakeBuildArtifact(),
        MakeFailingFirstArtifact(),
        forbidden_runtime_report,
        forbidden_checkpoint_report,
        BuildEmptyP0AssumptionLedgerDelta());

    DEC3D_CHECK(forbidden_summary.is_complete());
    DEC3D_CHECK(!forbidden_summary.success);
    DEC3D_CHECK(!forbidden_summary.runtime_contract_satisfied);
    DEC3D_CHECK(!forbidden_summary.no_physics_stage_registered);
    DEC3D_CHECK(forbidden_summary.runtime_ingested_stage_mask == 0u);
    DEC3D_CHECK(!forbidden_summary.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
