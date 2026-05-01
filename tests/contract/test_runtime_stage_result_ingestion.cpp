#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

dec3d::core::StageResult MakeStageResult(
    dec3d::core::AuthoritativeFieldMask updated_fields,
    const char* implementation_id) {
  dec3d::core::DiagnosticsPayload diagnostics;
  diagnostics.entries.push_back({"runtime.ingestion.contract", implementation_id});

  dec3d::core::ExecutionEvidence execution_evidence;
  execution_evidence.entered_stage = true;
  execution_evidence.implementation_id = implementation_id;
  execution_evidence.touched_cell_count = 4;

  return dec3d::core::StageResult::Successful(
      updated_fields,
      diagnostics,
      execution_evidence);
}

}  // namespace

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::Combine;
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::BuildStageReportAggregation;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::HasIngestedStageResult;
    using dec3d::runtime::IngestStageResult;
    using dec3d::runtime::RegisterStage;
    using dec3d::runtime::TryIngestStageResult;

    const auto hydro_stage_mask =
        static_cast<std::uint32_t>(1u << static_cast<std::uint32_t>(StageId::hydro));
    const auto expected_p2_stage_mask =
        static_cast<std::uint32_t>(
            (1u << static_cast<std::uint32_t>(StageId::hydro)) |
            (1u << static_cast<std::uint32_t>(StageId::thermal)) |
            (1u << static_cast<std::uint32_t>(StageId::equilibration)));

    auto p2_runtime = CreateRuntimeScaffold(PhaseId::p2, "p0-v0.1");
    const auto unregistered_ingestion = TryIngestStageResult(
        p2_runtime,
        StageId::hydro,
        MakeStageResult(Combine({AuthoritativeField::rho}), "p0.runtime.stage.hydro"));
    DEC3D_CHECK(unregistered_ingestion.is_complete());
    DEC3D_CHECK(!unregistered_ingestion.accepted);
    DEC3D_CHECK(!unregistered_ingestion.stage_registered);
    DEC3D_CHECK(!unregistered_ingestion.failure_reason.empty());

    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::hydro));
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::thermal));
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::equilibration));

    dec3d::core::StageResult incomplete_result;
    incomplete_result.success = true;
    incomplete_result.updated_fields = Combine({AuthoritativeField::rho});
    DEC3D_CHECK(!IngestStageResult(p2_runtime, StageId::hydro, incomplete_result));

    dec3d::core::StageResult empty_diagnostics_result;
    empty_diagnostics_result.success = true;
    empty_diagnostics_result.updated_fields = Combine({AuthoritativeField::rho});
    empty_diagnostics_result.execution_evidence = dec3d::core::ExecutionEvidence{
        true,
        "p0.runtime.stage.empty_diagnostics",
        4};
    const auto empty_diagnostics_ingestion = TryIngestStageResult(
        p2_runtime,
        StageId::hydro,
        empty_diagnostics_result);
    DEC3D_CHECK(empty_diagnostics_ingestion.is_complete());
    DEC3D_CHECK(!empty_diagnostics_ingestion.accepted);
    DEC3D_CHECK(!empty_diagnostics_ingestion.stage_result_complete);
    DEC3D_CHECK(!empty_diagnostics_ingestion.failure_reason.empty());

    const auto hydro_result = MakeStageResult(
        Combine({AuthoritativeField::rho, AuthoritativeField::e_electron}),
        "p0.runtime.stage.hydro");
    const auto thermal_result = MakeStageResult(
        Combine({AuthoritativeField::e_fluid_total}),
        "p0.runtime.stage.thermal");
    const auto equilibration_result = MakeStageResult(
        Combine({AuthoritativeField::e_electron}),
        "p0.runtime.stage.equilibration");

    DEC3D_CHECK(IngestStageResult(p2_runtime, StageId::hydro, hydro_result));
    DEC3D_CHECK(HasIngestedStageResult(p2_runtime, StageId::hydro));
    DEC3D_CHECK(!IngestStageResult(p2_runtime, StageId::hydro, hydro_result));
    DEC3D_CHECK(!IngestStageResult(p2_runtime, StageId::radiation, hydro_result));

    auto forbidden_p0_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    DEC3D_CHECK(RegisterStage(forbidden_p0_runtime, StageId::hydro));
    const auto forbidden_ingestion = TryIngestStageResult(
        forbidden_p0_runtime,
        StageId::hydro,
        hydro_result);
    DEC3D_CHECK(forbidden_ingestion.is_complete());
    DEC3D_CHECK(!forbidden_ingestion.accepted);
    DEC3D_CHECK(forbidden_ingestion.stage_registered);
    DEC3D_CHECK(!forbidden_ingestion.stage_required_for_phase);
    DEC3D_CHECK(!forbidden_ingestion.failure_reason.empty());

    const auto partial_aggregation = BuildStageReportAggregation(p2_runtime);
    DEC3D_CHECK(partial_aggregation.is_complete());
    DEC3D_CHECK(!partial_aggregation.success);
    DEC3D_CHECK(!partial_aggregation.stage_result_ingestion_complete);
    DEC3D_CHECK_EQ(partial_aggregation.ingested_stage_result_count, static_cast<std::size_t>(1));
    DEC3D_CHECK_EQ(partial_aggregation.ingested_stage_mask, hydro_stage_mask);
    DEC3D_CHECK(partial_aggregation.missing_stage_report_mask != 0u);
    DEC3D_CHECK(!partial_aggregation.ingestion_failure_reason.empty());

    DEC3D_CHECK(IngestStageResult(p2_runtime, StageId::thermal, thermal_result));
    DEC3D_CHECK(IngestStageResult(p2_runtime, StageId::equilibration, equilibration_result));

    const auto complete_aggregation = BuildStageReportAggregation(p2_runtime);
    DEC3D_CHECK(complete_aggregation.is_complete());
    DEC3D_CHECK(complete_aggregation.success);
    DEC3D_CHECK(complete_aggregation.stage_result_ingestion_complete);
    DEC3D_CHECK(complete_aggregation.execution_evidence_complete);
    DEC3D_CHECK_EQ(complete_aggregation.ingested_stage_result_count, static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(complete_aggregation.ingested_execution_evidence_count, static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(complete_aggregation.registered_stage_digest, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(complete_aggregation.validated_required_stage_digest, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(complete_aggregation.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(complete_aggregation.ingested_stage_mask, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(complete_aggregation.missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(
        complete_aggregation.authoritative_write_set_digest,
        Combine({
            AuthoritativeField::rho,
            AuthoritativeField::e_electron,
            AuthoritativeField::e_fluid_total}));

    const auto p2_report = BuildRuntimeSubstrateReport(p2_runtime);
    DEC3D_CHECK(p2_report.is_complete());
    DEC3D_CHECK(p2_report.diagnostics_sink_initialized);
    DEC3D_CHECK(p2_report.stage_result_ingestion_complete);
    DEC3D_CHECK(p2_report.execution_evidence_complete);
    DEC3D_CHECK_EQ(p2_report.ingested_stage_result_count, static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(p2_report.ingested_execution_evidence_count, static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(p2_report.registered_stage_digest, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(p2_report.validated_required_stage_digest, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(p2_report.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p2_report.ingested_stage_mask, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(p2_report.registered_stage_mask, expected_p2_stage_mask);
    DEC3D_CHECK_EQ(p2_report.missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(
        p2_report.authoritative_write_set_digest,
        Combine({
            AuthoritativeField::rho,
            AuthoritativeField::e_electron,
            AuthoritativeField::e_fluid_total}));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
