#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::runtime::BuildRuntimeAcceptanceSnapshot;
    using dec3d::runtime::BuildStageReportAggregation;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::CreateRuntimeScaffold;

    const auto p0_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    const auto p0_aggregation = BuildStageReportAggregation(p0_runtime);
    DEC3D_CHECK(p0_aggregation.is_complete());
    DEC3D_CHECK(p0_aggregation.success);
    DEC3D_CHECK(p0_aggregation.validation_succeeded);
    DEC3D_CHECK(p0_aggregation.stage_machine_initialized);
    DEC3D_CHECK(p0_aggregation.diagnostics_sink_initialized);
    DEC3D_CHECK(p0_aggregation.stage_result_ingestion_complete);
    DEC3D_CHECK(p0_aggregation.execution_evidence_complete);
    DEC3D_CHECK(p0_aggregation.diagnostics_complete);
    DEC3D_CHECK_EQ(p0_aggregation.registered_stage_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.ingested_stage_result_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.stage_report_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.registered_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.validated_required_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.ingested_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_aggregation.authoritative_write_set_digest, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(p0_aggregation.execution_evidence_summary.is_present());
    DEC3D_CHECK(p0_aggregation.diagnostics.has_entries());
    DEC3D_CHECK(p0_aggregation.validation_failure_reason.empty());
    DEC3D_CHECK(p0_aggregation.stage_machine_failure_reason.empty());
    DEC3D_CHECK(p0_aggregation.ingestion_failure_reason.empty());
    DEC3D_CHECK(p0_aggregation.report_line.find("diagnostics_complete=true") != std::string::npos);

    const auto p0_report = BuildRuntimeSubstrateReport(p0_runtime);
    DEC3D_CHECK(p0_report.stage_report_aggregation_performed);
    DEC3D_CHECK(p0_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(p0_report.diagnostics_sink_initialized);
    DEC3D_CHECK(p0_report.stage_result_ingestion_complete);
    DEC3D_CHECK(p0_report.execution_evidence_complete);
    DEC3D_CHECK(p0_report.diagnostics_complete);
    DEC3D_CHECK(p0_report.diagnostics_entry_count >= static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(p0_report.registered_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_report.validated_required_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_report.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(p0_report.validation_failure_reason.empty());
    DEC3D_CHECK(p0_report.stage_machine_failure_reason.empty());
    DEC3D_CHECK(p0_report.ingestion_failure_reason.empty());
    DEC3D_CHECK(p0_report.stage_report_failure_reason.empty());

    const auto p0_snapshot = BuildRuntimeAcceptanceSnapshot(p0_runtime);
    DEC3D_CHECK(p0_snapshot.is_complete());
    DEC3D_CHECK(p0_snapshot.runtime_report_complete);
    DEC3D_CHECK(p0_snapshot.execution_evidence_present);
    DEC3D_CHECK(p0_snapshot.diagnostics_complete);
    DEC3D_CHECK_EQ(p0_snapshot.registered_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_snapshot.validated_required_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_snapshot.ingested_stage_digest, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_snapshot.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(p0_snapshot.failure_reason.empty());

    const auto p1_runtime = CreateRuntimeScaffold(PhaseId::p1, "p0-v0.1");
    const auto p1_aggregation = BuildStageReportAggregation(p1_runtime);
    DEC3D_CHECK(p1_aggregation.is_complete());
    DEC3D_CHECK(!p1_aggregation.success);
    DEC3D_CHECK(!p1_aggregation.validation_succeeded);
    DEC3D_CHECK(!p1_aggregation.stage_machine_initialized);
    DEC3D_CHECK(p1_aggregation.diagnostics_sink_initialized);
    DEC3D_CHECK(p1_aggregation.stage_result_ingestion_complete);
    DEC3D_CHECK(p1_aggregation.diagnostics_complete);
    DEC3D_CHECK_EQ(
        p1_aggregation.validated_required_stage_digest,
        static_cast<std::uint32_t>(1u << static_cast<std::uint32_t>(dec3d::core::StageId::hydro)));
    DEC3D_CHECK_EQ(p1_aggregation.ingested_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p1_aggregation.registered_but_not_ingested_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p1_aggregation.missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(!p1_aggregation.validation_failure_reason.empty());
    DEC3D_CHECK(!p1_aggregation.stage_machine_failure_reason.empty());
    DEC3D_CHECK(!p1_aggregation.failure_reason.empty());

    const auto p1_report = BuildRuntimeSubstrateReport(p1_runtime);
    DEC3D_CHECK(p1_report.is_complete());
    DEC3D_CHECK(p1_report.stage_report_aggregation_performed);
    DEC3D_CHECK(!p1_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(p1_report.diagnostics_sink_initialized);
    DEC3D_CHECK(p1_report.stage_result_ingestion_complete);
    DEC3D_CHECK(p1_report.diagnostics_complete);
    DEC3D_CHECK(!p1_report.validation_failure_reason.empty());
    DEC3D_CHECK(!p1_report.stage_machine_failure_reason.empty());
    DEC3D_CHECK(p1_report.ingestion_failure_reason.empty());
    DEC3D_CHECK(!p1_report.stage_report_failure_reason.empty());

    auto broken_sink_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    broken_sink_runtime.diagnostics_sink.initialized = false;
    broken_sink_runtime.diagnostics_sink.sink_id.clear();
    const auto broken_sink_aggregation = BuildStageReportAggregation(broken_sink_runtime);
    DEC3D_CHECK(!broken_sink_aggregation.success);
    DEC3D_CHECK(!broken_sink_aggregation.is_complete());
    DEC3D_CHECK(!broken_sink_aggregation.diagnostics_sink_initialized);
    DEC3D_CHECK(!broken_sink_aggregation.failure_reason.empty());

    const auto broken_sink_report = BuildRuntimeSubstrateReport(broken_sink_runtime);
    DEC3D_CHECK(!broken_sink_report.is_complete());
    DEC3D_CHECK(!broken_sink_report.diagnostics_sink_initialized);
    DEC3D_CHECK(!broken_sink_report.stage_report_aggregation_succeeded);

    const auto broken_sink_snapshot = BuildRuntimeAcceptanceSnapshot(broken_sink_runtime);
    DEC3D_CHECK(broken_sink_snapshot.is_complete());
    DEC3D_CHECK(!broken_sink_snapshot.runtime_report_complete);
    DEC3D_CHECK(!broken_sink_snapshot.execution_evidence_present);
    DEC3D_CHECK(!broken_sink_snapshot.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
