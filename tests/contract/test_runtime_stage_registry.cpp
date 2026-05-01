#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::IsStageRegistered;
    using dec3d::runtime::RegisterStage;
    using dec3d::runtime::ValidateStageRegistry;

    auto p0_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    DEC3D_CHECK(!IsStageRegistered(p0_runtime, StageId::hydro));

    const auto p0_validation = ValidateStageRegistry(p0_runtime);
    DEC3D_CHECK(p0_validation.success);
    DEC3D_CHECK_EQ(p0_validation.required_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_validation.registered_stage_mask, static_cast<std::uint32_t>(0));

    const auto p0_report = BuildRuntimeSubstrateReport(p0_runtime);
    DEC3D_CHECK(p0_report.is_complete());
    DEC3D_CHECK(p0_report.zero_physics_mode);
    DEC3D_CHECK(p0_report.phase_id_valid);
    DEC3D_CHECK(p0_report.validation_performed);
    DEC3D_CHECK_EQ(p0_report.registered_stage_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_report.required_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_report.registered_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(p0_report.required_stage_contract_satisfied);
    DEC3D_CHECK(p0_report.execution_evidence.is_present());
    DEC3D_CHECK(p0_report.stage_report_aggregation_performed);
    DEC3D_CHECK(p0_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(p0_report.diagnostics_sink_initialized);
    DEC3D_CHECK(p0_report.stage_result_ingestion_complete);
    DEC3D_CHECK(p0_report.execution_evidence_complete);
    DEC3D_CHECK(p0_report.diagnostics_complete);
    DEC3D_CHECK_EQ(p0_report.ingested_stage_result_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_report.authoritative_write_set_digest, static_cast<dec3d::core::AuthoritativeFieldMask>(0));

    DEC3D_CHECK(RegisterStage(p0_runtime, StageId::hydro));
    DEC3D_CHECK(IsStageRegistered(p0_runtime, StageId::hydro));
    DEC3D_CHECK(!RegisterStage(p0_runtime, StageId::hydro));

    const auto forbidden_validation = ValidateStageRegistry(p0_runtime);
    DEC3D_CHECK(!forbidden_validation.success);
    DEC3D_CHECK_EQ(forbidden_validation.registered_stage_mask, static_cast<std::uint32_t>(1u));
    DEC3D_CHECK(forbidden_validation.forbidden_registered_stage_mask != 0u);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
