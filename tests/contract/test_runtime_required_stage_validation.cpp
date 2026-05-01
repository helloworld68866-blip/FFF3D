#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::RegisterStage;
    using dec3d::runtime::ValidateStageRegistry;

    const auto p1_runtime = CreateRuntimeScaffold(PhaseId::p1, "p0-v0.1");
    const auto p1_validation = ValidateStageRegistry(p1_runtime);
    DEC3D_CHECK(!p1_validation.success);
    DEC3D_CHECK_EQ(p1_validation.registered_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(p1_validation.missing_required_stage_mask != 0u);
    DEC3D_CHECK(p1_validation.required_stage_mask != 0u);

    auto p2_runtime = CreateRuntimeScaffold(PhaseId::p2, "p0-v0.1");
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::hydro));
    const auto incomplete_p2_validation = ValidateStageRegistry(p2_runtime);
    DEC3D_CHECK(!incomplete_p2_validation.success);
    DEC3D_CHECK(incomplete_p2_validation.registered_stage_mask != 0u);
    DEC3D_CHECK(incomplete_p2_validation.missing_required_stage_mask != 0u);

    auto complete_p2_runtime = CreateRuntimeScaffold(PhaseId::p2, "p0-v0.1");
    DEC3D_CHECK(RegisterStage(complete_p2_runtime, StageId::hydro));
    DEC3D_CHECK(RegisterStage(complete_p2_runtime, StageId::thermal));
    DEC3D_CHECK(RegisterStage(complete_p2_runtime, StageId::equilibration));
    const auto complete_p2_validation = ValidateStageRegistry(complete_p2_runtime);
    DEC3D_CHECK(complete_p2_validation.success);
    DEC3D_CHECK_EQ(
        complete_p2_validation.required_stage_mask,
        complete_p2_validation.registered_stage_mask);

    const auto complete_p2_report = dec3d::runtime::BuildRuntimeSubstrateReport(complete_p2_runtime);
    DEC3D_CHECK(complete_p2_report.phase_id_valid);
    DEC3D_CHECK(complete_p2_report.validation_performed);
    DEC3D_CHECK(complete_p2_report.execution_evidence.is_present());

    const auto invalid_runtime = CreateRuntimeScaffold(static_cast<PhaseId>(255), "p0-v0.1");
    const auto invalid_validation = ValidateStageRegistry(invalid_runtime);
    DEC3D_CHECK(!invalid_validation.success);
    DEC3D_CHECK(!invalid_validation.failure_reason.empty());
    const auto invalid_report = dec3d::runtime::BuildRuntimeSubstrateReport(invalid_runtime);
    DEC3D_CHECK(!invalid_report.is_complete());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
