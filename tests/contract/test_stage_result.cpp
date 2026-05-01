#include "core/diagnostics/stage_contracts.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using namespace dec3d::core;

    DiagnosticsPayload diagnostics;
    diagnostics.entries.push_back({"stage.result.contract", "stage result carries diagnostics"});

    ExecutionEvidence execution_evidence;
    execution_evidence.entered_stage = true;
    execution_evidence.implementation_id = "p0.contract.stage_result";
    execution_evidence.touched_cell_count = 8;

    const auto updated_fields = Combine({AuthoritativeField::rho, AuthoritativeField::e_electron});
    const StageResult result = StageResult::Successful(updated_fields, diagnostics, execution_evidence);

    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_semantically_complete());
    DEC3D_CHECK(result.diagnostics.has_entries());
    DEC3D_CHECK(result.execution_evidence.has_value());
    DEC3D_CHECK(result.execution_evidence->is_present());
    DEC3D_CHECK_EQ(result.updated_fields, updated_fields);

    const StageResult failed = StageResult::Failed(
        "stage contract failed for testing",
        diagnostics,
        execution_evidence);
    DEC3D_CHECK(!failed.success);
    DEC3D_CHECK(!failed.failure_reason.empty());
    DEC3D_CHECK(failed.diagnostics.has_entries());
    DEC3D_CHECK(failed.execution_evidence.has_value());
    DEC3D_CHECK_EQ(failed.updated_fields, static_cast<AuthoritativeFieldMask>(0));
    DEC3D_CHECK(!failed.is_semantically_complete());

    const StageContext context{
        0.0,
        1.0e-9,
        1,
        PhaseId::p0,
        "p0-v0.1",
        "mesh:initial",
        "ownership:rank0",
        "diagnostics:memory"};
    DEC3D_CHECK(context.is_complete());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
