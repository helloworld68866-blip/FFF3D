#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::InitializeStageMachine;
    using dec3d::runtime::RegisterStage;

    const auto p0_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    const auto p0_stage_machine = InitializeStageMachine(p0_runtime);
    DEC3D_CHECK(p0_stage_machine.is_complete());
    DEC3D_CHECK(p0_stage_machine.success);
    DEC3D_CHECK(p0_stage_machine.initialized);
    DEC3D_CHECK(p0_stage_machine.zero_physics_mode);
    DEC3D_CHECK_EQ(p0_stage_machine.registry_entry_count, static_cast<std::size_t>(0));
    DEC3D_CHECK_EQ(p0_stage_machine.registered_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(p0_stage_machine.required_stage_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK(p0_stage_machine.execution_evidence.is_present());

    auto p1_runtime = CreateRuntimeScaffold(PhaseId::p1, "p0-v0.1");
    const auto p1_stage_machine = InitializeStageMachine(p1_runtime);
    DEC3D_CHECK(p1_stage_machine.is_complete());
    DEC3D_CHECK(!p1_stage_machine.success);
    DEC3D_CHECK(!p1_stage_machine.initialized);
    DEC3D_CHECK(!p1_stage_machine.zero_physics_mode);
    DEC3D_CHECK(p1_stage_machine.missing_required_stage_mask != 0u);
    DEC3D_CHECK(!p1_stage_machine.failure_reason.empty());

    auto p2_runtime = CreateRuntimeScaffold(PhaseId::p2, "p0-v0.1");
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::hydro));
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::thermal));
    DEC3D_CHECK(RegisterStage(p2_runtime, StageId::equilibration));
    const auto p2_stage_machine = InitializeStageMachine(p2_runtime);
    DEC3D_CHECK(p2_stage_machine.is_complete());
    DEC3D_CHECK(p2_stage_machine.success);
    DEC3D_CHECK(p2_stage_machine.initialized);
    DEC3D_CHECK_EQ(p2_stage_machine.registry_entry_count, static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(
        p2_stage_machine.required_stage_mask,
        p2_stage_machine.registered_stage_mask);
    DEC3D_CHECK_EQ(
        p2_stage_machine.registry_entries.size(),
        static_cast<std::size_t>(3));
    for (const auto& entry : p2_stage_machine.registry_entries) {
      DEC3D_CHECK(entry.is_complete());
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
