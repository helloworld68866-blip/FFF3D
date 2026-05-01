#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/checkpoint_contract.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::core::Combine;
    using dec3d::core::PhaseId;
    using dec3d::runtime::BuildCheckpointContractReport;
    using dec3d::runtime::BuildCheckpointPayload;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::RestoreCanonicalStateFromCheckpoint;
    using dec3d::runtime::ValidateRestartFromCheckpoint;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    state.rho(0, 0, 0) = 1.25;
    state.e_fluid_total(0, 0, 0) = 7.5;
    state.e_electron(0, 0, 0) = 2.5;
    state.radiation_groups[0](0, 0, 0) = 0.75;
    state.alpha_state.storage(0, 0, 0) = 0.5;
    state.MarkCachedFieldValid(CachedField::electron_temperature);
    state.MarkCachedFieldValid(CachedField::chi_e);

    auto runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");

    const auto payload = BuildCheckpointPayload(state, runtime);
    DEC3D_CHECK(payload.is_complete());
    DEC3D_CHECK(payload.runtime_snapshot.is_complete());
    DEC3D_CHECK(payload.diagnostics.has_entries());
    DEC3D_CHECK_EQ(
        payload.authoritative_field_mask,
        Combine({
            AuthoritativeField::rho,
            AuthoritativeField::mom_r,
            AuthoritativeField::mom_theta,
            AuthoritativeField::mom_phi,
            AuthoritativeField::e_fluid_total,
            AuthoritativeField::e_electron,
            AuthoritativeField::radiation_groups,
            AuthoritativeField::alpha_state}));
    DEC3D_CHECK_EQ(payload.serialized_cached_field_mask, static_cast<dec3d::core::CachedFieldMask>(0));
    DEC3D_CHECK_EQ(payload.rho(0, 0, 0), 1.25);
    DEC3D_CHECK_EQ(payload.e_electron(0, 0, 0), 2.5);
    DEC3D_CHECK_EQ(payload.radiation_groups[0](0, 0, 0), 0.75);
    DEC3D_CHECK_EQ(payload.alpha_state.storage(0, 0, 0), 0.5);

    const auto restart_validation = ValidateRestartFromCheckpoint(payload);
    DEC3D_CHECK(restart_validation.is_complete());
    DEC3D_CHECK(restart_validation.success);
    DEC3D_CHECK(restart_validation.payload_complete);
    DEC3D_CHECK(restart_validation.runtime_snapshot_complete);
    DEC3D_CHECK(restart_validation.cached_restart_truth_rejected);
    DEC3D_CHECK(restart_validation.cached_fields_invalid_after_restart);

    const auto restored_state = RestoreCanonicalStateFromCheckpoint(payload);
    DEC3D_CHECK_EQ(restored_state.rho(0, 0, 0), 1.25);
    DEC3D_CHECK_EQ(restored_state.e_electron(0, 0, 0), 2.5);
    DEC3D_CHECK_EQ(restored_state.radiation_groups[0](0, 0, 0), 0.75);
    DEC3D_CHECK_EQ(restored_state.alpha_state.storage(0, 0, 0), 0.5);
    DEC3D_CHECK(!restored_state.IsCachedFieldValid(CachedField::electron_temperature));
    DEC3D_CHECK(!restored_state.IsCachedFieldValid(CachedField::chi_e));

    const auto report = BuildCheckpointContractReport(payload);
    DEC3D_CHECK(report.is_complete());
    DEC3D_CHECK(report.success);
    DEC3D_CHECK(report.payload_complete);
    DEC3D_CHECK(report.restart_validation_succeeded);
    DEC3D_CHECK(report.runtime_snapshot_complete);
    DEC3D_CHECK(report.diagnostics_complete);
    DEC3D_CHECK_EQ(report.serialized_cached_field_mask, static_cast<dec3d::core::CachedFieldMask>(0));
    DEC3D_CHECK(report.diagnostics.has_entries());
    DEC3D_CHECK(report.failure_reason.empty());

    auto corrupted_payload = payload;
    corrupted_payload.serialized_cached_field_mask =
        dec3d::core::ToMask(CachedField::chi_e);

    const auto corrupted_validation = ValidateRestartFromCheckpoint(corrupted_payload);
    DEC3D_CHECK(corrupted_validation.is_complete());
    DEC3D_CHECK(!corrupted_validation.success);
    DEC3D_CHECK(!corrupted_validation.cached_restart_truth_rejected);
    DEC3D_CHECK(!corrupted_validation.failure_reason.empty());

    const auto corrupted_report = BuildCheckpointContractReport(corrupted_payload);
    DEC3D_CHECK(corrupted_report.is_complete());
    DEC3D_CHECK(!corrupted_report.success);
    DEC3D_CHECK(!corrupted_report.restart_validation_succeeded);
    DEC3D_CHECK(corrupted_report.diagnostics_complete);
    DEC3D_CHECK(!corrupted_report.failure_reason.empty());

    auto broken_runtime = CreateRuntimeScaffold(PhaseId::p0, "p0-v0.1");
    broken_runtime.diagnostics_sink.initialized = false;
    broken_runtime.diagnostics_sink.sink_id.clear();

    const auto broken_runtime_payload = BuildCheckpointPayload(state, broken_runtime);
    DEC3D_CHECK(!broken_runtime_payload.is_complete());

    const auto broken_runtime_validation = ValidateRestartFromCheckpoint(broken_runtime_payload);
    DEC3D_CHECK(broken_runtime_validation.is_complete());
    DEC3D_CHECK(!broken_runtime_validation.success);
    DEC3D_CHECK(!broken_runtime_validation.payload_complete);
    DEC3D_CHECK(!broken_runtime_validation.runtime_snapshot_complete);
    DEC3D_CHECK(!broken_runtime_validation.failure_reason.empty());

    const auto broken_runtime_report = BuildCheckpointContractReport(broken_runtime_payload);
    DEC3D_CHECK(broken_runtime_report.is_complete());
    DEC3D_CHECK(!broken_runtime_report.success);
    DEC3D_CHECK(!broken_runtime_report.payload_complete);
    DEC3D_CHECK(!broken_runtime_report.runtime_snapshot_complete);
    DEC3D_CHECK_EQ(
        broken_runtime_report.restart_invalidated_cached_mask,
        dec3d::core::ToMask(CachedField::electron_temperature) |
            dec3d::core::ToMask(CachedField::ion_temperature) |
            dec3d::core::ToMask(CachedField::electron_pressure) |
            dec3d::core::ToMask(CachedField::ion_pressure) |
            dec3d::core::ToMask(CachedField::chi_e) |
            dec3d::core::ToMask(CachedField::conductivity) |
            dec3d::core::ToMask(CachedField::opacity) |
            dec3d::core::ToMask(CachedField::sound_speed) |
            dec3d::core::ToMask(CachedField::velocity));
    DEC3D_CHECK(!broken_runtime_report.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
