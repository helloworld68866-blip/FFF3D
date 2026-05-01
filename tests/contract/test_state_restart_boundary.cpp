#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/canonical_state/state_utilities.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::InvalidateRecoveredAndCachedFieldsForRestart;
    using dec3d::state::MakeAuthoritativeSnapshotView;
    using dec3d::state::RestartSerializationPolicy;
    using dec3d::state::SerializationPolicyForAuthoritativeField;
    using dec3d::state::SerializationPolicyForCachedField;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    state.MarkCachedFieldValid(CachedField::electron_temperature);
    state.MarkCachedFieldValid(CachedField::chi_e);

    DEC3D_CHECK_EQ(
        static_cast<int>(SerializationPolicyForAuthoritativeField(AuthoritativeField::rho)),
        static_cast<int>(RestartSerializationPolicy::authoritative_restart_truth));
    DEC3D_CHECK_EQ(
        static_cast<int>(SerializationPolicyForCachedField(CachedField::chi_e)),
        static_cast<int>(RestartSerializationPolicy::invalidate_on_restart));

    const auto snapshot = MakeAuthoritativeSnapshotView(state);
    DEC3D_CHECK(snapshot.authoritative_only);
    DEC3D_CHECK(snapshot.is_complete());
    DEC3D_CHECK(snapshot.rho != nullptr);
    DEC3D_CHECK(snapshot.e_electron != nullptr);
    DEC3D_CHECK(snapshot.radiation_groups != nullptr);
    DEC3D_CHECK(snapshot.alpha_state != nullptr);

    InvalidateRecoveredAndCachedFieldsForRestart(state);
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::electron_temperature));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::chi_e));

    auto zero_radiation_state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 0});
    const auto zero_radiation_snapshot = MakeAuthoritativeSnapshotView(zero_radiation_state);
    DEC3D_CHECK(zero_radiation_snapshot.is_complete());
    DEC3D_CHECK(zero_radiation_snapshot.radiation_groups != nullptr);
    DEC3D_CHECK(zero_radiation_snapshot.radiation_groups->empty());

    auto mismatched_state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    mismatched_state.alpha_state.storage = dec3d::core::Array3D<double>(1, 2, 2, 0.0);
    const auto mismatched_snapshot = MakeAuthoritativeSnapshotView(mismatched_state);
    DEC3D_CHECK(!mismatched_snapshot.is_complete());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
