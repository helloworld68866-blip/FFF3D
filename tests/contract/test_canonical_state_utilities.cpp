#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/canonical_state/state_utilities.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::core::Combine;
    using dec3d::core::MaskContains;
    using dec3d::state::BuildCanonicalStateSummary;
    using dec3d::state::BuildFieldUpdateSummary;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::QueryInvalidatedCachedFields;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    state.MarkCachedFieldValid(CachedField::electron_temperature);
    state.MarkCachedFieldValid(CachedField::velocity);
    const auto write_mask = Combine({AuthoritativeField::e_electron, AuthoritativeField::mom_r});
    state.ApplyAuthoritativeWrite(write_mask);

    const auto invalidated = QueryInvalidatedCachedFields(
        state,
        write_mask);
    DEC3D_CHECK(MaskContains(invalidated, CachedField::electron_temperature));
    DEC3D_CHECK(MaskContains(invalidated, CachedField::chi_e));
    DEC3D_CHECK(MaskContains(invalidated, CachedField::velocity));

    const auto summary = BuildCanonicalStateSummary(state);
    DEC3D_CHECK(summary.allocated);
    DEC3D_CHECK_EQ(summary.total_cell_count, static_cast<std::size_t>(8));
    DEC3D_CHECK_EQ(summary.radiation_group_count, static_cast<std::size_t>(1));
    DEC3D_CHECK(MaskContains(summary.authoritative_storage_mask, AuthoritativeField::rho));
    DEC3D_CHECK(MaskContains(summary.authoritative_storage_mask, AuthoritativeField::radiation_groups));
    DEC3D_CHECK(MaskContains(summary.restart_authoritative_mask, AuthoritativeField::rho));
    DEC3D_CHECK(MaskContains(summary.cached_supported_mask, CachedField::electron_temperature));
    DEC3D_CHECK(MaskContains(summary.restart_cached_policy_mask, CachedField::chi_e));
    DEC3D_CHECK_EQ(summary.last_authoritative_write_mask, write_mask);
    DEC3D_CHECK_EQ(summary.last_invalidated_cached_mask, invalidated);
    DEC3D_CHECK(!summary.report_line.empty());
    DEC3D_CHECK(summary.report_line.find("last_write_mask=") != std::string::npos);

    const auto field_update_summary = BuildFieldUpdateSummary(write_mask);
    DEC3D_CHECK_EQ(field_update_summary.mask, write_mask);
    DEC3D_CHECK(field_update_summary.field_names.size() >= static_cast<std::size_t>(2));
    DEC3D_CHECK(field_update_summary.report_line.find("e_electron") != std::string::npos);
    DEC3D_CHECK(field_update_summary.report_line.find("mom_r") != std::string::npos);

    const auto empty_update_summary = BuildFieldUpdateSummary(0);
    DEC3D_CHECK_EQ(empty_update_summary.mask, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(empty_update_summary.field_names.empty());
    DEC3D_CHECK(empty_update_summary.report_line.find("none") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
