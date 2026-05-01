#pragma once

#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <string>

namespace dec3d::state {

enum class RestartSerializationPolicy {
  authoritative_restart_truth,
  invalidate_on_restart,
  never_serialize,
};

struct CanonicalStateLayoutValidationResult {
  bool success{false};
  std::string failure_reason;
};

struct CanonicalStateSummary {
  bool allocated{false};
  std::size_t total_cell_count{0};
  std::size_t radiation_group_count{0};
  dec3d::core::AuthoritativeFieldMask authoritative_storage_mask{0};
  dec3d::core::AuthoritativeFieldMask restart_authoritative_mask{0};
  dec3d::core::CachedFieldMask cached_supported_mask{0};
  dec3d::core::CachedFieldMask cached_valid_mask{0};
  dec3d::core::CachedFieldMask restart_cached_policy_mask{0};
  dec3d::core::AuthoritativeFieldMask last_authoritative_write_mask{0};
  dec3d::core::CachedFieldMask last_invalidated_cached_mask{0};
  std::string report_line;
};

struct FieldUpdateSummary {
  dec3d::core::AuthoritativeFieldMask mask{0};
  std::vector<std::string> field_names;
  std::string report_line;
};

struct AuthoritativeSnapshotView {
  const dec3d::core::Array3D<double>* rho{nullptr};
  const dec3d::core::Array3D<double>* mom_r{nullptr};
  const dec3d::core::Array3D<double>* mom_theta{nullptr};
  const dec3d::core::Array3D<double>* mom_phi{nullptr};
  const dec3d::core::Array3D<double>* e_fluid_total{nullptr};
  const dec3d::core::Array3D<double>* e_electron{nullptr};
  const std::vector<dec3d::core::Array3D<double>>* radiation_groups{nullptr};
  const AuthoritativeAlphaState* alpha_state{nullptr};
  bool authoritative_only{true};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] CanonicalStateLayoutValidationResult ValidateCanonicalStateLayout(
    const CanonicalStateLayout& layout) noexcept;

[[nodiscard]] dec3d::core::CachedFieldMask QueryInvalidatedCachedFields(
    const CanonicalState& state,
    dec3d::core::AuthoritativeFieldMask updated_fields) noexcept;

[[nodiscard]] CanonicalStateSummary BuildCanonicalStateSummary(
    const CanonicalState& state) noexcept;

[[nodiscard]] FieldUpdateSummary BuildFieldUpdateSummary(
    dec3d::core::AuthoritativeFieldMask mask) noexcept;

[[nodiscard]] dec3d::core::CachedFieldMask RestartInvalidatedCachedMask(
    const CanonicalState& state) noexcept;

[[nodiscard]] RestartSerializationPolicy SerializationPolicyForAuthoritativeField(
    dec3d::core::AuthoritativeField field) noexcept;

[[nodiscard]] RestartSerializationPolicy SerializationPolicyForCachedField(
    dec3d::core::CachedField field) noexcept;

void InvalidateRecoveredAndCachedFieldsForRestart(CanonicalState& state) noexcept;

[[nodiscard]] AuthoritativeSnapshotView MakeAuthoritativeSnapshotView(
    const CanonicalState& state) noexcept;

}  // namespace dec3d::state
