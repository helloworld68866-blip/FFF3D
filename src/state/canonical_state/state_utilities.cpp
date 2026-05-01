#include "state/canonical_state/state_utilities.hpp"

#include <sstream>

namespace dec3d::state {

bool AuthoritativeSnapshotView::is_complete() const noexcept {
  if (!(authoritative_only &&
        rho != nullptr &&
        mom_r != nullptr &&
        mom_theta != nullptr &&
        mom_phi != nullptr &&
        e_fluid_total != nullptr &&
        e_electron != nullptr &&
        radiation_groups != nullptr &&
        alpha_state != nullptr)) {
    return false;
  }

  const auto matches_shape = [this](const dec3d::core::Array3D<double>& field) {
    return field.extent_r() == rho->extent_r() &&
           field.extent_theta() == rho->extent_theta() &&
           field.extent_phi() == rho->extent_phi();
  };

  if (!(matches_shape(*mom_r) &&
        matches_shape(*mom_theta) &&
        matches_shape(*mom_phi) &&
        matches_shape(*e_fluid_total) &&
        matches_shape(*e_electron) &&
        matches_shape(alpha_state->storage))) {
    return false;
  }

  for (const auto& radiation_group : *radiation_groups) {
    if (!matches_shape(radiation_group)) {
      return false;
    }
  }

  return true;
}

CanonicalStateLayoutValidationResult ValidateCanonicalStateLayout(
    const CanonicalStateLayout& layout) noexcept {
  CanonicalStateLayoutValidationResult result;
  result.success = layout.is_valid();

  if (result.success) {
    return result;
  }

  if (layout.radial_cells == 0) {
    result.failure_reason = "radial_cells must be positive";
    return result;
  }

  if (layout.theta_cells == 0) {
    result.failure_reason = "theta_cells must be positive";
    return result;
  }

  if (layout.phi_cells == 0) {
    result.failure_reason = "phi_cells must be positive";
    return result;
  }

  result.failure_reason = "layout is invalid";
  return result;
}

dec3d::core::CachedFieldMask QueryInvalidatedCachedFields(
    const CanonicalState& state,
    dec3d::core::AuthoritativeFieldMask updated_fields) noexcept {
  if (!state.layout.is_valid()) {
    return dec3d::core::ToMask(dec3d::core::CachedField::none);
  }

  return dec3d::core::InvalidatedCachesForWriteMask(updated_fields);
}

CanonicalStateSummary BuildCanonicalStateSummary(const CanonicalState& state) noexcept {
  CanonicalStateSummary summary;
  summary.allocated = state.layout.is_valid() && !state.rho.empty();
  summary.total_cell_count = state.rho.size();
  summary.radiation_group_count = state.radiation_groups.size();
  summary.cached_supported_mask = state.cached_fields.supported_mask();
  summary.cached_valid_mask = state.cached_fields.valid_mask();
  summary.last_authoritative_write_mask = state.last_authoritative_write_mask;
  summary.last_invalidated_cached_mask = state.last_invalidated_cached_mask;

  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::rho)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::rho);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_r)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_r);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_theta)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_theta);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_phi)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::mom_phi);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_fluid_total)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_electron)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::radiation_groups)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
  }
  if (state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state)) {
    summary.authoritative_storage_mask |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
  }
  summary.restart_authoritative_mask = summary.authoritative_storage_mask;
  summary.restart_cached_policy_mask = state.cached_fields.supported_mask();

  std::ostringstream report;
  report << "cells=" << summary.total_cell_count
         << "; radiation_groups=" << summary.radiation_group_count
         << "; authoritative_mask=" << summary.authoritative_storage_mask
         << "; restart_authoritative_mask=" << summary.restart_authoritative_mask
         << "; cached_valid_mask=" << summary.cached_valid_mask
         << "; restart_cached_policy_mask=" << summary.restart_cached_policy_mask
         << "; last_write_mask=" << summary.last_authoritative_write_mask
         << "; last_invalidated_cached_mask=" << summary.last_invalidated_cached_mask;
  summary.report_line = report.str();

  return summary;
}

FieldUpdateSummary BuildFieldUpdateSummary(dec3d::core::AuthoritativeFieldMask mask) noexcept {
  FieldUpdateSummary summary;
  summary.mask = mask;
  auto append_field = [&summary, mask](dec3d::core::AuthoritativeField field, const char* name) {
    if (dec3d::core::MaskContains(mask, field)) {
      summary.field_names.emplace_back(name);
    }
  };

  append_field(dec3d::core::AuthoritativeField::rho, "rho");
  append_field(dec3d::core::AuthoritativeField::mom_r, "mom_r");
  append_field(dec3d::core::AuthoritativeField::mom_theta, "mom_theta");
  append_field(dec3d::core::AuthoritativeField::mom_phi, "mom_phi");
  append_field(dec3d::core::AuthoritativeField::e_fluid_total, "e_fluid_total");
  append_field(dec3d::core::AuthoritativeField::e_electron, "e_electron");
  append_field(dec3d::core::AuthoritativeField::radiation_groups, "radiation_groups");
  append_field(dec3d::core::AuthoritativeField::alpha_state, "alpha_state");

  std::ostringstream report;
  report << "updated_fields=";
  if (summary.field_names.empty()) {
    report << "none";
  } else {
    for (std::size_t index = 0; index < summary.field_names.size(); ++index) {
      if (index != 0) {
        report << ",";
      }
      report << summary.field_names[index];
    }
  }
  summary.report_line = report.str();
  return summary;
}

dec3d::core::CachedFieldMask RestartInvalidatedCachedMask(
    const CanonicalState& state) noexcept {
  return state.cached_fields.supported_mask();
}

RestartSerializationPolicy SerializationPolicyForAuthoritativeField(
    dec3d::core::AuthoritativeField field) noexcept {
  if (field == dec3d::core::AuthoritativeField::none) {
    return RestartSerializationPolicy::never_serialize;
  }

  return RestartSerializationPolicy::authoritative_restart_truth;
}

RestartSerializationPolicy SerializationPolicyForCachedField(
    dec3d::core::CachedField field) noexcept {
  if (field == dec3d::core::CachedField::none) {
    return RestartSerializationPolicy::never_serialize;
  }

  return RestartSerializationPolicy::invalidate_on_restart;
}

void InvalidateRecoveredAndCachedFieldsForRestart(CanonicalState& state) noexcept {
  state.cached_fields.Invalidate(state.cached_fields.supported_mask());
}

AuthoritativeSnapshotView MakeAuthoritativeSnapshotView(const CanonicalState& state) noexcept {
  AuthoritativeSnapshotView view;
  view.rho = &state.rho;
  view.mom_r = &state.mom_r;
  view.mom_theta = &state.mom_theta;
  view.mom_phi = &state.mom_phi;
  view.e_fluid_total = &state.e_fluid_total;
  view.e_electron = &state.e_electron;
  view.radiation_groups = &state.radiation_groups;
  view.alpha_state = &state.alpha_state;
  view.authoritative_only = true;
  return view;
}

}  // namespace dec3d::state
