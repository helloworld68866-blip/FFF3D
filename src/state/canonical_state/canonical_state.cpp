#include "state/canonical_state/canonical_state.hpp"

#include <algorithm>
#include <limits>
#include <sstream>

namespace dec3d::state {

bool CachedFieldRegistry::HasMetadata(dec3d::core::CachedField field) const noexcept {
  return dec3d::core::MaskContains(supported_mask_, field);
}

bool CachedFieldRegistry::HasDedicatedStorage(dec3d::core::CachedField field) const noexcept {
  (void)field;
  return false;
}

bool CachedFieldRegistry::IsValid(dec3d::core::CachedField field) const noexcept {
  return dec3d::core::MaskContains(valid_mask_, field);
}

void CachedFieldRegistry::MarkValid(dec3d::core::CachedField field) noexcept {
  if (HasMetadata(field)) {
    valid_mask_ |= dec3d::core::ToMask(field);
  }
}

void CachedFieldRegistry::Invalidate(dec3d::core::CachedFieldMask mask) noexcept {
  valid_mask_ &= ~mask;
}

CanonicalState CanonicalState::Create(const CanonicalStateLayout& layout) {
  CanonicalState state;
  state.layout = layout;

  if (!layout.is_valid()) {
    return state;
  }

  state.rho = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  state.mom_r = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  state.mom_theta = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  state.mom_phi = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  state.e_fluid_total = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  state.e_electron = dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);

  state.radiation_groups.reserve(layout.radiation_group_count);
  for (std::size_t group_index = 0; group_index < layout.radiation_group_count; ++group_index) {
    state.radiation_groups.emplace_back(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  }

  state.alpha_state.storage =
      dec3d::core::Array3D<double>(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);

  return state;
}

bool CanonicalState::HasAuthoritativeStorage(dec3d::core::AuthoritativeField field) const noexcept {
  switch (field) {
    case dec3d::core::AuthoritativeField::rho:
      return !rho.empty();
    case dec3d::core::AuthoritativeField::mom_r:
      return !mom_r.empty();
    case dec3d::core::AuthoritativeField::mom_theta:
      return !mom_theta.empty();
    case dec3d::core::AuthoritativeField::mom_phi:
      return !mom_phi.empty();
    case dec3d::core::AuthoritativeField::e_fluid_total:
      return !e_fluid_total.empty();
    case dec3d::core::AuthoritativeField::e_electron:
      return !e_electron.empty();
    case dec3d::core::AuthoritativeField::radiation_groups:
      return !radiation_groups.empty();
    case dec3d::core::AuthoritativeField::alpha_state:
      return alpha_state.has_storage();
    case dec3d::core::AuthoritativeField::none:
      return false;
  }

  return false;
}

bool CanonicalState::HasCachedFieldMetadata(dec3d::core::CachedField field) const noexcept {
  return cached_fields.HasMetadata(field);
}

bool CanonicalState::HasCachedFieldStorage(dec3d::core::CachedField field) const noexcept {
  return cached_fields.HasDedicatedStorage(field);
}

bool CanonicalState::IsCachedFieldValid(dec3d::core::CachedField field) const noexcept {
  return cached_fields.IsValid(field);
}

void CanonicalState::MarkCachedFieldValid(dec3d::core::CachedField field) noexcept {
  cached_fields.MarkValid(field);
}

void CanonicalState::ApplyAuthoritativeWrite(dec3d::core::AuthoritativeFieldMask updated_fields) noexcept {
  last_authoritative_write_mask = updated_fields;
  last_invalidated_cached_mask = dec3d::core::InvalidatedCachesForWriteMask(updated_fields);
  cached_fields.Invalidate(last_invalidated_cached_mask);
}

bool CanonicalState::TryCommitCachedField(dec3d::core::CachedField field) noexcept {
  cached_fields.Invalidate(dec3d::core::ToMask(field));
  return false;
}

ThermodynamicFeasibilityResult ValidateThermodynamicFeasibility(
    const CanonicalState& state,
    const ThermodynamicFeasibilityThresholds& thresholds) noexcept {
  ThermodynamicFeasibilityResult result;
  result.success = true;
  result.minimum_electron_energy = std::numeric_limits<double>::infinity();
  result.minimum_ion_energy = std::numeric_limits<double>::infinity();

  if (state.rho.empty()) {
    result.success = false;
    result.minimum_electron_energy = 0.0;
    result.minimum_ion_energy = 0.0;
    result.failure_reason = "canonical state is not allocated";
    return result;
  }

  const auto& rho_storage = state.rho.storage();
  const auto& mom_r_storage = state.mom_r.storage();
  const auto& mom_theta_storage = state.mom_theta.storage();
  const auto& mom_phi_storage = state.mom_phi.storage();
  const auto& fluid_storage = state.e_fluid_total.storage();
  const auto& electron_storage = state.e_electron.storage();

  for (std::size_t index = 0; index < rho_storage.size(); ++index) {
    const double rho = rho_storage[index];
    if (rho <= 0.0) {
      result.success = false;
      result.minimum_electron_energy = 0.0;
      result.minimum_ion_energy = 0.0;
      std::ostringstream reason;
      reason << "rho must remain positive"
             << "; failing_linear_index=" << index
             << "; rho=" << rho;
      result.failure_reason = reason.str();
      return result;
    }

    const double kinetic =
        0.5 * ((mom_r_storage[index] * mom_r_storage[index]) +
               (mom_theta_storage[index] * mom_theta_storage[index]) +
               (mom_phi_storage[index] * mom_phi_storage[index])) /
        rho;
    const double electron = electron_storage[index];
    const double ion = fluid_storage[index] - kinetic - electron;

    result.minimum_electron_energy = std::min(result.minimum_electron_energy, electron);
    result.minimum_ion_energy = std::min(result.minimum_ion_energy, ion);

    if (electron < thresholds.electron_min) {
      result.success = false;
      std::ostringstream reason;
      reason << "electron energy fell below minimum threshold"
             << "; failing_linear_index=" << index
             << "; electron_energy=" << electron
             << "; electron_min=" << thresholds.electron_min
             << "; rho=" << rho
             << "; kinetic=" << kinetic
             << "; fluid_energy=" << fluid_storage[index]
             << "; ion_energy=" << ion
             << "; mom_r=" << mom_r_storage[index]
             << "; mom_theta=" << mom_theta_storage[index]
             << "; mom_phi=" << mom_phi_storage[index];
      result.failure_reason = reason.str();
      return result;
    }

    if (ion < thresholds.ion_min) {
      result.success = false;
      std::ostringstream reason;
      reason << "ion energy fell below minimum threshold"
             << "; failing_linear_index=" << index
             << "; ion_energy=" << ion
             << "; ion_min=" << thresholds.ion_min
             << "; rho=" << rho
             << "; kinetic=" << kinetic
             << "; fluid_energy=" << fluid_storage[index]
             << "; electron_energy=" << electron
             << "; mom_r=" << mom_r_storage[index]
             << "; mom_theta=" << mom_theta_storage[index]
             << "; mom_phi=" << mom_phi_storage[index];
      result.failure_reason = reason.str();
      return result;
    }
  }

  return result;
}

}  // namespace dec3d::state
