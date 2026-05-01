#pragma once

#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::state {

struct CanonicalStateLayout {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
  std::size_t radiation_group_count{0};

  [[nodiscard]] bool is_valid() const noexcept {
    return radial_cells > 0 && theta_cells > 0 && phi_cells > 0;
  }
};

struct AuthoritativeAlphaState {
  dec3d::core::Array3D<double> storage;

  [[nodiscard]] bool has_storage() const noexcept { return !storage.empty(); }
};

class CachedFieldRegistry {
 public:
  CachedFieldRegistry() = default;

  [[nodiscard]] bool HasMetadata(dec3d::core::CachedField field) const noexcept;
  [[nodiscard]] bool HasDedicatedStorage(dec3d::core::CachedField field) const noexcept;
  [[nodiscard]] bool IsValid(dec3d::core::CachedField field) const noexcept;

  void MarkValid(dec3d::core::CachedField field) noexcept;
  void Invalidate(dec3d::core::CachedFieldMask mask) noexcept;

  [[nodiscard]] dec3d::core::CachedFieldMask supported_mask() const noexcept { return supported_mask_; }
  [[nodiscard]] dec3d::core::CachedFieldMask valid_mask() const noexcept { return valid_mask_; }

 private:
  static constexpr dec3d::core::CachedFieldMask supported_mask_ =
      dec3d::core::ToMask(dec3d::core::CachedField::electron_temperature) |
      dec3d::core::ToMask(dec3d::core::CachedField::ion_temperature) |
      dec3d::core::ToMask(dec3d::core::CachedField::electron_pressure) |
      dec3d::core::ToMask(dec3d::core::CachedField::ion_pressure) |
      dec3d::core::ToMask(dec3d::core::CachedField::chi_e) |
      dec3d::core::ToMask(dec3d::core::CachedField::conductivity) |
      dec3d::core::ToMask(dec3d::core::CachedField::opacity) |
      dec3d::core::ToMask(dec3d::core::CachedField::sound_speed) |
      dec3d::core::ToMask(dec3d::core::CachedField::velocity);

  dec3d::core::CachedFieldMask valid_mask_{0};
};

struct CanonicalState {
  CanonicalStateLayout layout;
  dec3d::core::Array3D<double> rho;
  dec3d::core::Array3D<double> mom_r;
  dec3d::core::Array3D<double> mom_theta;
  dec3d::core::Array3D<double> mom_phi;
  dec3d::core::Array3D<double> e_fluid_total;
  dec3d::core::Array3D<double> e_electron;
  std::vector<dec3d::core::Array3D<double>> radiation_groups;
  AuthoritativeAlphaState alpha_state;
  CachedFieldRegistry cached_fields;
  dec3d::core::AuthoritativeFieldMask last_authoritative_write_mask{0};
  dec3d::core::CachedFieldMask last_invalidated_cached_mask{0};

  [[nodiscard]] static CanonicalState Create(const CanonicalStateLayout& layout);

  [[nodiscard]] bool HasAuthoritativeStorage(dec3d::core::AuthoritativeField field) const noexcept;
  [[nodiscard]] bool HasCachedFieldMetadata(dec3d::core::CachedField field) const noexcept;
  [[nodiscard]] bool HasCachedFieldStorage(dec3d::core::CachedField field) const noexcept;
  [[nodiscard]] bool IsCachedFieldValid(dec3d::core::CachedField field) const noexcept;

  void MarkCachedFieldValid(dec3d::core::CachedField field) noexcept;
  void ApplyAuthoritativeWrite(dec3d::core::AuthoritativeFieldMask updated_fields) noexcept;

  [[nodiscard]] bool TryCommitCachedField(dec3d::core::CachedField field) noexcept;
};

struct ThermodynamicFeasibilityThresholds {
  double electron_min{0.0};
  double ion_min{0.0};
};

struct ThermodynamicFeasibilityResult {
  bool success{false};
  double minimum_electron_energy{0.0};
  double minimum_ion_energy{0.0};
  std::string failure_reason;
};

[[nodiscard]] ThermodynamicFeasibilityResult ValidateThermodynamicFeasibility(
    const CanonicalState& state,
    const ThermodynamicFeasibilityThresholds& thresholds) noexcept;

}  // namespace dec3d::state
