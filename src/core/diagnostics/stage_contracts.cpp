#include "core/diagnostics/stage_contracts.hpp"

namespace dec3d::core {

AuthoritativeFieldMask Combine(std::initializer_list<AuthoritativeField> fields) noexcept {
  AuthoritativeFieldMask mask = 0;

  for (const auto field : fields) {
    mask |= ToMask(field);
  }

  return mask;
}

CachedFieldMask InvalidatedCachesForWriteMask(AuthoritativeFieldMask updated_fields) noexcept {
  constexpr AuthoritativeFieldMask kStateThermodynamicWrites =
      ToMask(AuthoritativeField::rho) |
      ToMask(AuthoritativeField::e_fluid_total) |
      ToMask(AuthoritativeField::e_electron);

  constexpr AuthoritativeFieldMask kMomentumWrites =
      ToMask(AuthoritativeField::mom_r) |
      ToMask(AuthoritativeField::mom_theta) |
      ToMask(AuthoritativeField::mom_phi);

  constexpr CachedFieldMask kThermodynamicCaches =
      ToMask(CachedField::electron_temperature) |
      ToMask(CachedField::ion_temperature) |
      ToMask(CachedField::electron_pressure) |
      ToMask(CachedField::ion_pressure) |
      ToMask(CachedField::chi_e) |
      ToMask(CachedField::conductivity) |
      ToMask(CachedField::opacity) |
      ToMask(CachedField::sound_speed);

  constexpr CachedFieldMask kVelocityCaches =
      ToMask(CachedField::velocity) |
      ToMask(CachedField::sound_speed);

  CachedFieldMask invalidated = ToMask(CachedField::none);

  if ((updated_fields & kStateThermodynamicWrites) != 0) {
    invalidated |= kThermodynamicCaches;
  }

  if ((updated_fields & kMomentumWrites) != 0) {
    invalidated |= kThermodynamicCaches;
    invalidated |= kVelocityCaches;
  }

  if ((updated_fields & ToMask(AuthoritativeField::alpha_state)) != 0) {
    invalidated |= kThermodynamicCaches;
  }

  return invalidated;
}

bool ExecutionEvidence::is_present() const noexcept {
  return entered_stage && !implementation_id.empty();
}

bool MeshUpdateProposal::is_complete(std::size_t expected_radial_face_count) const noexcept {
  return requested &&
         radial_ale &&
         dt_s > 0.0 &&
         expected_radial_face_count > 0u &&
         radial_face_velocities.size() == expected_radial_face_count &&
         proposed_radial_faces.size() == expected_radial_face_count &&
         !implementation_id.empty() &&
         !summary.empty();
}

bool MeshUpdateProposal::is_complete_global(
    std::size_t expected_global_radial_face_count) const noexcept {
  return requested &&
         radial_ale &&
         radial_face_indexing_is_global &&
         dt_s > 0.0 &&
         expected_global_radial_face_count > 0u &&
         global_radial_face_count == expected_global_radial_face_count &&
         radial_face_velocities.size() == expected_global_radial_face_count &&
         proposed_radial_faces.size() == expected_global_radial_face_count &&
         !implementation_id.empty() &&
         !summary.empty();
}

bool MeshUpdateProposal::has_radial_face_window(
    std::size_t global_face_begin,
    std::size_t local_face_count) const noexcept {
  if (!is_complete_global(global_radial_face_count) || local_face_count == 0u) {
    return false;
  }

  if (global_face_begin > global_radial_face_count) {
    return false;
  }

  return local_face_count <= global_radial_face_count - global_face_begin;
}

bool StageContext::is_complete() const noexcept {
  return dt_s >= 0.0 &&
         !contract_version.empty() &&
         !mesh_snapshot_handle.empty() &&
         !ownership_handle.empty() &&
         !diagnostics_sink_handle.empty();
}

bool DtAdvice::is_complete() const noexcept {
  return hard_cap_dt > 0.0 && soft_advice_dt > 0.0 &&
         !reason.empty() &&
         !evidence.empty();
}

bool StageResult::is_semantically_complete() const noexcept {
  return success &&
         updated_fields != 0 &&
         diagnostics.has_entries() &&
         execution_evidence.has_value() &&
         execution_evidence->is_present();
}

StageResult StageResult::Successful(
    AuthoritativeFieldMask updated_fields,
    DiagnosticsPayload diagnostics,
    ExecutionEvidence execution_evidence) {
  StageResult result;
  result.success = true;
  result.updated_fields = updated_fields;
  result.diagnostics = std::move(diagnostics);
  result.execution_evidence = std::move(execution_evidence);
  return result;
}

StageResult StageResult::Failed(
    std::string failure_reason,
    DiagnosticsPayload diagnostics,
    ExecutionEvidence execution_evidence) {
  StageResult result;
  result.success = false;
  result.updated_fields = 0;
  result.failure_reason = std::move(failure_reason);
  result.diagnostics = std::move(diagnostics);
  result.execution_evidence = std::move(execution_evidence);
  return result;
}

}  // namespace dec3d::core
