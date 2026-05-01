#pragma once

#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace dec3d::core {

enum class PhaseId : std::uint8_t {
  p0,
  p1,
  p2,
  p3,
  p4,
  p5,
};

enum class StageId : std::uint8_t {
  hydro,
  thermal,
  equilibration,
  radiation,
  alpha,
};

enum class AuthoritativeField : std::uint64_t {
  none = 0,
  rho = 1ull << 0,
  mom_r = 1ull << 1,
  mom_theta = 1ull << 2,
  mom_phi = 1ull << 3,
  e_fluid_total = 1ull << 4,
  e_electron = 1ull << 5,
  radiation_groups = 1ull << 6,
  alpha_state = 1ull << 7,
};

enum class CachedField : std::uint32_t {
  none = 0,
  electron_temperature = 1u << 0,
  ion_temperature = 1u << 1,
  electron_pressure = 1u << 2,
  ion_pressure = 1u << 3,
  chi_e = 1u << 4,
  conductivity = 1u << 5,
  opacity = 1u << 6,
  sound_speed = 1u << 7,
  velocity = 1u << 8,
};

using AuthoritativeFieldMask = std::uint64_t;
using CachedFieldMask = std::uint32_t;

[[nodiscard]] constexpr AuthoritativeFieldMask ToMask(AuthoritativeField field) noexcept {
  return static_cast<AuthoritativeFieldMask>(field);
}

[[nodiscard]] constexpr CachedFieldMask ToMask(CachedField field) noexcept {
  return static_cast<CachedFieldMask>(field);
}

[[nodiscard]] constexpr bool MaskContains(AuthoritativeFieldMask mask, AuthoritativeField field) noexcept {
  return (mask & ToMask(field)) != 0;
}

[[nodiscard]] constexpr bool MaskContains(CachedFieldMask mask, CachedField field) noexcept {
  return (mask & ToMask(field)) != 0;
}

[[nodiscard]] AuthoritativeFieldMask Combine(std::initializer_list<AuthoritativeField> fields) noexcept;
[[nodiscard]] CachedFieldMask InvalidatedCachesForWriteMask(AuthoritativeFieldMask updated_fields) noexcept;

struct DiagnosticMessage {
  std::string code;
  std::string message;
};

struct DiagnosticsPayload {
  std::vector<DiagnosticMessage> entries;

  [[nodiscard]] bool has_entries() const noexcept { return !entries.empty(); }
};

struct ExecutionEvidence {
  bool entered_stage{false};
  std::string implementation_id;
  std::uint64_t touched_cell_count{0};

  [[nodiscard]] bool is_present() const noexcept;
};

struct MeshUpdateProposal {
  bool requested{false};
  bool radial_ale{false};
  bool radial_face_indexing_is_global{false};
  std::size_t global_radial_face_count{0u};
  double dt_s{0.0};
  std::vector<double> radial_face_velocities;
  std::vector<double> proposed_radial_faces;
  std::string implementation_id;
  std::string summary;

  [[nodiscard]] bool is_complete(std::size_t expected_radial_face_count) const noexcept;
  [[nodiscard]] bool is_complete_global(std::size_t expected_global_radial_face_count) const noexcept;
  [[nodiscard]] bool has_radial_face_window(
      std::size_t global_face_begin,
      std::size_t local_face_count) const noexcept;
};

struct StageContext {
  double time_s{0.0};
  double dt_s{0.0};
  std::uint64_t step{0};
  PhaseId phase_id{PhaseId::p0};
  std::string contract_version;
  std::string mesh_snapshot_handle;
  std::string ownership_handle;
  std::string diagnostics_sink_handle;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct DtAdvice {
  double hard_cap_dt{std::numeric_limits<double>::infinity()};
  double soft_advice_dt{std::numeric_limits<double>::infinity()};
  std::string reason;
  std::string evidence;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct StageResult {
  bool success{false};
  AuthoritativeFieldMask updated_fields{0};
  std::string failure_reason;
  DiagnosticsPayload diagnostics;
  std::optional<ExecutionEvidence> execution_evidence;
  std::optional<MeshUpdateProposal> mesh_update_proposal;

  [[nodiscard]] bool is_semantically_complete() const noexcept;

  static StageResult Successful(
      AuthoritativeFieldMask updated_fields,
      DiagnosticsPayload diagnostics,
      ExecutionEvidence execution_evidence);

  static StageResult Failed(
      std::string failure_reason,
      DiagnosticsPayload diagnostics,
      ExecutionEvidence execution_evidence);
};

}  // namespace dec3d::core
