#pragma once

#include "runtime/runtime_scaffold.hpp"
#include "state/canonical_state/state_utilities.hpp"

#include <string>

namespace dec3d::runtime {

struct CheckpointPayload {
  bool serialized{false};
  dec3d::core::PhaseId phase_id{dec3d::core::PhaseId::p0};
  std::string contract_version;
  dec3d::state::CanonicalStateLayout layout;
  dec3d::core::AuthoritativeFieldMask authoritative_field_mask{0};
  dec3d::core::CachedFieldMask serialized_cached_field_mask{0};
  dec3d::core::Array3D<double> rho;
  dec3d::core::Array3D<double> mom_r;
  dec3d::core::Array3D<double> mom_theta;
  dec3d::core::Array3D<double> mom_phi;
  dec3d::core::Array3D<double> e_fluid_total;
  dec3d::core::Array3D<double> e_electron;
  std::vector<dec3d::core::Array3D<double>> radiation_groups;
  dec3d::state::AuthoritativeAlphaState alpha_state;
  RuntimeAcceptanceSnapshot runtime_snapshot;
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RestartValidationResult {
  bool validation_performed{false};
  bool success{false};
  bool payload_complete{false};
  bool runtime_snapshot_complete{false};
  bool authoritative_state_restored{false};
  bool authoritative_mask_matches_payload{false};
  bool cached_restart_truth_rejected{false};
  bool cached_fields_invalid_after_restart{false};
  dec3d::core::CachedFieldMask restart_invalidated_cached_mask{0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct CheckpointContractReport {
  bool built{false};
  bool success{false};
  bool payload_complete{false};
  bool restart_validation_succeeded{false};
  bool runtime_snapshot_complete{false};
  bool diagnostics_complete{false};
  dec3d::core::AuthoritativeFieldMask authoritative_field_mask{0};
  dec3d::core::CachedFieldMask serialized_cached_field_mask{0};
  dec3d::core::CachedFieldMask restart_invalidated_cached_mask{0};
  std::uint32_t runtime_registered_stage_digest{0};
  std::uint32_t runtime_validated_required_stage_digest{0};
  std::size_t runtime_stage_report_count{0};
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] CheckpointPayload BuildCheckpointPayload(
    const dec3d::state::CanonicalState& state,
    const RuntimeScaffold& runtime);

[[nodiscard]] dec3d::state::CanonicalState RestoreCanonicalStateFromCheckpoint(
    const CheckpointPayload& payload);

[[nodiscard]] RestartValidationResult ValidateRestartFromCheckpoint(
    const CheckpointPayload& payload);

[[nodiscard]] CheckpointContractReport BuildCheckpointContractReport(
    const CheckpointPayload& payload);

}  // namespace dec3d::runtime
