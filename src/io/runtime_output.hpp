#pragma once

#include "io/input_deck.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <filesystem>
#include <string>

namespace dec3d::io {

struct CheckpointTriggerConfig {
  int every_steps{0};
  double interval_s{0.0};
};

[[nodiscard]] bool ShouldWriteCheckpoint(
    const CheckpointTriggerConfig& config,
    std::size_t step,
    double time_s,
    double last_written_time_s) noexcept;

[[nodiscard]] std::string FormatStepSnapshotName(
    const std::string& prefix,
    std::size_t step);

[[nodiscard]] CheckpointTriggerConfig FieldCheckpointTrigger(
    const OutputConfig& output) noexcept;

[[nodiscard]] CheckpointTriggerConfig RestartCheckpointTrigger(
    const OutputConfig& output) noexcept;

[[nodiscard]] CheckpointTriggerConfig HistoryProfileTrigger(
    const OutputConfig& output) noexcept;

struct RuntimeOutputWriteResult {
  bool success{false};
  bool dec3d_out_written{false};
  bool field_checkpoint_written{false};
  bool restart_checkpoint_written{false};
  bool history_profile_written{false};
  double dec3d_out_write_s{0.0};
  double field_checkpoint_write_s{0.0};
  double restart_checkpoint_write_s{0.0};
  double history_profile_write_s{0.0};
  double runtime_output_total_s{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct RuntimeOutputStepState {
  double last_field_checkpoint_time_s{0.0};
  double last_restart_checkpoint_time_s{0.0};
  double last_history_profile_time_s{0.0};
};

struct RuntimeOutputExtrema {
  double rho_min{0.0};
  double rho_max{0.0};
  double te_min_keV{0.0};
  double te_max_keV{0.0};
  double ti_min_keV{0.0};
  double ti_max_keV{0.0};
};

struct RuntimeOutputExtremaResult {
  bool success{false};
  RuntimeOutputExtrema extrema;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] RuntimeOutputExtremaResult ComputeRuntimeOutputExtrema(
    const dec3d::state::CanonicalState& state) noexcept;

[[nodiscard]] RuntimeOutputWriteResult WriteRuntimeStepDec3DOut(
    const InputDeckConfig& config,
    std::size_t step,
    double time_s,
    double dt_s,
    double wall_elapsed_s,
    const std::string& status,
    const RuntimeOutputExtrema& extrema) noexcept;

[[nodiscard]] RuntimeOutputWriteResult WriteInitialRuntimeOutputs(
    const InputDeckConfig& config,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path) noexcept;

[[nodiscard]] RuntimeOutputWriteResult WriteRuntimeStepOutputs(
    const InputDeckConfig& config,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s,
    double dt_s,
    double wall_elapsed_s,
    const std::string& status,
    RuntimeOutputStepState& output_state) noexcept;

}  // namespace dec3d::io
