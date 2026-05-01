#include "io/runtime_output.hpp"

#include <cassert>

int main() {
  using dec3d::io::CheckpointTriggerConfig;
  using dec3d::io::ShouldWriteCheckpoint;

  CheckpointTriggerConfig disabled;
  disabled.every_steps = 0;
  disabled.interval_s = 0.0;
  assert(!ShouldWriteCheckpoint(disabled, 10, 1.0e-9, 0.0));

  CheckpointTriggerConfig step_only;
  step_only.every_steps = 5;
  step_only.interval_s = 0.0;
  assert(ShouldWriteCheckpoint(step_only, 10, 0.0, 0.0));
  assert(!ShouldWriteCheckpoint(step_only, 11, 0.0, 0.0));

  CheckpointTriggerConfig time_only;
  time_only.every_steps = 0;
  time_only.interval_s = 2.0e-12;
  assert(!ShouldWriteCheckpoint(time_only, 7, 1.0e-12, 0.0));
  assert(ShouldWriteCheckpoint(time_only, 7, 2.0e-12, 0.0));
  assert(ShouldWriteCheckpoint(time_only, 7, 4.1e-12, 2.0e-12));

  CheckpointTriggerConfig either;
  either.every_steps = 10;
  either.interval_s = 1.0e-12;
  assert(ShouldWriteCheckpoint(either, 20, 0.0, 0.0));
  assert(ShouldWriteCheckpoint(either, 21, 1.0e-12, 0.0));

  assert(dec3d::io::FormatStepSnapshotName("fields", 10) == "fields_000010.snap");
  assert(dec3d::io::FormatStepSnapshotName("restart", 100) == "restart_000100.snap");
  return 0;
}
