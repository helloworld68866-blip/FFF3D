#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"

#include <functional>
#include <string>

namespace dec3d::runtime {

struct P4ProductionOptions {
  double dt_s{0.0};
};

struct P4ProductionStageResult {
  bool success{false};
  std::string stage_id;
  double dt_s{0.0};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  bool ale_geometry_committed{false};
  bool canonical_state_mutated{false};
  bool metadata_written{false};
  std::string geometry_epoch;
  std::string state_epoch;
  std::string radiation_groups_epoch;
  std::string electron_thermal_state_epoch;
  std::string radiation_boundary_model;
  std::string radiation_flux_limiter_model;
  std::string alpha_state_epoch;
  std::string alpha_geometry_epoch;
  std::string alpha_thermal_state_epoch;
  std::string report_line;
  std::string failure_diagnostics;
};

struct P4ProductionOperatorHooks {
  std::function<P4ProductionStageResult(double)> hydro;
  std::function<P4ProductionStageResult(double)> thermal;
  std::function<P4ProductionStageResult(double)> equilibration;
  std::function<P4ProductionStageResult(double)> radiation;
  std::function<P4ProductionStageResult(double)> alpha;
};

struct P4ProductionStepResult {
  bool success{false};
  double dt_s{0.0};
  bool advance_to_next_timestep{false};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  P4ProductionStageResult hydro;
  P4ProductionStageResult thermal;
  P4ProductionStageResult equilibration;
  P4ProductionStageResult radiation;
  P4ProductionStageResult alpha;
};

[[nodiscard]] P4ProductionStepResult ExecuteP4ProductionStep(
    const P4ProductionOptions& options,
    const P4ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] P4ProductionStepResult ExecuteP4ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P4ProductionOptions& options,
    const P4ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] bool ValidateP4ProductionStepDiagnostics(
    const P4ProductionStepResult& result) noexcept;

[[nodiscard]] P4ProductionStageResult MakeP4ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success);

}  // namespace dec3d::runtime
