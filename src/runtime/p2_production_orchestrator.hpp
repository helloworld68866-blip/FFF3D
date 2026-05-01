#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"

#include <functional>
#include <string>

namespace dec3d::runtime {

struct P2ProductionOptions {
  double dt_s{0.0};
  bool register_r_stage{false};
  bool register_a_stage{false};
};

struct P2ProductionStageResult {
  bool success{false};
  std::string stage_id;
  double dt_s{0.0};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  bool ale_geometry_committed{false};
  std::string geometry_epoch;
  std::string report_line;
  std::string failure_diagnostics;
};

struct P2ProductionOperatorHooks {
  std::function<P2ProductionStageResult(double)> hydro;
  std::function<P2ProductionStageResult(double)> thermal;
  std::function<P2ProductionStageResult(double)> equilibration;
};

struct P2ProductionStepResult {
  bool success{false};
  double dt_s{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  P2ProductionStageResult hydro;
  P2ProductionStageResult thermal;
  P2ProductionStageResult equilibration;
};

[[nodiscard]] P2ProductionStepResult ExecuteP2ProductionStep(
    const P2ProductionOptions& options,
    const P2ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] P2ProductionStepResult ExecuteP2ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P2ProductionOptions& options,
    const P2ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] bool ValidateP2ProductionStepDiagnostics(
    const P2ProductionStepResult& result) noexcept;

[[nodiscard]] P2ProductionStageResult MakeP2ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success);

}  // namespace dec3d::runtime
