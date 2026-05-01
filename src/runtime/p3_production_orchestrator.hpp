#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "runtime/runtime_scaffold.hpp"

#include <functional>
#include <string>

namespace dec3d::runtime {

struct P3ProductionOptions {
  double dt_s{0.0};
  bool register_a_stage{false};
};

struct P3ProductionStageResult {
  bool success{false};
  std::string stage_id;
  double dt_s{0.0};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  bool ale_geometry_committed{false};
  bool canonical_state_mutated{false};
  std::string geometry_epoch;
  std::string state_epoch;
  std::string radiation_groups_epoch;
  std::string electron_thermal_state_epoch;
  std::string radiation_boundary_model;
  std::string radiation_flux_limiter_model;
  std::string report_line;
  std::string failure_diagnostics;
};

struct P3ProductionOperatorHooks {
  std::function<P3ProductionStageResult(double)> hydro;
  std::function<P3ProductionStageResult(double)> thermal;
  std::function<P3ProductionStageResult(double)> equilibration;
  std::function<P3ProductionStageResult(double)> radiation;
};

struct P3ProductionStepResult {
  bool success{false};
  double dt_s{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  P3ProductionStageResult hydro;
  P3ProductionStageResult thermal;
  P3ProductionStageResult equilibration;
  P3ProductionStageResult radiation;
};

[[nodiscard]] P3ProductionStepResult ExecuteP3ProductionStep(
    const P3ProductionOptions& options,
    const P3ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] P3ProductionStepResult ExecuteP3ProductionStepIntoRuntime(
    RuntimeScaffold& runtime,
    const P3ProductionOptions& options,
    const P3ProductionOperatorHooks& hooks) noexcept;

[[nodiscard]] bool ValidateP3ProductionStepDiagnostics(
    const P3ProductionStepResult& result) noexcept;

[[nodiscard]] P3ProductionStageResult MakeP3ProductionStageResultForTest(
    const char* stage_id,
    double dt_s,
    bool success);

}  // namespace dec3d::runtime
