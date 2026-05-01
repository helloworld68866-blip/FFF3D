#include "runtime/p3_production_orchestrator.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool Contains(const std::string& text, const char* token) {
  return text.find(token) != std::string::npos;
}

dec3d::runtime::P3ProductionStageResult MakeStage(
    const char* stage_id,
    double dt_s,
    bool success) {
  return dec3d::runtime::MakeP3ProductionStageResultForTest(stage_id, dt_s, success);
}

dec3d::runtime::P3ProductionStageResult MakeRadiationStage(
    double dt_s,
    bool success) {
  auto stage = MakeStage("R", dt_s, success);
  stage.geometry_epoch = "post_H_committed_ALE_geometry";
  stage.state_epoch = "post_E_committed";
  stage.radiation_groups_epoch = "post_H_committed";
  stage.electron_thermal_state_epoch = "post_E_committed";
  if (success) {
    stage.report_line =
        "diagnostic_id=p3.production.stage.R"
        "; stage_id=R"
        "; success=true"
        "; radiation_stage_report_present=true"
        "; radiation_groups_epoch=post_H_committed"
        "; electron_thermal_state_epoch=post_E_committed"
        "; radiation_geometry_epoch=post_H_committed_ALE_geometry"
        "; radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry"
        "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
        "; radiation_boundary_model=thesis_marshak_vacuum"
        "; radiation_flux_limiter_model=harmonic"
        "; opacity_source_thesis_exact_match=false"
        "; parity_claim_allowed=false";
  } else {
    stage.failure_diagnostics =
        "diagnostic_id=p3.production.stage.R.failure"
        "; stage_id=R"
        "; success=false"
        "; canonical_state_mutated=false";
  }
  return stage;
}

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::ExecuteP3ProductionStep;
    using dec3d::runtime::ExecuteP3ProductionStepIntoRuntime;
    using dec3d::runtime::HasIngestedStageResult;
    using dec3d::runtime::IsStageRegistered;
    using dec3d::runtime::P3ProductionOperatorHooks;
    using dec3d::runtime::P3ProductionOptions;
    using dec3d::runtime::ValidateP3ProductionStepDiagnostics;

    constexpr double kDt = 2.5e-12;

    std::vector<std::string> order;
    P3ProductionOperatorHooks hooks;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.ale_geometry_committed = true;
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      auto stage = MakeStage("E", dt_s, true);
      stage.state_epoch = "post_E_committed";
      return stage;
    };
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      return MakeRadiationStage(dt_s, true);
    };

    const auto success = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    if (!success.success) {
      std::cerr << success.failure_diagnostics << '\n';
    }
    DEC3D_CHECK(success.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{4});
    DEC3D_CHECK(order[0] == "H");
    DEC3D_CHECK(order[1] == "T");
    DEC3D_CHECK(order[2] == "E");
    DEC3D_CHECK(order[3] == "R");
    DEC3D_CHECK(Contains(success.report_line, "diagnostic_id=p3.production.step"));
    DEC3D_CHECK(Contains(success.report_line, "stage_order=H,T,E,R"));
    DEC3D_CHECK(Contains(success.report_line, "same_dt_for_H_T_E_R=true"));
    DEC3D_CHECK(Contains(success.report_line, "stage_local_atomic=true"));
    DEC3D_CHECK(Contains(success.report_line, "step_atomic_across_H_T_E_R=false"));
    DEC3D_CHECK(Contains(success.report_line, "r_stage_registered=true"));
    DEC3D_CHECK(Contains(success.report_line, "a_stage_registered=false"));
    DEC3D_CHECK(Contains(success.report_line, "radiation_stage_report_present=true"));
    DEC3D_CHECK(Contains(success.report_line, "radiation_groups_epoch=post_H_committed"));
    DEC3D_CHECK(Contains(success.report_line, "electron_thermal_state_epoch=post_E_committed"));
    DEC3D_CHECK(Contains(
        success.report_line,
        "radiation_geometry_epoch=post_H_committed_ALE_geometry"));
    DEC3D_CHECK(Contains(
        success.report_line,
        "radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry"));
    DEC3D_CHECK(Contains(success.report_line, "opacity_source_thesis_exact_match=false"));
    DEC3D_CHECK(Contains(success.report_line, "parity_claim_allowed=false"));
    DEC3D_CHECK(ValidateP3ProductionStepDiagnostics(success));

    order.clear();
    bool t_called = false;
    bool e_called = false;
    bool r_called = false;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      return MakeStage("H", dt_s, false);
    };
    hooks.thermal = [&](double dt_s) {
      t_called = true;
      order.push_back("T");
      return MakeStage("T", dt_s, true);
    };
    hooks.equilibration = [&](double dt_s) {
      e_called = true;
      order.push_back("E");
      return MakeStage("E", dt_s, true);
    };
    hooks.radiation = [&](double dt_s) {
      r_called = true;
      order.push_back("R");
      return MakeRadiationStage(dt_s, true);
    };
    const auto h_failure = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!h_failure.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{1});
    DEC3D_CHECK(!t_called);
    DEC3D_CHECK(!e_called);
    DEC3D_CHECK(!r_called);
    DEC3D_CHECK(Contains(h_failure.failure_diagnostics, "failure_stage=H"));
    DEC3D_CHECK(Contains(h_failure.failure_diagnostics, "t_stage_executed=false"));
    DEC3D_CHECK(Contains(h_failure.failure_diagnostics, "e_stage_executed=false"));
    DEC3D_CHECK(Contains(h_failure.failure_diagnostics, "r_stage_executed=false"));

    order.clear();
    e_called = false;
    r_called = false;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      return MakeStage("T", dt_s, false);
    };
    hooks.equilibration = [&](double dt_s) {
      e_called = true;
      order.push_back("E");
      return MakeStage("E", dt_s, true);
    };
    hooks.radiation = [&](double dt_s) {
      r_called = true;
      order.push_back("R");
      return MakeRadiationStage(dt_s, true);
    };
    const auto t_failure = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!t_failure.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{2});
    DEC3D_CHECK(!e_called);
    DEC3D_CHECK(!r_called);
    DEC3D_CHECK(Contains(t_failure.failure_diagnostics, "diagnostic_id=p3.production.step.failure"));
    DEC3D_CHECK(Contains(t_failure.failure_diagnostics, "failure_stage=T"));
    DEC3D_CHECK(Contains(t_failure.failure_diagnostics, "r_stage_executed=false"));
    DEC3D_CHECK(Contains(t_failure.failure_diagnostics, "later_stages_executed=false"));

    order.clear();
    r_called = false;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      return MakeStage("E", dt_s, false);
    };
    hooks.radiation = [&](double dt_s) {
      r_called = true;
      order.push_back("R");
      return MakeRadiationStage(dt_s, true);
    };
    const auto e_failure = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!e_failure.success);
    DEC3D_CHECK(!r_called);
    DEC3D_CHECK(Contains(e_failure.failure_diagnostics, "failure_stage=E"));
    DEC3D_CHECK(Contains(e_failure.failure_diagnostics, "committed_prior_stages=H,T"));

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.ale_geometry_committed = true;
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      auto stage = MakeStage("E", dt_s, true);
      stage.state_epoch = "post_E_committed";
      return stage;
    };
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      return MakeRadiationStage(dt_s, false);
    };
    const auto r_failure = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!r_failure.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{4});
    DEC3D_CHECK(Contains(r_failure.failure_diagnostics, "failure_stage=R"));
    DEC3D_CHECK(Contains(r_failure.failure_diagnostics, "committed_prior_stages=H,T,E"));
    DEC3D_CHECK(Contains(r_failure.failure_diagnostics, "r_stage_published=false"));
    DEC3D_CHECK(Contains(
        r_failure.failure_diagnostics,
        "step_atomic_across_H_T_E_R=false"));

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.ale_geometry_committed = true;
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      auto stage = MakeStage("E", dt_s, true);
      stage.state_epoch = "post_E_committed";
      return stage;
    };
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      auto stage = MakeRadiationStage(dt_s, true);
      stage.geometry_epoch = "stale_geometry";
      stage.report_line =
          "diagnostic_id=p3.production.stage.R"
          "; stage_id=R"
          "; success=true"
          "; radiation_stage_report_present=true"
          "; radiation_groups_epoch=post_H_committed"
          "; electron_thermal_state_epoch=post_E_committed"
          "; radiation_geometry_epoch=stale_geometry"
          "; radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry"
          "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
          "; radiation_boundary_model=thesis_marshak_vacuum"
          "; radiation_flux_limiter_model=harmonic"
          "; opacity_source_thesis_exact_match=false"
          "; parity_claim_allowed=false";
      return stage;
    };
    const auto stale_r_geometry = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!stale_r_geometry.success);
    DEC3D_CHECK(Contains(stale_r_geometry.failure_diagnostics, "failure_stage=R"));
    DEC3D_CHECK(Contains(
        stale_r_geometry.failure_diagnostics,
        "radiation stage did not consume committed ALE geometry"));

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      auto stage = MakeStage("E", dt_s, true);
      stage.state_epoch = "post_E_committed";
      return stage;
    };
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      auto stage = MakeRadiationStage(dt_s, true);
      stage.electron_thermal_state_epoch = "pre_E_state";
      stage.report_line =
          "diagnostic_id=p3.production.stage.R"
          "; stage_id=R"
          "; success=true"
          "; radiation_stage_report_present=true"
          "; radiation_groups_epoch=post_H_committed"
          "; electron_thermal_state_epoch=pre_E_state"
          "; radiation_geometry_epoch=post_H_committed_ALE_geometry"
          "; radiation_hydro_coupling_mode=frozen_hydro_post_H_geometry"
          "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
          "; radiation_boundary_model=thesis_marshak_vacuum"
          "; radiation_flux_limiter_model=harmonic"
          "; opacity_source_thesis_exact_match=false"
          "; parity_claim_allowed=false";
      return stage;
    };
    const auto stale_r_state = ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!stale_r_state.success);
    DEC3D_CHECK(Contains(stale_r_state.failure_diagnostics, "failure_stage=R"));
    DEC3D_CHECK(Contains(stale_r_state.failure_diagnostics, "post-E electron thermal state epoch"));

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = MakeStage("H", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = MakeStage("T", dt_s, true);
      stage.geometry_epoch = "post_H_committed_ALE_geometry";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      auto stage = MakeStage("E", dt_s, true);
      stage.state_epoch = "post_E_committed";
      return stage;
    };
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      auto stage = MakeRadiationStage(dt_s, true);
      stage.report_line =
          "diagnostic_id=p3.production.stage.R"
          "; stage_id=R"
          "; success=true";
      return stage;
    };
    const auto missing_r_diagnostics =
        ExecuteP3ProductionStep(P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!missing_r_diagnostics.success);
    DEC3D_CHECK(Contains(missing_r_diagnostics.failure_diagnostics, "failure_stage=R"));
    DEC3D_CHECK(Contains(
        missing_r_diagnostics.failure_diagnostics,
        "radiation stage diagnostics are incomplete"));

    auto runtime = CreateRuntimeScaffold(PhaseId::p3, "p3-v0.7");
    hooks.radiation = [&](double dt_s) {
      order.push_back("R");
      return MakeRadiationStage(dt_s, true);
    };
    order.clear();
    const auto runtime_success =
        ExecuteP3ProductionStepIntoRuntime(runtime, P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(runtime_success.success);
    DEC3D_CHECK(IsStageRegistered(runtime, StageId::hydro));
    DEC3D_CHECK(IsStageRegistered(runtime, StageId::thermal));
    DEC3D_CHECK(IsStageRegistered(runtime, StageId::equilibration));
    DEC3D_CHECK(IsStageRegistered(runtime, StageId::radiation));
    DEC3D_CHECK(!IsStageRegistered(runtime, StageId::alpha));
    DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::hydro));
    DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::thermal));
    DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::equilibration));
    DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::radiation));
    DEC3D_CHECK(Contains(runtime_success.report_line, "runtime_scaffold_ingested=true"));

    auto alpha_runtime = CreateRuntimeScaffold(PhaseId::p3, "p3-v0.7");
    const auto alpha_rejected =
        ExecuteP3ProductionStepIntoRuntime(alpha_runtime, P3ProductionOptions{kDt, true}, hooks);
    DEC3D_CHECK(!alpha_rejected.success);
    DEC3D_CHECK(Contains(alpha_rejected.failure_diagnostics, "forbidden alpha registration"));

    auto p2_runtime = CreateRuntimeScaffold(PhaseId::p2, "p3-v0.7");
    const auto wrong_phase =
        ExecuteP3ProductionStepIntoRuntime(p2_runtime, P3ProductionOptions{kDt}, hooks);
    DEC3D_CHECK(!wrong_phase.success);
    DEC3D_CHECK(Contains(
        wrong_phase.failure_diagnostics,
        "runtime scaffold must be initialized for P3"));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
