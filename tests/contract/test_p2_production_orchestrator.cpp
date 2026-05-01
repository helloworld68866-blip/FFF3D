#include "runtime/p2_production_orchestrator.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "test_assert.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <vector>

int main() {
  try {
    using dec3d::runtime::ExecuteP2ProductionStep;
    using dec3d::runtime::P2ProductionOperatorHooks;
    using dec3d::runtime::P2ProductionOptions;

    std::vector<std::string> order;
    P2ProductionOperatorHooks hooks;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("H", dt_s, true);
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("T", dt_s, true);
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("E", dt_s, true);
    };

    const auto result = ExecuteP2ProductionStep(P2ProductionOptions{1.0e-12}, hooks);
    if (!result.success) {
      std::cerr << result.failure_diagnostics << '\n';
    }
    DEC3D_CHECK(result.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{3});
    DEC3D_CHECK(order[0] == "H");
    DEC3D_CHECK(order[1] == "T");
    DEC3D_CHECK(order[2] == "E");
    DEC3D_CHECK(result.report_line.find("stage_order=H,T,E") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("r_stage_registered=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("a_stage_registered=false") != std::string::npos);
    DEC3D_CHECK(dec3d::runtime::ValidateP2ProductionStepDiagnostics(result));

    order.clear();
    bool equilibration_called = false;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("H", dt_s, true);
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("T", dt_s, false);
    };
    hooks.equilibration = [&](double dt_s) {
      equilibration_called = true;
      order.push_back("E");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("E", dt_s, true);
    };
    const auto failed_thermal =
        ExecuteP2ProductionStep(P2ProductionOptions{1.0e-12}, hooks);
    DEC3D_CHECK(!failed_thermal.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{2});
    DEC3D_CHECK(order[0] == "H");
    DEC3D_CHECK(order[1] == "T");
    DEC3D_CHECK(!equilibration_called);
    DEC3D_CHECK(failed_thermal.failure_diagnostics.find("stage_id=T") != std::string::npos);
    DEC3D_CHECK(failed_thermal.failure_diagnostics.find("e_stage_executed=false") !=
                std::string::npos);

    order.clear();
    equilibration_called = false;
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = dec3d::runtime::MakeP2ProductionStageResultForTest("H", dt_s, true);
      stage.ale_geometry_committed = true;
      stage.geometry_epoch = "geometry_n_plus_1";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = dec3d::runtime::MakeP2ProductionStageResultForTest("T", dt_s, true);
      stage.geometry_epoch = "geometry_n";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      equilibration_called = true;
      order.push_back("E");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("E", dt_s, true);
    };
    const auto stale_thermal_geometry =
        ExecuteP2ProductionStep(P2ProductionOptions{1.0e-12}, hooks);
    DEC3D_CHECK(!stale_thermal_geometry.success);
    DEC3D_CHECK_EQ(order.size(), std::size_t{2});
    DEC3D_CHECK(!equilibration_called);
    DEC3D_CHECK(stale_thermal_geometry.failure_diagnostics.find("stage_id=T") !=
                std::string::npos);
    DEC3D_CHECK(stale_thermal_geometry.failure_diagnostics.find(
                    "failure_reason=thermal stage did not consume committed ALE geometry") !=
                std::string::npos);

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      auto stage = dec3d::runtime::MakeP2ProductionStageResultForTest("H", dt_s, true);
      stage.ale_geometry_committed = true;
      stage.geometry_epoch = "geometry_n_plus_1";
      return stage;
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      auto stage = dec3d::runtime::MakeP2ProductionStageResultForTest("T", dt_s, true);
      stage.geometry_epoch = "geometry_n_plus_1";
      return stage;
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("E", dt_s, true);
    };
    const auto ale_geometry_result =
        ExecuteP2ProductionStep(P2ProductionOptions{1.0e-12}, hooks);
    DEC3D_CHECK(ale_geometry_result.success);
    DEC3D_CHECK(ale_geometry_result.report_line.find("ale_geometry_committed=true") !=
                std::string::npos);
    DEC3D_CHECK(ale_geometry_result.report_line.find(
                    "thermal_geometry_epoch=geometry_n_plus_1") != std::string::npos);

    order.clear();
    hooks.hydro = [&](double dt_s) {
      order.push_back("H");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("H", dt_s, true);
    };
    hooks.thermal = [&](double dt_s) {
      order.push_back("T");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("T", dt_s, true);
    };
    hooks.equilibration = [&](double dt_s) {
      order.push_back("E");
      return dec3d::runtime::MakeP2ProductionStageResultForTest("E", dt_s, true);
    };

    auto runtime = dec3d::runtime::CreateRuntimeScaffold(
        dec3d::core::PhaseId::p2,
        "p2-v0.7");
    const auto runtime_result = dec3d::runtime::ExecuteP2ProductionStepIntoRuntime(
        runtime,
        P2ProductionOptions{1.0e-12},
        hooks);
    DEC3D_CHECK(runtime_result.success);
    DEC3D_CHECK(dec3d::runtime::IsStageRegistered(runtime, dec3d::core::StageId::hydro));
    DEC3D_CHECK(dec3d::runtime::IsStageRegistered(runtime, dec3d::core::StageId::thermal));
    DEC3D_CHECK(dec3d::runtime::IsStageRegistered(runtime, dec3d::core::StageId::equilibration));
    DEC3D_CHECK(!dec3d::runtime::IsStageRegistered(runtime, dec3d::core::StageId::radiation));
    DEC3D_CHECK(!dec3d::runtime::IsStageRegistered(runtime, dec3d::core::StageId::alpha));
    DEC3D_CHECK(dec3d::runtime::HasIngestedStageResult(runtime, dec3d::core::StageId::hydro));
    DEC3D_CHECK(dec3d::runtime::HasIngestedStageResult(runtime, dec3d::core::StageId::thermal));
    DEC3D_CHECK(dec3d::runtime::HasIngestedStageResult(runtime, dec3d::core::StageId::equilibration));
    const auto runtime_report = dec3d::runtime::BuildRuntimeSubstrateReport(runtime);
    DEC3D_CHECK(runtime_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(runtime_result.report_line.find("runtime_scaffold_ingested=true") !=
                std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
