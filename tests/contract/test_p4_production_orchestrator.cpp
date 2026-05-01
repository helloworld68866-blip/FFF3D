#include "runtime/p3_production_orchestrator.hpp"
#include "runtime/p4_production_orchestrator.hpp"
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

dec3d::runtime::P4ProductionStageResult MakeStage(
    const char* stage_id,
    double dt_s,
    bool success) {
  return dec3d::runtime::MakeP4ProductionStageResultForTest(stage_id, dt_s, success);
}

dec3d::runtime::P4ProductionStageResult MakeRadiationStage(
    double dt_s,
    bool success) {
  auto stage = MakeStage("R", dt_s, success);
  stage.geometry_epoch = "post_H_committed_ALE_geometry";
  stage.state_epoch = "post_R_committed";
  stage.radiation_groups_epoch = "post_H_committed";
  stage.electron_thermal_state_epoch = "post_E_committed";
  stage.radiation_boundary_model = "thesis_marshak_vacuum";
  stage.radiation_flux_limiter_model = "harmonic";
  if (success) {
    stage.report_line =
        "diagnostic_id=p4.production.stage.R"
        "; stage_id=R"
        "; success=true"
        "; radiation_stage_report_present=true"
        "; radiation_hydro_terms_report_present=true"
        "; radiation_groups_epoch=post_H_committed"
        "; electron_thermal_state_epoch=post_E_committed"
        "; radiation_geometry_epoch=post_H_committed_ALE_geometry"
        "; radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"
        "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
        "; radiation_boundary_model=thesis_marshak_vacuum"
        "; radiation_flux_limiter_model=harmonic"
        "; opacity_source_thesis_exact_match=false"
        "; parity_claim_allowed=false";
  } else {
    stage.failure_diagnostics =
        "diagnostic_id=p4.production.stage.R.failure"
        "; stage_id=R"
        "; success=false"
        "; canonical_state_mutated=false";
  }
  return stage;
}

dec3d::runtime::P4ProductionStageResult MakeAlphaStage(
    double dt_s,
    bool success) {
  auto stage = MakeStage("A", dt_s, success);
  stage.geometry_epoch = "post_H_committed_ALE_geometry";
  stage.state_epoch = "post_R_committed";
  stage.alpha_geometry_epoch = "post_H_committed_ALE_geometry";
  stage.alpha_state_epoch = "post_H_committed";
  stage.alpha_thermal_state_epoch = "post_R_committed";
  if (success) {
    stage.report_line =
        "diagnostic_id=p4.production.stage.A"
        "; stage_id=A"
        "; success=true"
        "; alpha_stage_report_present=true"
        "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
        "; alpha_state_epoch=post_H_committed"
        "; alpha_thermal_state_epoch=post_R_committed"
        "; alpha_transport_model=atzeni_one_group"
        "; reactivity_model_executed=bosch_hale_dt"
        "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
        "; fuel_depletion_enabled=false"
        "; separate_dt_species_authoritative=false"
        "; alpha_updated_fields=alpha_state,e_electron,e_fluid_total"
        "; alpha_backend_executed=hypre_parcsr_gmres_boomeramg"
        "; owned_slab_writeback_only=true"
        "; rank0_gather_solve_used=false"
        "; serial_dense_fallback_used=false"
        "; parity_claim_allowed=false";
  } else {
    stage.failure_diagnostics =
        "diagnostic_id=p4.production.stage.A.failure"
        "; stage_id=A"
        "; success=false"
        "; canonical_state_mutated=false"
        "; failed_stage_published=false"
        "; a_stage_published=false";
  }
  return stage;
}

dec3d::runtime::P4ProductionOperatorHooks MakeSuccessfulHooks(
    std::vector<std::string>& order,
    double dt_s) {
  dec3d::runtime::P4ProductionOperatorHooks hooks;
  hooks.hydro = [&](double step_dt) {
    order.push_back("H");
    auto stage = MakeStage("H", step_dt, true);
    stage.ale_geometry_committed = true;
    stage.geometry_epoch = "post_H_committed_ALE_geometry";
    stage.alpha_state_epoch = "post_H_committed";
    stage.report_line +=
        "; alpha_hydro_terms_report_present=true"
        "; alpha_hydro_terms_enabled=true"
        "; alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure";
    return stage;
  };
  hooks.thermal = [&](double step_dt) {
    order.push_back("T");
    auto stage = MakeStage("T", step_dt, true);
    stage.geometry_epoch = "post_H_committed_ALE_geometry";
    return stage;
  };
  hooks.equilibration = [&](double step_dt) {
    order.push_back("E");
    auto stage = MakeStage("E", step_dt, true);
    stage.state_epoch = "post_E_committed";
    return stage;
  };
  hooks.radiation = [&](double step_dt) {
    order.push_back("R");
    return MakeRadiationStage(step_dt, true);
  };
  hooks.alpha = [&](double step_dt) {
    order.push_back("A");
    return MakeAlphaStage(step_dt, true);
  };
  DEC3D_CHECK(dt_s >= 0.0);
  return hooks;
}

}  // namespace

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::PhaseId;
    using dec3d::core::StageId;
    using dec3d::core::ToMask;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::ExecuteP4ProductionStep;
    using dec3d::runtime::ExecuteP4ProductionStepIntoRuntime;
    using dec3d::runtime::HasIngestedStageResult;
    using dec3d::runtime::IsStageRegistered;
    using dec3d::runtime::P4ProductionOptions;
    using dec3d::runtime::ValidateP4ProductionStepDiagnostics;

    constexpr double kDt = 2.5e-12;

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      if (!result.success) {
        std::cerr << result.failure_diagnostics << '\n';
      }
      DEC3D_CHECK(result.success);
      DEC3D_CHECK_EQ(order.size(), std::size_t{5});
      DEC3D_CHECK(order[0] == "H");
      DEC3D_CHECK(order[1] == "T");
      DEC3D_CHECK(order[2] == "E");
      DEC3D_CHECK(order[3] == "R");
      DEC3D_CHECK(order[4] == "A");
      DEC3D_CHECK(Contains(result.report_line, "diagnostic_id=p4.production.step"));
      DEC3D_CHECK(Contains(result.report_line, "stage_order=H,T,E,R,A"));
      DEC3D_CHECK(Contains(result.report_line, "same_dt_for_H_T_E_R_A=true"));
      DEC3D_CHECK(Contains(result.report_line, "stage_local_atomic=true"));
      DEC3D_CHECK(Contains(result.report_line, "step_atomic_across_H_T_E_R_A=false"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_hydro_terms_report_present=true"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_hydro_terms_enabled=true"));
      DEC3D_CHECK(Contains(
          result.report_line,
          "alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure"));
      DEC3D_CHECK(Contains(result.report_line, "radiation_hydro_terms_report_present=true"));
      DEC3D_CHECK(Contains(
          result.report_line,
          "radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_stage_report_present=true"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_geometry_epoch=post_H_committed_ALE_geometry"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_state_epoch=post_H_committed"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_thermal_state_epoch=post_R_committed"));
      DEC3D_CHECK(Contains(result.report_line, "coefficient_time_level=old_time_lagged_to_A_stage_entry"));
      DEC3D_CHECK(Contains(result.report_line, "alpha_updated_fields=alpha_state,e_electron,e_fluid_total"));
      DEC3D_CHECK(Contains(result.report_line, "fuel_depletion_enabled=false"));
      DEC3D_CHECK(Contains(result.report_line, "advance_to_next_timestep=true"));
      DEC3D_CHECK(ValidateP4ProductionStepDiagnostics(result));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        return MakeAlphaStage(step_dt, false);
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(order.size(), std::size_t{5});
      DEC3D_CHECK(Contains(result.failure_diagnostics, "diagnostic_id=p4.production.step.failure"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=A"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "committed_prior_stages=H,T,E,R"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "a_stage_executed=true"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "a_stage_published=false"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "nested_failure_diagnostics={"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "diagnostic_id=p4.production.stage.A.failure"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "advance_to_next_timestep=false"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "step_atomic_across_H_T_E_R_A=false"));
      DEC3D_CHECK(ValidateP4ProductionStepDiagnostics(result));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.hydro = [&](double step_dt) {
        order.push_back("H");
        return MakeStage("H", step_dt, false);
      };
      bool t_called = false;
      bool e_called = false;
      bool r_called = false;
      bool a_called = false;
      hooks.thermal = [&](double step_dt) {
        t_called = true;
        order.push_back("T");
        return MakeStage("T", step_dt, true);
      };
      hooks.equilibration = [&](double step_dt) {
        e_called = true;
        order.push_back("E");
        return MakeStage("E", step_dt, true);
      };
      hooks.radiation = [&](double step_dt) {
        r_called = true;
        order.push_back("R");
        return MakeRadiationStage(step_dt, true);
      };
      hooks.alpha = [&](double step_dt) {
        a_called = true;
        order.push_back("A");
        return MakeAlphaStage(step_dt, true);
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(order.size(), std::size_t{1});
      DEC3D_CHECK(!t_called);
      DEC3D_CHECK(!e_called);
      DEC3D_CHECK(!r_called);
      DEC3D_CHECK(!a_called);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=H"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "t_stage_executed=false"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "a_stage_executed=false"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.hydro = [&](double step_dt) {
        order.push_back("H");
        auto stage = MakeStage("H", step_dt, true);
        stage.ale_geometry_committed = true;
        stage.geometry_epoch = "post_H_committed_ALE_geometry";
        stage.alpha_state_epoch = "post_H_committed";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=H"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "hydro stage did not report alpha hydro terms"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.radiation = [&](double step_dt) {
        order.push_back("R");
        auto stage = MakeRadiationStage(step_dt, true);
        stage.report_line =
            "diagnostic_id=p4.production.stage.R"
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
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=R"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "radiation stage diagnostics are incomplete"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.thermal = [&](double step_dt) {
        order.push_back("T");
        return MakeStage("T", step_dt, false);
      };
      bool e_called = false;
      bool r_called = false;
      bool a_called = false;
      hooks.equilibration = [&](double step_dt) {
        e_called = true;
        order.push_back("E");
        return MakeStage("E", step_dt, true);
      };
      hooks.radiation = [&](double step_dt) {
        r_called = true;
        order.push_back("R");
        return MakeRadiationStage(step_dt, true);
      };
      hooks.alpha = [&](double step_dt) {
        a_called = true;
        order.push_back("A");
        return MakeAlphaStage(step_dt, true);
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!e_called);
      DEC3D_CHECK(!r_called);
      DEC3D_CHECK(!a_called);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=T"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "committed_prior_stages=H"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.equilibration = [&](double step_dt) {
        order.push_back("E");
        return MakeStage("E", step_dt, false);
      };
      bool r_called = false;
      bool a_called = false;
      hooks.radiation = [&](double step_dt) {
        r_called = true;
        order.push_back("R");
        return MakeRadiationStage(step_dt, true);
      };
      hooks.alpha = [&](double step_dt) {
        a_called = true;
        order.push_back("A");
        return MakeAlphaStage(step_dt, true);
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!r_called);
      DEC3D_CHECK(!a_called);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=E"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "committed_prior_stages=H,T"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.radiation = [&](double step_dt) {
        order.push_back("R");
        return MakeRadiationStage(step_dt, false);
      };
      bool a_called = false;
      hooks.alpha = [&](double step_dt) {
        a_called = true;
        order.push_back("A");
        return MakeAlphaStage(step_dt, true);
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!a_called);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=R"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "committed_prior_stages=H,T,E"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "a_stage_executed=false"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = {};
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=preflight"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "missing A hook"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "a_stage_executed=false"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        auto stage = MakeAlphaStage(step_dt, true);
        stage.alpha_geometry_epoch = "stale_geometry";
        stage.report_line =
            "diagnostic_id=p4.production.stage.A"
            "; stage_id=A"
            "; success=true"
            "; alpha_stage_report_present=true"
            "; alpha_geometry_epoch=stale_geometry"
            "; alpha_state_epoch=post_H_committed"
            "; alpha_thermal_state_epoch=post_R_committed"
            "; alpha_transport_model=atzeni_one_group"
            "; reactivity_model_executed=bosch_hale_dt"
            "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
            "; fuel_depletion_enabled=false"
            "; separate_dt_species_authoritative=false"
            "; alpha_updated_fields=alpha_state,e_electron,e_fluid_total"
            "; alpha_backend_executed=hypre_parcsr_gmres_boomeramg"
            "; owned_slab_writeback_only=true"
            "; rank0_gather_solve_used=false"
            "; serial_dense_fallback_used=false"
            "; parity_claim_allowed=false";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "failure_stage=A"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "alpha stage did not consume committed ALE geometry"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        auto stage = MakeAlphaStage(step_dt, true);
        stage.alpha_state_epoch = "pre_H_alpha";
        stage.report_line =
            "diagnostic_id=p4.production.stage.A"
            "; stage_id=A"
            "; success=true"
            "; alpha_stage_report_present=true"
            "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
            "; alpha_state_epoch=pre_H_alpha"
            "; alpha_thermal_state_epoch=post_R_committed"
            "; alpha_transport_model=atzeni_one_group"
            "; reactivity_model_executed=bosch_hale_dt"
            "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
            "; fuel_depletion_enabled=false"
            "; separate_dt_species_authoritative=false"
            "; alpha_updated_fields=alpha_state,e_electron,e_fluid_total"
            "; alpha_backend_executed=hypre_parcsr_gmres_boomeramg"
            "; owned_slab_writeback_only=true"
            "; rank0_gather_solve_used=false"
            "; serial_dense_fallback_used=false"
            "; parity_claim_allowed=false";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "alpha stage did not report post-H alpha state epoch"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        auto stage = MakeAlphaStage(step_dt, true);
        stage.alpha_thermal_state_epoch = "pre_R_thermal";
        stage.report_line =
            "diagnostic_id=p4.production.stage.A"
            "; stage_id=A"
            "; success=true"
            "; alpha_stage_report_present=true"
            "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
            "; alpha_state_epoch=post_H_committed"
            "; alpha_thermal_state_epoch=pre_R_thermal"
            "; alpha_transport_model=atzeni_one_group"
            "; reactivity_model_executed=bosch_hale_dt"
            "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
            "; fuel_depletion_enabled=false"
            "; separate_dt_species_authoritative=false"
            "; alpha_updated_fields=alpha_state,e_electron,e_fluid_total"
            "; alpha_backend_executed=hypre_parcsr_gmres_boomeramg"
            "; owned_slab_writeback_only=true"
            "; rank0_gather_solve_used=false"
            "; serial_dense_fallback_used=false"
            "; parity_claim_allowed=false";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "alpha stage did not report post-R thermal state epoch"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        auto stage = MakeAlphaStage(step_dt, true);
        stage.updated_fields = ToMask(AuthoritativeField::alpha_state) |
                               ToMask(AuthoritativeField::e_electron) |
                               ToMask(AuthoritativeField::radiation_groups);
        stage.report_line =
            "diagnostic_id=p4.production.stage.A"
            "; stage_id=A"
            "; success=true"
            "; alpha_stage_report_present=true"
            "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
            "; alpha_state_epoch=post_H_committed"
            "; alpha_thermal_state_epoch=post_R_committed"
            "; alpha_transport_model=atzeni_one_group"
            "; reactivity_model_executed=bosch_hale_dt"
            "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
            "; fuel_depletion_enabled=false"
            "; separate_dt_species_authoritative=false"
            "; alpha_updated_fields=alpha_state,e_electron,radiation_groups"
            "; alpha_backend_executed=hypre_parcsr_gmres_boomeramg"
            "; owned_slab_writeback_only=true"
            "; rank0_gather_solve_used=false"
            "; serial_dense_fallback_used=false"
            "; parity_claim_allowed=false";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "alpha stage write set is invalid"));
    }

    {
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, 0.0);
      hooks.alpha = [&](double step_dt) {
        order.push_back("A");
        auto stage = MakeAlphaStage(step_dt, true);
        stage.updated_fields = 0u;
        stage.canonical_state_mutated = false;
        stage.report_line =
            "diagnostic_id=p4.production.stage.A"
            "; stage_id=A"
            "; success=true"
            "; alpha_stage_report_present=true"
            "; alpha_geometry_epoch=post_H_committed_ALE_geometry"
            "; alpha_state_epoch=post_H_committed"
            "; alpha_thermal_state_epoch=post_R_committed"
            "; alpha_transport_model=atzeni_one_group"
            "; reactivity_model_executed=bosch_hale_dt"
            "; coefficient_time_level=old_time_lagged_to_A_stage_entry"
            "; fuel_depletion_enabled=false"
            "; separate_dt_species_authoritative=false"
            "; alpha_updated_fields=none"
            "; canonical_state_mutated=false"
            "; metadata_written=false"
            "; alpha_backend_executed=none"
            "; owned_slab_writeback_only=true"
            "; rank0_gather_solve_used=false"
            "; serial_dense_fallback_used=false"
            "; parity_claim_allowed=false";
        return stage;
      };
      const auto result = ExecuteP4ProductionStep(P4ProductionOptions{0.0}, hooks);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(Contains(result.report_line, "alpha_updated_fields=none"));
      DEC3D_CHECK(Contains(result.report_line, "advance_to_next_timestep=true"));
      DEC3D_CHECK(ValidateP4ProductionStepDiagnostics(result));
    }

    {
      auto runtime = CreateRuntimeScaffold(PhaseId::p4, "p4-v0.4");
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      const auto result = ExecuteP4ProductionStepIntoRuntime(runtime, P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(IsStageRegistered(runtime, StageId::hydro));
      DEC3D_CHECK(IsStageRegistered(runtime, StageId::thermal));
      DEC3D_CHECK(IsStageRegistered(runtime, StageId::equilibration));
      DEC3D_CHECK(IsStageRegistered(runtime, StageId::radiation));
      DEC3D_CHECK(IsStageRegistered(runtime, StageId::alpha));
      DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::hydro));
      DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::thermal));
      DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::equilibration));
      DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::radiation));
      DEC3D_CHECK(HasIngestedStageResult(runtime, StageId::alpha));
      DEC3D_CHECK(Contains(result.report_line, "runtime_scaffold_registered=true"));
      DEC3D_CHECK(Contains(result.report_line, "runtime_scaffold_ingested=true"));
    }

    {
      auto runtime = CreateRuntimeScaffold(PhaseId::p3, "p3-v0.4");
      std::vector<std::string> order;
      auto hooks = MakeSuccessfulHooks(order, kDt);
      const auto result = ExecuteP4ProductionStepIntoRuntime(runtime, P4ProductionOptions{kDt}, hooks);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "runtime scaffold must be initialized for P4"));
    }

    {
      dec3d::runtime::P3ProductionOptions p3_options;
      p3_options.dt_s = kDt;
      p3_options.register_a_stage = true;
      dec3d::runtime::P3ProductionOperatorHooks p3_hooks;
      p3_hooks.hydro = [&](double step_dt) {
        auto stage = dec3d::runtime::MakeP3ProductionStageResultForTest("H", step_dt, true);
        stage.geometry_epoch = "post_H_committed_ALE_geometry";
        return stage;
      };
      p3_hooks.thermal = [&](double step_dt) {
        auto stage = dec3d::runtime::MakeP3ProductionStageResultForTest("T", step_dt, true);
        stage.geometry_epoch = "post_H_committed_ALE_geometry";
        return stage;
      };
      p3_hooks.equilibration = [&](double step_dt) {
        return dec3d::runtime::MakeP3ProductionStageResultForTest("E", step_dt, true);
      };
      p3_hooks.radiation = [&](double step_dt) {
        auto stage = dec3d::runtime::MakeP3ProductionStageResultForTest("R", step_dt, true);
        stage.geometry_epoch = "post_H_committed_ALE_geometry";
        stage.radiation_groups_epoch = "post_H_committed";
        stage.electron_thermal_state_epoch = "post_E_committed";
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
        return stage;
      };
      const auto p3_result = dec3d::runtime::ExecuteP3ProductionStep(p3_options, p3_hooks);
      DEC3D_CHECK(!p3_result.success);
      DEC3D_CHECK(Contains(p3_result.failure_diagnostics, "forbidden alpha registration"));
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
