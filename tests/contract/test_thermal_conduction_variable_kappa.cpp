#include "core/diagnostics/stage_contracts.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "test_assert.hpp"
#include "thermal_conduction_test_support.hpp"
#include "transport/thermal/thermal_conduction.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace {

bool Contains(const std::string& haystack, const std::string& needle) {
  return haystack.find(needle) != std::string::npos;
}

dec3d::transport::ThermalConductionOptions VariableOptions(
    dec3d::transport::ThermalConductionKappaModel model,
    double dt_s) {
  dec3d::transport::ThermalConductionOptions options;
  options.dt_s = dt_s;
  options.kappa_model = model;
  options.backend = dec3d::transport::ThermalConductionBackend::serial_dense_reference;
  options.boundary_policy = dec3d::test_support::PatchBoundary();
  options.serial_reference_options.residual_tolerance = 1.0e-2;
  return options;
}

}  // namespace

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::MaskContains;
    using dec3d::transport::ApplyThermalConduction;
    using dec3d::transport::ThermalConductionKappaModel;
    using dec3d::transport::ValidateThermalConductionDiagnostics;
    using dec3d::test_support::BuildPatchGeometry;
    using dec3d::test_support::BuildThermalState;
    using dec3d::test_support::CheckNear;
    using dec3d::test_support::Layout;
    using dec3d::test_support::VolumeWeightedThermalTotal;

    const auto layout = Layout(2, 1, 1);
    const auto geometry = BuildPatchGeometry(layout, 1.0, 3.0, 0.5, 1.2);

    {
      auto state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      const auto before = state;
      const auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::spitzer_no_degeneracy, 0.0));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.is_complete());
      DEC3D_CHECK(ValidateThermalConductionDiagnostics(result));
      DEC3D_CHECK(Contains(result.report_line, "diagnostic_id=p2.thermal_conduction.variable_kappa"));
      DEC3D_CHECK(Contains(result.report_line, "coefficient_time_level=old_time_lagged"));
      DEC3D_CHECK(Contains(result.report_line, "provider_failure_count=0"));
      DEC3D_CHECK(Contains(result.report_line, "kappa_model_executed=spitzer_no_degeneracy"));
      DEC3D_CHECK(result.min_kappa_e_cm_inv_s > 0.0);
      DEC3D_CHECK(result.max_kappa_e_cm_inv_s >= result.min_kappa_e_cm_inv_s);
      DEC3D_CHECK(result.min_kappa_i_cm_inv_s > 0.0);
      DEC3D_CHECK(result.max_kappa_i_cm_inv_s >= result.min_kappa_i_cm_inv_s);
      CheckNear(result.max_abs_delta_Te_erg, 0.0, 0.0, "zero-dt variable Te delta");
      CheckNear(result.max_abs_delta_Ti_erg, 0.0, 0.0, "zero-dt variable Ti delta");
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0));
      DEC3D_CHECK_EQ(state.last_authoritative_write_mask, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
      DEC3D_CHECK(Contains(result.report_line, "canonical_state_mutated=false"));
    }

    {
      auto state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      const double before_total = VolumeWeightedThermalTotal(state, geometry);
      const auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::spitzer_no_degeneracy, 1.0e-16));
      if (!result.success) {
        std::cerr << result.failure_diagnostics << '\n';
      }
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(Contains(result.report_line, "fallback_used=false"));
      DEC3D_CHECK(MaskContains(result.updated_fields, AuthoritativeField::e_electron));
      DEC3D_CHECK(MaskContains(result.updated_fields, AuthoritativeField::e_fluid_total));
      DEC3D_CHECK(result.max_abs_delta_Te_erg > 0.0);
      DEC3D_CHECK(result.max_abs_delta_Ti_erg > 0.0);
      CheckNear(
          VolumeWeightedThermalTotal(state, geometry),
          before_total,
          1.0e-18,
          "variable zero-flux global thermal conservation");
    }

    {
      auto state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      const auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::lee_more_with_degeneracy, 1.0e-16));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(Contains(result.report_line, "kappa_model_executed=lee_more_with_degeneracy"));
      DEC3D_CHECK(Contains(result.report_line, "provider_model_executed=lee_more_with_degeneracy"));
      DEC3D_CHECK(result.min_provider_lnLambda >= 2.0);
      DEC3D_CHECK(result.min_provider_f_LM >= 1.0);
    }

    {
      auto state = BuildThermalState(layout, 0.001, 0.001, 0.001, 0.001);
      const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
      const double n = 1.0e40;
      const double rho = n * dec3d::physics::DefaultMeanDTIonMassG();
      const double te = dec3d::physics::ErgFromKeV(0.001);
      for (std::size_t r = 0; r < layout.radial_cells; ++r) {
        state.rho(r, 0, 0) = rho;
        state.e_electron(r, 0, 0) = n * te / gamma_minus_one;
        state.e_fluid_total(r, 0, 0) = 2.0 * n * te / gamma_minus_one;
      }
      const auto before = state;
      const auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::spitzer_no_degeneracy, 1.0e-16));
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(Contains(result.failure_diagnostics, "canonical_state_mutated=false"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "diagnostic_id=p2.thermal_conduction.variable_kappa.failure"));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "provider_failure_count="));
      DEC3D_CHECK(Contains(result.failure_diagnostics, "first_provider_failure_cell="));
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_fluid_total(1, 0, 0), before.e_fluid_total(1, 0, 0));
      DEC3D_CHECK_EQ(state.last_authoritative_write_mask, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    }

    {
      const auto one_cell_layout = Layout(1, 1, 1);
      const auto one_cell_geometry = BuildPatchGeometry(one_cell_layout, 1.0, 2.0, 0.5, 1.2);
      auto state = BuildThermalState(one_cell_layout, 0.01, 0.01, 0.008, 0.008);
      const auto before = state;
      const auto result = ApplyThermalConduction(
          state,
          one_cell_geometry,
          VariableOptions(ThermalConductionKappaModel::spitzer_no_degeneracy, 1.0e-16));
      DEC3D_CHECK(result.success);
      CheckNear(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0), 1.0e-20, "one-cell variable e_e unchanged");
      CheckNear(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0), 1.0e-20, "one-cell variable e_total unchanged");
      CheckNear(result.max_abs_delta_Te_erg, 0.0, 1.0e-24, "one-cell variable Te no-change");
      CheckNear(result.max_abs_delta_Ti_erg, 0.0, 1.0e-24, "one-cell variable Ti no-change");
    }

    {
      using dec3d::transport::ApplyConstantKappaThermalConduction;
      using dec3d::transport::ThermalConductionBackend;
      using dec3d::transport::ThermalConductionOptions;

      auto via_general = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      auto via_wrapper = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      ThermalConductionOptions options;
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = dec3d::test_support::PatchBoundary();
      const auto general = ApplyThermalConduction(via_general, geometry, options);
      const auto wrapper = ApplyConstantKappaThermalConduction(via_wrapper, geometry, options);
      DEC3D_CHECK(general.success);
      DEC3D_CHECK(wrapper.success);
      for (std::size_t r = 0; r < layout.radial_cells; ++r) {
        CheckNear(via_general.e_electron(r, 0, 0), via_wrapper.e_electron(r, 0, 0), 0.0, "constant wrapper e_e parity");
        CheckNear(via_general.e_fluid_total(r, 0, 0), via_wrapper.e_fluid_total(r, 0, 0), 0.0, "constant wrapper e_total parity");
      }
    }

    {
      auto state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      const auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::lee_more_with_degeneracy, 1.0e-16));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK_EQ(result.kappa_model_requested, std::string("lee_more_with_degeneracy"));
      DEC3D_CHECK_EQ(result.kappa_model_executed, std::string("lee_more_with_degeneracy"));
      DEC3D_CHECK_EQ(result.provider_model_requested, std::string("lee_more_with_degeneracy"));
      DEC3D_CHECK_EQ(result.provider_model_executed, std::string("lee_more_with_degeneracy"));
      DEC3D_CHECK(Contains(result.report_line, "fallback_used=false"));
      DEC3D_CHECK(!Contains(result.report_line, "provider_model_executed=spitzer_no_degeneracy"));
    }

    {
      auto state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      auto result = ApplyThermalConduction(
          state,
          geometry,
          VariableOptions(ThermalConductionKappaModel::spitzer_no_degeneracy, 1.0e-16));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(ValidateThermalConductionDiagnostics(result));
      const auto token_pos = result.report_line.find("provider_failure_count=");
      DEC3D_CHECK(token_pos != std::string::npos);
      result.report_line.erase(token_pos, std::string("provider_failure_count=").size());
      DEC3D_CHECK(!ValidateThermalConductionDiagnostics(result));
    }

    std::cout << "P2-5 variable-kappa thermal conduction contract passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
