#include "core/diagnostics/stage_contracts.hpp"
#include "test_assert.hpp"
#include "thermal_conduction_test_support.hpp"
#include "transport/thermal/thermal_conduction.hpp"

#include <iostream>
#include <string>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::MaskContains;
    using dec3d::transport::ApplyConstantKappaThermalConduction;
    using dec3d::transport::ThermalConductionBackend;
    using dec3d::transport::ThermalConductionKappaModel;
    using dec3d::transport::ThermalConductionOptions;
    using dec3d::transport::ValidateThermalConductionDiagnostics;

    using dec3d::test_support::BuildPatchGeometry;
    using dec3d::test_support::BuildThermalState;
    using dec3d::test_support::CheckNear;
    using dec3d::test_support::Layout;
    using dec3d::test_support::PatchBoundary;
    using dec3d::test_support::VolumeWeightedThermalTotal;

    const auto layout = Layout(2, 1, 1);
    const auto geometry = BuildPatchGeometry(layout, 1.0, 3.0, 0.5, 1.2);

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 3.0, 2.0);
      const auto before = state;
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.0;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 1.0;
      options.kappa_i_cm_inv_s = 1.0;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.is_complete());
      DEC3D_CHECK(ValidateThermalConductionDiagnostics(result));
      CheckNear(result.max_abs_delta_Te_erg, 0.0, 0.0, "zero dt Te delta");
      CheckNear(result.max_abs_delta_Ti_erg, 0.0, 0.0, "zero dt Ti delta");
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0));
      DEC3D_CHECK_EQ(state.last_authoritative_write_mask, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
      DEC3D_CHECK(result.report_line.find("metadata_written=false") != std::string::npos);
    }

    {
      const auto one_cell_layout = Layout(1, 1, 1);
      const auto one_cell_geometry = BuildPatchGeometry(one_cell_layout, 1.0, 2.0, 0.5, 1.2);
      auto state = BuildThermalState(one_cell_layout, 5.0, 5.0, 2.0, 2.0);
      const auto before = state;
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.7;
      options.kappa_i_cm_inv_s = 0.0;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, one_cell_geometry, options);
      DEC3D_CHECK(result.success);
      CheckNear(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0), 1.0e-20, "one-cell electron energy unchanged");
      CheckNear(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0), 1.0e-20, "one-cell electron total energy unchanged");
      CheckNear(result.max_abs_delta_Te_erg, 0.0, 1.0e-30, "one-cell electron no-source Te delta");
    }

    {
      const auto one_cell_layout = Layout(1, 1, 1);
      const auto one_cell_geometry = BuildPatchGeometry(one_cell_layout, 1.0, 2.0, 0.5, 1.2);
      auto state = BuildThermalState(one_cell_layout, 2.0, 2.0, 5.0, 5.0);
      const auto before = state;
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.0;
      options.kappa_i_cm_inv_s = 0.7;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, one_cell_geometry, options);
      DEC3D_CHECK(result.success);
      CheckNear(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0), 1.0e-20, "one-cell ion electron energy unchanged");
      CheckNear(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0), 1.0e-20, "one-cell ion total energy unchanged");
      CheckNear(result.max_abs_delta_Ti_erg, 0.0, 1.0e-30, "one-cell ion no-source Ti delta");
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 3.0, 3.0);
      const double before_total = VolumeWeightedThermalTotal(state, geometry);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.0;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(result.success);
      const auto recovered = dec3d::state::RecoverThermodynamicState(state);
      DEC3D_CHECK(recovered.success);
      DEC3D_CHECK(recovered.cells(0, 0, 0).t_e_erg_per_particle <
                  dec3d::physics::ErgFromKeV(5.0));
      DEC3D_CHECK(recovered.cells(1, 0, 0).t_e_erg_per_particle >
                  dec3d::physics::ErgFromKeV(1.0));
      CheckNear(VolumeWeightedThermalTotal(state, geometry), before_total, 1.0e-20, "electron-only global thermal conservation");
      DEC3D_CHECK(result.report_line.find("electron_equation_executed=true") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("ion_equation_executed=true") != std::string::npos);
    }

    {
      auto state = BuildThermalState(layout, 3.0, 3.0, 5.0, 1.0);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.0;
      options.kappa_i_cm_inv_s = 0.2;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(result.success);
      const auto recovered = dec3d::state::RecoverThermodynamicState(state);
      DEC3D_CHECK(recovered.success);
      DEC3D_CHECK(recovered.cells(0, 0, 0).t_i_erg_per_particle <
                  dec3d::physics::ErgFromKeV(5.0));
      DEC3D_CHECK(recovered.cells(1, 0, 0).t_i_erg_per_particle >
                  dec3d::physics::ErgFromKeV(1.0));
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      const double before_total = VolumeWeightedThermalTotal(state, geometry);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK_EQ(
          result.updated_fields,
          dec3d::core::ToMask(AuthoritativeField::e_electron) |
              dec3d::core::ToMask(AuthoritativeField::e_fluid_total));
      DEC3D_CHECK(MaskContains(result.updated_fields, AuthoritativeField::e_electron));
      DEC3D_CHECK(MaskContains(result.updated_fields, AuthoritativeField::e_fluid_total));
      DEC3D_CHECK_EQ(state.last_authoritative_write_mask, result.updated_fields);
      DEC3D_CHECK(result.report_line.find("rho_momentum_unchanged=true") != std::string::npos);
      CheckNear(VolumeWeightedThermalTotal(state, geometry), before_total, 1.0e-20, "dual global thermal conservation");
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      const auto before = state;
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = -1.0;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0));
      DEC3D_CHECK(result.failure_diagnostics.find("kappa_e_cm_inv_s must be non-negative") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") != std::string::npos);
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(ValidateThermalConductionDiagnostics(result));
      const auto backend_requested_pos = result.report_line.find("backend_requested=");
      DEC3D_CHECK(backend_requested_pos != std::string::npos);
      result.report_line.erase(backend_requested_pos, std::string("backend_requested=").size());
      DEC3D_CHECK(!ValidateThermalConductionDiagnostics(result));
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::unsupported;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("unsupported thermal conduction backend") != std::string::npos);
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      state.rho(0, 0, 0) = 0.0;
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::serial_dense_reference;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("P2-0 thermodynamic recovery failed") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("rho must be positive") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") != std::string::npos);
    }

    {
      auto state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
      auto options = ThermalConductionOptions{};
      options.dt_s = 0.01;
      options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
      options.kappa_e_cm_inv_s = 0.2;
      options.kappa_i_cm_inv_s = 0.3;
      options.backend = ThermalConductionBackend::hypre_parcsr_gmres_boomeramg;
      options.boundary_policy = PatchBoundary();
      const auto result = ApplyConstantKappaThermalConduction(state, geometry, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("requested backend is unavailable without a solve hook") != std::string::npos);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
