#include "core/diagnostics/stage_contracts.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/electron_ion_equilibration.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

void CheckNear(double lhs, double rhs, double tolerance, const std::string& label) {
  if (std::abs(lhs - rhs) > tolerance) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("near equality", __FILE__, __LINE__, detail.str());
  }
}

dec3d::state::CanonicalState BuildOneCellThermalState(double te_kev, double ti_kev) {
  using dec3d::physics::DefaultMeanDTIonMassG;
  using dec3d::physics::ErgFromKeV;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;

  auto state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 0});
  const double rho = 2.0 * DefaultMeanDTIonMassG();
  const double n_e = 2.0;
  const double n_i = 2.0;
  const double te = ErgFromKeV(te_kev);
  const double ti = ErgFromKeV(ti_kev);
  const double pe = n_e * te;
  const double pi = n_i * ti;
  const double e_e = pe / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double e_i = pi / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double mom_r = rho * 1.0e5;
  const double mom_theta = rho * 2.0e5;
  const double mom_phi = rho * 3.0e5;
  const double kinetic =
      0.5 * (mom_r * mom_r + mom_theta * mom_theta + mom_phi * mom_phi) / rho;

  state.rho(0, 0, 0) = rho;
  state.mom_r(0, 0, 0) = mom_r;
  state.mom_theta(0, 0, 0) = mom_theta;
  state.mom_phi(0, 0, 0) = mom_phi;
  state.e_electron(0, 0, 0) = e_e;
  state.e_fluid_total(0, 0, 0) = kinetic + e_e + e_i;
  return state;
}

double ThermalTotal(const dec3d::state::CanonicalState& state) {
  const auto recovered = dec3d::state::RecoverThermodynamicState(state);
  DEC3D_CHECK(recovered.success);
  const auto& cell = recovered.cells(0, 0, 0);
  return cell.e_electron_erg_per_cm3 + cell.e_ion_erg_per_cm3;
}

}  // namespace

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::MaskContains;
    using dec3d::physics::ErgFromKeV;
    using dec3d::state::ApplyLocalElectronIonEquilibration;
    using dec3d::state::ElectronIonEquilibrationOptions;
    using dec3d::state::ElectronIonTauModel;
    using dec3d::state::RecoverThermodynamicState;
    using dec3d::state::ValidateElectronIonEquilibrationDiagnostics;

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const auto before_e_e = state.e_electron(0, 0, 0);
      const auto before_e_total = state.e_fluid_total(0, 0, 0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.0, ElectronIonTauModel::constant_user_supplied, 2.0});
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.is_complete());
      DEC3D_CHECK(ValidateElectronIonEquilibrationDiagnostics(result));
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before_e_e);
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before_e_total);
      CheckNear(result.max_dt_over_tau, 0.0, 0.0, "zero dt/tau");
      CheckNear(result.max_expected_decay_residual_erg, 0.0, 0.0, "zero-dt residual");
      DEC3D_CHECK(result.report_line.find("e_fluid_total_unchanged=true") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const double before_e_total = state.e_fluid_total(0, 0, 0);
      const double before_thermal_total = ThermalTotal(state);
      const double dt_s = 0.25;
      const double tau_s = 2.0;
      const double decay = std::exp(-2.0 * dt_s / tau_s);
      const double te_before = ErgFromKeV(1.0);
      const double ti_before = ErgFromKeV(3.0);
      const double teq = 0.5 * (te_before + ti_before);
      const double expected_te = teq + 0.5 * (te_before - ti_before) * decay;

      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{dt_s, ElectronIonTauModel::constant_user_supplied, tau_s});
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(MaskContains(result.updated_fields, AuthoritativeField::e_electron));
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::core::ToMask(AuthoritativeField::e_electron));
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before_e_total);

      const auto recovered = RecoverThermodynamicState(state);
      DEC3D_CHECK(recovered.success);
      const auto& cell = recovered.cells(0, 0, 0);
      CheckNear(cell.t_e_erg_per_particle, expected_te, 1.0e-24, "electron temperature after");
      CheckNear(
          cell.t_e_erg_per_particle - cell.t_i_erg_per_particle,
          (te_before - ti_before) * decay,
          1.0e-24,
          "delta T decay");
      CheckNear(ThermalTotal(state), before_thermal_total, 1.0e-24, "thermal conservation");
      CheckNear(result.max_dt_over_tau, dt_s / tau_s, 1.0e-15, "dt/tau diagnostic");
      CheckNear(result.max_abs_ne_minus_ni, 0.0, 1.0e-14, "ne-ni diagnostic");
      DEC3D_CHECK(result.report_line.find("tau_model=constant_user_supplied") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("updated_fields=e_electron") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("max_abs_ne_minus_ni=") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("max_dt_over_tau=") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const double before_e_e = state.e_electron(0, 0, 0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.1, ElectronIonTauModel::constant_user_supplied, -1.0});
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before_e_e);
      DEC3D_CHECK(result.failure_diagnostics.find("tau_ei_s must be positive") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const double before_e_e = state.e_electron(0, 0, 0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{-0.1, ElectronIonTauModel::constant_user_supplied, 1.0});
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before_e_e);
      DEC3D_CHECK(result.failure_diagnostics.find("dt_s must be non-negative") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const double before_e_e = state.e_electron(0, 0, 0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{
              0.1,
              ElectronIonTauModel::unsupported,
              1.0});
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before_e_e);
      DEC3D_CHECK(result.failure_diagnostics.find("unsupported tau_model") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      state.rho(0, 0, 0) = 0.0;
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.1, ElectronIonTauModel::constant_user_supplied, 1.0});
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("p2_0_recovery_failed=true") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("P2-0 thermodynamic recovery failed") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("rho must be positive") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("first_bad_cell=0,0,0") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{
              1.0e-12,
              ElectronIonTauModel::thesis_spitzer_eq_5_241,
              0.0});
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.report_line.find("diagnostic_id=p2.production.equilibration_stage") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("stage_id=E") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("tau_model_requested=thesis_spitzer_eq_5_241") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("tau_model_executed=thesis_spitzer_eq_5_241") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("min_tau_ei_s=") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("min_lnLambda_ei_spitzer=") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("max_lnLambda_ei_spitzer=") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("double_kB_guard_passed=true") != std::string::npos);
      DEC3D_CHECK(result.updated_fields == dec3d::core::ToMask(AuthoritativeField::e_electron));
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{
              0.1,
              ElectronIonTauModel::constant_user_supplied,
              std::numeric_limits<double>::infinity()});
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("tau_ei_s must be finite") != std::string::npos);
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.1, ElectronIonTauModel::constant_user_supplied, 1.0});
      DEC3D_CHECK(result.success);
      result.report_line.clear();
      DEC3D_CHECK(!result.is_complete());
      DEC3D_CHECK(!ValidateElectronIonEquilibrationDiagnostics(result));
    }

    {
      auto state = BuildOneCellThermalState(1.0, 3.0);
      auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.1, ElectronIonTauModel::constant_user_supplied, 1.0});
      DEC3D_CHECK(result.success);
      const auto ne_ni_pos = result.report_line.find("max_abs_ne_minus_ni=");
      DEC3D_CHECK(ne_ni_pos != std::string::npos);
      result.report_line.erase(ne_ni_pos, std::string("max_abs_ne_minus_ni=").size());
      DEC3D_CHECK(!ValidateElectronIonEquilibrationDiagnostics(result));

      result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.1, ElectronIonTauModel::constant_user_supplied, 1.0});
      DEC3D_CHECK(result.success);
      const auto dt_tau_pos = result.report_line.find("max_dt_over_tau=");
      DEC3D_CHECK(dt_tau_pos != std::string::npos);
      result.report_line.erase(dt_tau_pos, std::string("max_dt_over_tau=").size());
      DEC3D_CHECK(!ValidateElectronIonEquilibrationDiagnostics(result));
    }

    {
      auto state = dec3d::state::CanonicalState::Create(
          dec3d::state::CanonicalStateLayout{2, 1, 1, 0});
      auto cell_a = BuildOneCellThermalState(1.0, 3.0);
      auto cell_b = BuildOneCellThermalState(4.0, 2.0);
      state.rho(0, 0, 0) = cell_a.rho(0, 0, 0);
      state.mom_r(0, 0, 0) = cell_a.mom_r(0, 0, 0);
      state.mom_theta(0, 0, 0) = cell_a.mom_theta(0, 0, 0);
      state.mom_phi(0, 0, 0) = cell_a.mom_phi(0, 0, 0);
      state.e_electron(0, 0, 0) = cell_a.e_electron(0, 0, 0);
      state.e_fluid_total(0, 0, 0) = cell_a.e_fluid_total(0, 0, 0);
      state.rho(1, 0, 0) = cell_b.rho(0, 0, 0);
      state.mom_r(1, 0, 0) = cell_b.mom_r(0, 0, 0);
      state.mom_theta(1, 0, 0) = cell_b.mom_theta(0, 0, 0);
      state.mom_phi(1, 0, 0) = cell_b.mom_phi(0, 0, 0);
      state.e_electron(1, 0, 0) = cell_b.e_electron(0, 0, 0);
      state.e_fluid_total(1, 0, 0) = cell_b.e_fluid_total(0, 0, 0);

      const auto result = ApplyLocalElectronIonEquilibration(
          state,
          ElectronIonEquilibrationOptions{0.5, ElectronIonTauModel::constant_user_supplied, 2.0});
      DEC3D_CHECK(result.success);
      DEC3D_CHECK_EQ(result.cell_count, static_cast<std::size_t>(2));
      DEC3D_CHECK(result.max_deltaT_before_erg > 0.0);
      DEC3D_CHECK(result.max_deltaT_after_erg > 0.0);
      CheckNear(result.max_abs_ne_minus_ni, 0.0, 1.0e-14, "multi-cell ne-ni");
      DEC3D_CHECK(result.report_line.find("cell_count=2") != std::string::npos);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
