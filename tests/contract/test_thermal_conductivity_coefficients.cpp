#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "test_assert.hpp"
#include "thermal_conduction_test_support.hpp"
#include "transport/thermal/thermal_conductivity.hpp"

#include <cmath>
#include <iostream>
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

dec3d::transport::ThermalConductivityInput BuildHotSpotInput(
    dec3d::transport::ThermalConductivityModel model) {
  using dec3d::physics::DefaultMeanDTIonMassG;
  using dec3d::physics::DefaultZbar;
  using dec3d::physics::ErgFromKeV;

  dec3d::transport::ThermalConductivityInput input;
  input.Te_erg = ErgFromKeV(5.0);
  input.Ti_erg = ErgFromKeV(4.0);
  input.ne_cm3 = 1.0e25;
  input.ni_cm3 = 1.0e25;
  input.zbar = DefaultZbar();
  input.mean_ion_mass_g = DefaultMeanDTIonMassG();
  input.model = model;
  return input;
}

double SpitzerInternalOracle(const dec3d::transport::ThermalConductivityResult& result) {
  using dec3d::physics::PhysicsConstantsCGS;
  const double c0 = 20.0 * std::pow(2.0 / PhysicsConstantsCGS::pi, 1.5);
  return c0 *
         std::pow(result.Te_erg, 2.5) /
         (std::sqrt(PhysicsConstantsCGS::electron_mass_g) *
          result.zbar *
          std::pow(PhysicsConstantsCGS::elementary_charge_statcoulomb, 4.0) *
          result.lnLambda) *
         result.delta_prefactor *
         result.f_LM;
}

double WrongDoubleKBOracle(const dec3d::transport::ThermalConductivityResult& result) {
  using dec3d::physics::PhysicsConstantsCGS;
  const double c0 = 20.0 * std::pow(2.0 / PhysicsConstantsCGS::pi, 1.5);
  return c0 *
         std::pow(PhysicsConstantsCGS::boltzmann_erg_per_k * result.Te_erg, 2.5) /
         (std::sqrt(PhysicsConstantsCGS::electron_mass_g) *
          result.zbar *
          std::pow(PhysicsConstantsCGS::elementary_charge_statcoulomb, 4.0) *
          result.lnLambda) *
         result.delta_prefactor *
         result.f_LM;
}

}  // namespace

int main() {
  try {
    using dec3d::physics::PhysicsConstantsCGS;
    using dec3d::transport::ComputeThermalConductivity;
    using dec3d::transport::ThermalConductivityModel;
    using dec3d::transport::ValidateThermalConductivityDiagnostics;

    {
      const auto result = ComputeThermalConductivity(
          BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy));
      const auto values = dec3d::transport::ComputeThermalConductivityValues(
          BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(values.success);
      DEC3D_CHECK(result.is_complete());
      DEC3D_CHECK(ValidateThermalConductivityDiagnostics(result));
      DEC3D_CHECK(result.kappa_e_cm_inv_s > 0.0);
      DEC3D_CHECK(result.kappa_i_cm_inv_s > 0.0);
      CheckNear(
          values.kappa_e_cm_inv_s,
          result.kappa_e_cm_inv_s,
          result.kappa_e_cm_inv_s * 1.0e-15,
          "Spitzer value-only electron kappa parity");
      CheckNear(
          values.kappa_i_cm_inv_s,
          result.kappa_i_cm_inv_s,
          result.kappa_i_cm_inv_s * 1.0e-15,
          "Spitzer value-only ion kappa parity");
      DEC3D_CHECK(result.lnLambda > 0.0);
      DEC3D_CHECK(result.report_line.find("model_executed=spitzer_no_degeneracy") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("fallback_used=false") != std::string::npos);
      CheckNear(result.f_LM, 1.0, 0.0, "Spitzer f_LM");
      CheckNear(
          result.kappa_i_cm_inv_s,
          result.kappa_e_cm_inv_s *
              std::sqrt(PhysicsConstantsCGS::electron_mass_g / result.mean_ion_mass_g),
          result.kappa_i_cm_inv_s * 1.0e-14,
          "Spitzer ion mass-ratio rule");
      CheckNear(
          result.kappa_e_cm_inv_s,
          SpitzerInternalOracle(result),
          result.kappa_e_cm_inv_s * 1.0e-14,
          "Spitzer internal kappa oracle");
      const double wrong = WrongDoubleKBOracle(result);
      DEC3D_CHECK(std::abs(result.kappa_e_cm_inv_s - wrong) >
                  std::abs(result.kappa_e_cm_inv_s) * 0.999);
    }

    {
      const auto result = ComputeThermalConductivity(
          BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy));
      const auto values = dec3d::transport::ComputeThermalConductivityValues(
          BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(values.success);
      DEC3D_CHECK(result.is_complete());
      DEC3D_CHECK(ValidateThermalConductivityDiagnostics(result));
      DEC3D_CHECK(result.kappa_e_cm_inv_s > 0.0);
      DEC3D_CHECK(result.kappa_i_cm_inv_s > 0.0);
      CheckNear(
          values.kappa_e_cm_inv_s,
          result.kappa_e_cm_inv_s,
          result.kappa_e_cm_inv_s * 1.0e-15,
          "Lee-More value-only electron kappa parity");
      CheckNear(
          values.kappa_i_cm_inv_s,
          result.kappa_i_cm_inv_s,
          result.kappa_i_cm_inv_s * 1.0e-15,
          "Lee-More value-only ion kappa parity");
      DEC3D_CHECK(result.lnLambda >= 2.0);
      DEC3D_CHECK(result.f_LM >= 1.0);
      DEC3D_CHECK(result.report_line.find("model_executed=lee_more_with_degeneracy") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("fallback_used=false") != std::string::npos);
      CheckNear(
          result.kappa_i_cm_inv_s,
          result.kappa_e_cm_inv_s *
              std::sqrt(PhysicsConstantsCGS::electron_mass_g / result.mean_ion_mass_g),
          result.kappa_i_cm_inv_s * 1.0e-14,
          "Lee-More ion mass-ratio rule");
    }

    {
      const auto result = ComputeThermalConductivity(
          BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy));
      const double expected_TF =
          (PhysicsConstantsCGS::hbar_erg_s * PhysicsConstantsCGS::hbar_erg_s) /
          (2.0 * PhysicsConstantsCGS::electron_mass_g) *
          std::pow(3.0 * PhysicsConstantsCGS::pi * PhysicsConstantsCGS::pi *
                       result.ne_cm3,
                   2.0 / 3.0);
      CheckNear(result.TF_erg, expected_TF, expected_TF * 1.0e-14, "Lee-More Fermi temperature");
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy);
      input.Te_erg = dec3d::physics::ErgFromKeV(1000.0);
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(result.success);
      CheckNear(result.f_LM, 1.0, 1.0e-7, "Lee-More hot-limit f_LM");
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy);
      input.ne_cm3 = 1.0e31;
      input.ni_cm3 = 1.0e31;
      input.Te_erg = dec3d::physics::ErgFromKeV(0.05);
      input.Ti_erg = dec3d::physics::ErgFromKeV(0.05);
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(result.lnLambda_raw < 2.0);
      CheckNear(result.lnLambda, 2.0, 0.0, "Lee-More Coulomb floor");
      DEC3D_CHECK(result.lee_more_floor_active);
      DEC3D_CHECK(result.report_line.find("lee_more_floor_active=true") != std::string::npos);
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy);
      input.ne_cm3 = 1.0e27;
      input.ni_cm3 = 1.0e27;
      input.Te_erg = dec3d::physics::ErgFromKeV(0.5);
      input.Ti_erg = dec3d::physics::ErgFromKeV(0.5);
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(result.success);
      const double expected_bmax =
          1.0 /
          std::sqrt(
              4.0 * PhysicsConstantsCGS::pi * input.ne_cm3 *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb /
                  input.Te_erg +
              4.0 * PhysicsConstantsCGS::pi * input.ni_cm3 * input.zbar * input.zbar *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb /
                  input.Ti_erg);
      const double wrongly_degenerate_bmax =
          1.0 /
          std::sqrt(
              4.0 * PhysicsConstantsCGS::pi * input.ne_cm3 *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb /
                  std::sqrt(input.Te_erg * input.Te_erg + result.TF_erg * result.TF_erg) +
              4.0 * PhysicsConstantsCGS::pi * input.ni_cm3 * input.zbar * input.zbar *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb *
                  PhysicsConstantsCGS::elementary_charge_statcoulomb /
                  input.Ti_erg);
      CheckNear(result.bmax_cm, expected_bmax, expected_bmax * 1.0e-14, "Spitzer nondegenerate bmax");
      DEC3D_CHECK(std::abs(result.bmax_cm - wrongly_degenerate_bmax) >
                  std::abs(result.bmax_cm) * 1.0e-3);
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy);
      input.Te_erg = 0.0;
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("Te_erg must be positive") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("fallback_used=false") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") != std::string::npos);
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy);
      input.ne_cm3 = 0.0;
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("ne_cm3 must be positive") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("fallback_used=false") != std::string::npos);
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::unsupported);
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("unsupported conductivity model") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("fallback_used=false") != std::string::npos);
    }

    {
      auto input = BuildHotSpotInput(ThermalConductivityModel::spitzer_no_degeneracy);
      input.ne_cm3 = 1.0e40;
      input.ni_cm3 = 1.0e40;
      input.Te_erg = dec3d::physics::ErgFromKeV(0.001);
      input.Ti_erg = dec3d::physics::ErgFromKeV(0.001);
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.failure_diagnostics.find("Spitzer Coulomb logarithm requires bmax greater than bmin") != std::string::npos ||
                  result.failure_diagnostics.find("lnLambda must be positive and finite") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("fallback_used=false") != std::string::npos);
    }

    {
      auto result = ComputeThermalConductivity(
          BuildHotSpotInput(ThermalConductivityModel::lee_more_with_degeneracy));
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(ValidateThermalConductivityDiagnostics(result));
      const auto pos = result.report_line.find("ion_mass_ratio_rule=");
      DEC3D_CHECK(pos != std::string::npos);
      result.report_line.erase(pos, std::string("ion_mass_ratio_rule=").size());
      DEC3D_CHECK(!ValidateThermalConductivityDiagnostics(result));
    }

    {
      using dec3d::test_support::BuildThermalState;
      using dec3d::test_support::Layout;
      const auto layout = Layout(1, 1, 1);
      auto state = BuildThermalState(layout, 5.0, 5.0, 4.0, 4.0);
      const auto before = state;
      const auto recovered = dec3d::state::RecoverThermodynamicState(state);
      DEC3D_CHECK(recovered.success);
      const auto& cell = recovered.cells(0, 0, 0);
      dec3d::transport::ThermalConductivityInput input;
      input.Te_erg = cell.t_e_erg_per_particle;
      input.Ti_erg = cell.t_i_erg_per_particle;
      input.ne_cm3 = cell.n_e_cm3;
      input.ni_cm3 = cell.n_i_cm3;
      input.zbar = cell.zbar;
      input.mean_ion_mass_g = cell.mean_ion_mass_g;
      input.model = ThermalConductivityModel::spitzer_no_degeneracy;
      const auto result = ComputeThermalConductivity(input);
      DEC3D_CHECK(result.success);
      DEC3D_CHECK_EQ(state.rho(0, 0, 0), before.rho(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before.e_electron(0, 0, 0));
      DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before.e_fluid_total(0, 0, 0));
      DEC3D_CHECK_EQ(state.last_authoritative_write_mask, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
      DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
    }

    std::cout << "P2-4 thermal conductivity coefficient contract passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
