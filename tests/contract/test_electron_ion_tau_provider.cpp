#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/electron_ion_tau.hpp"
#include "test_assert.hpp"

#include <cmath>
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

dec3d::state::ElectronIonTauInput HotSpotInput() {
  return dec3d::state::ElectronIonTauInput{
      dec3d::state::ElectronIonTauModel::thesis_spitzer_eq_5_241,
      dec3d::physics::ErgFromKeV(5.0),
      dec3d::physics::ErgFromKeV(0.5),
      1.2e25,
      1.2e25,
      dec3d::physics::DefaultZbar(),
      dec3d::physics::DefaultMeanDTIonMassG(),
      0.0};
}

}  // namespace

int main() {
  using dec3d::state::ComputeElectronIonTau;
  using dec3d::state::ElectronIonTauModel;
  using dec3d::state::ValidateElectronIonTauDiagnostics;

  {
    auto input = HotSpotInput();
    const auto result = ComputeElectronIonTau(input);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.tau_ei_s > 0.0);
    DEC3D_CHECK(result.lnLambda_ei_spitzer > 0.0);
    DEC3D_CHECK(result.double_kB_guard_passed);
    DEC3D_CHECK(result.report_line.find("tau_model_executed=thesis_spitzer_eq_5_241") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("temperature_internal_unit=erg_per_particle") !=
                std::string::npos);
    DEC3D_CHECK(ValidateElectronIonTauDiagnostics(result));
  }

  {
    auto input = HotSpotInput();
    const auto correct = ComputeElectronIonTau(input);
    DEC3D_CHECK(correct.success);
    input.Te_erg_per_particle *= dec3d::physics::PhysicsConstantsCGS::boltzmann_erg_per_k;
    const auto wrong_double_kb = ComputeElectronIonTau(input);
    DEC3D_CHECK(wrong_double_kb.success);
    DEC3D_CHECK(std::abs(wrong_double_kb.tau_ei_s - correct.tau_ei_s) >
                correct.tau_ei_s * 1.0e-3);
  }

  {
    auto input = HotSpotInput();
    input.model = ElectronIonTauModel::constant_user_supplied;
    input.constant_tau_ei_s = 3.0e-11;
    const auto result = ComputeElectronIonTau(input);
    DEC3D_CHECK(result.success);
    CheckNear(result.tau_ei_s, 3.0e-11, 0.0, "constant tau");
    DEC3D_CHECK(result.report_line.find("thesis_spitzer_tau_used=false") != std::string::npos);
  }

  {
    auto input = HotSpotInput();
    input.ne_cm3 = 0.0;
    const auto result = ComputeElectronIonTau(input);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_diagnostics.find("electron density must be positive") !=
                std::string::npos);
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);
  }

  {
    auto input = HotSpotInput();
    input.model = ElectronIonTauModel::unsupported;
    const auto result = ComputeElectronIonTau(input);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_diagnostics.find("unsupported tau_model") != std::string::npos);
  }

  return 0;
}
