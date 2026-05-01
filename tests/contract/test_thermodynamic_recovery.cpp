#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "test_assert.hpp"

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

dec3d::state::CanonicalState BuildOneCellState() {
  using dec3d::physics::DefaultMeanDTIonMassG;
  using dec3d::physics::ErgFromKeV;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;

  auto state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 0});
  const double rho = 2.0 * DefaultMeanDTIonMassG();
  const double ne = 2.0;
  const double ni = 2.0;
  const double te = ErgFromKeV(1.0);
  const double ti = ErgFromKeV(2.0);
  const double pe = ne * te;
  const double pi = ni * ti;
  const double e_e = pe / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double e_i = pi / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double mom_r = rho * 3.0;
  const double mom_theta = rho * 4.0;
  const double mom_phi = rho * 12.0;
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

}  // namespace

int main() {
  try {
    using dec3d::physics::DefaultMeanDTIonMassG;
    using dec3d::physics::KeVFromErg;
    using dec3d::state::RecoverThermodynamicState;
    using dec3d::state::ThermodynamicRecoveryOptions;

    auto state = BuildOneCellState();
    const auto before_rho = state.rho(0, 0, 0);
    const auto before_e_total = state.e_fluid_total(0, 0, 0);
    const auto before_e_electron = state.e_electron(0, 0, 0);

    const auto result = RecoverThermodynamicState(state, ThermodynamicRecoveryOptions{});
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK_EQ(result.cell_count, static_cast<std::size_t>(1));
    DEC3D_CHECK_EQ(result.cells.extent_r(), static_cast<std::size_t>(1));

    const auto& cell = result.cells(0, 0, 0);
    CheckNear(cell.mean_ion_mass_g, DefaultMeanDTIonMassG(), 1.0e-36, "mean ion mass");
    CheckNear(cell.zbar, 1.0, 0.0, "zbar");
    CheckNear(cell.n_i_cm3, 2.0, 1.0e-14, "ion number density");
    CheckNear(cell.n_e_cm3, 2.0, 1.0e-14, "electron number density");
    CheckNear(KeVFromErg(cell.t_e_erg_per_particle), 1.0, 1.0e-13, "Te keV");
    CheckNear(KeVFromErg(cell.t_i_erg_per_particle), 2.0, 1.0e-13, "Ti keV");
    CheckNear(cell.p_e_erg_per_cm3, cell.n_e_cm3 * cell.t_e_erg_per_particle, 1.0e-20, "Pe=nT");
    CheckNear(cell.p_i_erg_per_cm3, cell.n_i_cm3 * cell.t_i_erg_per_particle, 1.0e-20, "Pi=nT");
    CheckNear(
        cell.e_electron_erg_per_cm3,
        cell.p_e_erg_per_cm3 / (dec3d::state::HydroIdealGasGamma() - 1.0),
        1.0e-20,
        "Ee=Pe/(gamma-1)");
    CheckNear(
        cell.e_ion_erg_per_cm3,
        cell.p_i_erg_per_cm3 / (dec3d::state::HydroIdealGasGamma() - 1.0),
        1.0e-20,
        "Ei=Pi/(gamma-1)");
    CheckNear(cell.energy_partition_residual, 0.0, 1.0e-30, "partition residual");

    DEC3D_CHECK_EQ(state.rho(0, 0, 0), before_rho);
    DEC3D_CHECK_EQ(state.e_fluid_total(0, 0, 0), before_e_total);
    DEC3D_CHECK_EQ(state.e_electron(0, 0, 0), before_e_electron);

    DEC3D_CHECK(result.unit_diagnostics.find("diagnostic_id=p2.units.contract") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("diagnostic_id=p2.thermo.recovery") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("implementation_id=p2.thermo.recovery.cgs_v1") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("temperature_internal_unit=erg_per_particle") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("min_ne_cm3=") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("max_ne_cm3=") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("min_ni_cm3=") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("max_ni_cm3=") != std::string::npos);
    DEC3D_CHECK(result.recovery_diagnostics.find("max_energy_partition_residual=") != std::string::npos);
    auto diagnostics_missing = result;
    diagnostics_missing.recovery_diagnostics.clear();
    DEC3D_CHECK(!diagnostics_missing.is_complete());

    diagnostics_missing = result;
    diagnostics_missing.unit_diagnostics.clear();
    DEC3D_CHECK(!diagnostics_missing.is_complete());

    DEC3D_CHECK(dec3d::state::ValidateThermodynamicRecoveryDiagnostics(result));
    DEC3D_CHECK(!dec3d::state::ValidateThermodynamicRecoveryDiagnostics(diagnostics_missing));

    auto bad_rho = BuildOneCellState();
    bad_rho.rho(0, 0, 0) = 0.0;
    const auto bad_rho_result = RecoverThermodynamicState(bad_rho, ThermodynamicRecoveryOptions{});
    DEC3D_CHECK(!bad_rho_result.success);
    DEC3D_CHECK(bad_rho_result.failure_diagnostics.find("failure_reason=rho must be positive") != std::string::npos);
    DEC3D_CHECK(bad_rho_result.failure_diagnostics.find("first_bad_cell=0,0,0") != std::string::npos);

    auto bad_electron = BuildOneCellState();
    bad_electron.e_electron(0, 0, 0) = -1.0e-20;
    const auto bad_electron_result =
        RecoverThermodynamicState(bad_electron, ThermodynamicRecoveryOptions{});
    DEC3D_CHECK(!bad_electron_result.success);
    DEC3D_CHECK(bad_electron_result.failure_diagnostics.find("electron energy fell below minimum threshold") != std::string::npos);

    auto bad_ion = BuildOneCellState();
    bad_ion.e_fluid_total(0, 0, 0) = bad_ion.e_electron(0, 0, 0);
    const auto bad_ion_result = RecoverThermodynamicState(bad_ion, ThermodynamicRecoveryOptions{});
    DEC3D_CHECK(!bad_ion_result.success);
    DEC3D_CHECK(bad_ion_result.failure_diagnostics.find("ion energy fell below minimum threshold") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
