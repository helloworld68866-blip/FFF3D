#include "physics/units/physical_constants.hpp"
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

void CheckGreater(double lhs, double rhs, const std::string& label) {
  if (!(lhs > rhs)) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("greater-than", __FILE__, __LINE__, detail.str());
  }
}

}  // namespace

int main() {
  try {
    using dec3d::physics::DefaultMeanDTIonMassG;
    using dec3d::physics::DefaultZbar;
    using dec3d::physics::ErgFromEV;
    using dec3d::physics::ErgFromKelvin;
    using dec3d::physics::ErgFromKeV;
    using dec3d::physics::KeVFromErg;
    using dec3d::physics::PhysicsConstantsCGS;

    CheckGreater(PhysicsConstantsCGS::pi, 3.14159, "pi");
    CheckGreater(PhysicsConstantsCGS::speed_of_light_cm_per_s, 2.9e10, "speed of light");
    CheckGreater(PhysicsConstantsCGS::boltzmann_erg_per_k, 0.0, "Boltzmann constant");
    CheckGreater(PhysicsConstantsCGS::hbar_erg_s, 0.0, "hbar");
    CheckGreater(PhysicsConstantsCGS::electron_mass_g, 0.0, "electron mass");
    CheckGreater(PhysicsConstantsCGS::proton_mass_g, 0.0, "proton mass");
    CheckGreater(
        PhysicsConstantsCGS::deuteron_mass_g,
        PhysicsConstantsCGS::proton_mass_g,
        "deuteron mass");
    CheckGreater(
        PhysicsConstantsCGS::triton_mass_g,
        PhysicsConstantsCGS::deuteron_mass_g,
        "triton mass");
    CheckGreater(
        PhysicsConstantsCGS::elementary_charge_statcoulomb,
        0.0,
        "elementary charge");
    CheckGreater(PhysicsConstantsCGS::erg_per_ev, 0.0, "erg per eV");
    CheckGreater(PhysicsConstantsCGS::erg_per_kev, 0.0, "erg per keV");

    CheckNear(
        PhysicsConstantsCGS::erg_per_kev,
        1000.0 * PhysicsConstantsCGS::erg_per_ev,
        1.0e-24,
        "keV/eV conversion");
    CheckNear(ErgFromEV(1.0), PhysicsConstantsCGS::erg_per_ev, 1.0e-30, "1 eV");
    CheckNear(ErgFromKeV(1.0), PhysicsConstantsCGS::erg_per_kev, 1.0e-27, "1 keV");
    CheckNear(
        ErgFromKelvin(1.0),
        PhysicsConstantsCGS::boltzmann_erg_per_k,
        1.0e-32,
        "1 K");
    CheckNear(KeVFromErg(PhysicsConstantsCGS::erg_per_kev), 1.0, 1.0e-15, "erg to keV");
    CheckNear(
        DefaultMeanDTIonMassG(),
        0.5 * (PhysicsConstantsCGS::deuteron_mass_g + PhysicsConstantsCGS::triton_mass_g),
        1.0e-36,
        "mean DT ion mass");
    CheckNear(DefaultZbar(), 1.0, 0.0, "default Zbar");

    const std::string diagnostics = dec3d::physics::BuildUnitsDiagnosticsLine();
    DEC3D_CHECK(diagnostics.find("diagnostic_id=p2.units.contract") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("unit_system=cgs") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("temperature_internal_unit=erg_per_particle") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("temperature_output_unit=keV") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("constants_source=PhysicsConstantsCGS") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("speed_of_light_cm_per_s=") != std::string::npos);
    DEC3D_CHECK(diagnostics.find("mean_dt_ion_mass_g=") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
