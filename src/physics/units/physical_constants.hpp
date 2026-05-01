#pragma once

#include <iomanip>
#include <sstream>
#include <string>

namespace dec3d::physics {

struct PhysicsConstantsCGS {
  static constexpr double pi = 3.141592653589793238462643383279502884;
  static constexpr double speed_of_light_cm_per_s = 2.99792458e10;
  static constexpr double boltzmann_erg_per_k = 1.380649e-16;
  static constexpr double hbar_erg_s = 1.054571817e-27;
  static constexpr double electron_mass_g = 9.1093837015e-28;
  static constexpr double proton_mass_g = 1.67262192369e-24;
  static constexpr double deuteron_mass_g = 3.3435837724e-24;
  static constexpr double triton_mass_g = 5.0073567446e-24;
  static constexpr double alpha_mass_g = 6.6446573357e-24;
  static constexpr double elementary_charge_statcoulomb = 4.803204712570263e-10;
  static constexpr double erg_per_ev = 1.602176634e-12;
  static constexpr double erg_per_kev = 1000.0 * erg_per_ev;
  static constexpr double alpha_charge_number = 2.0;
  static constexpr double alpha_birth_energy_erg = 3.5e6 * erg_per_ev;
};

[[nodiscard]] constexpr double ErgFromEV(double electron_volts) noexcept {
  return electron_volts * PhysicsConstantsCGS::erg_per_ev;
}

[[nodiscard]] constexpr double ErgFromKeV(double kilo_electron_volts) noexcept {
  return kilo_electron_volts * PhysicsConstantsCGS::erg_per_kev;
}

[[nodiscard]] constexpr double ErgFromKelvin(double kelvin) noexcept {
  return kelvin * PhysicsConstantsCGS::boltzmann_erg_per_k;
}

[[nodiscard]] constexpr double KeVFromErg(double erg_per_particle) noexcept {
  return erg_per_particle / PhysicsConstantsCGS::erg_per_kev;
}

[[nodiscard]] constexpr double DefaultMeanDTIonMassG() noexcept {
  return 0.5 * (PhysicsConstantsCGS::deuteron_mass_g + PhysicsConstantsCGS::triton_mass_g);
}

[[nodiscard]] constexpr double DefaultZbar() noexcept {
  return 1.0;
}

[[nodiscard]] inline std::string BuildUnitsDiagnosticsLine() {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p2.units.contract"
         << "; unit_system=cgs"
         << "; length_unit=cm"
         << "; mass_unit=g"
         << "; time_unit=s"
         << "; energy_unit=erg"
         << "; density_unit=g_per_cm3"
         << "; pressure_unit=erg_per_cm3"
         << "; number_density_unit=cm^-3"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; constants_source=PhysicsConstantsCGS"
         << "; pi=" << PhysicsConstantsCGS::pi
         << "; speed_of_light_cm_per_s=" << PhysicsConstantsCGS::speed_of_light_cm_per_s
         << "; k_B_erg_per_K=" << PhysicsConstantsCGS::boltzmann_erg_per_k
         << "; hbar_erg_s=" << PhysicsConstantsCGS::hbar_erg_s
         << "; electron_mass_g=" << PhysicsConstantsCGS::electron_mass_g
         << "; proton_mass_g=" << PhysicsConstantsCGS::proton_mass_g
         << "; deuteron_mass_g=" << PhysicsConstantsCGS::deuteron_mass_g
         << "; triton_mass_g=" << PhysicsConstantsCGS::triton_mass_g
         << "; alpha_mass_g=" << PhysicsConstantsCGS::alpha_mass_g
         << "; elementary_charge_statcoulomb="
         << PhysicsConstantsCGS::elementary_charge_statcoulomb
         << "; erg_per_eV=" << PhysicsConstantsCGS::erg_per_ev
         << "; erg_per_keV=" << PhysicsConstantsCGS::erg_per_kev
         << "; mean_dt_ion_mass_g=" << DefaultMeanDTIonMassG()
         << "; zbar=" << DefaultZbar()
         << "; alpha_charge_number=" << PhysicsConstantsCGS::alpha_charge_number
         << "; alpha_birth_energy_erg=" << PhysicsConstantsCGS::alpha_birth_energy_erg;
  return report.str();
}

}  // namespace dec3d::physics
