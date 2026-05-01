#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::state {

namespace {

struct RecoverySummary {
  double min_e_electron{std::numeric_limits<double>::infinity()};
  double min_e_ion{std::numeric_limits<double>::infinity()};
  double min_ne{std::numeric_limits<double>::infinity()};
  double max_ne{-std::numeric_limits<double>::infinity()};
  double min_ni{std::numeric_limits<double>::infinity()};
  double max_ni{-std::numeric_limits<double>::infinity()};
  double min_te_kev{std::numeric_limits<double>::infinity()};
  double max_te_kev{-std::numeric_limits<double>::infinity()};
  double min_ti_kev{std::numeric_limits<double>::infinity()};
  double max_ti_kev{-std::numeric_limits<double>::infinity()};
  double max_partition_residual{0.0};
};

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

void Fail(
    ThermodynamicRecoveryResult& result,
    const std::string& reason,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "diagnostic_id=p2.thermo.recovery.failure"
         << "; failure_reason=" << reason
         << "; first_bad_cell=" << radial << "," << theta << "," << phi;
  result.failure_diagnostics = report.str();
}

[[nodiscard]] std::string BuildRecoveryDiagnostics(
    const ThermodynamicRecoveryOptions& options,
    std::size_t cell_count,
    const RecoverySummary& summary) {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p2.thermo.recovery"
         << "; implementation_id=p2.thermo.recovery.cgs_v1"
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; gamma=" << HydroIdealGasGamma()
         << "; zbar=" << options.zbar
         << "; mean_ion_mass_g=" << options.mean_ion_mass_g
         << "; electron_energy_floor=" << options.electron_energy_floor
         << "; ion_energy_floor=" << options.ion_energy_floor
         << "; cell_count=" << cell_count
         << "; min_e_electron=" << summary.min_e_electron
         << "; min_e_ion=" << summary.min_e_ion
         << "; min_ne_cm3=" << summary.min_ne
         << "; max_ne_cm3=" << summary.max_ne
         << "; min_ni_cm3=" << summary.min_ni
         << "; max_ni_cm3=" << summary.max_ni
         << "; min_Te_keV=" << summary.min_te_kev
         << "; max_Te_keV=" << summary.max_te_kev
         << "; min_Ti_keV=" << summary.min_ti_kev
         << "; max_Ti_keV=" << summary.max_ti_kev
         << "; max_energy_partition_residual=" << summary.max_partition_residual
         << "; invalid_cell_count=0";
  return report.str();
}

}  // namespace

bool ThermodynamicRecoveryResult::is_complete() const noexcept {
  return success &&
         cell_count > 0 &&
         !cells.empty() &&
         !unit_diagnostics.empty() &&
         unit_diagnostics.find("diagnostic_id=p2.units.contract") != std::string::npos &&
         !recovery_diagnostics.empty() &&
         recovery_diagnostics.find("diagnostic_id=p2.thermo.recovery") != std::string::npos &&
         recovery_diagnostics.find("min_ne_cm3=") != std::string::npos &&
         recovery_diagnostics.find("max_ne_cm3=") != std::string::npos &&
         recovery_diagnostics.find("min_ni_cm3=") != std::string::npos &&
         recovery_diagnostics.find("max_ni_cm3=") != std::string::npos;
}

ThermodynamicRecoveryResult RecoverThermodynamicState(
    const CanonicalState& state,
    ThermodynamicRecoveryOptions options) noexcept {
  ThermodynamicRecoveryResult result;
  result.unit_diagnostics = dec3d::physics::BuildUnitsDiagnosticsLine();

  if (options.mean_ion_mass_g == 0.0) {
    options.mean_ion_mass_g = dec3d::physics::DefaultMeanDTIonMassG();
  }
  if (options.zbar == 0.0) {
    options.zbar = dec3d::physics::DefaultZbar();
  }

  if (!(state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::rho) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_r) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_theta) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_phi) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_fluid_total) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_electron))) {
    Fail(result, "canonical state is missing authoritative hydro storage", 0, 0, 0);
    return result;
  }

  if (!(Finite(options.mean_ion_mass_g) && options.mean_ion_mass_g > 0.0 &&
        Finite(options.zbar) && options.zbar > 0.0 &&
        Finite(options.electron_energy_floor) &&
        Finite(options.ion_energy_floor))) {
    Fail(result, "thermodynamic recovery options are not physical", 0, 0, 0);
    return result;
  }

  result.cells = dec3d::core::Array3D<ThermodynamicCellState>(
      state.rho.extent_r(),
      state.rho.extent_theta(),
      state.rho.extent_phi(),
      ThermodynamicCellState{});

  RecoverySummary summary;
  const double gamma_minus_one = HydroIdealGasGamma() - 1.0;

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double rho = state.rho(radial, theta, phi);
        const double mom_r = state.mom_r(radial, theta, phi);
        const double mom_theta = state.mom_theta(radial, theta, phi);
        const double mom_phi = state.mom_phi(radial, theta, phi);
        const double e_total = state.e_fluid_total(radial, theta, phi);
        const double e_electron = state.e_electron(radial, theta, phi);

        if (!(Finite(rho) && rho > 0.0)) {
          Fail(result, "rho must be positive", radial, theta, phi);
          return result;
        }
        if (!(Finite(mom_r) && Finite(mom_theta) && Finite(mom_phi) &&
              Finite(e_total) && Finite(e_electron))) {
          Fail(result, "authoritative state contains non-finite values", radial, theta, phi);
          return result;
        }

        const double kinetic =
            0.5 * (mom_r * mom_r + mom_theta * mom_theta + mom_phi * mom_phi) / rho;
        const double e_ion = e_total - kinetic - e_electron;
        if (e_electron < options.electron_energy_floor) {
          Fail(result, "electron energy fell below minimum threshold", radial, theta, phi);
          return result;
        }
        if (e_ion < options.ion_energy_floor) {
          Fail(result, "ion energy fell below minimum threshold", radial, theta, phi);
          return result;
        }

        const double n_i = rho / options.mean_ion_mass_g;
        const double n_e = options.zbar * n_i;
        const double p_e = gamma_minus_one * e_electron;
        const double p_i = gamma_minus_one * e_ion;
        const double t_e = p_e / n_e;
        const double t_i = p_i / n_i;

        if (!(Finite(n_i) && n_i > 0.0 && Finite(n_e) && n_e > 0.0 &&
              Finite(p_e) && Finite(p_i) && Finite(t_e) && Finite(t_i))) {
          Fail(result, "derived thermodynamic state is not finite", radial, theta, phi);
          return result;
        }

        ThermodynamicCellState cell;
        cell.rho_g_per_cm3 = rho;
        cell.kinetic_energy_density_erg_per_cm3 = kinetic;
        cell.e_electron_erg_per_cm3 = e_electron;
        cell.e_ion_erg_per_cm3 = e_ion;
        cell.p_e_erg_per_cm3 = p_e;
        cell.p_i_erg_per_cm3 = p_i;
        cell.n_e_cm3 = n_e;
        cell.n_i_cm3 = n_i;
        cell.t_e_erg_per_particle = t_e;
        cell.t_i_erg_per_particle = t_i;
        cell.t_e_keV = dec3d::physics::KeVFromErg(t_e);
        cell.t_i_keV = dec3d::physics::KeVFromErg(t_i);
        cell.mean_ion_mass_g = options.mean_ion_mass_g;
        cell.zbar = options.zbar;
        cell.energy_partition_residual = e_total - kinetic - e_electron - e_ion;
        result.cells(radial, theta, phi) = cell;

        summary.min_e_electron = std::min(summary.min_e_electron, e_electron);
        summary.min_e_ion = std::min(summary.min_e_ion, e_ion);
        summary.min_ne = std::min(summary.min_ne, n_e);
        summary.max_ne = std::max(summary.max_ne, n_e);
        summary.min_ni = std::min(summary.min_ni, n_i);
        summary.max_ni = std::max(summary.max_ni, n_i);
        summary.min_te_kev = std::min(summary.min_te_kev, cell.t_e_keV);
        summary.max_te_kev = std::max(summary.max_te_kev, cell.t_e_keV);
        summary.min_ti_kev = std::min(summary.min_ti_kev, cell.t_i_keV);
        summary.max_ti_kev = std::max(summary.max_ti_kev, cell.t_i_keV);
        summary.max_partition_residual =
            std::max(summary.max_partition_residual, std::abs(cell.energy_partition_residual));
        ++result.cell_count;
      }
    }
  }

  result.success = true;
  result.recovery_diagnostics = BuildRecoveryDiagnostics(options, result.cell_count, summary);
  return result;
}

bool ValidateThermodynamicRecoveryDiagnostics(
    const ThermodynamicRecoveryResult& result) noexcept {
  return result.is_complete() &&
         result.unit_diagnostics.find("unit_system=cgs") != std::string::npos &&
         result.unit_diagnostics.find("constants_source=PhysicsConstantsCGS") != std::string::npos &&
         result.recovery_diagnostics.find("implementation_id=p2.thermo.recovery.cgs_v1") != std::string::npos &&
         result.recovery_diagnostics.find("temperature_internal_unit=erg_per_particle") != std::string::npos &&
         result.recovery_diagnostics.find("temperature_output_unit=keV") != std::string::npos &&
         result.recovery_diagnostics.find("min_ne_cm3=") != std::string::npos &&
         result.recovery_diagnostics.find("max_ne_cm3=") != std::string::npos &&
         result.recovery_diagnostics.find("min_ni_cm3=") != std::string::npos &&
         result.recovery_diagnostics.find("max_ni_cm3=") != std::string::npos &&
         result.recovery_diagnostics.find("max_energy_partition_residual=") != std::string::npos;
}

}  // namespace dec3d::state
