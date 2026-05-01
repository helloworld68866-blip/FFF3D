#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/electron_ion_equilibration.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

[[nodiscard]] dec3d::state::CanonicalState BuildOneCellState(
    double te_keV,
    double ti_keV) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 0u});
  const double rho = 2.0 * dec3d::physics::DefaultMeanDTIonMassG();
  const double ne = 2.0;
  const double ni = 2.0;
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double te = dec3d::physics::ErgFromKeV(te_keV);
  const double ti = dec3d::physics::ErgFromKeV(ti_keV);
  const double e_electron = ne * te / gamma_minus_one;
  const double e_ion = ni * ti / gamma_minus_one;

  state.rho(0u, 0u, 0u) = rho;
  state.mom_r(0u, 0u, 0u) = 0.0;
  state.mom_theta(0u, 0u, 0u) = 0.0;
  state.mom_phi(0u, 0u, 0u) = 0.0;
  state.e_electron(0u, 0u, 0u) = e_electron;
  state.e_fluid_total(0u, 0u, 0u) = e_electron + e_ion;
  return state;
}

[[nodiscard]] double KeVFromErg(double value) {
  return dec3d::physics::KeVFromErg(value);
}

struct Temperatures {
  double te_keV{0.0};
  double ti_keV{0.0};
};

[[nodiscard]] Temperatures RecoverTemperatures(
    const dec3d::state::CanonicalState& state) {
  const auto recovered = dec3d::state::RecoverThermodynamicState(state);
  if (!recovered.success) {
    throw std::runtime_error(recovered.failure_reason);
  }
  const auto& cell = recovered.cells(0u, 0u, 0u);
  return {
      KeVFromErg(cell.t_e_erg_per_particle),
      KeVFromErg(cell.t_i_erg_per_particle)};
}

[[nodiscard]] Temperatures ExactTemperatures(
    double te0_keV,
    double ti0_keV,
    double tau_s,
    double time_s) {
  const double equilibrium = 0.5 * (te0_keV + ti0_keV);
  const double delta = (te0_keV - ti0_keV) * std::exp(-2.0 * time_s / tau_s);
  return {equilibrium + 0.5 * delta, equilibrium - 0.5 * delta};
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const std::filesystem::path output_dir =
        argc > 1
            ? std::filesystem::path(argv[1])
            : std::filesystem::path(
                  "F:/dec3d/analysis/output/p2_ei_equilibration_runtime_exact");
    const double te0_keV = argc > 2 ? std::stod(argv[2]) : 5.0;
    const double ti0_keV = argc > 3 ? std::stod(argv[3]) : 0.5;
    const double tau_s = argc > 4 ? std::stod(argv[4]) : 1.0e-12;
    const double final_time_s = argc > 5 ? std::stod(argv[5]) : 3.0e-12;
    const std::size_t samples =
        argc > 6 ? static_cast<std::size_t>(std::stoull(argv[6])) : 41u;

    if (!(tau_s > 0.0) || !(final_time_s >= 0.0) || samples < 2u) {
      throw std::runtime_error("invalid tau/final_time/samples");
    }

    std::filesystem::create_directories(output_dir);
    const auto csv_path = output_dir / "ei_runtime_vs_exact.csv";
    std::ofstream csv(csv_path);
    csv << std::setprecision(17)
        << "time_s,Te_numeric_keV,Ti_numeric_keV,"
        << "Te_exact_keV,Ti_exact_keV,"
        << "Te_abs_error_keV,Ti_abs_error_keV,"
        << "max_expected_decay_residual_erg,"
        << "max_total_thermal_energy_residual,"
        << "report_line\n";

    double max_error = 0.0;
    for (std::size_t sample = 0u; sample < samples; ++sample) {
      const double time_s =
          final_time_s * static_cast<double>(sample) /
          static_cast<double>(samples - 1u);
      auto state = BuildOneCellState(te0_keV, ti0_keV);
      const auto result = dec3d::state::ApplyLocalElectronIonEquilibration(
          state,
          dec3d::state::ElectronIonEquilibrationOptions{
              time_s,
              dec3d::state::ElectronIonTauModel::constant_user_supplied,
              tau_s,
              0.0,
              0.0});
      if (!result.success || !result.is_complete()) {
        throw std::runtime_error(
            result.failure_diagnostics.empty()
                ? "electron-ion equilibration failed"
                : result.failure_diagnostics);
      }
      const auto numeric = RecoverTemperatures(state);
      const auto exact = ExactTemperatures(te0_keV, ti0_keV, tau_s, time_s);
      const double te_error = std::abs(numeric.te_keV - exact.te_keV);
      const double ti_error = std::abs(numeric.ti_keV - exact.ti_keV);
      max_error = std::max(max_error, std::max(te_error, ti_error));
      std::string report = result.report_line;
      for (char& ch : report) {
        if (ch == ',') {
          ch = ' ';
        }
      }
      csv << time_s << ','
          << numeric.te_keV << ','
          << numeric.ti_keV << ','
          << exact.te_keV << ','
          << exact.ti_keV << ','
          << te_error << ','
          << ti_error << ','
          << result.max_expected_decay_residual_erg << ','
          << result.max_total_thermal_energy_residual << ','
          << '"' << report << '"' << '\n';
    }

    const auto summary_path = output_dir / "ei_runtime_vs_exact_summary.txt";
    std::ofstream summary(summary_path);
    summary << std::setprecision(17)
            << "diagnostic_id=p2.ei.runtime_exact_benchmark"
            << "; numerical_source=ApplyLocalElectronIonEquilibration"
            << "; tau_model=constant_user_supplied"
            << "; te0_keV=" << te0_keV
            << "; ti0_keV=" << ti0_keV
            << "; tau_ei_s=" << tau_s
            << "; final_time_s=" << final_time_s
            << "; samples=" << samples
            << "; max_abs_error_keV=" << max_error
            << "; csv=" << csv_path.generic_string()
            << '\n';

    std::cout << summary_path.generic_string() << '\n'
              << csv_path.generic_string() << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
