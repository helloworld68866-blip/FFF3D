#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <memory>
#include <string>
#include <vector>

namespace dec3d::state {

[[nodiscard]] constexpr double HydroIdealGasGamma() noexcept {
  return 5.0 / 3.0;
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask BuildHydroAuthorizedWriteMask() noexcept;
[[nodiscard]] dec3d::core::AuthoritativeFieldMask BuildHydroElectronWritebackMask() noexcept;

[[nodiscard]] double ElectronPressureFromElectronEnergyDensity(double electron_energy_density) noexcept;
[[nodiscard]] double ElectronEnergyDensityFromPressure(double electron_pressure) noexcept;
[[nodiscard]] double ChiEFromElectronPressure(double electron_pressure) noexcept;
[[nodiscard]] double ElectronPressureFromChiE(double chi_e) noexcept;

struct HydroStateView {
  dec3d::core::Array3D<double>* rho{nullptr};
  dec3d::core::Array3D<double>* mom_r{nullptr};
  dec3d::core::Array3D<double>* mom_theta{nullptr};
  dec3d::core::Array3D<double>* mom_phi{nullptr};
  dec3d::core::Array3D<double>* e_fluid_total{nullptr};
  dec3d::core::Array3D<double>* e_electron{nullptr};
  std::unique_ptr<dec3d::core::Array3D<double>> owned_rho;
  std::unique_ptr<dec3d::core::Array3D<double>> owned_mom_r;
  std::unique_ptr<dec3d::core::Array3D<double>> owned_mom_theta;
  std::unique_ptr<dec3d::core::Array3D<double>> owned_mom_phi;
  std::unique_ptr<dec3d::core::Array3D<double>> owned_e_fluid_total;
  std::unique_ptr<dec3d::core::Array3D<double>> owned_e_electron;
  dec3d::core::Array3D<double> chi_e;
  dec3d::core::Array3D<double> alpha_chi;
  std::vector<dec3d::core::Array3D<double>> radiation_chi;
  bool operator_local_chi_e{true};
  bool operator_local_alpha_chi{false};
  bool operator_local_radiation_chi{true};
  bool transactional_work_copy{false};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct HydroAlphaViewOptions {
  bool enabled{false};
  double alpha_energy_floor_erg_cm3{0.0};
};

struct HydroWritebackResult {
  bool success{false};
  bool electron_channel_updated{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] HydroStateView BuildHydroStateView(CanonicalState& state) noexcept;
[[nodiscard]] HydroStateView BuildHydroStateView(
    CanonicalState& state,
    HydroAlphaViewOptions alpha_options) noexcept;
[[nodiscard]] HydroStateView BuildHydroWorkView(const CanonicalState& state) noexcept;
[[nodiscard]] HydroStateView BuildHydroWorkView(
    const CanonicalState& state,
    HydroAlphaViewOptions alpha_options) noexcept;

[[nodiscard]] HydroWritebackResult CommitHydroWriteback(
    CanonicalState& state,
    const HydroStateView& hydro_view,
    dec3d::core::AuthoritativeFieldMask updated_fields = BuildHydroElectronWritebackMask()) noexcept;

}  // namespace dec3d::state
