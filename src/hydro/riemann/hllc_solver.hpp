#pragma once

#include <string>
#include <vector>

namespace dec3d::hydro {

struct HydroPrimitiveState {
  double rho{0.0};
  double v_r{0.0};
  double v_theta{0.0};
  double v_phi{0.0};
  double pressure{0.0};
  double chi_e{0.0};
  double alpha_chi{0.0};
  std::vector<double> radiation_chi;

  [[nodiscard]] bool is_physical() const noexcept;
};

struct HydroConservativeState {
  double rho{0.0};
  double mom_r{0.0};
  double mom_theta{0.0};
  double mom_phi{0.0};
  double e_fluid_total{0.0};
  double chi_e{0.0};
  double alpha_chi{0.0};
  std::vector<double> radiation_chi;

  [[nodiscard]] bool is_finite() const noexcept;
};

struct HllcWaveStructure {
  double s_left{0.0};
  double s_star{0.0};
  double s_right{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

enum class HllcActiveRegion {
  invalid,
  left_flux,
  left_star,
  right_star,
  right_flux,
};

struct MovingInterfaceBranch {
  HllcActiveRegion active_region{HllcActiveRegion::invalid};
  HydroConservativeState upwind_state;
  HydroConservativeState upwind_flux;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct HllcResult {
  bool success{false};
  HllcWaveStructure waves;
  HydroConservativeState left_star_state;
  HydroConservativeState right_star_state;
  HydroConservativeState interface_flux;
  HllcActiveRegion active_region{HllcActiveRegion::invalid};
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MovingInterfaceHllcFluxResult {
  bool success{false};
  HllcWaveStructure waves;
  HllcActiveRegion active_region{HllcActiveRegion::invalid};
  HllcActiveRegion static_zero_region{HllcActiveRegion::invalid};
  HydroConservativeState selected_state;
  HydroConservativeState selected_flux;
  HydroConservativeState flux;
  double face_speed{0.0};
  std::string mode;
  std::string remap_order;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] HydroConservativeState MakeConservativeState(
    const HydroPrimitiveState& primitive) noexcept;

[[nodiscard]] HydroPrimitiveState RecoverPrimitiveState(
    const HydroConservativeState& conservative) noexcept;

[[nodiscard]] HydroConservativeState ComputePhysicalFlux(
    const HydroConservativeState& conservative) noexcept;

[[nodiscard]] HllcResult SolveHllcRiemann(
    const HydroConservativeState& left,
    const HydroConservativeState& right) noexcept;

[[nodiscard]] MovingInterfaceBranch SelectMovingInterfaceBranch(
    const HllcResult& result,
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept;

[[nodiscard]] HydroConservativeState ComputeAleCorrectedFlux(
    const HllcResult& result,
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept;

[[nodiscard]] MovingInterfaceHllcFluxResult SolveMovingInterfaceHllcFlux(
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept;

}  // namespace dec3d::hydro
