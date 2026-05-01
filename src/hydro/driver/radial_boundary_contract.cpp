#include "hydro/driver/radial_boundary_contract.hpp"

#include <sstream>

namespace dec3d::hydro {

std::size_t MapThetaAcrossOrigin(std::size_t theta, std::size_t theta_cells) noexcept {
  if (theta_cells == 0u || theta >= theta_cells) {
    return theta_cells;
  }

  return theta_cells - 1u - theta;
}

std::size_t MapPhiAcrossOrigin(std::size_t phi, std::size_t phi_cells) noexcept {
  return MapPhiAcrossPole(phi, phi_cells);
}

HydroConservativeState BuildOriginRadialGhostState(
    const HydroConservativeState& mapped_interior_state) noexcept {
  return {
      mapped_interior_state.rho,
      -mapped_interior_state.mom_r,
      mapped_interior_state.mom_theta,
      -mapped_interior_state.mom_phi,
      mapped_interior_state.e_fluid_total,
      mapped_interior_state.chi_e,
      mapped_interior_state.alpha_chi,
      mapped_interior_state.radiation_chi};
}

PreparedBoundaryStates PrepareRadialBoundaryStates(
    BoundaryFace face,
    const HydroConservativeState& boundary_interior_state) noexcept {
  PreparedBoundaryStates prepared;

  const auto primitive = RecoverPrimitiveState(boundary_interior_state);
  if (!boundary_interior_state.is_finite() || !primitive.is_physical()) {
    prepared.failure_reason = "radial boundary interior state is not physical";
  } else {
    HydroPrimitiveState ghost_primitive = primitive;
    bool zero_inflow_clamped = false;

    if (face == BoundaryFace::upper && ghost_primitive.v_r < 0.0) {
      ghost_primitive.v_r = 0.0;
      zero_inflow_clamped = true;
    }

    const auto ghost_state = MakeConservativeState(ghost_primitive);

    switch (face) {
      case BoundaryFace::lower:
        prepared.left_state = ghost_state;
        prepared.right_state = boundary_interior_state;
        prepared.success = true;
        break;
      case BoundaryFace::upper:
        prepared.left_state = boundary_interior_state;
        prepared.right_state = ghost_state;
        prepared.success = true;
        break;
    }

    std::ostringstream report;
    report << "radial_boundary_success=" << (prepared.success ? "true" : "false")
           << "; face=" << (face == BoundaryFace::lower ? "lower" : "upper")
           << "; zero_inflow_clamped=" << (zero_inflow_clamped ? "true" : "false");
    prepared.report_line = report.str();
  }

  if (!prepared.success && prepared.report_line.empty()) {
    std::ostringstream report;
    report << "radial_boundary_success=false";
    if (!prepared.failure_reason.empty()) {
      report << "; failure_reason=" << prepared.failure_reason;
    }
    prepared.report_line = report.str();
  }

  return prepared;
}

}  // namespace dec3d::hydro
