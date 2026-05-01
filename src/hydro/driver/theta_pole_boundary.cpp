#include "hydro/driver/theta_pole_boundary.hpp"

#include <sstream>

namespace dec3d::hydro {

bool PreparedBoundaryStates::is_complete() const noexcept {
  return success && left_state.is_finite() && right_state.is_finite() && !report_line.empty();
}

std::size_t MapPhiAcrossPole(std::size_t phi, std::size_t phi_cells) noexcept {
  if (phi_cells == 0u || (phi_cells % 2u) != 0u || phi >= phi_cells) {
    return phi_cells;
  }

  return (phi + (phi_cells / 2u)) % phi_cells;
}

HydroConservativeState BuildThetaPoleGhostState(
    const HydroConservativeState& mapped_interior_state) noexcept {
  return {
      mapped_interior_state.rho,
      mapped_interior_state.mom_r,
      -mapped_interior_state.mom_theta,
      -mapped_interior_state.mom_phi,
      mapped_interior_state.e_fluid_total,
      mapped_interior_state.chi_e,
      mapped_interior_state.alpha_chi,
      mapped_interior_state.radiation_chi};
}

PreparedBoundaryStates PrepareThetaPoleBoundaryStates(
    BoundaryFace face,
    const HydroConservativeState& boundary_interior_state,
    const HydroConservativeState& mapped_interior_state) noexcept {
  PreparedBoundaryStates prepared;
  const auto ghost_state = BuildThetaPoleGhostState(mapped_interior_state);

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

  if (!prepared.success) {
    prepared.failure_reason = "theta pole boundary face is invalid";
  }

  std::ostringstream report;
  report << "theta_pole_boundary_success=" << (prepared.success ? "true" : "false");
  if (!prepared.failure_reason.empty()) {
    report << "; failure_reason=" << prepared.failure_reason;
  }
  prepared.report_line = report.str();
  return prepared;
}

PreparedBoundaryStates PreparePhiPeriodicBoundaryStates(
    const HydroConservativeState& low_edge_state,
    const HydroConservativeState& high_edge_state) noexcept {
  PreparedBoundaryStates prepared;
  prepared.success = true;
  prepared.left_state = high_edge_state;
  prepared.right_state = low_edge_state;

  std::ostringstream report;
  report << "phi_periodic_boundary_success=true";
  prepared.report_line = report.str();
  return prepared;
}

}  // namespace dec3d::hydro
