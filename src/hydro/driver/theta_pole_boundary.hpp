#pragma once

#include "hydro/riemann/hllc_solver.hpp"

#include <cstddef>
#include <string>

namespace dec3d::hydro {

enum class BoundaryFace {
  lower,
  upper,
};

struct PreparedBoundaryStates {
  bool success{false};
  HydroConservativeState left_state;
  HydroConservativeState right_state;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] std::size_t MapPhiAcrossPole(std::size_t phi, std::size_t phi_cells) noexcept;

[[nodiscard]] HydroConservativeState BuildThetaPoleGhostState(
    const HydroConservativeState& mapped_interior_state) noexcept;

[[nodiscard]] PreparedBoundaryStates PrepareThetaPoleBoundaryStates(
    BoundaryFace face,
    const HydroConservativeState& boundary_interior_state,
    const HydroConservativeState& mapped_interior_state) noexcept;

[[nodiscard]] PreparedBoundaryStates PreparePhiPeriodicBoundaryStates(
    const HydroConservativeState& low_edge_state,
    const HydroConservativeState& high_edge_state) noexcept;

}  // namespace dec3d::hydro
