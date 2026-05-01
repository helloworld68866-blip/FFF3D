#pragma once

#include "hydro/driver/theta_pole_boundary.hpp"

namespace dec3d::hydro {

[[nodiscard]] std::size_t MapThetaAcrossOrigin(
    std::size_t theta,
    std::size_t theta_cells) noexcept;

[[nodiscard]] std::size_t MapPhiAcrossOrigin(
    std::size_t phi,
    std::size_t phi_cells) noexcept;

[[nodiscard]] HydroConservativeState BuildOriginRadialGhostState(
    const HydroConservativeState& mapped_interior_state) noexcept;

[[nodiscard]] PreparedBoundaryStates PrepareRadialBoundaryStates(
    BoundaryFace face,
    const HydroConservativeState& boundary_interior_state) noexcept;

}  // namespace dec3d::hydro
