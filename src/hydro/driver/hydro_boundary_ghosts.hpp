#pragma once

#include "core/array/array3d.hpp"
#include "hydro/riemann/hllc_solver.hpp"

#include <cstddef>
#include <string>

namespace dec3d::hydro {

enum class HydroDirection {
  radial,
  theta,
  phi,
};

enum class RadialLowerBoundaryMode {
  origin_remap,
  boundary_contract,
};

struct DirectionalGhostPreparation {
  bool success{false};
  HydroDirection direction{HydroDirection::radial};
  std::size_t ghost_layers{0};
  dec3d::core::Array3D<HydroConservativeState> ghosted_states;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete(
      std::size_t interior_radial,
      std::size_t interior_theta,
      std::size_t interior_phi,
      std::size_t required_ghost_layers) const noexcept;
};

struct RadialGhostOverride {
  bool has_inner_neighbor{false};
  bool has_outer_neighbor{false};
  bool inner_ghost_reuses_boundary_partition{false};
  bool outer_ghost_reuses_boundary_partition{false};
  std::size_t ghost_layers{0};
  dec3d::core::Array3D<HydroConservativeState> inner_ghost_states;
  dec3d::core::Array3D<HydroConservativeState> outer_ghost_states;
  std::string report_line;

  [[nodiscard]] bool is_complete(
      std::size_t theta_cells,
      std::size_t phi_cells,
      std::size_t required_ghost_layers) const noexcept;
};

[[nodiscard]] DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept;

[[nodiscard]] DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    RadialLowerBoundaryMode lower_boundary_mode) noexcept;

[[nodiscard]] DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    const RadialGhostOverride& ghost_override) noexcept;

[[nodiscard]] DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    const RadialGhostOverride& ghost_override,
    RadialLowerBoundaryMode lower_boundary_mode) noexcept;

[[nodiscard]] DirectionalGhostPreparation FillHydroThetaPoleGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept;

[[nodiscard]] DirectionalGhostPreparation FillHydroPhiPeriodicGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept;

[[nodiscard]] DirectionalGhostPreparation PrepareDirectionalGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    HydroDirection direction,
    std::size_t ghost_layers) noexcept;

}  // namespace dec3d::hydro
