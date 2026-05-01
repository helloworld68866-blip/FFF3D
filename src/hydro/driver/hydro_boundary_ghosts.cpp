#include "hydro/driver/hydro_boundary_ghosts.hpp"

#include "hydro/driver/radial_boundary_contract.hpp"
#include "hydro/driver/theta_pole_boundary.hpp"

#include <sstream>

namespace dec3d::hydro {

namespace {

[[nodiscard]] const char* DirectionName(HydroDirection direction) noexcept {
  switch (direction) {
    case HydroDirection::radial:
      return "radial";
    case HydroDirection::theta:
      return "theta";
    case HydroDirection::phi:
      return "phi";
  }

  return "invalid";
}

[[nodiscard]] DirectionalGhostPreparation FailedPreparation(
    HydroDirection direction,
    std::size_t ghost_layers,
    std::string failure_reason) noexcept {
  DirectionalGhostPreparation preparation;
  preparation.direction = direction;
  preparation.ghost_layers = ghost_layers;
  preparation.failure_reason = std::move(failure_reason);

  std::ostringstream report;
  report << "directional_ghost_success=false; direction=" << DirectionName(direction)
         << "; ghost_layers=" << ghost_layers;
  if (!preparation.failure_reason.empty()) {
    report << "; failure_reason=" << preparation.failure_reason;
  }
  preparation.report_line = report.str();
  return preparation;
}

[[nodiscard]] bool HasInteriorCells(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states) noexcept {
  return interior_states.extent_r() > 0u &&
         interior_states.extent_theta() > 0u &&
         interior_states.extent_phi() > 0u;
}

void CopyInteriorWithOffset(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    dec3d::core::Array3D<HydroConservativeState>& ghosted_states,
    HydroDirection direction,
    std::size_t ghost_layers) {
  for (std::size_t radial = 0; radial < interior_states.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < interior_states.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < interior_states.extent_phi(); ++phi) {
        switch (direction) {
          case HydroDirection::radial:
            ghosted_states(radial + ghost_layers, theta, phi) =
                interior_states(radial, theta, phi);
            break;
          case HydroDirection::theta:
            ghosted_states(radial, theta + ghost_layers, phi) =
                interior_states(radial, theta, phi);
            break;
          case HydroDirection::phi:
            ghosted_states(radial, theta, phi + ghost_layers) =
                interior_states(radial, theta, phi);
            break;
        }
      }
    }
  }
}

void MarkPreparationSuccess(DirectionalGhostPreparation& preparation) {
  preparation.success = true;
  std::ostringstream report;
  report << "directional_ghost_success=true; direction=" << DirectionName(preparation.direction)
         << "; ghost_layers=" << preparation.ghost_layers;
  preparation.report_line = report.str();
}

}  // namespace

bool DirectionalGhostPreparation::is_complete(
    std::size_t interior_radial,
    std::size_t interior_theta,
    std::size_t interior_phi,
    std::size_t required_ghost_layers) const noexcept {
  if (!success || ghost_layers != required_ghost_layers || report_line.empty()) {
    return false;
  }

  switch (direction) {
    case HydroDirection::radial:
      return ghosted_states.extent_r() == interior_radial + (2u * ghost_layers) &&
             ghosted_states.extent_theta() == interior_theta &&
             ghosted_states.extent_phi() == interior_phi;
    case HydroDirection::theta:
      return ghosted_states.extent_r() == interior_radial &&
             ghosted_states.extent_theta() == interior_theta + (2u * ghost_layers) &&
             ghosted_states.extent_phi() == interior_phi;
    case HydroDirection::phi:
      return ghosted_states.extent_r() == interior_radial &&
             ghosted_states.extent_theta() == interior_theta &&
             ghosted_states.extent_phi() == interior_phi + (2u * ghost_layers);
  }

  return false;
}

bool RadialGhostOverride::is_complete(
    std::size_t theta_cells,
    std::size_t phi_cells,
    std::size_t required_ghost_layers) const noexcept {
  if (ghost_layers != required_ghost_layers) {
    return false;
  }

  const auto shape_matches =
      [theta_cells, phi_cells, required_ghost_layers](
          const dec3d::core::Array3D<HydroConservativeState>& ghost_states) noexcept {
        return ghost_states.extent_r() == required_ghost_layers &&
               ghost_states.extent_theta() == theta_cells &&
               ghost_states.extent_phi() == phi_cells;
      };

  if (has_inner_neighbor && !shape_matches(inner_ghost_states)) {
    return false;
  }
  if (has_outer_neighbor && !shape_matches(outer_ghost_states)) {
    return false;
  }

  return !report_line.empty();
}

DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept {
  return FillHydroRadialGhosts(
      interior_states,
      ghost_layers,
      RadialLowerBoundaryMode::origin_remap);
}

DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    RadialLowerBoundaryMode lower_boundary_mode) noexcept {
  DirectionalGhostPreparation preparation;
  preparation.direction = HydroDirection::radial;
  preparation.ghost_layers = ghost_layers;

  if (!HasInteriorCells(interior_states)) {
    return FailedPreparation(
        HydroDirection::radial,
        ghost_layers,
        "radial ghost fill requires non-empty interior hydro states");
  }
  if (ghost_layers == 0u) {
    return FailedPreparation(
        HydroDirection::radial,
        ghost_layers,
        "radial ghost fill requires at least one ghost layer");
  }

  preparation.ghosted_states = dec3d::core::Array3D<HydroConservativeState>(
      interior_states.extent_r() + (2u * ghost_layers),
      interior_states.extent_theta(),
      interior_states.extent_phi());
  CopyInteriorWithOffset(
      interior_states,
      preparation.ghosted_states,
      HydroDirection::radial,
      ghost_layers);

  const std::size_t outer_radial = interior_states.extent_r() - 1u;
  const bool use_origin_remap =
      lower_boundary_mode == RadialLowerBoundaryMode::origin_remap;
  if (use_origin_remap && (interior_states.extent_phi() % 2u) != 0u) {
    return FailedPreparation(
        HydroDirection::radial,
        ghost_layers,
        "radial origin ghost remap requires an even phi cell count");
  }
  if (use_origin_remap && interior_states.extent_r() < ghost_layers) {
    return FailedPreparation(
        HydroDirection::radial,
        ghost_layers,
        "radial origin ghost remap requires at least ghost_layers interior radial cells");
  }

  for (std::size_t theta = 0; theta < interior_states.extent_theta(); ++theta) {
    for (std::size_t phi = 0; phi < interior_states.extent_phi(); ++phi) {
      PreparedBoundaryStates lower;
      const auto upper = PrepareRadialBoundaryStates(
          BoundaryFace::upper,
          interior_states(outer_radial, theta, phi));
      if (!upper.is_complete()) {
        return FailedPreparation(
            HydroDirection::radial,
            ghost_layers,
            upper.failure_reason.empty() ? "radial upper ghost fill failed"
                                         : upper.failure_reason);
      }
      if (!use_origin_remap) {
        lower = PrepareRadialBoundaryStates(
            BoundaryFace::lower,
            interior_states(0u, theta, phi));
        if (!lower.is_complete()) {
          return FailedPreparation(
              HydroDirection::radial,
              ghost_layers,
              lower.failure_reason.empty() ? "radial lower ghost fill failed"
                                           : lower.failure_reason);
        }
      }

      std::size_t mapped_theta = 0u;
      std::size_t mapped_phi = 0u;
      if (use_origin_remap) {
        mapped_theta = MapThetaAcrossOrigin(theta, interior_states.extent_theta());
        mapped_phi = MapPhiAcrossOrigin(phi, interior_states.extent_phi());
        if (mapped_theta >= interior_states.extent_theta() ||
            mapped_phi >= interior_states.extent_phi()) {
          return FailedPreparation(
              HydroDirection::radial,
              ghost_layers,
              "radial origin ghost remap could not map antipodal angular cell");
        }
      }

      for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
        if (use_origin_remap) {
          const std::size_t requested_radial = ghost_layers - 1u - ghost;
          preparation.ghosted_states(ghost, theta, phi) = BuildOriginRadialGhostState(
              interior_states(requested_radial, mapped_theta, mapped_phi));
        } else {
          preparation.ghosted_states(ghost, theta, phi) = lower.left_state;
        }
        preparation.ghosted_states(ghost_layers + interior_states.extent_r() + ghost, theta, phi) =
            upper.right_state;
      }
    }
  }

  MarkPreparationSuccess(preparation);
  preparation.report_line += use_origin_remap
                                 ? "; radial_lower_boundary=origin_remap"
                                 : "; radial_lower_boundary=boundary_contract";
  return preparation;
}

DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    const RadialGhostOverride& ghost_override) noexcept {
  return FillHydroRadialGhosts(
      interior_states,
      ghost_layers,
      ghost_override,
      RadialLowerBoundaryMode::origin_remap);
}

DirectionalGhostPreparation FillHydroRadialGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers,
    const RadialGhostOverride& ghost_override,
    RadialLowerBoundaryMode lower_boundary_mode) noexcept {
  if (!ghost_override.is_complete(
          interior_states.extent_theta(),
          interior_states.extent_phi(),
          ghost_layers)) {
    return FailedPreparation(
        HydroDirection::radial,
        ghost_layers,
        "radial ghost override is incomplete");
  }

  auto preparation = FillHydroRadialGhosts(
      interior_states,
      ghost_layers,
      lower_boundary_mode);
  if (!preparation.success) {
    return preparation;
  }

  for (std::size_t theta = 0; theta < interior_states.extent_theta(); ++theta) {
    for (std::size_t phi = 0; phi < interior_states.extent_phi(); ++phi) {
      for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
        if (ghost_override.has_inner_neighbor) {
          preparation.ghosted_states(ghost, theta, phi) =
              ghost_override.inner_ghost_states(ghost, theta, phi);
        }
        if (ghost_override.has_outer_neighbor) {
          preparation.ghosted_states(
              ghost_layers + interior_states.extent_r() + ghost,
              theta,
              phi) = ghost_override.outer_ghost_states(ghost, theta, phi);
        }
      }
    }
  }

  std::ostringstream report;
  report << "directional_ghost_success=true; direction=radial"
         << "; ghost_layers=" << ghost_layers
         << "; inner_neighbor=" << (ghost_override.has_inner_neighbor ? "true" : "false")
         << "; outer_neighbor=" << (ghost_override.has_outer_neighbor ? "true" : "false");
  preparation.report_line = report.str();
  return preparation;
}

DirectionalGhostPreparation FillHydroThetaPoleGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept {
  DirectionalGhostPreparation preparation;
  preparation.direction = HydroDirection::theta;
  preparation.ghost_layers = ghost_layers;

  if (!HasInteriorCells(interior_states)) {
    return FailedPreparation(
        HydroDirection::theta,
        ghost_layers,
        "theta ghost fill requires non-empty interior hydro states");
  }
  if (ghost_layers == 0u) {
    return FailedPreparation(
        HydroDirection::theta,
        ghost_layers,
        "theta ghost fill requires at least one ghost layer");
  }
  if (interior_states.extent_theta() < ghost_layers) {
    return FailedPreparation(
        HydroDirection::theta,
        ghost_layers,
        "theta ghost fill requires at least ng interior theta cells");
  }
  if ((interior_states.extent_phi() % 2u) != 0u) {
    return FailedPreparation(
        HydroDirection::theta,
        ghost_layers,
        "theta pole ghost fill requires an even phi cell count");
  }

  preparation.ghosted_states = dec3d::core::Array3D<HydroConservativeState>(
      interior_states.extent_r(),
      interior_states.extent_theta() + (2u * ghost_layers),
      interior_states.extent_phi());
  CopyInteriorWithOffset(
      interior_states,
      preparation.ghosted_states,
      HydroDirection::theta,
      ghost_layers);

  for (std::size_t radial = 0; radial < interior_states.extent_r(); ++radial) {
    for (std::size_t phi = 0; phi < interior_states.extent_phi(); ++phi) {
      const auto mapped_phi = MapPhiAcrossPole(phi, interior_states.extent_phi());
      if (mapped_phi >= interior_states.extent_phi()) {
        return FailedPreparation(
            HydroDirection::theta,
            ghost_layers,
            "theta pole ghost fill produced an invalid mapped phi index");
      }

      for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
        preparation.ghosted_states(radial, ghost_layers - ghost, phi) =
            BuildThetaPoleGhostState(interior_states(radial, ghost - 1u, mapped_phi));
        preparation.ghosted_states(
            radial,
            ghost_layers + interior_states.extent_theta() + ghost - 1u,
            phi) = BuildThetaPoleGhostState(
            interior_states(radial, interior_states.extent_theta() - ghost, mapped_phi));
      }
    }
  }

  MarkPreparationSuccess(preparation);
  return preparation;
}

DirectionalGhostPreparation FillHydroPhiPeriodicGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    std::size_t ghost_layers) noexcept {
  DirectionalGhostPreparation preparation;
  preparation.direction = HydroDirection::phi;
  preparation.ghost_layers = ghost_layers;

  if (!HasInteriorCells(interior_states)) {
    return FailedPreparation(
        HydroDirection::phi,
        ghost_layers,
        "phi ghost fill requires non-empty interior hydro states");
  }
  if (ghost_layers == 0u) {
    return FailedPreparation(
        HydroDirection::phi,
        ghost_layers,
        "phi ghost fill requires at least one ghost layer");
  }

  const std::size_t phi_cells = interior_states.extent_phi();
  preparation.ghosted_states = dec3d::core::Array3D<HydroConservativeState>(
      interior_states.extent_r(),
      interior_states.extent_theta(),
      phi_cells + (2u * ghost_layers));
  CopyInteriorWithOffset(
      interior_states,
      preparation.ghosted_states,
      HydroDirection::phi,
      ghost_layers);

  for (std::size_t radial = 0; radial < interior_states.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < interior_states.extent_theta(); ++theta) {
      for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
        const auto lower_source_phi =
            (phi_cells - (ghost % phi_cells)) % phi_cells;
        const auto upper_source_phi = (ghost - 1u) % phi_cells;
        preparation.ghosted_states(radial, theta, ghost_layers - ghost) =
            interior_states(radial, theta, lower_source_phi);
        preparation.ghosted_states(
            radial,
            theta,
            ghost_layers + phi_cells + ghost - 1u) =
            interior_states(radial, theta, upper_source_phi);
      }
    }
  }

  MarkPreparationSuccess(preparation);
  return preparation;
}

DirectionalGhostPreparation PrepareDirectionalGhosts(
    const dec3d::core::Array3D<HydroConservativeState>& interior_states,
    HydroDirection direction,
    std::size_t ghost_layers) noexcept {
  switch (direction) {
    case HydroDirection::radial:
      return FillHydroRadialGhosts(interior_states, ghost_layers);
    case HydroDirection::theta:
      return FillHydroThetaPoleGhosts(interior_states, ghost_layers);
    case HydroDirection::phi:
      return FillHydroPhiPeriodicGhosts(interior_states, ghost_layers);
  }

  return FailedPreparation(direction, ghost_layers, "invalid hydro ghost direction");
}

}  // namespace dec3d::hydro
