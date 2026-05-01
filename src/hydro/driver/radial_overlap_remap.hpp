#pragma once

#include "core/array/array3d.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <string>
#include <vector>

namespace dec3d::hydro {

struct RadialOverlapRemapDiagnostics {
  bool executed{false};
  std::string remap_order;
  double mass_residual{0.0};
  double mom_r_residual{0.0};
  double mom_theta_residual{0.0};
  double mom_phi_residual{0.0};
  double e_fluid_total_residual{0.0};
  double chi_e_residual{0.0};
  double max_volume_coverage_error{0.0};
  double max_density_angular_spread_before_remap{0.0};
  double max_density_angular_spread_after_remap{0.0};
  double max_tangential_momentum_ratio_after_remap{0.0};
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RadialOverlapRemapResult {
  bool success{false};
  dec3d::core::Array3D<HydroConservativeState> remapped_cells;
  RadialOverlapRemapDiagnostics diagnostics;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] RadialOverlapRemapResult RemapHydroStateRadiallyConservative(
    const dec3d::core::Array3D<HydroConservativeState>& old_cell_averages,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry,
    const std::vector<double>& proposal_mapped_source_faces) noexcept;

}  // namespace dec3d::hydro
