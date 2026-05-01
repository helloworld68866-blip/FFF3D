#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::hydro {

struct GeometricSourceSnapshot {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
  std::vector<HydroConservativeState> states;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct GeometricSourceStepResult {
  bool success{false};
  std::size_t updated_cell_count{0};
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct GeometricSourceTermsResult {
  bool success{false};
  std::size_t updated_cell_count{0};
  std::vector<HydroConservativeState> source_terms;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete(std::size_t expected_cell_count) const noexcept;
};

[[nodiscard]] GeometricSourceSnapshot CaptureGeometricSourceSnapshot(
    const dec3d::state::HydroStateView& hydro_view) noexcept;

[[nodiscard]] GeometricSourceTermsResult ComputeGeometricSourceTerms(
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept;

[[nodiscard]] GeometricSourceStepResult ApplyGeometricSourceStep(
    dec3d::state::HydroStateView& hydro_view,
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) noexcept;

}  // namespace dec3d::hydro
