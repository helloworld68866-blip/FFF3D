#pragma once

#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>

namespace dec3d::state::diagnostics {

struct IntegratedCartesianMomentum {
  double px{0.0};
  double py{0.0};
  double pz{0.0};
  std::string report_line;
};

[[nodiscard]] inline double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::state::CanonicalStateLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return geometry.cell_volumes[(radial * layout.theta_cells + theta) * layout.phi_cells + phi];
}

[[nodiscard]] inline IntegratedCartesianMomentum
ComputeAxisymmetricCartesianMomentum(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  IntegratedCartesianMomentum out;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      const double theta = 0.5 * (geometry.theta_faces[t] + geometry.theta_faces[t + 1u]);
      const double volume = CellVolume(geometry, state.layout, r, t, 0u);
      out.pz +=
          (state.mom_r(r, t, 0u) * std::cos(theta) -
           state.mom_theta(r, t, 0u) * std::sin(theta)) *
          volume;
    }
  }
  std::ostringstream report;
  report << "diagnostic_id=axisymmetric.cartesian_momentum"
         << "; cartesian_momentum_phi_average=analytic"
         << "; px=" << out.px << "; py=" << out.py << "; pz=" << out.pz;
  out.report_line = report.str();
  return out;
}

[[nodiscard]] inline IntegratedCartesianMomentum ComputeFull3DCartesianMomentum(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  IntegratedCartesianMomentum out;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      const double theta = 0.5 * (geometry.theta_faces[t] + geometry.theta_faces[t + 1u]);
      const double st = std::sin(theta);
      const double ct = std::cos(theta);
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double phi = 0.5 * (geometry.phi_faces[p] + geometry.phi_faces[p + 1u]);
        const double cp = std::cos(phi);
        const double sp = std::sin(phi);
        const double volume = CellVolume(geometry, state.layout, r, t, p);
        const double mr = state.mom_r(r, t, p);
        const double mt = state.mom_theta(r, t, p);
        const double mp = state.mom_phi(r, t, p);
        out.px += (mr * st * cp + mt * ct * cp - mp * sp) * volume;
        out.py += (mr * st * sp + mt * ct * sp + mp * cp) * volume;
        out.pz += (mr * ct - mt * st) * volume;
      }
    }
  }
  std::ostringstream report;
  report << "diagnostic_id=full3d.cartesian_momentum"
         << "; cartesian_momentum_phi_average=numeric"
         << "; px=" << out.px << "; py=" << out.py << "; pz=" << out.pz;
  out.report_line = report.str();
  return out;
}

}  // namespace dec3d::state::diagnostics
