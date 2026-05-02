#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/diagnostics/axisymmetric_vector_diagnostics.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <sstream>
#include <string>

namespace {

void CheckNear(double lhs, double rhs, double tol, const char* label) {
  if (std::abs(lhs - rhs) > tol) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("near equality", __FILE__, __LINE__, detail.str());
  }
}

double AveragePhi(const dec3d::core::Array3D<double>& field,
                  std::size_t radial,
                  std::size_t theta) {
  double sum = 0.0;
  for (std::size_t p = 0; p < field.extent_phi(); ++p) {
    sum += field(radial, theta, p);
  }
  return sum / static_cast<double>(field.extent_phi());
}

}  // namespace

int main() {
  {
    dec3d::state::CanonicalStateLayout layout;
    layout.radial_cells = 2u;
    layout.theta_cells = 4u;
    layout.phi_cells = 1u;
    layout.radiation_group_count = 1u;
    auto state = dec3d::state::CanonicalState::Create(layout);
    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{2u, 4u, 1u, 1.0, 2.0});
    DEC3D_CHECK(geometry.valid);
    for (std::size_t r = 0; r < layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < layout.theta_cells; ++t) {
        state.rho(r, t, 0u) = 2.0;
        state.mom_r(r, t, 0u) = 6.0;
        state.mom_theta(r, t, 0u) = 2.0;
        state.mom_phi(r, t, 0u) = 0.0;
      }
    }
    const auto momentum =
        dec3d::state::diagnostics::ComputeAxisymmetricCartesianMomentum(state, geometry);
    CheckNear(momentum.px, 0.0, 0.0, "axisymmetric px");
    CheckNear(momentum.py, 0.0, 0.0, "axisymmetric py");
    DEC3D_CHECK(std::isfinite(momentum.pz));
    DEC3D_CHECK(momentum.report_line.find("cartesian_momentum_phi_average=analytic") !=
                std::string::npos);
  }

  dec3d::state::CanonicalStateLayout axisym_layout;
  axisym_layout.radial_cells = 3u;
  axisym_layout.theta_cells = 4u;
  axisym_layout.phi_cells = 1u;
  axisym_layout.radiation_group_count = 1u;
  dec3d::state::CanonicalStateLayout full_layout = axisym_layout;
  full_layout.phi_cells = 8u;
  auto axisym = dec3d::state::CanonicalState::Create(axisym_layout);
  auto full3d = dec3d::state::CanonicalState::Create(full_layout);
  for (std::size_t r = 0; r < axisym_layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < axisym_layout.theta_cells; ++t) {
      const double rho = 1.0 + 0.1 * static_cast<double>(r) +
                         0.01 * static_cast<double>(t);
      const double mom_r = -2.0 + 0.2 * static_cast<double>(r);
      const double mom_theta = 0.05 * static_cast<double>(t);
      const double e_electron = 10.0 + static_cast<double>(r + t);
      const double e_total = e_electron + 4.0;
      axisym.rho(r, t, 0u) = rho;
      axisym.mom_r(r, t, 0u) = mom_r;
      axisym.mom_theta(r, t, 0u) = mom_theta;
      axisym.mom_phi(r, t, 0u) = 0.0;
      axisym.e_electron(r, t, 0u) = e_electron;
      axisym.e_fluid_total(r, t, 0u) = e_total;
      for (std::size_t p = 0; p < full_layout.phi_cells; ++p) {
        full3d.rho(r, t, p) = rho;
        full3d.mom_r(r, t, p) = mom_r;
        full3d.mom_theta(r, t, p) = mom_theta;
        full3d.mom_phi(r, t, p) = 0.0;
        full3d.e_electron(r, t, p) = e_electron;
        full3d.e_fluid_total(r, t, p) = e_total;
      }
      CheckNear(axisym.rho(r, t, 0u), AveragePhi(full3d.rho, r, t),
                1.0e-14, "rho phi average");
      CheckNear(axisym.mom_r(r, t, 0u), AveragePhi(full3d.mom_r, r, t),
                1.0e-14, "mom_r phi average");
      CheckNear(axisym.mom_theta(r, t, 0u), AveragePhi(full3d.mom_theta, r, t),
                1.0e-14, "mom_theta phi average");
      CheckNear(axisym.e_electron(r, t, 0u), AveragePhi(full3d.e_electron, r, t),
                1.0e-14, "e_electron phi average");
      CheckNear(axisym.e_fluid_total(r, t, 0u),
                AveragePhi(full3d.e_fluid_total, r, t),
                1.0e-14, "e_fluid_total phi average");
    }
  }

  const auto axisym_geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{3u, 4u, 1u, 1.0, 2.0});
  const auto full_geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{3u, 4u, 8u, 1.0, 2.0});
  DEC3D_CHECK(axisym_geometry.valid);
  DEC3D_CHECK(full_geometry.valid);
  const auto axisym_momentum =
      dec3d::state::diagnostics::ComputeAxisymmetricCartesianMomentum(axisym, axisym_geometry);
  const auto full_momentum =
      dec3d::state::diagnostics::ComputeFull3DCartesianMomentum(full3d, full_geometry);
  CheckNear(axisym_momentum.px, 0.0, 0.0, "axisym px");
  CheckNear(axisym_momentum.py, 0.0, 0.0, "axisym py");
  CheckNear(full_momentum.px, 0.0, 1.0e-12, "full3d replicated px");
  CheckNear(full_momentum.py, 0.0, 1.0e-12, "full3d replicated py");
  CheckNear(axisym_momentum.pz, full_momentum.pz, 1.0e-12, "axisym full3d pz");
  return 0;
}
