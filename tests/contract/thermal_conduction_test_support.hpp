#pragma once

#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>

namespace dec3d::test_support {

constexpr double kPi = 3.141592653589793238462643383279502884;

void CheckNear(double lhs, double rhs, double tolerance, const std::string& label) {
  if (std::abs(lhs - rhs) > tolerance) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("near equality", __FILE__, __LINE__, detail.str());
  }
}

dec3d::transport::DiffusionGridLayout Layout(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return dec3d::transport::DiffusionGridLayout{radial, theta, phi};
}

dec3d::mesh::SphericalGeometryMetadata BuildPatchGeometry(
    const dec3d::transport::DiffusionGridLayout& layout,
    double r_inner,
    double r_outer,
    double theta_lower,
    double theta_upper) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(layout.radial_cells + 1u);
  geometry.theta_faces.resize(layout.theta_cells + 1u);
  geometry.phi_faces.resize(layout.phi_cells + 1u);
  geometry.cell_volumes.resize(layout.cell_count(), 0.0);

  const double dr = (r_outer - r_inner) / static_cast<double>(layout.radial_cells);
  const double dtheta = (theta_upper - theta_lower) / static_cast<double>(layout.theta_cells);
  const double dphi = 2.0 * kPi / static_cast<double>(layout.phi_cells);

  for (std::size_t r = 0; r <= layout.radial_cells; ++r) {
    geometry.radial_faces[r] = r_inner + dr * static_cast<double>(r);
  }
  for (std::size_t t = 0; t <= layout.theta_cells; ++t) {
    geometry.theta_faces[t] = theta_lower + dtheta * static_cast<double>(t);
  }
  for (std::size_t p = 0; p <= layout.phi_cells; ++p) {
    geometry.phi_faces[p] = dphi * static_cast<double>(p);
  }

  geometry.global_volume = 0.0;
  std::size_t linear = 0;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) / 3.0;
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      const double polar_factor =
          std::cos(geometry.theta_faces[t]) - std::cos(geometry.theta_faces[t + 1u]);
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const double azimuthal_factor = geometry.phi_faces[p + 1u] - geometry.phi_faces[p];
        geometry.cell_volumes[linear] = radial_factor * polar_factor * azimuthal_factor;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }

  return geometry;
}

dec3d::transport::GenericDiffusionBoundaryPolicy PatchBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::interior_patch_no_origin,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

std::size_t FlatIndex(
    const dec3d::transport::DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return (radial * layout.theta_cells + theta) * layout.phi_cells + phi;
}

dec3d::state::CanonicalState BuildThermalState(
    const dec3d::transport::DiffusionGridLayout& layout,
    double te_left_keV,
    double te_right_keV,
    double ti_left_keV,
    double ti_right_keV) {
  using dec3d::physics::DefaultMeanDTIonMassG;
  using dec3d::physics::ErgFromKeV;
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;

  auto state = CanonicalState::Create(
      CanonicalStateLayout{layout.radial_cells, layout.theta_cells, layout.phi_cells, 0});
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double n = 2.0;
  const double rho = n * DefaultMeanDTIonMassG();

  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    const double te = ErgFromKeV(r == 0 ? te_left_keV : te_right_keV);
    const double ti = ErgFromKeV(r == 0 ? ti_left_keV : ti_right_keV);
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        state.rho(r, t, p) = rho;
        state.mom_r(r, t, p) = 0.0;
        state.mom_theta(r, t, p) = 0.0;
        state.mom_phi(r, t, p) = 0.0;
        state.e_electron(r, t, p) = n * te / gamma_minus_one;
        const double e_ion = n * ti / gamma_minus_one;
        state.e_fluid_total(r, t, p) = state.e_electron(r, t, p) + e_ion;
      }
    }
  }

  return state;
}

double VolumeWeightedThermalTotal(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const auto recovered = dec3d::state::RecoverThermodynamicState(state);
  DEC3D_CHECK(recovered.success);
  const dec3d::transport::DiffusionGridLayout layout{
      state.layout.radial_cells, state.layout.theta_cells, state.layout.phi_cells};
  double total = 0.0;

  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const auto& cell = recovered.cells(r, t, p);
        total += geometry.cell_volumes[FlatIndex(layout, r, t, p)] *
                 (cell.e_electron_erg_per_cm3 + cell.e_ion_erg_per_cm3);
      }
    }
  }

  return total;
}

}  // namespace dec3d::test_support
