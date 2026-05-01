#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/radial_overlap_remap.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] std::size_t Linear(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  return (radial * theta_cells + theta) * phi_cells + phi;
}

[[nodiscard]] double TotalMass(
    const dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>& cells,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double total = 0.0;
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < cells.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < cells.extent_phi(); ++phi) {
        total += cells(radial, theta, phi).rho *
                 geometry.cell_volumes[Linear(
                     radial,
                     theta,
                     phi,
                     cells.extent_theta(),
                     cells.extent_phi())];
      }
    }
  }
  return total;
}

[[nodiscard]] double ColumnMass(
    const dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>& cells,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta,
    std::size_t phi) {
  double total = 0.0;
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    total += cells(radial, theta, phi).rho *
             geometry.cell_volumes[Linear(
                 radial,
                 theta,
                 phi,
                 cells.extent_theta(),
                 cells.extent_phi())];
  }
  return total;
}

}  // namespace

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::hydro::HydroConservativeState;
    using dec3d::hydro::RemapHydroStateRadiallyConservative;
    using dec3d::mesh::ApplyRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;

    constexpr std::size_t kRadial = 4u;
    constexpr std::size_t kTheta = 3u;
    constexpr std::size_t kPhi = 2u;
    const auto old_geometry =
        BuildSphericalGeometry(SphericalMeshDescriptor{kRadial, kTheta, kPhi, 0.25, 1.25});
    auto new_geometry = old_geometry;

    dec3d::core::MeshUpdateProposal proposal;
    proposal.requested = true;
    proposal.radial_ale = true;
    proposal.radial_face_indexing_is_global = true;
    proposal.global_radial_face_count = old_geometry.radial_faces.size();
    proposal.dt_s = 1.0e-3;
    proposal.radial_face_velocities.assign(old_geometry.radial_faces.size(), 0.0);
    proposal.proposed_radial_faces = {0.25, 0.47, 0.73, 1.02, 1.25};
    proposal.implementation_id = "p1.mesh.ale.radial_proposal";
    proposal.summary = "radial_face_indexing=global; test_nonuniform_faces=true";
    DEC3D_CHECK(ApplyRadialAleMeshUpdateProposal(proposal, new_geometry).success);
    DEC3D_CHECK(old_geometry.is_valid());
    DEC3D_CHECK(new_geometry.is_valid());

    Array3D<HydroConservativeState> cells(kRadial, kTheta, kPhi);
    for (std::size_t radial = 0; radial < kRadial; ++radial) {
      for (std::size_t theta = 0; theta < kTheta; ++theta) {
        for (std::size_t phi = 0; phi < kPhi; ++phi) {
          const double base = 1.0 + 0.2 * static_cast<double>(radial) +
                              0.03 * static_cast<double>(theta) +
                              0.01 * static_cast<double>(phi);
          cells(radial, theta, phi) = HydroConservativeState{
              base,
              0.1 * base,
              -0.02 * base,
              0.03 * base,
              2.0 + 0.4 * base,
              0.5 + 0.05 * base};
        }
      }
    }

    const double old_mass = TotalMass(cells, old_geometry);
    const double old_column_mass = ColumnMass(cells, old_geometry, 1u, 1u);
    const auto remap = RemapHydroStateRadiallyConservative(
        cells,
        old_geometry,
        new_geometry,
        proposal.proposed_radial_faces);

    DEC3D_CHECK(remap.success);
    DEC3D_CHECK(remap.diagnostics.is_complete());
    DEC3D_CHECK(remap.diagnostics.remap_order == "first_order_proposal_mapped_overlap");
    DEC3D_CHECK(std::abs(TotalMass(remap.remapped_cells, new_geometry) - old_mass) < 1.0e-11);
    DEC3D_CHECK(
        std::abs(ColumnMass(remap.remapped_cells, new_geometry, 1u, 1u) - old_column_mass) <
        1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.mass_residual) < 1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.mom_r_residual) < 1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.mom_theta_residual) < 1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.mom_phi_residual) < 1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.e_fluid_total_residual) < 1.0e-11);
    DEC3D_CHECK(std::abs(remap.diagnostics.chi_e_residual) < 1.0e-11);
    DEC3D_CHECK(remap.diagnostics.max_volume_coverage_error < 1.0e-12);
    DEC3D_CHECK(!remap.diagnostics.report_line.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
