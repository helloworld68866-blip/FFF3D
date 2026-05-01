#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace {

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

  std::size_t linear = 0;
  geometry.global_volume = 0.0;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) - std::pow(geometry.radial_faces[r], 3)) / 3.0;
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

std::size_t TestLinearIndex(
    const dec3d::transport::DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return ((radial * layout.theta_cells) + theta) * layout.phi_cells + phi;
}

dec3d::mesh::SphericalGeometryMetadata BuildFullSphereGeometry(
    const dec3d::transport::DiffusionGridLayout& layout,
    double r_outer) {
  return BuildPatchGeometry(layout, 0.0, r_outer, 0.0, kPi);
}

dec3d::core::Array3D<double> Filled(
    const dec3d::transport::DiffusionGridLayout& layout,
    double value) {
  return dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, value);
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

dec3d::transport::GenericDiffusionBoundaryPolicy MarshakOuterBoundary() {
  auto policy = PatchBoundary();
  policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::radiation_marshak_vacuum;
  return policy;
}

dec3d::transport::GenericDiffusionBoundaryPolicy FullSphereBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::periodic};
}

dec3d::transport::GenericDiffusionProblem BuildPatchProblem(
    const dec3d::transport::DiffusionGridLayout& layout,
    double dt_s,
    double a,
    double d,
    double c,
    double b,
    double scalar_old) {
  return dec3d::transport::GenericDiffusionProblem{
      layout,
      BuildPatchGeometry(layout, 1.0, 3.0, 0.5, 1.2),
      dt_s,
      Filled(layout, a),
      Filled(layout, d),
      Filled(layout, c),
      Filled(layout, b),
      Filled(layout, scalar_old),
      PatchBoundary()};
}

dec3d::transport::GenericDiffusionProblem BuildFullSphereProblem(
    const dec3d::transport::DiffusionGridLayout& layout,
    double dt_s,
    double a,
    double d,
    double c,
    double b,
    double scalar_old) {
  return dec3d::transport::GenericDiffusionProblem{
      layout,
      BuildFullSphereGeometry(layout, 2.0),
      dt_s,
      Filled(layout, a),
      Filled(layout, d),
      Filled(layout, c),
      Filled(layout, b),
      Filled(layout, scalar_old),
      FullSphereBoundary()};
}

double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::transport::DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return geometry.cell_volumes[TestLinearIndex(layout, radial, theta, phi)];
}

double RadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) {
  const double radius = geometry.radial_faces[radial_face];
  const double polar_factor =
      std::cos(geometry.theta_faces[theta]) - std::cos(geometry.theta_faces[theta + 1u]);
  const double dphi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return radius * radius * polar_factor * dphi;
}

double RadialCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

double MarshakGFace(double d_face, double center_to_boundary_distance_cm) {
  const double h = 0.5 * dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double d_over_distance = d_face / center_to_boundary_distance_cm;
  return (h * d_over_distance) / (h + d_over_distance);
}

double ThetaCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

double PhiCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t phi) {
  return 0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
}

std::vector<double> L1M1ScalarAtCenters(
    const dec3d::transport::DiffusionGridLayout& layout,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  std::vector<double> values(layout.cell_count(), 0.0);
  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      const double polar = ThetaCenter(geometry, theta);
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const double azimuth = PhiCenter(geometry, phi);
        values[TestLinearIndex(layout, radial, theta, phi)] =
            std::sin(polar) * std::cos(azimuth);
      }
    }
  }
  return values;
}

std::vector<double> CartesianXScalarAtCenters(
    const dec3d::transport::DiffusionGridLayout& layout,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  std::vector<double> values(layout.cell_count(), 0.0);
  for (std::size_t radial = 0; radial < layout.radial_cells; ++radial) {
    const double radius =
        0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
    for (std::size_t theta = 0; theta < layout.theta_cells; ++theta) {
      const double polar = ThetaCenter(geometry, theta);
      for (std::size_t phi = 0; phi < layout.phi_cells; ++phi) {
        const double azimuth = PhiCenter(geometry, phi);
        values[TestLinearIndex(layout, radial, theta, phi)] =
            radius * std::sin(polar) * std::cos(azimuth);
      }
    }
  }
  return values;
}

double DiffusionOperatorFromCsr(
    const dec3d::transport::GenericDiffusionAssemblyResult& assembly,
    const dec3d::transport::GenericDiffusionProblem& problem,
    const std::vector<double>& scalar,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  const std::size_t row = TestLinearIndex(problem.layout, radial, theta, phi);
  const double volume = CellVolume(problem.geometry, problem.layout, radial, theta, phi);
  double conductance_matrix_times_scalar = 0.0;
  for (std::size_t slot = assembly.matrix.row_offsets[row];
       slot < assembly.matrix.row_offsets[row + 1u];
       ++slot) {
    const std::size_t column = assembly.matrix.column_indices[slot];
    double value = assembly.matrix.values[slot];
    if (column == row) {
      value -= volume * problem.coefficient_A(radial, theta, phi) / problem.dt_s;
    }
    conductance_matrix_times_scalar += value * scalar[column];
  }
  return -conductance_matrix_times_scalar / volume;
}

double ExactL1M1FiniteVolumeLaplacian(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  const double r0 = geometry.radial_faces[radial];
  const double r1 = geometry.radial_faces[radial + 1u];
  const double radial_ratio =
      (r1 - r0) / ((std::pow(r1, 3) - std::pow(r0, 3)) / 3.0);
  const double theta0 = geometry.theta_faces[theta];
  const double theta1 = geometry.theta_faces[theta + 1u];
  const double phi0 = geometry.phi_faces[phi];
  const double phi1 = geometry.phi_faces[phi + 1u];
  const double integral_sin_squared =
      0.5 * (theta1 - theta0) -
      0.25 * (std::sin(2.0 * theta1) - std::sin(2.0 * theta0));
  const double integral_cos_phi = std::sin(phi1) - std::sin(phi0);
  const double solid_angle =
      (std::cos(theta0) - std::cos(theta1)) * (phi1 - phi0);
  return radial_ratio * (-2.0 * integral_sin_squared * integral_cos_phi / solid_angle);
}

}  // namespace

int main() {
  try {
    using dec3d::transport::AssembleGenericImplicitDiffusionSystem;
    using dec3d::transport::DiffusionBoundaryKind;
    using dec3d::transport::GenericDiffusionBoundaryPolicy;
    using dec3d::transport::GenericDiffusionReferenceSolveOptions;
    using dec3d::transport::GetCsrValue;
    using dec3d::transport::SolveGenericDiffusionReference;
    using dec3d::transport::ValidateGenericDiffusionAssemblyDiagnostics;
    using dec3d::transport::ValidateGenericDiffusionSolveDiagnostics;

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.0, 2.0, 1.0, 0.25, 3.0, 10.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(assembly.is_complete());
      DEC3D_CHECK(ValidateGenericDiffusionAssemblyDiagnostics(assembly));
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);
      DEC3D_CHECK(solve.is_complete());
      DEC3D_CHECK(ValidateGenericDiffusionSolveDiagnostics(solve));
      DEC3D_CHECK_EQ(solve.scalar_new.size(), std::size_t{1});
      CheckNear(solve.scalar_new[0], 10.0, 0.0, "zero dt returns old scalar");
      DEC3D_CHECK(assembly.report_line.find("dt_s=0") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("canonical_state_mutated=false") != std::string::npos);
    }

    {
      const double dt_s = 0.5;
      const double a = 2.0;
      const double c = 0.25;
      const double b = 3.0;
      const double old_scalar = 10.0;
      auto problem = BuildPatchProblem(Layout(1, 1, 1), dt_s, a, 0.0, c, b, old_scalar);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);
      const double expected = ((a / dt_s) * old_scalar + b) / ((a / dt_s) - c);
      CheckNear(solve.scalar_new[0], expected, 1.0e-12, "one-cell B/C sign oracle");
      DEC3D_CHECK(assembly.report_line.find("interface_diffusion_mean=arithmetic") != std::string::npos);
      DEC3D_CHECK(GetCsrValue(assembly.matrix, 0, 0) > 0.0);
    }

    {
      auto zero_flux =
          BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 2.0, 0.0, 0.0, 8.0);
      const auto zero_flux_assembly = AssembleGenericImplicitDiffusionSystem(zero_flux);
      DEC3D_CHECK(zero_flux_assembly.success);

      auto marshak = zero_flux;
      marshak.boundary_policy = MarshakOuterBoundary();
      const auto marshak_assembly = AssembleGenericImplicitDiffusionSystem(marshak);
      DEC3D_CHECK(marshak_assembly.success);
      DEC3D_CHECK(marshak_assembly.marshak_boundary_used);
      DEC3D_CHECK(marshak_assembly.marshak_outer_face_count > 0u);
      DEC3D_CHECK(marshak_assembly.marshak_nonzero_diagonal_loss_count > 0u);
      DEC3D_CHECK(marshak_assembly.report_line.find(
                      "outer_boundary_policy=radiation_marshak_vacuum") != std::string::npos);
      DEC3D_CHECK(marshak_assembly.report_line.find("marshak_boundary_used=true") !=
                  std::string::npos);
      DEC3D_CHECK(
          GetCsrValue(marshak_assembly.matrix, 0u, 0u) >
          GetCsrValue(zero_flux_assembly.matrix, 0u, 0u));
    }

    {
      const double dt_s = 0.25;
      const double d_face = 3.0;
      const double old_scalar = 12.0;
      auto problem =
          BuildPatchProblem(Layout(1, 1, 1), dt_s, 1.0, d_face, 0.0, 0.0, old_scalar);
      problem.boundary_policy = MarshakOuterBoundary();
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);

      const double outer_area = RadialFaceArea(problem.geometry, 1u, 0u, 0u);
      const double volume = CellVolume(problem.geometry, problem.layout, 0u, 0u, 0u);
      const double distance =
          problem.geometry.radial_faces.back() - RadialCenter(problem.geometry, 0u);
      const double g_face = MarshakGFace(d_face, distance);
      const double expected =
          old_scalar / (1.0 + dt_s * (outer_area / volume) * g_face);
      CheckNear(solve.scalar_new[0], expected, 1.0e-12, "one-cell Marshak loss oracle");
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.0, 1.0, 2.0, 0.0, 0.0, 6.0);
      problem.boundary_policy = MarshakOuterBoundary();
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(assembly.is_complete());
      CheckNear(GetCsrValue(assembly.matrix, 0u, 0u), 1.0, 0.0,
                "dt zero Marshak identity matrix");
      CheckNear(assembly.rhs[0], 6.0, 0.0, "dt zero Marshak rhs");
      DEC3D_CHECK(!assembly.marshak_boundary_used);
      DEC3D_CHECK_EQ(assembly.marshak_outer_face_count, std::size_t{0});
      DEC3D_CHECK_EQ(assembly.marshak_nonzero_diagonal_loss_count, std::size_t{0});
      DEC3D_CHECK(assembly.report_line.find(
                      "outer_boundary_policy=radiation_marshak_vacuum") != std::string::npos);
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);
      CheckNear(solve.scalar_new[0], 6.0, 0.0, "dt zero Marshak solution");
    }

    {
      auto inner = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 2.0, 0.0, 0.0, 1.0);
      inner.boundary_policy.inner_radial =
          DiffusionBoundaryKind::radiation_marshak_vacuum;
      const auto inner_assembly = AssembleGenericImplicitDiffusionSystem(inner);
      DEC3D_CHECK(!inner_assembly.success);
      DEC3D_CHECK(inner_assembly.failure_diagnostics.find(
                      "inner radial boundary policy is unsupported") != std::string::npos);

      auto theta = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 2.0, 0.0, 0.0, 1.0);
      theta.boundary_policy.theta_lower =
          DiffusionBoundaryKind::radiation_marshak_vacuum;
      const auto theta_assembly = AssembleGenericImplicitDiffusionSystem(theta);
      DEC3D_CHECK(!theta_assembly.success);
      DEC3D_CHECK(theta_assembly.failure_diagnostics.find(
                      "theta boundary policy is unsupported") != std::string::npos);

      auto phi = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 2.0, 0.0, 0.0, 1.0);
      phi.boundary_policy.phi = DiffusionBoundaryKind::radiation_marshak_vacuum;
      const auto phi_assembly = AssembleGenericImplicitDiffusionSystem(phi);
      DEC3D_CHECK(!phi_assembly.success);
      DEC3D_CHECK(phi_assembly.failure_diagnostics.find(
                      "phi boundary policy must be periodic") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(2, 1, 1), 0.1, 2.0, 4.0, 0.0, 0.0, 7.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK_EQ(assembly.matrix.row_count, std::size_t{2});
      DEC3D_CHECK_EQ(assembly.matrix.column_count, std::size_t{2});
      DEC3D_CHECK_EQ(assembly.matrix.row_offsets.size(), std::size_t{3});
      DEC3D_CHECK(GetCsrValue(assembly.matrix, 0, 1) < 0.0);
      CheckNear(
          GetCsrValue(assembly.matrix, 0, 1),
          GetCsrValue(assembly.matrix, 1, 0),
          1.0e-12,
          "shared radial face symmetric conductance");
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);
      CheckNear(solve.scalar_new[0], 7.0, 1.0e-12, "constant field cell 0");
      CheckNear(solve.scalar_new[1], 7.0, 1.0e-12, "constant field cell 1");
      DEC3D_CHECK(!assembly.marshak_boundary_used);
      DEC3D_CHECK_EQ(assembly.marshak_outer_face_count, std::size_t{0});
      DEC3D_CHECK(assembly.report_line.find("marshak_boundary_used=false") !=
                  std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(2, 2, 2), 0.1, 1.0, 0.25, 0.0, 0.0, 2.0);
      const auto baseline = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(baseline.success);
      const auto baseline_solution =
          SolveGenericDiffusionReference(baseline, GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(baseline_solution.success);

      auto disabled = problem;
      disabled.face_effective_coefficients.enabled = false;
      const auto with_disabled = AssembleGenericImplicitDiffusionSystem(disabled);
      DEC3D_CHECK(with_disabled.success);
      const auto disabled_solution =
          SolveGenericDiffusionReference(with_disabled, GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(disabled_solution.success);

      DEC3D_CHECK_EQ(with_disabled.matrix.row_offsets.size(), baseline.matrix.row_offsets.size());
      DEC3D_CHECK_EQ(with_disabled.matrix.column_indices.size(), baseline.matrix.column_indices.size());
      DEC3D_CHECK_EQ(with_disabled.matrix.values.size(), baseline.matrix.values.size());
      for (std::size_t i = 0; i < baseline.matrix.row_offsets.size(); ++i) {
        DEC3D_CHECK_EQ(with_disabled.matrix.row_offsets[i], baseline.matrix.row_offsets[i]);
      }
      for (std::size_t i = 0; i < baseline.matrix.column_indices.size(); ++i) {
        DEC3D_CHECK_EQ(with_disabled.matrix.column_indices[i], baseline.matrix.column_indices[i]);
        CheckNear(with_disabled.matrix.values[i], baseline.matrix.values[i], 0.0,
                  "disabled limiter matrix value parity");
      }
      for (std::size_t i = 0; i < baseline.rhs.size(); ++i) {
        CheckNear(with_disabled.rhs[i], baseline.rhs[i], 0.0,
                  "disabled limiter rhs parity");
        CheckNear(disabled_solution.scalar_new[i], baseline_solution.scalar_new[i], 0.0,
                  "disabled limiter solution parity");
      }
    }

    {
      auto problem = BuildPatchProblem(Layout(2, 2, 2), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK_EQ(assembly.row_count, std::size_t{8});
      DEC3D_CHECK_EQ(assembly.matrix.row_offsets.size(), std::size_t{9});
      DEC3D_CHECK(assembly.nonzero_count >= assembly.row_count);
      DEC3D_CHECK(assembly.report_line.find("direction_splitting=false") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.coefficient_A(0, 0, 0) = -1.0;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(assembly.failure_diagnostics.find("coefficient_A must be positive") != std::string::npos);
      DEC3D_CHECK(assembly.failure_diagnostics.find("first_bad_cell=0,0,0") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.boundary_policy = GenericDiffusionBoundaryPolicy{};
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(assembly.failure_diagnostics.find("boundary policy is missing") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.boundary_policy.inner_radial = DiffusionBoundaryKind::neumann_zero_flux;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(
          assembly.failure_diagnostics.find("inner radial boundary policy is unsupported") !=
          std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.boundary_policy.theta_lower = DiffusionBoundaryKind::periodic;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(
          assembly.failure_diagnostics.find("theta boundary policy is unsupported") !=
          std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.geometry = BuildPatchGeometry(Layout(1, 1, 1), 0.0, 1.0, 0.5, 1.2);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(
          assembly.failure_diagnostics.find("patch boundary policy touches origin") !=
          std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      problem.geometry = BuildPatchGeometry(Layout(1, 1, 1), 1.0, 2.0, 0.0, kPi);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(
          assembly.failure_diagnostics.find("patch boundary policy touches pole") !=
          std::string::npos);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(assembly.report_line.find("origin_remap_used=true") != std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("pole_remap_used=true") != std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("singular_boundary_face_conductance_used=false") !=
                  std::string::npos);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      problem.boundary_policy.outer_radial = DiffusionBoundaryKind::radiation_marshak_vacuum;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(assembly.origin_remap_used);
      DEC3D_CHECK(assembly.pole_remap_used);
      DEC3D_CHECK(assembly.marshak_boundary_used);
      DEC3D_CHECK(assembly.report_line.find("origin_remap_used=true") != std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("pole_remap_used=true") != std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("marshak_boundary_used=true") != std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("singular_boundary_face_conductance_used=false") !=
                  std::string::npos);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(
          assembly.report_line.find("pole_metric_mode=axis_regular_polar_phi") !=
          std::string::npos);
      const auto origin_lhs = TestLinearIndex(problem.layout, 0, 1, 2);
      const auto origin_rhs = TestLinearIndex(problem.layout, 0, 2, 6);
      const double origin_value = GetCsrValue(assembly.matrix, origin_lhs, origin_rhs);
      const double origin_symmetric_value = GetCsrValue(assembly.matrix, origin_rhs, origin_lhs);
      DEC3D_CHECK(std::isfinite(origin_value));
      CheckNear(origin_value, 0.0, 0.0, "origin metric closure does not add artificial half-turn edge");
      CheckNear(
          origin_symmetric_value,
          origin_value,
          1.0e-14,
          "origin mapped CSR coupling is symmetric");
      DEC3D_CHECK(
          assembly.report_line.find("origin_remap_coupling_mode=metric_regular_origin") !=
          std::string::npos);

      const auto pole_lhs = TestLinearIndex(problem.layout, 1, 0, 3);
      const auto pole_rhs = TestLinearIndex(problem.layout, 1, 0, 7);
      const double pole_value = GetCsrValue(assembly.matrix, pole_lhs, pole_rhs);
      const double pole_symmetric_value = GetCsrValue(assembly.matrix, pole_rhs, pole_lhs);
      DEC3D_CHECK(std::isfinite(pole_value));
      CheckNear(pole_value, 0.0, 0.0, "pole metric closure does not add artificial opposite-phi edge");
      CheckNear(
          pole_symmetric_value,
          pole_value,
          1.0e-14,
          "pole mapped CSR coupling is symmetric");
      DEC3D_CHECK(
          assembly.report_line.find("pole_remap_coupling_mode=metric_regular_axis") !=
          std::string::npos);
    }

    {
      auto axis_regular_problem =
          BuildFullSphereProblem(Layout(1, 16, 32), 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
      axis_regular_problem.geometry =
          BuildPatchGeometry(axis_regular_problem.layout, 1.0, 2.0, 0.0, kPi);
      axis_regular_problem.boundary_policy.inner_radial =
          DiffusionBoundaryKind::interior_patch_no_origin;
      const auto axis_regular = AssembleGenericImplicitDiffusionSystem(axis_regular_problem);
      DEC3D_CHECK(axis_regular.success);
      DEC3D_CHECK(
          axis_regular.report_line.find("pole_metric_mode=axis_regular_polar_phi") !=
          std::string::npos);

      auto legacy_problem = axis_regular_problem;
      legacy_problem.pole_metric_mode =
          dec3d::transport::DiffusionPoleMetricMode::cell_centered_spherical;
      const auto legacy = AssembleGenericImplicitDiffusionSystem(legacy_problem);
      DEC3D_CHECK(legacy.success);
      DEC3D_CHECK(
          legacy.report_line.find("pole_metric_mode=cell_centered_spherical") !=
          std::string::npos);

      const auto phi_lhs = TestLinearIndex(axis_regular_problem.layout, 0, 0, 0);
      const auto phi_rhs = TestLinearIndex(axis_regular_problem.layout, 0, 0, 1);
      const double legacy_polar_phi = GetCsrValue(legacy.matrix, phi_lhs, phi_rhs);
      const double axis_regular_polar_phi =
          GetCsrValue(axis_regular.matrix, phi_lhs, phi_rhs);
      DEC3D_CHECK(legacy_polar_phi < 0.0);
      DEC3D_CHECK(axis_regular_polar_phi < 0.0);
      DEC3D_CHECK(std::abs(axis_regular_polar_phi) < std::abs(legacy_polar_phi));
    }

    {
      auto problem = BuildFullSphereProblem(Layout(8, 16, 32), 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto scalar = CartesianXScalarAtCenters(problem.layout, problem.geometry);
      double max_origin_error = 0.0;
      std::size_t max_origin_theta = 0;
      std::size_t max_origin_phi = 0;
      for (std::size_t theta = 1; theta + 1u < problem.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
          const double numeric =
              DiffusionOperatorFromCsr(assembly, problem, scalar, 0u, theta, phi);
          if (std::abs(numeric) > max_origin_error) {
            max_origin_error = std::abs(numeric);
            max_origin_theta = theta;
            max_origin_phi = phi;
          }
        }
      }
      if (!(max_origin_error < 0.1)) {
        std::cerr << "origin harmonic max error=" << max_origin_error
                  << " theta=" << max_origin_theta << " phi=" << max_origin_phi
                  << '\n';
      }
      DEC3D_CHECK(max_origin_error < 0.1);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(1, 16, 32), 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
      problem.geometry = BuildPatchGeometry(problem.layout, 1.0, 2.0, 0.0, kPi);
      problem.boundary_policy.inner_radial = DiffusionBoundaryKind::interior_patch_no_origin;
      problem.pole_metric_mode =
          dec3d::transport::DiffusionPoleMetricMode::cell_centered_spherical;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto scalar = L1M1ScalarAtCenters(problem.layout, problem.geometry);
      double max_pole_error = 0.0;
      double max_second_band_error = 0.0;
      for (std::size_t phi = 0; phi < problem.layout.phi_cells; ++phi) {
        const double pole_numeric =
            DiffusionOperatorFromCsr(assembly, problem, scalar, 0u, 0u, phi);
        const double pole_exact =
            ExactL1M1FiniteVolumeLaplacian(problem.geometry, 0u, 0u, phi);
        max_pole_error = std::max(max_pole_error, std::abs(pole_numeric - pole_exact));

        const double second_numeric =
            DiffusionOperatorFromCsr(assembly, problem, scalar, 0u, 1u, phi);
        const double second_exact =
            ExactL1M1FiniteVolumeLaplacian(problem.geometry, 0u, 1u, phi);
        max_second_band_error =
            std::max(max_second_band_error, std::abs(second_numeric - second_exact));
      }
      DEC3D_CHECK(max_pole_error < 0.02);
      DEC3D_CHECK(max_second_band_error < 0.01);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 5.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-11});
      DEC3D_CHECK(solve.success);
      for (double value : solve.scalar_new) {
        CheckNear(value, 5.0, 1.0e-10, "full-sphere constant scalar remains constant");
      }
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 7), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(assembly.failure_diagnostics.find("even phi") != std::string::npos);
      DEC3D_CHECK(assembly.failure_diagnostics.find("canonical_state_mutated=false") !=
                  std::string::npos);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      problem.boundary_policy.inner_radial = DiffusionBoundaryKind::interior_patch_no_origin;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(assembly.failure_diagnostics.find("patch boundary policy touches origin") !=
                  std::string::npos);
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      problem.boundary_policy.theta_lower = DiffusionBoundaryKind::interior_patch_no_pole;
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(!assembly.success);
      DEC3D_CHECK(assembly.failure_diagnostics.find("patch boundary policy touches pole") !=
                  std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(5, 5, 3), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(!solve.success);
      DEC3D_CHECK(solve.failure_diagnostics.find("reference backend row_count exceeds limit") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string token = "marshak_boundary_used=false";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string token = "interface_diffusion_mean=arithmetic";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string token = "direction_splitting=false";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string token = "origin_remap_nonzero_offdiag_count=";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string pole_metric_token = "pole_metric_mode=";
      const auto pole_metric_pos = assembly.report_line.find(pole_metric_token);
      DEC3D_CHECK(pole_metric_pos != std::string::npos);
      assembly.report_line.erase(pole_metric_pos, pole_metric_token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildFullSphereProblem(Layout(2, 4, 8), 0.1, 1.0, 2.0, 0.0, 0.0, 3.0);
      auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const std::string token = "singular_boundary_face_conductance_used=false";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      auto solve = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(solve.success);
      const std::string token = "residual_linf=";
      const auto pos = solve.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      solve.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionSolveDiagnostics(solve));
    }

    std::cout << "P2-2 generic implicit diffusion substrate contract checks passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
