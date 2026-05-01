#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/generic_diffusion.hpp"
#include "transport/diffusion/hypre_diffusion_solver.hpp"

#include <mpi.h>

#include <cmath>
#include <exception>
#include <iostream>
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

void CheckHypreSuccess(
    const dec3d::transport::GenericDiffusionHypreSolveResult& result,
    const char* label,
    int line) {
  if (!result.success) {
    std::ostringstream detail;
    detail << label << " failure_diagnostics=" << result.failure_diagnostics;
    dec3d::test::Fail("HYPRE solve success", __FILE__, line, detail.str());
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

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    using dec3d::transport::AssembleGenericImplicitDiffusionSystem;
    using dec3d::transport::GenericDiffusionHypreSolveOptions;
    using dec3d::transport::GenericDiffusionReferenceSolveOptions;
    using dec3d::transport::SolveGenericDiffusionHypre;
    using dec3d::transport::SolveGenericDiffusionReference;
    using dec3d::transport::ValidateGenericDiffusionHypreSolveDiagnostics;

    const bool world_reject_mode = argc > 1 && std::string(argv[1]) == "--expect-world-reject";
    if (world_reject_mode) {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      auto options = GenericDiffusionHypreSolveOptions{};
      options.communicator = MPI_COMM_WORLD;
      const auto hypre = SolveGenericDiffusionHypre(assembly, options);
      int world_size = 1;
      MPI_Comm_size(MPI_COMM_WORLD, &world_size);
      if (world_size > 1) {
        DEC3D_CHECK(!hypre.success);
        DEC3D_CHECK(
            hypre.failure_diagnostics.find(
                "distributed HYPRE row ownership is out of scope for P2-2b") !=
            std::string::npos);
      }
      MPI_Finalize();
      return 0;
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
      const auto solve = SolveGenericDiffusionHypre(assembly, GenericDiffusionHypreSolveOptions{});
      CheckHypreSuccess(solve, "one-cell oracle", __LINE__);
      DEC3D_CHECK(solve.is_complete());
      DEC3D_CHECK(ValidateGenericDiffusionHypreSolveDiagnostics(solve));
      const double expected = ((a / dt_s) * old_scalar + b) / ((a / dt_s) - c);
      CheckNear(solve.scalar_new[0], expected, 1.0e-10, "HYPRE one-cell B/C oracle");
      DEC3D_CHECK(solve.report_line.find("backend=hypre_parcsr_gmres_boomeramg") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("hypre_enabled=true") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("ij_matrix_created=true") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("ij_vector_created=true") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("parcsr_matrix_created=true") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("solver=ParCSRGMRES") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("preconditioner=BoomerAMG") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find(
                      "boomeramg_parameters=preconditioner_max_iter_1_tol_0") !=
                  std::string::npos);
      DEC3D_CHECK(solve.report_line.find("boomeramg_max_iterations=1") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("boomeramg_tolerance=0") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("rhs_linf=") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("residual_linf_relative=") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("fallback_used=false") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(2, 1, 1), 0.1, 2.0, 4.0, 0.0, 0.0, 7.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto reference = SolveGenericDiffusionReference(
          assembly,
          GenericDiffusionReferenceSolveOptions{64, 1.0e-12});
      DEC3D_CHECK(reference.success);
      const auto hypre = SolveGenericDiffusionHypre(assembly, GenericDiffusionHypreSolveOptions{});
      CheckHypreSuccess(hypre, "two-cell reference parity", __LINE__);
      DEC3D_CHECK_EQ(hypre.scalar_new.size(), reference.scalar_new.size());
      for (std::size_t i = 0; i < hypre.scalar_new.size(); ++i) {
        CheckNear(hypre.scalar_new[i], reference.scalar_new[i], 1.0e-9, "HYPRE/reference parity");
      }
      DEC3D_CHECK(hypre.report_line.find("row_ownership=whole_matrix_single_rank") != std::string::npos);
      DEC3D_CHECK(hypre.report_line.find("mpi_rank_count=1") != std::string::npos);
      DEC3D_CHECK(hypre.report_line.find("local_row_begin=0") != std::string::npos);
      DEC3D_CHECK(hypre.report_line.find("local_row_end=2") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(2, 2, 2), 0.1, 1.0, 1.0, 0.0, 0.0, 4.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      DEC3D_CHECK(assembly.success);
      const auto hypre = SolveGenericDiffusionHypre(assembly, GenericDiffusionHypreSolveOptions{});
      CheckHypreSuccess(hypre, "constant field", __LINE__);
      for (double value : hypre.scalar_new) {
        CheckNear(value, 4.0, 1.0e-9, "constant field no-change");
      }
    }

    {
      auto failed_assembly = dec3d::transport::GenericDiffusionAssemblyResult{};
      const auto hypre = SolveGenericDiffusionHypre(failed_assembly, GenericDiffusionHypreSolveOptions{});
      DEC3D_CHECK(!hypre.success);
      DEC3D_CHECK(hypre.failure_diagnostics.find("assembly is incomplete or failed") != std::string::npos);
      DEC3D_CHECK(hypre.failure_diagnostics.find("fallback_used=false") != std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      auto options = GenericDiffusionHypreSolveOptions{};
      options.max_iterations = 0;
      const auto hypre = SolveGenericDiffusionHypre(assembly, options);
      DEC3D_CHECK(!hypre.success);
      DEC3D_CHECK(hypre.failure_diagnostics.find("HYPRE solver options are invalid") != std::string::npos);
      DEC3D_CHECK(hypre.failure_diagnostics.find("backend=serial_dense_reference") == std::string::npos);
    }

    {
      auto problem = BuildPatchProblem(Layout(1, 1, 1), 0.1, 1.0, 0.0, 0.0, 0.0, 1.0);
      const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
      auto hypre = SolveGenericDiffusionHypre(assembly, GenericDiffusionHypreSolveOptions{});
      CheckHypreSuccess(hypre, "missing diagnostics guard", __LINE__);
      const std::string token = "hypre_enabled=true";
      const auto pos = hypre.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      hypre.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateGenericDiffusionHypreSolveDiagnostics(hypre));
    }

    std::cout << "P2-2b HYPRE diffusion solver contract passed\n";
    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << "\n";
    MPI_Finalize();
    return 1;
  }
}
