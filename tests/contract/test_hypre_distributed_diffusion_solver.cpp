#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/generic_diffusion.hpp"
#include "transport/diffusion/hypre_diffusion_solver.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void CheckNear(double lhs, double rhs, double tolerance, const std::string& label) {
  if (std::abs(lhs - rhs) > tolerance) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("near equality", __FILE__, __LINE__, detail.str());
  }
}

double CsrValue(
    const dec3d::transport::DistributedLocalCsrMatrix& matrix,
    std::size_t local_row,
    std::size_t global_column) {
  for (std::size_t slot = matrix.row_offsets[local_row];
       slot < matrix.row_offsets[local_row + 1u];
       ++slot) {
    if (matrix.column_indices[slot] == global_column) {
      return matrix.values[slot];
    }
  }
  return 0.0;
}

dec3d::transport::DiffusionGridLayout Layout(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return dec3d::transport::DiffusionGridLayout{radial, theta, phi};
}

dec3d::mesh::SphericalGeometryMetadata BuildGeometry(
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

dec3d::transport::GenericDiffusionBoundaryPolicy PatchBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::interior_patch_no_origin,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
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

dec3d::core::Array3D<double> Filled(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    double value) {
  return dec3d::core::Array3D<double>(radial, theta, phi, value);
}

dec3d::transport::GenericDiffusionProblem BuildWholePatchProblem(
    const dec3d::transport::DiffusionGridLayout& layout,
    double dt_s,
    double scalar_old) {
  return dec3d::transport::GenericDiffusionProblem{
      layout,
      BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2),
      dt_s,
      Filled(layout.radial_cells, layout.theta_cells, layout.phi_cells, 1.0),
      Filled(layout.radial_cells, layout.theta_cells, layout.phi_cells, 2.0),
      Filled(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0),
      Filled(layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0),
      Filled(layout.radial_cells, layout.theta_cells, layout.phi_cells, scalar_old),
      PatchBoundary()};
}

dec3d::transport::DistributedGenericDiffusionProblem BuildDistributedProblem(
    MPI_Comm communicator,
    const dec3d::transport::DiffusionGridLayout& global_layout,
    const dec3d::mesh::SphericalGeometryMetadata& global_geometry,
    double dt_s,
    double scalar_old,
    dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy) {
  const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
      communicator,
      global_layout.radial_cells,
      global_layout.theta_cells,
      global_layout.phi_cells);
  const std::size_t local_radial = ownership.global_radial_end - ownership.global_radial_begin;
  return dec3d::transport::DistributedGenericDiffusionProblem{
      ownership,
      global_geometry,
      dt_s,
      Filled(local_radial, global_layout.theta_cells, global_layout.phi_cells, 1.0),
      Filled(local_radial, global_layout.theta_cells, global_layout.phi_cells, 2.0),
      Filled(local_radial, global_layout.theta_cells, global_layout.phi_cells, 0.0),
      Filled(local_radial, global_layout.theta_cells, global_layout.phi_cells, 0.0),
      Filled(local_radial, global_layout.theta_cells, global_layout.phi_cells, scalar_old),
      boundary_policy};
}

void RequireTwoRanks() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    dec3d::test::Fail("rank count", __FILE__, __LINE__, "P2-6b contract requires mpiexec -n 2");
  }
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    RequireTwoRanks();

    using dec3d::transport::AssembleDistributedGenericDiffusionSystem;
    using dec3d::transport::AssembleGenericImplicitDiffusionSystem;
    using dec3d::transport::GatherDistributedScalarForTest;
    using dec3d::transport::SolveDistributedGenericDiffusionHypre;
    using dec3d::transport::SolveGenericDiffusionHypre;
    using dec3d::transport::ValidateDistributedGenericDiffusionAssemblyDiagnostics;
    using dec3d::transport::ValidateDistributedGenericDiffusionHypreSolveDiagnostics;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    {
      const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
          MPI_COMM_WORLD,
          4,
          2,
          2);
      DEC3D_CHECK(ownership.success);
      DEC3D_CHECK_EQ(ownership.rank_count, 2);
      DEC3D_CHECK_EQ(ownership.global_row_count, std::size_t{16});
      DEC3D_CHECK_EQ(ownership.local_row_count, std::size_t{8});
      DEC3D_CHECK(ownership.report_line.find("row_ownership=distributed_radial_slab") !=
                  std::string::npos);
      DEC3D_CHECK(ownership.report_line.find("collective_row_ownership_valid=true") !=
                  std::string::npos);
      DEC3D_CHECK(ownership.report_line.find("global_layout_consistent=true") !=
                  std::string::npos);
    }

    {
      // distributed_row_ownership_rejects_inconsistent_global_radial_layout
      const std::size_t rank_local_radial = rank == 0 ? 4u : 5u;
      const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
          MPI_COMM_WORLD,
          rank_local_radial,
          2,
          2);
      DEC3D_CHECK(!ownership.success);
      DEC3D_CHECK(ownership.failure_diagnostics.find("inconsistent global layout") !=
                  std::string::npos);
      DEC3D_CHECK(ownership.failure_diagnostics.find("collective_row_ownership_valid=false") !=
                  std::string::npos);
    }

    {
      // distributed_row_ownership_rejects_inconsistent_global_angular_layout
      const std::size_t rank_local_theta = rank == 0 ? 2u : 3u;
      const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
          MPI_COMM_WORLD,
          4,
          rank_local_theta,
          2);
      DEC3D_CHECK(!ownership.success);
      DEC3D_CHECK(ownership.failure_diagnostics.find("inconsistent global layout") !=
                  std::string::npos);
      DEC3D_CHECK(ownership.failure_diagnostics.find("collective_row_ownership_valid=false") !=
                  std::string::npos);
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, PatchBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(ValidateDistributedGenericDiffusionAssemblyDiagnostics(assembly));
      DEC3D_CHECK(assembly.global_radial_seam_coupling_count > 0);
      DEC3D_CHECK(assembly.global_off_rank_column_count > 0);
      DEC3D_CHECK(assembly.report_line.find("local_off_rank_column_count=") != std::string::npos);
      DEC3D_CHECK(
          assembly.report_line.find("pole_metric_mode=axis_regular_polar_phi") !=
          std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("coefficient_D_halo_lower_received=") !=
                  std::string::npos);
      DEC3D_CHECK(assembly.report_line.find("coefficient_D_halo_upper_received=") !=
                  std::string::npos);
    }

    {
      // distributed_marshak_only_on_true_outer_global_rank
      using dec3d::transport::DiffusionBoundaryKind;
      const auto layout = Layout(4, 1, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      auto zero_boundary = PatchBoundary();
      zero_boundary.outer_radial = DiffusionBoundaryKind::neumann_zero_flux;
      auto marshak_boundary = PatchBoundary();
      marshak_boundary.outer_radial = DiffusionBoundaryKind::radiation_marshak_vacuum;

      const auto zero_problem =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, zero_boundary);
      const auto marshak_problem =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, marshak_boundary);
      const auto zero_assembly = AssembleDistributedGenericDiffusionSystem(zero_problem);
      const auto marshak_assembly = AssembleDistributedGenericDiffusionSystem(marshak_problem);
      DEC3D_CHECK(zero_assembly.success);
      DEC3D_CHECK(marshak_assembly.success);

      if (rank == 0) {
        for (std::size_t row = 0; row < zero_assembly.ownership.local_row_count; ++row) {
          const std::size_t global_row = zero_assembly.ownership.local_row_begin + row;
          CheckNear(
              CsrValue(marshak_assembly.local_matrix, row, global_row),
              CsrValue(zero_assembly.local_matrix, row, global_row),
              0.0,
              "Marshak must not alter non-outer rank diagonal");
        }
      } else {
        const std::size_t local_radial =
            marshak_assembly.ownership.global_radial_end -
            marshak_assembly.ownership.global_radial_begin;
        const std::size_t outer_local_row =
            (local_radial - 1u) * layout.theta_cells * layout.phi_cells;
        const std::size_t outer_global_row =
            marshak_assembly.ownership.local_row_begin + outer_local_row;
        const double zero_diag =
            CsrValue(zero_assembly.local_matrix, outer_local_row, outer_global_row);
        const double marshak_diag =
            CsrValue(marshak_assembly.local_matrix, outer_local_row, outer_global_row);
        DEC3D_CHECK(marshak_diag > zero_diag);
      }
    }

    {
      // distributed_marshak_uses_outer_radial_face_effective_D
      using dec3d::transport::DiffusionBoundaryKind;
      const auto layout = Layout(4, 1, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      auto boundary = PatchBoundary();
      boundary.outer_radial = DiffusionBoundaryKind::radiation_marshak_vacuum;
      const auto unlimited =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, boundary);
      auto limited =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, boundary);
      const std::size_t local_radial =
          limited.ownership.global_radial_end - limited.ownership.global_radial_begin;
      limited.face_effective_coefficients.enabled = true;
      limited.face_effective_coefficients.radial_face_D =
          dec3d::core::Array3D<double>(
              local_radial + 1u,
              limited.ownership.global_theta_cells,
              limited.ownership.global_phi_cells,
              2.0);
      limited.face_effective_coefficients.theta_face_D =
          dec3d::core::Array3D<double>(
              local_radial,
              limited.ownership.global_theta_cells + 1u,
              limited.ownership.global_phi_cells,
              2.0);
      limited.face_effective_coefficients.phi_face_D =
          dec3d::core::Array3D<double>(
              local_radial,
              limited.ownership.global_theta_cells,
              limited.ownership.global_phi_cells,
              2.0);
      if (limited.ownership.global_radial_end ==
          limited.ownership.global_radial_cells) {
        for (std::size_t t = 0; t < limited.ownership.global_theta_cells; ++t) {
          for (std::size_t p = 0; p < limited.ownership.global_phi_cells; ++p) {
            limited.face_effective_coefficients.radial_face_D(local_radial, t, p) = 0.25;
          }
        }
      }

      const auto unlimited_assembly = AssembleDistributedGenericDiffusionSystem(unlimited);
      const auto limited_assembly = AssembleDistributedGenericDiffusionSystem(limited);
      DEC3D_CHECK(unlimited_assembly.success);
      DEC3D_CHECK(limited_assembly.success);
      if (limited.ownership.global_radial_end ==
          limited.ownership.global_radial_cells) {
        const std::size_t local_outer = local_radial - 1u;
        const std::size_t row =
            (local_outer * limited.ownership.global_theta_cells) *
            limited.ownership.global_phi_cells;
        const std::size_t global_row = limited.ownership.local_row_begin + row;
        const double unlimited_diag =
            CsrValue(unlimited_assembly.local_matrix, row, global_row);
        const double limited_diag =
            CsrValue(limited_assembly.local_matrix, row, global_row);
        DEC3D_CHECK(limited_diag < unlimited_diag);
      }
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto whole = BuildWholePatchProblem(layout, 0.1, 4.0);
      const auto whole_assembly = AssembleGenericImplicitDiffusionSystem(whole);
      DEC3D_CHECK(whole_assembly.success);
      const auto whole_solve = SolveGenericDiffusionHypre(
          whole_assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{});
      DEC3D_CHECK(whole_solve.success);

      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 4.0, PatchBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      const auto solve = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{});
      DEC3D_CHECK(solve.success);
      DEC3D_CHECK(ValidateDistributedGenericDiffusionHypreSolveDiagnostics(solve));
      DEC3D_CHECK(solve.report_line.find("row_ownership=distributed_radial_slab") !=
                  std::string::npos);
      DEC3D_CHECK(solve.report_line.find("gathered_solve_used=false") != std::string::npos);
      DEC3D_CHECK(solve.report_line.find("fallback_used=false") != std::string::npos);

      const auto gathered = GatherDistributedScalarForTest(solve, MPI_COMM_WORLD);
      if (rank == 0) {
        DEC3D_CHECK_EQ(gathered.size(), whole_solve.scalar_new.size());
        for (std::size_t i = 0; i < gathered.size(); ++i) {
          CheckNear(gathered[i], whole_solve.scalar_new[i], 1.0e-9, "distributed/single parity");
        }
      }
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 4.0, PatchBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      dec3d::transport::DistributedLaggedBoomerAmgCache cache;
      auto lagged = dec3d::transport::DistributedLaggedBoomerAmgSolveOptions{};
      lagged.enabled = true;
      lagged.group_index = 0u;
      lagged.rebuild_every = 4;
      lagged.max_matrix_relative_change = 0.1;
      lagged.max_iteration_growth = 2.0;
      const auto first = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{},
          &cache,
          lagged);
      DEC3D_CHECK(first.success);
      DEC3D_CHECK(first.lagged_amg_rebuild_used);
      DEC3D_CHECK(!first.lagged_amg_reuse_accepted);
      const auto second = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{},
          &cache,
          lagged);
      DEC3D_CHECK(second.success);
      DEC3D_CHECK(second.lagged_amg_candidate);
      DEC3D_CHECK(second.lagged_amg_reuse_attempted);
      DEC3D_CHECK(second.lagged_amg_reuse_accepted);
      DEC3D_CHECK_EQ(second.lagged_amg_cached_reuse_count, 0);
      DEC3D_CHECK_EQ(second.lagged_amg_reuse_count_after, 1);
      DEC3D_CHECK(ValidateDistributedGenericDiffusionHypreSolveDiagnostics(second));
      lagged.rebuild_every = 2;
      const auto third = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{},
          &cache,
          lagged);
      DEC3D_CHECK(third.success);
      DEC3D_CHECK(!third.lagged_amg_candidate);
      DEC3D_CHECK(!third.lagged_amg_reuse_attempted);
      DEC3D_CHECK(third.lagged_amg_rebuild_used);
      DEC3D_CHECK_EQ(third.lagged_amg_cached_reuse_count, 1);
      DEC3D_CHECK_EQ(third.lagged_amg_reuse_count_after, 0);
      DEC3D_CHECK(ValidateDistributedGenericDiffusionHypreSolveDiagnostics(third));
      const auto first_gathered = GatherDistributedScalarForTest(first, MPI_COMM_WORLD);
      const auto second_gathered = GatherDistributedScalarForTest(second, MPI_COMM_WORLD);
      const auto third_gathered = GatherDistributedScalarForTest(third, MPI_COMM_WORLD);
      if (rank == 0) {
        DEC3D_CHECK_EQ(first_gathered.size(), second_gathered.size());
        DEC3D_CHECK_EQ(first_gathered.size(), third_gathered.size());
        for (std::size_t i = 0; i < first_gathered.size(); ++i) {
          CheckNear(first_gathered[i], second_gathered[i], 1.0e-9, "lagged AMG reuse parity");
          CheckNear(first_gathered[i], third_gathered[i], 1.0e-9, "lagged AMG rebuild parity");
        }
      }
    }

    {
      // distributed_full_sphere_two_rank_origin_remap_and_radial_seam_coexist
      const auto layout = Layout(4, 4, 4);
      const auto geometry = BuildGeometry(layout, 0.0, 2.0, 0.0, kPi);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 5.0, FullSphereBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      DEC3D_CHECK(assembly.global_radial_seam_coupling_count > 0);
      DEC3D_CHECK(assembly.global_off_rank_column_count > 0);
      DEC3D_CHECK(assembly.global_origin_remap_used);
      DEC3D_CHECK(assembly.global_pole_remap_used);
      const auto solve = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{});
      DEC3D_CHECK(solve.success);
      const auto gathered = GatherDistributedScalarForTest(solve, MPI_COMM_WORLD);
      if (rank == 0) {
        for (double value : gathered) {
          CheckNear(value, 5.0, 1.0e-9, "full-sphere constant scalar remains constant");
        }
      }
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, PatchBoundary());
      auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      const std::string token = "global_off_rank_column_count=";
      const auto pos = assembly.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      assembly.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateDistributedGenericDiffusionAssemblyDiagnostics(assembly));
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, PatchBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      auto solve = SolveDistributedGenericDiffusionHypre(
          assembly,
          dec3d::transport::GenericDiffusionHypreSolveOptions{});
      DEC3D_CHECK(solve.success);
      const std::string token = "gathered_solve_used=false";
      const auto pos = solve.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      solve.report_line.erase(pos, token.size());
      DEC3D_CHECK(!ValidateDistributedGenericDiffusionHypreSolveDiagnostics(solve));
    }

    {
      const auto layout = Layout(4, 2, 2);
      const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
      const auto distributed =
          BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 3.0, PatchBoundary());
      const auto assembly = AssembleDistributedGenericDiffusionSystem(distributed);
      DEC3D_CHECK(assembly.success);
      auto options = dec3d::transport::GenericDiffusionHypreSolveOptions{};
      options.max_iterations = 0;
      const auto solve = SolveDistributedGenericDiffusionHypre(assembly, options);
      DEC3D_CHECK(!solve.success);
      DEC3D_CHECK(solve.failure_diagnostics.find("HYPRE solver options are invalid") !=
                  std::string::npos);
      DEC3D_CHECK(solve.failure_diagnostics.find("fallback_used=false") != std::string::npos);
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
