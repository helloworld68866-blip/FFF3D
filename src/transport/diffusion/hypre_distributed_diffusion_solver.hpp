#pragma once

#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "transport/diffusion/generic_diffusion.hpp"
#include "transport/diffusion/hypre_diffusion_solver.hpp"

#include <mpi.h>

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::transport {

struct DistributedDiffusionRowOwnership {
  bool success{false};
  MPI_Comm communicator{MPI_COMM_WORLD};
  int rank{0};
  int rank_count{1};
  std::size_t global_radial_cells{0};
  std::size_t global_theta_cells{0};
  std::size_t global_phi_cells{0};
  std::size_t global_row_count{0};
  std::size_t global_radial_begin{0};
  std::size_t global_radial_end{0};
  std::size_t local_row_begin{0};
  std::size_t local_row_end{0};
  std::size_t local_row_count{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
};

struct DistributedLocalCsrMatrix {
  std::size_t global_row_count{0};
  std::size_t local_row_begin{0};
  std::size_t local_row_end{0};
  std::vector<std::size_t> row_offsets;
  std::vector<std::size_t> column_indices;
  std::vector<double> values;

  [[nodiscard]] bool is_shape_complete() const noexcept;
};

struct DistributedGenericDiffusionProblem {
  DistributedDiffusionRowOwnership ownership;
  dec3d::mesh::SphericalGeometryMetadata global_geometry;
  double dt_s{0.0};
  dec3d::core::Array3D<double> local_coefficient_A;
  dec3d::core::Array3D<double> local_coefficient_D;
  dec3d::core::Array3D<double> local_coefficient_C;
  dec3d::core::Array3D<double> local_coefficient_B;
  dec3d::core::Array3D<double> local_scalar_old;
  GenericDiffusionBoundaryPolicy boundary_policy;
  GenericDiffusionFaceEffectiveCoefficients face_effective_coefficients;
  DiffusionPoleMetricMode pole_metric_mode{
      DiffusionPoleMetricMode::axis_regular_polar_phi};
};

struct DistributedGenericDiffusionAssemblyResult {
  bool success{false};
  DistributedLocalCsrMatrix local_matrix;
  std::vector<double> local_rhs;
  std::vector<double> local_scalar_old_flat;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  DistributedDiffusionRowOwnership ownership;
  std::size_t local_nonzero_count{0};
  std::size_t local_off_rank_column_count{0};
  std::size_t local_radial_seam_coupling_count{0};
  std::size_t global_nonzero_count{0};
  std::size_t global_off_rank_column_count{0};
  std::size_t global_radial_seam_coupling_count{0};
  bool coefficient_D_halo_lower_received{false};
  bool coefficient_D_halo_upper_received{false};
  bool local_origin_remap_used{false};
  bool local_pole_remap_used{false};
  bool global_origin_remap_used{false};
  bool global_pole_remap_used{false};
  bool singular_boundary_face_conductance_used{false};
  DiffusionPoleMetricMode pole_metric_mode{
      DiffusionPoleMetricMode::axis_regular_polar_phi};

  [[nodiscard]] bool is_complete() const noexcept;
};

struct DistributedGenericDiffusionHypreSolveResult {
  bool success{false};
  std::vector<double> local_scalar_new;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string backend;
  bool hypre_enabled{false};
  DistributedDiffusionRowOwnership ownership;
  std::size_t local_nonzero_count{0};
  std::size_t global_nonzero_count{0};
  std::size_t local_off_rank_column_count{0};
  std::size_t global_off_rank_column_count{0};
  int gmres_iterations{0};
  double gmres_final_relative_residual{0.0};
  double hypre_setup_wall_s{0.0};
  double hypre_solve_wall_s{0.0};
  double local_residual_l2{0.0};
  double global_residual_l2{0.0};
  double local_residual_linf{0.0};
  double global_residual_linf{0.0};
  double local_rhs_linf{0.0};
  double global_rhs_linf{0.0};
  double global_residual_linf_relative{0.0};
  double max_abs_delta_local{0.0};
  double max_abs_delta_global{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] DistributedDiffusionRowOwnership BuildDistributedDiffusionRowOwnership(
    MPI_Comm communicator,
    std::size_t global_radial_cells,
    std::size_t global_theta_cells,
    std::size_t global_phi_cells) noexcept;

[[nodiscard]] DistributedGenericDiffusionAssemblyResult AssembleDistributedGenericDiffusionSystem(
    const DistributedGenericDiffusionProblem& problem) noexcept;

[[nodiscard]] DistributedGenericDiffusionHypreSolveResult SolveDistributedGenericDiffusionHypre(
    const DistributedGenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionHypreSolveOptions& options) noexcept;

[[nodiscard]] std::vector<double> GatherDistributedScalarForTest(
    const DistributedGenericDiffusionHypreSolveResult& result,
    MPI_Comm communicator);

[[nodiscard]] bool ValidateDistributedGenericDiffusionAssemblyDiagnostics(
    const DistributedGenericDiffusionAssemblyResult& result) noexcept;

[[nodiscard]] bool ValidateDistributedGenericDiffusionHypreSolveDiagnostics(
    const DistributedGenericDiffusionHypreSolveResult& result) noexcept;

}  // namespace dec3d::transport
