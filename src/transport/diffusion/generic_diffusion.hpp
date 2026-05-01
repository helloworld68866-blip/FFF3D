#pragma once

#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::transport {

struct DiffusionGridLayout {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};

  [[nodiscard]] bool is_valid() const noexcept {
    return radial_cells > 0 && theta_cells > 0 && phi_cells > 0;
  }

  [[nodiscard]] std::size_t cell_count() const noexcept {
    return radial_cells * theta_cells * phi_cells;
  }
};

enum class DiffusionBoundaryKind {
  missing,
  neumann_zero_flux,
  periodic,
  scalar_origin_remap_required,
  scalar_pole_remap_required,
  interior_patch_no_origin,
  interior_patch_no_pole,
  radiation_marshak_vacuum
};

enum class DiffusionPoleMetricMode {
  cell_centered_spherical,
  axis_regular_polar_phi
};

struct GenericDiffusionBoundaryPolicy {
  DiffusionBoundaryKind inner_radial{DiffusionBoundaryKind::missing};
  DiffusionBoundaryKind outer_radial{DiffusionBoundaryKind::missing};
  DiffusionBoundaryKind theta_lower{DiffusionBoundaryKind::missing};
  DiffusionBoundaryKind theta_upper{DiffusionBoundaryKind::missing};
  DiffusionBoundaryKind phi{DiffusionBoundaryKind::missing};
};

struct GenericDiffusionFaceEffectiveCoefficients {
  bool enabled{false};
  dec3d::core::Array3D<double> radial_face_D;
  dec3d::core::Array3D<double> theta_face_D;
  dec3d::core::Array3D<double> phi_face_D;
};

struct SparseMatrixCsr {
  std::size_t row_count{0};
  std::size_t column_count{0};
  std::vector<std::size_t> row_offsets;
  std::vector<std::size_t> column_indices;
  std::vector<double> values;

  [[nodiscard]] bool is_shape_complete() const noexcept;
};

struct GenericDiffusionProblem {
  DiffusionGridLayout layout;
  dec3d::mesh::SphericalGeometryMetadata geometry;
  double dt_s{0.0};
  dec3d::core::Array3D<double> coefficient_A;
  dec3d::core::Array3D<double> coefficient_D;
  dec3d::core::Array3D<double> coefficient_C;
  dec3d::core::Array3D<double> coefficient_B;
  dec3d::core::Array3D<double> scalar_old;
  GenericDiffusionBoundaryPolicy boundary_policy;
  GenericDiffusionFaceEffectiveCoefficients face_effective_coefficients;
  DiffusionPoleMetricMode pole_metric_mode{
      DiffusionPoleMetricMode::axis_regular_polar_phi};
};

struct GenericDiffusionAssemblyResult {
  bool success{false};
  SparseMatrixCsr matrix;
  std::vector<double> rhs;
  std::vector<double> scalar_old_flat;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::size_t row_count{0};
  std::size_t nonzero_count{0};
  std::size_t boundary_row_count{0};
  double min_diagonal{0.0};
  double max_abs_offdiagonal_row_sum{0.0};
  double max_rhs_abs{0.0};
  bool scalar_origin_remap_available{false};
  bool scalar_pole_remap_available{false};
  bool origin_remap_used{false};
  bool pole_remap_used{false};
  bool singular_boundary_face_conductance_used{false};
  std::size_t origin_remap_unique_pair_count{0};
  std::size_t origin_remap_row_coupling_count{0};
  std::size_t pole_remap_unique_pair_count{0};
  std::size_t pole_remap_row_coupling_count{0};
  std::size_t origin_remap_nonzero_offdiag_count{0};
  std::size_t pole_remap_nonzero_offdiag_count{0};
  double origin_remap_min_abs_offdiag{0.0};
  double pole_remap_min_abs_offdiag{0.0};
  bool marshak_boundary_used{false};
  std::size_t marshak_outer_face_count{0};
  std::size_t marshak_nonzero_diagonal_loss_count{0};
  double marshak_min_G_face{0.0};
  double marshak_max_G_face{0.0};
  double marshak_min_center_to_boundary_distance_cm{0.0};
  double marshak_max_center_to_boundary_distance_cm{0.0};
  double marshak_min_abs_diagonal_loss{0.0};
  double marshak_max_abs_diagonal_loss{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

struct GenericDiffusionReferenceSolveOptions {
  std::size_t row_limit{64};
  double residual_tolerance{1.0e-10};
};

struct GenericDiffusionSolveResult {
  bool success{false};
  std::vector<double> scalar_new;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string backend;
  std::size_t row_count{0};
  double residual_l2{0.0};
  double residual_linf{0.0};
  double rhs_linf{0.0};
  double residual_linf_relative{0.0};
  double max_abs_delta{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] GenericDiffusionAssemblyResult AssembleGenericImplicitDiffusionSystem(
    const GenericDiffusionProblem& problem) noexcept;

[[nodiscard]] GenericDiffusionSolveResult SolveGenericDiffusionReference(
    const GenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionReferenceSolveOptions& options) noexcept;

[[nodiscard]] double GetCsrValue(
    const SparseMatrixCsr& matrix,
    std::size_t row,
    std::size_t column) noexcept;

[[nodiscard]] bool ValidateGenericDiffusionAssemblyDiagnostics(
    const GenericDiffusionAssemblyResult& result) noexcept;

[[nodiscard]] bool ValidateGenericDiffusionSolveDiagnostics(
    const GenericDiffusionSolveResult& result) noexcept;

}  // namespace dec3d::transport
