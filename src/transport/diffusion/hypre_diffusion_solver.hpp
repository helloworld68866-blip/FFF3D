#pragma once

#include "transport/diffusion/generic_diffusion.hpp"

#include <mpi.h>

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::transport {

struct GenericDiffusionHypreSolveOptions {
  double relative_tolerance{1.0e-10};
  int max_iterations{200};
  int krylov_dimension{30};
  bool logging_enabled{true};
  int print_level{0};
  bool require_boomeramg_preconditioner{true};
  MPI_Comm communicator{MPI_COMM_SELF};
};

struct GenericDiffusionHypreSolveResult {
  bool success{false};
  std::vector<double> scalar_new;
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string backend;
  bool hypre_enabled{false};
  bool ij_matrix_created{false};
  bool ij_vector_created{false};
  bool parcsr_matrix_created{false};
  std::size_t row_count{0};
  std::size_t nonzero_count{0};
  std::size_t global_row_count{0};
  std::size_t local_row_begin{0};
  std::size_t local_row_end{0};
  std::size_t local_row_count{0};
  std::size_t matrix_nnz{0};
  int mpi_rank_count{1};
  int gmres_iterations{0};
  double gmres_final_relative_residual{0.0};
  double residual_l2{0.0};
  double residual_linf{0.0};
  double rhs_linf{0.0};
  double residual_linf_relative{0.0};
  double max_abs_delta{0.0};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] GenericDiffusionHypreSolveResult SolveGenericDiffusionHypre(
    const GenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionHypreSolveOptions& options) noexcept;

[[nodiscard]] bool ValidateGenericDiffusionHypreSolveDiagnostics(
    const GenericDiffusionHypreSolveResult& result) noexcept;

}  // namespace dec3d::transport
