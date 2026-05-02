#include "transport/diffusion/hypre_diffusion_solver.hpp"

#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_utilities.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace dec3d::transport {

namespace {

constexpr const char* kHypreBackend = "hypre_parcsr_gmres_boomeramg";
constexpr int kBoomerAmgPreconditionerMaxIterations = 1;
constexpr double kBoomerAmgPreconditionerTolerance = 0.0;
constexpr int kBoomerAmgPreconditionerPrintLevel = 0;
constexpr double kBoomerAmgStrongThreshold = 0.5;

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] bool OptionsAreValid(const GenericDiffusionHypreSolveOptions& options) noexcept {
  return std::isfinite(options.relative_tolerance) &&
         options.relative_tolerance > 0.0 &&
         options.max_iterations > 0 &&
         options.krylov_dimension > 0 &&
         options.print_level >= 0 &&
         options.require_boomeramg_preconditioner;
}

[[nodiscard]] bool FitsHypreInt(std::size_t value) noexcept {
  return value <= static_cast<std::size_t>(std::numeric_limits<HYPRE_Int>::max());
}

[[nodiscard]] GenericDiffusionHypreSolveResult HypreFailure(
    const std::string& reason,
    const char* phase = "none",
    int error_code = 0) {
  GenericDiffusionHypreSolveResult result;
  result.backend = kHypreBackend;
  result.hypre_enabled = true;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p2.diffusion.hypre_failure"
      << "; backend=" << kHypreBackend
      << "; hypre_enabled=true"
      << "; failure_reason=" << reason
      << "; hypre_phase=" << phase
      << "; hypre_error_code=" << error_code
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  result.failure_diagnostics = out.str();
  return result;
}

struct ResidualSummary {
  double l2{0.0};
  double linf{0.0};
  double rhs_linf{0.0};
  double residual_linf_relative{0.0};
  double max_abs_delta{0.0};
};

[[nodiscard]] ResidualSummary ComputeResidualSummary(
    const GenericDiffusionAssemblyResult& assembly,
    const std::vector<double>& solution) noexcept {
  ResidualSummary summary;
  double l2_accum = 0.0;
  for (std::size_t row = 0; row < assembly.matrix.row_count; ++row) {
    double ax = 0.0;
    for (std::size_t entry = assembly.matrix.row_offsets[row];
         entry < assembly.matrix.row_offsets[row + 1u];
         ++entry) {
      ax += assembly.matrix.values[entry] *
            solution[assembly.matrix.column_indices[entry]];
    }
    const double residual = ax - assembly.rhs[row];
    l2_accum += residual * residual;
    summary.linf = std::max(summary.linf, std::abs(residual));
    summary.rhs_linf = std::max(summary.rhs_linf, std::abs(assembly.rhs[row]));
    if (row < assembly.scalar_old_flat.size()) {
      summary.max_abs_delta =
          std::max(summary.max_abs_delta, std::abs(solution[row] - assembly.scalar_old_flat[row]));
    }
  }
  summary.l2 = std::sqrt(l2_accum);
  const double residual_scale = std::max(1.0, summary.rhs_linf);
  summary.residual_linf_relative = summary.linf / residual_scale;
  return summary;
}

[[nodiscard]] bool VectorIsFinite(const std::vector<double>& values) noexcept {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value);
  });
}

[[nodiscard]] int EnsureHypreInitialized() noexcept {
  const HYPRE_Int initialized = HYPRE_Initialized();
  if (initialized != 0) {
    return 0;
  }
  return HYPRE_Initialize();
}

struct HypreObjects {
  HYPRE_IJMatrix ij_matrix{nullptr};
  HYPRE_IJVector ij_rhs{nullptr};
  HYPRE_IJVector ij_solution{nullptr};
  HYPRE_Solver gmres{nullptr};
  HYPRE_Solver boomeramg{nullptr};

  ~HypreObjects() {
    if (gmres != nullptr) {
      HYPRE_ParCSRGMRESDestroy(gmres);
    }
    if (boomeramg != nullptr) {
      HYPRE_BoomerAMGDestroy(boomeramg);
    }
    if (ij_solution != nullptr) {
      HYPRE_IJVectorDestroy(ij_solution);
    }
    if (ij_rhs != nullptr) {
      HYPRE_IJVectorDestroy(ij_rhs);
    }
    if (ij_matrix != nullptr) {
      HYPRE_IJMatrixDestroy(ij_matrix);
    }
  }
};

[[nodiscard]] std::string BuildSuccessReport(
    const GenericDiffusionHypreSolveResult& result,
    const GenericDiffusionHypreSolveOptions& options) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p2.diffusion.hypre_solve"
      << "; implementation_id=p2.diffusion.hypre_parcsr_gmres_boomeramg_v1"
      << "; backend=" << kHypreBackend
      << "; hypre_enabled=true"
      << "; ij_matrix_created=true"
      << "; ij_vector_created=true"
      << "; parcsr_matrix_created=true"
      << "; unit_system=cgs"
      << "; unknown=generic_scalar"
      << "; matrix_source=p2.diffusion.generic_implicit_assembly_v1"
      << "; matrix_format=hypre_ij_parcsr"
      << "; direction_splitting=false"
      << "; row_ownership=whole_matrix_single_rank"
      << "; mpi_rank_count=" << result.mpi_rank_count
      << "; row_count=" << result.row_count
      << "; nonzero_count=" << result.nonzero_count
      << "; global_row_count=" << result.global_row_count
      << "; local_row_begin=" << result.local_row_begin
      << "; local_row_end=" << result.local_row_end
      << "; local_row_count=" << result.local_row_count
      << "; matrix_nnz=" << result.matrix_nnz
      << "; solver=ParCSRGMRES"
      << "; preconditioner=BoomerAMG"
      << "; boomeramg_parameters=preconditioner_max_iter_1_tol_0"
      << "; boomeramg_max_iterations=" << kBoomerAmgPreconditionerMaxIterations
      << "; boomeramg_tolerance=" << kBoomerAmgPreconditionerTolerance
      << "; boomeramg_print_level=" << kBoomerAmgPreconditionerPrintLevel
      << "; boomeramg_strong_threshold=" << kBoomerAmgStrongThreshold
      << "; gmres_relative_tolerance=" << options.relative_tolerance
      << "; gmres_max_iterations=" << options.max_iterations
      << "; gmres_krylov_dimension=" << options.krylov_dimension
      << "; gmres_iterations=" << result.gmres_iterations
      << "; gmres_final_relative_residual=" << result.gmres_final_relative_residual
      << "; final_residual_norm=" << result.gmres_final_relative_residual
      << "; residual_l2=" << result.residual_l2
      << "; residual_linf=" << result.residual_linf
      << "; rhs_linf=" << result.rhs_linf
      << "; residual_linf_relative=" << result.residual_linf_relative
      << "; max_abs_delta=" << result.max_abs_delta
      << "; canonical_state_mutated=false"
      << "; fallback_used=false";
  return out.str();
}

}  // namespace

bool GenericDiffusionHypreSolveResult::is_complete() const noexcept {
  return success && backend == kHypreBackend && hypre_enabled &&
         ij_matrix_created && ij_vector_created && parcsr_matrix_created &&
         !scalar_new.empty() && !report_line.empty() &&
         ValidateGenericDiffusionHypreSolveDiagnostics(*this);
}

GenericDiffusionHypreSolveResult SolveGenericDiffusionHypre(
    const GenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionHypreSolveOptions& options) noexcept {
  if (!assembly.success || !assembly.is_complete()) {
    return HypreFailure("assembly is incomplete or failed");
  }
  if (!assembly.matrix.is_shape_complete()) {
    return HypreFailure("assembly matrix shape is incomplete");
  }
  if (assembly.rhs.size() != assembly.row_count) {
    return HypreFailure("rhs size does not match row_count");
  }
  if (assembly.scalar_old_flat.size() != assembly.row_count) {
    return HypreFailure("scalar_old_flat size does not match row_count");
  }
  if (!FitsHypreInt(assembly.row_count)) {
    return HypreFailure("row_count exceeds HYPRE_Int range");
  }
  if (!OptionsAreValid(options)) {
    return HypreFailure("HYPRE solver options are invalid");
  }

  int communicator_size = 1;
  int status = MPI_Comm_size(options.communicator, &communicator_size);
  if (status != MPI_SUCCESS) {
    return HypreFailure("MPI_Comm_size failed for HYPRE communicator", "MPI_Comm_size", status);
  }
  if (communicator_size != 1) {
    return HypreFailure("distributed HYPRE row ownership is out of scope for P2-2b");
  }

  status = EnsureHypreInitialized();
  if (status != 0) {
    return HypreFailure("HYPRE library initialization failed", "HYPRE_Initialize", status);
  }

  const auto row_count = static_cast<HYPRE_BigInt>(assembly.row_count);
  const auto row_count_int = static_cast<HYPRE_Int>(assembly.row_count);
  HypreObjects objects;
  HYPRE_ParCSRMatrix parcsr_matrix = nullptr;
  HYPRE_ParVector par_rhs = nullptr;
  HYPRE_ParVector par_solution = nullptr;

  GenericDiffusionHypreSolveResult result;
  result.backend = kHypreBackend;
  result.hypre_enabled = true;
  result.row_count = assembly.row_count;
  result.nonzero_count = assembly.nonzero_count;
  result.global_row_count = assembly.row_count;
  result.local_row_begin = 0;
  result.local_row_end = assembly.row_count;
  result.local_row_count = assembly.row_count;
  result.matrix_nnz = assembly.nonzero_count;
  result.mpi_rank_count = communicator_size;

  status = HYPRE_IJMatrixCreate(
      options.communicator,
      0,
      row_count - 1,
      0,
      row_count - 1,
      &objects.ij_matrix);
  if (status != 0) {
    return HypreFailure("HYPRE matrix creation failed", "IJMatrixCreate", status);
  }
  status = HYPRE_IJMatrixSetObjectType(objects.ij_matrix, HYPRE_PARCSR);
  if (status != 0) {
    return HypreFailure("HYPRE matrix object type failed", "IJMatrixSetObjectType", status);
  }
  status = HYPRE_IJMatrixInitialize(objects.ij_matrix);
  if (status != 0) {
    return HypreFailure("HYPRE matrix initialize failed", "IJMatrixInitialize", status);
  }

  std::vector<HYPRE_BigInt> rows(assembly.row_count);
  std::vector<HYPRE_Int> row_nonzero_counts(assembly.row_count);
  std::vector<HYPRE_BigInt> matrix_columns;
  std::vector<HYPRE_Complex> matrix_values;
  matrix_columns.reserve(assembly.matrix.column_indices.size());
  matrix_values.reserve(assembly.matrix.values.size());
  for (std::size_t row = 0; row < assembly.matrix.row_count; ++row) {
    const auto begin = assembly.matrix.row_offsets[row];
    const auto end = assembly.matrix.row_offsets[row + 1u];
    if (!FitsHypreInt(end - begin)) {
      return HypreFailure("row nonzero count exceeds HYPRE_Int range");
    }
    rows[row] = static_cast<HYPRE_BigInt>(row);
    row_nonzero_counts[row] = static_cast<HYPRE_Int>(end - begin);
    for (std::size_t entry = begin; entry < end; ++entry) {
      matrix_columns.push_back(
          static_cast<HYPRE_BigInt>(assembly.matrix.column_indices[entry]));
      matrix_values.push_back(static_cast<HYPRE_Complex>(assembly.matrix.values[entry]));
    }
  }
  status = HYPRE_IJMatrixSetValues(
      objects.ij_matrix,
      row_count_int,
      row_nonzero_counts.data(),
      rows.data(),
      matrix_columns.data(),
      matrix_values.data());
  if (status != 0) {
    return HypreFailure("HYPRE matrix row insertion failed", "IJMatrixSetValues", status);
  }

  status = HYPRE_IJMatrixAssemble(objects.ij_matrix);
  if (status != 0) {
    return HypreFailure("HYPRE matrix assemble failed", "IJMatrixAssemble", status);
  }
  status = HYPRE_IJMatrixGetObject(objects.ij_matrix, reinterpret_cast<void**>(&parcsr_matrix));
  if (status != 0) {
    return HypreFailure("HYPRE matrix object extraction failed", "IJMatrixGetObject", status);
  }

  std::vector<HYPRE_Complex> rhs_values(assembly.row_count);
  std::vector<HYPRE_Complex> initial_values(assembly.row_count);
  for (std::size_t row = 0; row < assembly.row_count; ++row) {
    rhs_values[row] = static_cast<HYPRE_Complex>(assembly.rhs[row]);
    initial_values[row] = static_cast<HYPRE_Complex>(assembly.scalar_old_flat[row]);
  }

  status = HYPRE_IJVectorCreate(options.communicator, 0, row_count - 1, &objects.ij_rhs);
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector creation failed", "IJVectorCreate", status);
  }
  status = HYPRE_IJVectorSetObjectType(objects.ij_rhs, HYPRE_PARCSR);
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector object type failed", "IJVectorSetObjectType", status);
  }
  status = HYPRE_IJVectorInitialize(objects.ij_rhs);
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector initialize failed", "IJVectorInitialize", status);
  }
  status = HYPRE_IJVectorSetValues(objects.ij_rhs, row_count_int, rows.data(), rhs_values.data());
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector values failed", "IJVectorSetValues", status);
  }
  status = HYPRE_IJVectorAssemble(objects.ij_rhs);
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector assemble failed", "IJVectorAssemble", status);
  }
  status = HYPRE_IJVectorGetObject(objects.ij_rhs, reinterpret_cast<void**>(&par_rhs));
  if (status != 0) {
    return HypreFailure("HYPRE RHS vector object extraction failed", "IJVectorGetObject", status);
  }

  status = HYPRE_IJVectorCreate(options.communicator, 0, row_count - 1, &objects.ij_solution);
  if (status != 0) {
    return HypreFailure("HYPRE solution vector creation failed", "IJVectorCreate", status);
  }
  status = HYPRE_IJVectorSetObjectType(objects.ij_solution, HYPRE_PARCSR);
  if (status != 0) {
    return HypreFailure("HYPRE solution vector object type failed", "IJVectorSetObjectType", status);
  }
  status = HYPRE_IJVectorInitialize(objects.ij_solution);
  if (status != 0) {
    return HypreFailure("HYPRE solution vector initialize failed", "IJVectorInitialize", status);
  }
  status = HYPRE_IJVectorSetValues(
      objects.ij_solution,
      row_count_int,
      rows.data(),
      initial_values.data());
  if (status != 0) {
    return HypreFailure("HYPRE solution vector values failed", "IJVectorSetValues", status);
  }
  status = HYPRE_IJVectorAssemble(objects.ij_solution);
  if (status != 0) {
    return HypreFailure("HYPRE solution vector assemble failed", "IJVectorAssemble", status);
  }
  status = HYPRE_IJVectorGetObject(objects.ij_solution, reinterpret_cast<void**>(&par_solution));
  if (status != 0) {
    return HypreFailure("HYPRE solution vector object extraction failed", "IJVectorGetObject", status);
  }

  status = HYPRE_ParCSRGMRESCreate(options.communicator, &objects.gmres);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES creation failed", "ParCSRGMRESCreate", status);
  }
  status = HYPRE_BoomerAMGCreate(&objects.boomeramg);
  if (status != 0) {
    return HypreFailure("HYPRE BoomerAMG creation failed", "BoomerAMGCreate", status);
  }
  status = HYPRE_BoomerAMGSetPrintLevel(
      objects.boomeramg,
      kBoomerAmgPreconditionerPrintLevel);
  if (status != 0) {
    return HypreFailure("HYPRE BoomerAMG print-level setup failed",
                        "BoomerAMGSetPrintLevel",
                        status);
  }
  status = HYPRE_BoomerAMGSetMaxIter(
      objects.boomeramg,
      kBoomerAmgPreconditionerMaxIterations);
  if (status != 0) {
    return HypreFailure("HYPRE BoomerAMG max-iteration setup failed",
                        "BoomerAMGSetMaxIter",
                        status);
  }
  status = HYPRE_BoomerAMGSetTol(
      objects.boomeramg,
      kBoomerAmgPreconditionerTolerance);
  if (status != 0) {
    return HypreFailure("HYPRE BoomerAMG tolerance setup failed",
                        "BoomerAMGSetTol",
                        status);
  }
  status = HYPRE_BoomerAMGSetStrongThreshold(
      objects.boomeramg,
      kBoomerAmgStrongThreshold);
  if (status != 0) {
    return HypreFailure("HYPRE BoomerAMG strong-threshold setup failed",
                        "BoomerAMGSetStrongThreshold",
                        status);
  }
  status = HYPRE_ParCSRGMRESSetTol(objects.gmres, options.relative_tolerance);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES tolerance setup failed", "ParCSRGMRESSetTol", status);
  }
  status = HYPRE_ParCSRGMRESSetMaxIter(objects.gmres, options.max_iterations);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES max-iteration setup failed", "ParCSRGMRESSetMaxIter", status);
  }
  status = HYPRE_ParCSRGMRESSetKDim(objects.gmres, options.krylov_dimension);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES Krylov dimension setup failed", "ParCSRGMRESSetKDim", status);
  }
  status = HYPRE_ParCSRGMRESSetLogging(objects.gmres, options.logging_enabled ? 1 : 0);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES logging setup failed", "ParCSRGMRESSetLogging", status);
  }
  status = HYPRE_ParCSRGMRESSetPrintLevel(objects.gmres, options.print_level);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES print-level setup failed", "ParCSRGMRESSetPrintLevel", status);
  }
  status = HYPRE_ParCSRGMRESSetPrecond(
      objects.gmres,
      reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSolve),
      reinterpret_cast<HYPRE_PtrToParSolverFcn>(HYPRE_BoomerAMGSetup),
      objects.boomeramg);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES preconditioner setup failed", "ParCSRGMRESSetPrecond", status);
  }
  status = HYPRE_ParCSRGMRESSetup(objects.gmres, parcsr_matrix, par_rhs, par_solution);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES setup failed", "ParCSRGMRESSetup", status);
  }
  status = HYPRE_ParCSRGMRESSolve(objects.gmres, parcsr_matrix, par_rhs, par_solution);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES solve failed", "ParCSRGMRESSolve", status);
  }

  status = HYPRE_ParCSRGMRESGetNumIterations(objects.gmres, &result.gmres_iterations);
  if (status != 0) {
    return HypreFailure("HYPRE GMRES iteration query failed", "ParCSRGMRESGetNumIterations", status);
  }
  status = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(
      objects.gmres,
      &result.gmres_final_relative_residual);
  if (status != 0) {
    return HypreFailure(
        "HYPRE GMRES final residual query failed",
        "ParCSRGMRESGetFinalRelativeResidualNorm",
        status);
  }

  std::vector<HYPRE_Complex> solution_values(assembly.row_count);
  status = HYPRE_IJVectorGetValues(
      objects.ij_solution,
      row_count_int,
      rows.data(),
      solution_values.data());
  if (status != 0) {
    return HypreFailure("HYPRE solution extraction failed", "IJVectorGetValues", status);
  }
  result.scalar_new.resize(assembly.row_count);
  for (std::size_t row = 0; row < assembly.row_count; ++row) {
    result.scalar_new[row] = static_cast<double>(solution_values[row]);
  }

  const auto residual = ComputeResidualSummary(assembly, result.scalar_new);
  result.residual_l2 = residual.l2;
  result.residual_linf = residual.linf;
  result.rhs_linf = residual.rhs_linf;
  result.residual_linf_relative = residual.residual_linf_relative;
  result.max_abs_delta = residual.max_abs_delta;
  if (!VectorIsFinite(result.scalar_new)) {
    return HypreFailure("HYPRE solution contains non-finite values");
  }
  if (result.residual_linf_relative > options.relative_tolerance * 100.0) {
    return HypreFailure("HYPRE recomputed residual exceeds tolerance");
  }

  result.success = true;
  result.ij_matrix_created = true;
  result.ij_vector_created = true;
  result.parcsr_matrix_created = true;
  result.report_line = BuildSuccessReport(result, options);
  return result;
}

bool ValidateGenericDiffusionHypreSolveDiagnostics(
    const GenericDiffusionHypreSolveResult& result) noexcept {
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p2.diffusion.hypre_solve") &&
         Contains(line, "implementation_id=p2.diffusion.hypre_parcsr_gmres_boomeramg_v1") &&
         Contains(line, "backend=hypre_parcsr_gmres_boomeramg") &&
         Contains(line, "hypre_enabled=true") &&
         Contains(line, "ij_matrix_created=true") &&
         Contains(line, "ij_vector_created=true") &&
         Contains(line, "parcsr_matrix_created=true") &&
         Contains(line, "unit_system=cgs") &&
         Contains(line, "unknown=generic_scalar") &&
         Contains(line, "matrix_source=p2.diffusion.generic_implicit_assembly_v1") &&
         Contains(line, "matrix_format=hypre_ij_parcsr") &&
         Contains(line, "direction_splitting=false") &&
         Contains(line, "row_ownership=whole_matrix_single_rank") &&
         Contains(line, "mpi_rank_count=1") &&
         Contains(line, "global_row_count=") &&
         Contains(line, "local_row_begin=0") &&
         Contains(line, "local_row_end=") &&
         Contains(line, "local_row_count=") &&
         Contains(line, "matrix_nnz=") &&
         Contains(line, "solver=ParCSRGMRES") &&
         Contains(line, "preconditioner=BoomerAMG") &&
         Contains(line, "boomeramg_parameters=preconditioner_max_iter_1_tol_0") &&
         Contains(line, "boomeramg_max_iterations=1") &&
         Contains(line, "boomeramg_tolerance=0") &&
         Contains(line, "boomeramg_print_level=0") &&
         Contains(line, "boomeramg_strong_threshold=0.5") &&
         Contains(line, "gmres_relative_tolerance=") &&
         Contains(line, "gmres_max_iterations=") &&
         Contains(line, "gmres_krylov_dimension=") &&
         Contains(line, "gmres_iterations=") &&
         Contains(line, "gmres_final_relative_residual=") &&
         Contains(line, "final_residual_norm=") &&
         Contains(line, "residual_l2=") &&
         Contains(line, "residual_linf=") &&
         Contains(line, "rhs_linf=") &&
         Contains(line, "residual_linf_relative=") &&
         Contains(line, "max_abs_delta=") &&
         Contains(line, "canonical_state_mutated=false") &&
         Contains(line, "fallback_used=false");
}

}  // namespace dec3d::transport
