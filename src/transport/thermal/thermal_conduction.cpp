#include "transport/thermal/thermal_conduction.hpp"

#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/thermal/thermal_conductivity.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace dec3d::transport {

namespace {

constexpr double kConservationTolerance = 1.0e-18;

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] std::size_t FlatIndex(
    const DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return (radial * layout.theta_cells + theta) * layout.phi_cells + phi;
}

[[nodiscard]] DiffusionGridLayout LayoutFromState(
    const dec3d::state::CanonicalState& state) noexcept {
  return DiffusionGridLayout{
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells};
}

struct ThermalCoefficientArrays {
  bool success{false};
  std::string failure_reason;
  std::string failure_diagnostics;
  dec3d::core::Array3D<double> kappa_e;
  dec3d::core::Array3D<double> kappa_i;
  std::size_t provider_failure_count{0};
  std::string first_provider_failure_cell;
  double min_kappa_e{std::numeric_limits<double>::infinity()};
  double max_kappa_e{0.0};
  double min_kappa_i{std::numeric_limits<double>::infinity()};
  double max_kappa_i{0.0};
  double min_lnLambda{std::numeric_limits<double>::infinity()};
  double max_lnLambda{0.0};
  double min_f_LM{std::numeric_limits<double>::infinity()};
  double max_f_LM{0.0};
  std::string provider_model_requested;
  std::string provider_model_executed;
};

[[nodiscard]] ThermalConductivityModel ProviderModelFromKappaModel(
    ThermalConductionKappaModel model) noexcept {
  switch (model) {
    case ThermalConductionKappaModel::spitzer_no_degeneracy:
      return ThermalConductivityModel::spitzer_no_degeneracy;
    case ThermalConductionKappaModel::lee_more_with_degeneracy:
      return ThermalConductivityModel::lee_more_with_degeneracy;
    case ThermalConductionKappaModel::constant_user_supplied:
    case ThermalConductionKappaModel::unsupported:
      return ThermalConductivityModel::unsupported;
  }
  return ThermalConductivityModel::unsupported;
}

[[nodiscard]] std::string CellLabel(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  std::ostringstream label;
  label << radial << "," << theta << "," << phi;
  return label.str();
}

void Fail(
    ThermalConductionResult& result,
    const std::string& reason,
    const std::string& nested = {}) {
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p2.thermal_conduction.failure"
      << "; failure_reason=" << reason
      << "; canonical_state_mutated=false";
  if (!nested.empty()) {
    out << "; nested_diagnostics={" << nested << "}";
  }
  result.failure_diagnostics = out.str();
}

[[nodiscard]] bool IsBoundaryPolicyComplete(
    const GenericDiffusionBoundaryPolicy& policy) noexcept {
  return policy.inner_radial != DiffusionBoundaryKind::missing &&
         policy.outer_radial != DiffusionBoundaryKind::missing &&
         policy.theta_lower != DiffusionBoundaryKind::missing &&
         policy.theta_upper != DiffusionBoundaryKind::missing &&
         policy.phi != DiffusionBoundaryKind::missing;
}

[[nodiscard]] bool ValidateOptions(
    const ThermalConductionOptions& options,
    ThermalConductionResult& result) {
  if (!Finite(options.dt_s)) {
    Fail(result, "dt_s must be finite");
    return false;
  }
  if (options.dt_s < 0.0) {
    Fail(result, "dt_s must be non-negative");
    return false;
  }
  if (options.kappa_model == ThermalConductionKappaModel::unsupported) {
    Fail(result, "unsupported kappa_model");
    return false;
  }
  if (options.kappa_model == ThermalConductionKappaModel::constant_user_supplied) {
    if (!Finite(options.kappa_e_cm_inv_s)) {
      Fail(result, "kappa_e_cm_inv_s must be finite");
      return false;
    }
    if (options.kappa_e_cm_inv_s < 0.0) {
      Fail(result, "kappa_e_cm_inv_s must be non-negative");
      return false;
    }
    if (!Finite(options.kappa_i_cm_inv_s)) {
      Fail(result, "kappa_i_cm_inv_s must be finite");
      return false;
    }
    if (options.kappa_i_cm_inv_s < 0.0) {
      Fail(result, "kappa_i_cm_inv_s must be non-negative");
      return false;
    }
  }
  if (options.backend == ThermalConductionBackend::unsupported) {
    Fail(result, "unsupported thermal conduction backend");
    return false;
  }
  if (!IsBoundaryPolicyComplete(options.boundary_policy)) {
    Fail(result, "boundary policy is missing");
    return false;
  }
  return true;
}

[[nodiscard]] bool ValidateGeometryForLayout(
    const DiffusionGridLayout& layout,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    ThermalConductionResult& result) {
  if (!geometry.valid) {
    Fail(result, "thermal conduction geometry must be valid");
    return false;
  }
  if (geometry.cell_volumes.size() != layout.cell_count()) {
    Fail(result, "thermal conduction geometry cell volume count mismatch");
    return false;
  }
  for (const double volume : geometry.cell_volumes) {
    if (!Finite(volume) || volume <= 0.0) {
      Fail(result, "thermal conduction geometry cell volume must be positive");
      return false;
    }
  }
  return true;
}

[[nodiscard]] ThermalConductionSolveResult SerialSolve(
    const GenericDiffusionAssemblyResult& assembly,
    const GenericDiffusionReferenceSolveOptions& options) {
  ThermalConductionSolveResult result;
  const auto solve = SolveGenericDiffusionReference(assembly, options);
  result.success = solve.success;
  result.scalar_new = solve.scalar_new;
  result.report_line = solve.report_line;
  result.failure_diagnostics = solve.failure_diagnostics;
  result.failure_reason = solve.failure_reason;
  result.backend_executed = solve.backend;
  result.residual_l2 = solve.residual_l2;
  result.residual_linf = solve.residual_linf;
  result.max_abs_delta = solve.max_abs_delta;
  return result;
}

[[nodiscard]] ThermalCoefficientArrays BuildThermalCoefficientArrays(
    const dec3d::state::ThermodynamicRecoveryResult& recovered,
    const DiffusionGridLayout& layout,
    const ThermalConductionOptions& options) {
  ThermalCoefficientArrays out;
  out.kappa_e = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  out.kappa_i = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  out.provider_model_requested = "none";
  out.provider_model_executed = "none";

  if (options.kappa_model == ThermalConductionKappaModel::constant_user_supplied) {
    for (std::size_t r = 0; r < layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < layout.theta_cells; ++t) {
        for (std::size_t p = 0; p < layout.phi_cells; ++p) {
          out.kappa_e(r, t, p) = options.kappa_e_cm_inv_s;
          out.kappa_i(r, t, p) = options.kappa_i_cm_inv_s;
        }
      }
    }
    out.min_kappa_e = options.kappa_e_cm_inv_s;
    out.max_kappa_e = options.kappa_e_cm_inv_s;
    out.min_kappa_i = options.kappa_i_cm_inv_s;
    out.max_kappa_i = options.kappa_i_cm_inv_s;
    out.min_lnLambda = 0.0;
    out.max_lnLambda = 0.0;
    out.min_f_LM = 1.0;
    out.max_f_LM = 1.0;
    out.success = true;
    return out;
  }

  const auto provider_model = ProviderModelFromKappaModel(options.kappa_model);
  out.provider_model_requested = ThermalConductivityModelName(provider_model);
  if (provider_model == ThermalConductivityModel::unsupported) {
    out.failure_reason = "unsupported variable conductivity model";
    out.provider_failure_count = 1;
    out.failure_diagnostics =
        "diagnostic_id=p2.thermal_conduction.variable_kappa.failure"
        "; failure_reason=unsupported variable conductivity model"
        "; provider_failure_count=1"
        "; canonical_state_mutated=false";
    return out;
  }

  std::string nested_failure;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const auto& cell = recovered.cells(r, t, p);
        ThermalConductivityInput input;
        input.Te_erg = cell.t_e_erg_per_particle;
        input.Ti_erg = cell.t_i_erg_per_particle;
        input.ne_cm3 = cell.n_e_cm3;
        input.ni_cm3 = cell.n_i_cm3;
        input.zbar = cell.zbar;
        input.mean_ion_mass_g = cell.mean_ion_mass_g;
        input.model = provider_model;

        const auto coefficient = ComputeThermalConductivity(input);
        if (!coefficient.success) {
          ++out.provider_failure_count;
          if (out.first_provider_failure_cell.empty()) {
            out.first_provider_failure_cell = CellLabel(r, t, p);
            nested_failure = coefficient.failure_diagnostics;
          }
          continue;
        }
        if (coefficient.model_executed != out.provider_model_requested) {
          ++out.provider_failure_count;
          if (out.first_provider_failure_cell.empty()) {
            out.first_provider_failure_cell = CellLabel(r, t, p);
            nested_failure = "provider model mismatch";
          }
          continue;
        }

        out.provider_model_executed = coefficient.model_executed;
        out.kappa_e(r, t, p) = coefficient.kappa_e_cm_inv_s;
        out.kappa_i(r, t, p) = coefficient.kappa_i_cm_inv_s;
        out.min_kappa_e = std::min(out.min_kappa_e, coefficient.kappa_e_cm_inv_s);
        out.max_kappa_e = std::max(out.max_kappa_e, coefficient.kappa_e_cm_inv_s);
        out.min_kappa_i = std::min(out.min_kappa_i, coefficient.kappa_i_cm_inv_s);
        out.max_kappa_i = std::max(out.max_kappa_i, coefficient.kappa_i_cm_inv_s);
        out.min_lnLambda = std::min(out.min_lnLambda, coefficient.lnLambda);
        out.max_lnLambda = std::max(out.max_lnLambda, coefficient.lnLambda);
        out.min_f_LM = std::min(out.min_f_LM, coefficient.f_LM);
        out.max_f_LM = std::max(out.max_f_LM, coefficient.f_LM);
      }
    }
  }

  if (out.provider_failure_count > 0) {
    out.failure_reason = "thermal conductivity provider failed";
    std::ostringstream diag;
    diag << "diagnostic_id=p2.thermal_conduction.variable_kappa.failure"
         << "; failure_reason=thermal conductivity provider failed"
         << "; provider_failure_count=" << out.provider_failure_count
         << "; first_provider_failure_cell=" << out.first_provider_failure_cell;
    if (!nested_failure.empty()) {
      diag << "; nested_diagnostics={" << nested_failure << "}";
    }
    diag << "; canonical_state_mutated=false";
    out.failure_diagnostics = diag.str();
    return out;
  }
  if (out.provider_model_executed.empty() || out.provider_model_executed == "none") {
    out.failure_reason = "thermal conductivity provider did not execute";
    out.provider_failure_count = 1;
    out.failure_diagnostics =
        "diagnostic_id=p2.thermal_conduction.variable_kappa.failure"
        "; failure_reason=thermal conductivity provider did not execute"
        "; provider_failure_count=1"
        "; canonical_state_mutated=false";
    return out;
  }

  out.success = true;
  return out;
}

[[nodiscard]] GenericDiffusionProblem BuildThermalDiffusionProblem(
    const dec3d::state::ThermodynamicRecoveryResult& recovered,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const DiffusionGridLayout& layout,
    const GenericDiffusionBoundaryPolicy& boundary_policy,
    double dt_s,
    const dec3d::core::Array3D<double>& kappa,
    bool electron) {
  GenericDiffusionProblem problem;
  problem.layout = layout;
  problem.geometry = geometry;
  problem.dt_s = dt_s;
  problem.coefficient_A = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.coefficient_D = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.coefficient_C = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.coefficient_B = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.scalar_old = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.boundary_policy = boundary_policy;

  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const auto& cell = recovered.cells(r, t, p);
        problem.coefficient_A(r, t, p) =
            (electron ? cell.n_e_cm3 : cell.n_i_cm3) / gamma_minus_one;
        problem.coefficient_D(r, t, p) = kappa(r, t, p);
        problem.scalar_old(r, t, p) =
            electron ? cell.t_e_erg_per_particle : cell.t_i_erg_per_particle;
      }
    }
  }

  return problem;
}

[[nodiscard]] std::string BuildConstantSuccessReport(
    const ThermalConductionOptions& options,
    const ThermalConductionResult& result,
    bool metadata_written) {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p2.thermal_conduction.constant_kappa"
         << "; implementation_id=p2.thermal_conduction.constant_kappa_v1"
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; kappa_internal_unit=cm^-1_s^-1"
         << "; kappa_model=" << ThermalConductionKappaModelName(options.kappa_model)
         << "; kappa_e_cm_inv_s=" << options.kappa_e_cm_inv_s
         << "; kappa_i_cm_inv_s=" << options.kappa_i_cm_inv_s
         << "; dt_s=" << options.dt_s
         << "; backend_requested=" << result.backend_requested
         << "; backend_executed=" << result.backend_executed
         << "; backend=" << result.backend_executed
         << "; electron_equation_executed=true"
         << "; ion_equation_executed=true"
         << "; electron_assembly_diagnostic_id=p2.diffusion.assembly"
         << "; ion_assembly_diagnostic_id=p2.diffusion.assembly"
         << "; electron_solve_backend=" << result.backend_executed
         << "; ion_solve_backend=" << result.backend_executed
         << "; updated_fields=e_electron,e_fluid_total"
         << "; rho_momentum_unchanged=true"
         << "; canonical_state_mutated=" << (metadata_written ? "true" : "false")
         << "; metadata_written=" << (metadata_written ? "true" : "false")
         << "; fallback_used=false"
         << "; cell_count=" << result.cell_count
         << "; max_abs_delta_Te_erg=" << result.max_abs_delta_Te_erg
         << "; max_abs_delta_Ti_erg=" << result.max_abs_delta_Ti_erg
         << "; max_electron_energy_residual=" << result.max_electron_energy_residual
         << "; max_ion_energy_residual=" << result.max_ion_energy_residual
         << "; global_thermal_energy_residual=" << result.global_thermal_energy_residual
         << "; min_e_electron_after=" << result.min_e_electron_after
         << "; min_e_ion_after=" << result.min_e_ion_after;
  return report.str();
}

[[nodiscard]] std::string BuildVariableSuccessReport(
    const ThermalConductionOptions& options,
    const ThermalConductionResult& result,
    bool metadata_written) {
  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p2.thermal_conduction.variable_kappa"
         << "; implementation_id=p2.thermal_conduction.variable_kappa_v1"
         << "; unit_system=cgs"
         << "; temperature_internal_unit=erg_per_particle"
         << "; temperature_output_unit=keV"
         << "; kappa_internal_unit=cm^-1_s^-1"
         << "; kappa_model_requested=" << result.kappa_model_requested
         << "; kappa_model_executed=" << result.kappa_model_executed
         << "; provider_model_requested=" << result.provider_model_requested
         << "; provider_model_executed=" << result.provider_model_executed
         << "; coefficient_time_level=old_time_lagged"
         << "; provider_failure_count=" << result.provider_failure_count
         << "; first_provider_failure_cell=" << result.first_provider_failure_cell
         << "; min_kappa_e_cm_inv_s=" << result.min_kappa_e_cm_inv_s
         << "; max_kappa_e_cm_inv_s=" << result.max_kappa_e_cm_inv_s
         << "; min_kappa_i_cm_inv_s=" << result.min_kappa_i_cm_inv_s
         << "; max_kappa_i_cm_inv_s=" << result.max_kappa_i_cm_inv_s
         << "; min_provider_lnLambda=" << result.min_provider_lnLambda
         << "; max_provider_lnLambda=" << result.max_provider_lnLambda
         << "; min_provider_f_LM=" << result.min_provider_f_LM
         << "; max_provider_f_LM=" << result.max_provider_f_LM
         << "; dt_s=" << options.dt_s
         << "; backend_requested=" << result.backend_requested
         << "; backend_executed=" << result.backend_executed
         << "; backend=" << result.backend_executed
         << "; electron_equation_executed=true"
         << "; ion_equation_executed=true"
         << "; electron_assembly_diagnostic_id=p2.diffusion.assembly"
         << "; ion_assembly_diagnostic_id=p2.diffusion.assembly"
         << "; electron_solve_backend=" << result.backend_executed
         << "; ion_solve_backend=" << result.backend_executed
         << "; updated_fields=e_electron,e_fluid_total"
         << "; rho_momentum_unchanged=true"
         << "; canonical_state_mutated=" << (metadata_written ? "true" : "false")
         << "; metadata_written=" << (metadata_written ? "true" : "false")
         << "; fallback_used=false"
         << "; cell_count=" << result.cell_count
         << "; max_abs_delta_Te_erg=" << result.max_abs_delta_Te_erg
         << "; max_abs_delta_Ti_erg=" << result.max_abs_delta_Ti_erg
         << "; max_electron_energy_residual=" << result.max_electron_energy_residual
         << "; max_ion_energy_residual=" << result.max_ion_energy_residual
         << "; global_thermal_energy_residual=" << result.global_thermal_energy_residual
         << "; min_e_electron_after=" << result.min_e_electron_after
         << "; min_e_ion_after=" << result.min_e_ion_after;
  return report.str();
}

}  // namespace

const char* ThermalConductionBackendName(ThermalConductionBackend backend) noexcept {
  switch (backend) {
    case ThermalConductionBackend::serial_dense_reference:
      return "serial_dense_reference";
    case ThermalConductionBackend::hypre_parcsr_gmres_boomeramg:
      return "hypre_parcsr_gmres_boomeramg";
    case ThermalConductionBackend::unsupported:
      return "unsupported";
  }
  return "unknown";
}

const char* ThermalConductionKappaModelName(ThermalConductionKappaModel model) noexcept {
  switch (model) {
    case ThermalConductionKappaModel::constant_user_supplied:
      return "constant_user_supplied";
    case ThermalConductionKappaModel::spitzer_no_degeneracy:
      return "spitzer_no_degeneracy";
    case ThermalConductionKappaModel::lee_more_with_degeneracy:
      return "lee_more_with_degeneracy";
    case ThermalConductionKappaModel::unsupported:
      return "unsupported";
  }
  return "unknown";
}

bool ThermalConductionResult::is_complete() const noexcept {
  return success &&
         updated_fields ==
             (dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
              dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total)) &&
         cell_count > 0 &&
         !report_line.empty() &&
         (report_line.find("diagnostic_id=p2.thermal_conduction.constant_kappa") != std::string::npos ||
          report_line.find("diagnostic_id=p2.thermal_conduction.variable_kappa") != std::string::npos) &&
         report_line.find("backend_requested=") != std::string::npos &&
         report_line.find("backend_executed=") != std::string::npos &&
         report_line.find("updated_fields=e_electron,e_fluid_total") != std::string::npos;
}

ThermalConductionResult ApplyThermalConduction(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ThermalConductionOptions& options,
    const ThermalConductionSolveHook& solve_hook) noexcept {
  ThermalConductionResult result;
  result.dt_s = options.dt_s;
  result.kappa_e_cm_inv_s = options.kappa_e_cm_inv_s;
  result.kappa_i_cm_inv_s = options.kappa_i_cm_inv_s;
  result.backend_requested = ThermalConductionBackendName(options.backend);
  result.kappa_model_requested = ThermalConductionKappaModelName(options.kappa_model);
  result.kappa_model_executed = ThermalConductionKappaModelName(options.kappa_model);

  if (!ValidateOptions(options, result)) {
    return result;
  }

  const auto layout = LayoutFromState(state);
  if (!layout.is_valid()) {
    Fail(result, "thermal conduction state layout is invalid");
    return result;
  }
  if (!ValidateGeometryForLayout(layout, geometry, result)) {
    return result;
  }

  const auto recovered = dec3d::state::RecoverThermodynamicState(
      state,
      dec3d::state::ThermodynamicRecoveryOptions{
          options.electron_energy_floor,
          options.ion_energy_floor,
          1.0,
          0.0});
  if (!recovered.success) {
    Fail(result, "P2-0 thermodynamic recovery failed", recovered.failure_diagnostics);
    return result;
  }

  const auto coefficients = BuildThermalCoefficientArrays(recovered, layout, options);
  result.provider_failure_count = coefficients.provider_failure_count;
  result.first_provider_failure_cell = coefficients.first_provider_failure_cell;
  result.provider_model_requested = coefficients.provider_model_requested;
  result.provider_model_executed = coefficients.provider_model_executed;
  if (!coefficients.success) {
    result.provider_failure_count = coefficients.provider_failure_count;
    result.first_provider_failure_cell = coefficients.first_provider_failure_cell;
    Fail(result, coefficients.failure_reason, coefficients.failure_diagnostics);
    return result;
  }
  result.min_kappa_e_cm_inv_s = coefficients.min_kappa_e;
  result.max_kappa_e_cm_inv_s = coefficients.max_kappa_e;
  result.min_kappa_i_cm_inv_s = coefficients.min_kappa_i;
  result.max_kappa_i_cm_inv_s = coefficients.max_kappa_i;
  result.min_provider_lnLambda = coefficients.min_lnLambda;
  result.max_provider_lnLambda = coefficients.max_lnLambda;
  result.min_provider_f_LM = coefficients.min_f_LM;
  result.max_provider_f_LM = coefficients.max_f_LM;

  auto electron_problem = BuildThermalDiffusionProblem(
      recovered,
      geometry,
      layout,
      options.boundary_policy,
      options.dt_s,
      coefficients.kappa_e,
      true);
  auto ion_problem = BuildThermalDiffusionProblem(
      recovered,
      geometry,
      layout,
      options.boundary_policy,
      options.dt_s,
      coefficients.kappa_i,
      false);

  const auto electron_assembly = AssembleGenericImplicitDiffusionSystem(electron_problem);
  result.electron_assembly_report =
      electron_assembly.success ? electron_assembly.report_line : electron_assembly.failure_diagnostics;
  if (!electron_assembly.success) {
    Fail(result, "electron diffusion assembly failed", electron_assembly.failure_diagnostics);
    return result;
  }

  const auto ion_assembly = AssembleGenericImplicitDiffusionSystem(ion_problem);
  result.ion_assembly_report =
      ion_assembly.success ? ion_assembly.report_line : ion_assembly.failure_diagnostics;
  if (!ion_assembly.success) {
    Fail(result, "ion diffusion assembly failed", ion_assembly.failure_diagnostics);
    return result;
  }

  ThermalConductionSolveResult electron_solve;
  ThermalConductionSolveResult ion_solve;
  if (options.backend == ThermalConductionBackend::serial_dense_reference) {
    electron_solve = SerialSolve(electron_assembly, options.serial_reference_options);
    ion_solve = SerialSolve(ion_assembly, options.serial_reference_options);
  } else {
    if (!solve_hook) {
      Fail(result, "requested backend is unavailable without a solve hook");
      return result;
    }
    try {
      electron_solve = solve_hook(electron_assembly);
      ion_solve = solve_hook(ion_assembly);
    } catch (...) {
      Fail(result, "thermal conduction solve hook threw an exception");
      return result;
    }
  }

  result.electron_solve_report =
      electron_solve.success ? electron_solve.report_line : electron_solve.failure_diagnostics;
  result.ion_solve_report =
      ion_solve.success ? ion_solve.report_line : ion_solve.failure_diagnostics;
  if (!electron_solve.success) {
    Fail(result, "electron diffusion solve failed", electron_solve.failure_diagnostics);
    return result;
  }
  if (!ion_solve.success) {
    Fail(result, "ion diffusion solve failed", ion_solve.failure_diagnostics);
    return result;
  }
  if (electron_solve.scalar_new.size() != layout.cell_count() ||
      ion_solve.scalar_new.size() != layout.cell_count()) {
    Fail(result, "thermal conduction solve result size mismatch");
    return result;
  }

  result.backend_executed = electron_solve.backend_executed;
  if (result.backend_executed.empty() || result.backend_executed != ion_solve.backend_executed) {
    Fail(result, "electron and ion solve backends are inconsistent");
    return result;
  }
  if (result.backend_executed != result.backend_requested) {
    Fail(result, "backend fallback or mismatch detected");
    return result;
  }

  auto candidate = state;
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  result.min_e_electron_after = std::numeric_limits<double>::infinity();
  result.min_e_ion_after = std::numeric_limits<double>::infinity();
  double thermal_before_volume = 0.0;
  double thermal_after_volume = 0.0;
  double electron_before_volume = 0.0;
  double electron_after_volume = 0.0;
  double ion_before_volume = 0.0;
  double ion_after_volume = 0.0;

  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const std::size_t flat = FlatIndex(layout, r, t, p);
        const auto& before = recovered.cells(r, t, p);
        const double te_new = electron_solve.scalar_new[flat];
        const double ti_new = ion_solve.scalar_new[flat];
        const double e_electron_new = before.n_e_cm3 * te_new / gamma_minus_one;
        const double e_ion_new = before.n_i_cm3 * ti_new / gamma_minus_one;
        const double e_total_new =
            before.kinetic_energy_density_erg_per_cm3 + e_electron_new + e_ion_new;
        if (!(Finite(te_new) && Finite(ti_new) && Finite(e_electron_new) &&
              Finite(e_ion_new) && Finite(e_total_new))) {
          Fail(result, "thermal conduction updated state is not finite");
          return result;
        }
        if (e_electron_new < options.electron_energy_floor) {
          Fail(result, "electron energy fell below minimum threshold");
          return result;
        }
        if (e_ion_new < options.ion_energy_floor) {
          Fail(result, "ion energy fell below minimum threshold");
          return result;
        }

        result.max_abs_delta_Te_erg =
            std::max(result.max_abs_delta_Te_erg, std::abs(te_new - before.t_e_erg_per_particle));
        result.max_abs_delta_Ti_erg =
            std::max(result.max_abs_delta_Ti_erg, std::abs(ti_new - before.t_i_erg_per_particle));
        result.min_e_electron_after = std::min(result.min_e_electron_after, e_electron_new);
        result.min_e_ion_after = std::min(result.min_e_ion_after, e_ion_new);

        const double volume = geometry.cell_volumes[flat];
        electron_before_volume += volume * before.e_electron_erg_per_cm3;
        electron_after_volume += volume * e_electron_new;
        ion_before_volume += volume * before.e_ion_erg_per_cm3;
        ion_after_volume += volume * e_ion_new;
        thermal_before_volume += volume * (before.e_electron_erg_per_cm3 + before.e_ion_erg_per_cm3);
        thermal_after_volume += volume * (e_electron_new + e_ion_new);

        candidate.e_electron(r, t, p) = e_electron_new;
        candidate.e_fluid_total(r, t, p) = e_total_new;
        ++result.cell_count;
      }
    }
  }

  result.max_electron_energy_residual = std::abs(electron_after_volume - electron_before_volume);
  result.max_ion_energy_residual = std::abs(ion_after_volume - ion_before_volume);
  result.global_thermal_energy_residual = std::abs(thermal_after_volume - thermal_before_volume);
  if (result.global_thermal_energy_residual > kConservationTolerance) {
    Fail(result, "global thermal energy residual exceeds tolerance");
    return result;
  }

  result.updated_fields =
      dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
      dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);

  if (options.dt_s > 0.0) {
    state.e_electron = std::move(candidate.e_electron);
    state.e_fluid_total = std::move(candidate.e_fluid_total);
    state.ApplyAuthoritativeWrite(result.updated_fields);
  }

  result.success = true;
  if (options.kappa_model == ThermalConductionKappaModel::constant_user_supplied) {
    result.report_line = BuildConstantSuccessReport(options, result, options.dt_s > 0.0);
  } else {
    result.report_line = BuildVariableSuccessReport(options, result, options.dt_s > 0.0);
  }
  return result;
}

ThermalConductionResult ApplyConstantKappaThermalConduction(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const ThermalConductionOptions& options,
    const ThermalConductionSolveHook& solve_hook) noexcept {
  return ApplyThermalConduction(state, geometry, options, solve_hook);
}

bool ValidateThermalConductionDiagnostics(const ThermalConductionResult& result) noexcept {
  if (!result.is_complete()) {
    return false;
  }
  const auto& line = result.report_line;
  const bool variable =
      line.find("diagnostic_id=p2.thermal_conduction.variable_kappa") != std::string::npos;
  const bool base =
      line.find("unit_system=cgs") != std::string::npos &&
      line.find("temperature_internal_unit=erg_per_particle") != std::string::npos &&
      line.find("temperature_output_unit=keV") != std::string::npos &&
      line.find("kappa_internal_unit=cm^-1_s^-1") != std::string::npos &&
      line.find("backend_requested=") != std::string::npos &&
      line.find("backend_executed=") != std::string::npos &&
      line.find("electron_equation_executed=true") != std::string::npos &&
      line.find("ion_equation_executed=true") != std::string::npos &&
      line.find("electron_assembly_diagnostic_id=p2.diffusion.assembly") != std::string::npos &&
      line.find("ion_assembly_diagnostic_id=p2.diffusion.assembly") != std::string::npos &&
      line.find("updated_fields=e_electron,e_fluid_total") != std::string::npos &&
      line.find("rho_momentum_unchanged=true") != std::string::npos &&
      line.find("fallback_used=false") != std::string::npos &&
      line.find("max_abs_delta_Te_erg=") != std::string::npos &&
      line.find("max_abs_delta_Ti_erg=") != std::string::npos &&
      line.find("global_thermal_energy_residual=") != std::string::npos &&
      line.find("min_e_electron_after=") != std::string::npos &&
      line.find("min_e_ion_after=") != std::string::npos;
  if (!base) {
    return false;
  }
  if (variable) {
    return line.find("kappa_model_requested=") != std::string::npos &&
           line.find("kappa_model_executed=") != std::string::npos &&
           line.find("provider_model_requested=") != std::string::npos &&
           line.find("provider_model_executed=") != std::string::npos &&
           line.find("coefficient_time_level=old_time_lagged") != std::string::npos &&
           line.find("provider_failure_count=0") != std::string::npos &&
           line.find("min_kappa_e_cm_inv_s=") != std::string::npos &&
           line.find("max_kappa_e_cm_inv_s=") != std::string::npos &&
           line.find("min_kappa_i_cm_inv_s=") != std::string::npos &&
           line.find("max_kappa_i_cm_inv_s=") != std::string::npos &&
           line.find("min_provider_lnLambda=") != std::string::npos &&
           line.find("max_provider_lnLambda=") != std::string::npos &&
           line.find("min_provider_f_LM=") != std::string::npos &&
           line.find("max_provider_f_LM=") != std::string::npos;
  }
  return line.find("kappa_model=constant_user_supplied") != std::string::npos &&
         line.find("kappa_e_cm_inv_s=") != std::string::npos &&
         line.find("kappa_i_cm_inv_s=") != std::string::npos;
}

}  // namespace dec3d::transport
