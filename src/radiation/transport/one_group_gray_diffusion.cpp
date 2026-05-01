#include "radiation/transport/one_group_gray_diffusion.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

namespace dec3d::radiation {
namespace {

constexpr const char* kImplementationId = "p3.radiation.one_group_gray_diffusion_v1";

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::string BoolToken(bool value) {
  return value ? "true" : "false";
}

[[nodiscard]] std::string UpdatedFieldsName(dec3d::core::AuthoritativeFieldMask mask) {
  return mask == RadiationWriteMask() ? "radiation_groups" : "none";
}

[[nodiscard]] std::string UpdatedMatterFieldsName(
    dec3d::core::AuthoritativeFieldMask mask) {
  return mask == RadiationMatterWriteMask()
             ? "radiation_groups,e_electron,e_fluid_total"
             : "none";
}

[[nodiscard]] const char* ToString(OneGroupGrayRadiationBoundaryModel model) noexcept {
  switch (model) {
    case OneGroupGrayRadiationBoundaryModel::contract_zero_flux_or_scalar_remap:
      return "contract_zero_flux_or_scalar_remap";
    case OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum:
      return "thesis_marshak_vacuum";
  }
  return "unknown";
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy EffectiveBoundaryPolicy(
    const OneGroupGrayRadiationOptions& options) noexcept {
  auto policy = options.boundary_policy;
  if (options.radiation_boundary_model ==
      OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum) {
    policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::radiation_marshak_vacuum;
  }
  return policy;
}

[[nodiscard]] std::string RequestedRadiationBackendName(
    const OneGroupGrayRadiationOptions& options) {
  (void)options;
  return "serial_dense_reference";
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t linear) noexcept {
  return geometry.cell_volumes[linear];
}

[[nodiscard]] double BudgetScale(double value) noexcept {
  return std::max(1.0, std::abs(value));
}

void BuildReport(
    OneGroupGrayRadiationResult& result,
    bool canonical_state_mutated,
    const char* diagnostic_id = "p3.radiation.one_group_gray_diffusion") {
  const bool metadata_written = result.updated_fields != 0u;
  std::ostringstream out;
  out << "diagnostic_id=" << diagnostic_id
      << "; stage_id=R"
      << "; implementation_id=" << kImplementationId
      << "; radiation_group_count=1"
      << "; radiation_group_mode=gray_full_spectrum"
      << "; frequency_edges_used=false"
      << "; unknown=Ug"
      << "; Ug_unit=erg_per_cm3"
      << "; coefficient_time_level=old_time_lagged"
      << "; Dbar_cm2_per_s=" << result.Dbar_cm2_per_s
      << "; min_Dbar_cm2_per_s=" << result.Dbar_cm2_per_s
      << "; max_Dbar_cm2_per_s=" << result.Dbar_cm2_per_s
      << "; kappaP_cm_inv=" << result.kappaP_cm_inv
      << "; min_kappaP_cm_inv=" << result.kappaP_cm_inv
      << "; max_kappaP_cm_inv=" << result.kappaP_cm_inv
      << "; Bgray_erg_per_cm3=" << result.Bgray_erg_per_cm3
      << "; min_Bgray_erg_per_cm3=" << result.Bgray_erg_per_cm3
      << "; max_Bgray_erg_per_cm3=" << result.Bgray_erg_per_cm3
      << "; Bgray_source=fixed_user_supplied"
      << "; radiation_energy_floor_erg_per_cm3=" << result.radiation_energy_floor_erg_per_cm3
      << "; matter_coupling=fixed_source_no_electron_writeback"
      << "; electron_writeback=false"
      << "; radiation_boundary_model=" << result.radiation_boundary_model
      << "; marshak_enabled=" << BoolToken(result.marshak_enabled)
      << "; marshak_outer_face_count=" << result.marshak_outer_face_count
      << "; marshak_nonzero_diagonal_loss_count="
      << result.marshak_nonzero_diagonal_loss_count
      << "; parity_claim_allowed=false"
      << "; fallback_used=false"
      << "; speed_of_light_cm_per_s=" << dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s
      << "; backend_requested=" << result.backend_requested
      << "; backend_executed=" << result.backend_executed
      << "; matrix_report_present=" << BoolToken(!result.assembly_report.empty())
      << "; solve_report_present=" << BoolToken(!result.solve_report.empty())
      << "; cell_count=" << result.cell_count
      << "; dt_s=" << result.dt_s
      << "; min_Ug_before=" << result.min_Ug_before
      << "; max_Ug_before=" << result.max_Ug_before
      << "; min_Ug_after=" << result.min_Ug_after
      << "; max_Ug_after=" << result.max_Ug_after
      << "; updated_fields=" << UpdatedFieldsName(result.updated_fields)
      << "; metadata_written=" << BoolToken(metadata_written)
      << "; canonical_state_mutated=" << BoolToken(canonical_state_mutated);
  result.report_line = out.str();
}

void BuildMatterCouplingReport(
    OneGroupGrayRadiationMatterCouplingResult& result,
    bool canonical_state_mutated,
    const char* diagnostic_id =
        "p3.radiation.one_group_gray_matter_coupling") {
  std::ostringstream out;
  out << "diagnostic_id=" << diagnostic_id
      << "; stage_id=R"
      << "; implementation_id=p3.radiation.one_group_gray_matter_coupling_v1"
      << "; radiation_electron_coupling_enabled=true"
      << "; radiation_group_mode=gray_full_spectrum"
      << "; group_count=1"
      << "; boundary_model="
      << (result.marshak_enabled ? "thesis_marshak_vacuum"
                                 : "contract_zero_flux_or_scalar_remap")
      << "; marshak_enabled=" << BoolToken(result.marshak_enabled)
      << "; backend_requested=" << result.backend_requested
      << "; backend_executed=" << result.backend_executed
      << "; radiation_report_present=" << BoolToken(!result.radiation_report.empty())
      << "; thermodynamic_recovery_report_present="
      << BoolToken(!result.thermodynamic_recovery_report.empty())
      << "; electron_writeback=true"
      << "; updated_fields=" << UpdatedMatterFieldsName(result.updated_fields)
      << "; opacity_provider=none"
      << "; Bgray_source=fixed_user_supplied"
      << "; Bgray_feedback=disabled"
      << "; advection_enabled=false"
      << "; radiation_pressure_work_enabled=false"
      << "; parity_claim_allowed=false"
      << "; cell_count=" << result.cell_count
      << "; dt_s=" << result.dt_s
      << "; delta_radiation_total=" << result.delta_radiation_total
      << "; delta_electron_total=" << result.delta_electron_total
      << "; source_gain_radiation_total=" << result.source_gain_radiation_total
      << "; boundary_leak_total=" << result.boundary_leak_total
      << "; radiation_electron_exchange_residual="
      << result.radiation_electron_exchange_residual
      << "; global_radiation_plus_electron_residual="
      << result.global_radiation_plus_electron_residual
      << "; min_Ug_after=" << result.min_Ug_after
      << "; min_e_electron_after=" << result.min_e_electron_after
      << "; min_e_ion_after=" << result.min_e_ion_after
      << "; canonical_state_mutated=" << BoolToken(canonical_state_mutated);
  result.report_line = out.str();
}

[[nodiscard]] OneGroupGrayRadiationResult Fail(
    OneGroupGrayRadiationResult result,
    const std::string& reason) {
  const std::string nested_diagnostics = result.failure_diagnostics;
  result.success = false;
  result.updated_fields = 0u;
  result.failure_reason = reason;
  BuildReport(result, false, "p3.radiation.one_group_gray_diffusion.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  if (!nested_diagnostics.empty()) {
    result.failure_diagnostics += "; nested_failure_diagnostics={" + nested_diagnostics + "}";
  }
  return result;
}

[[nodiscard]] OneGroupGrayRadiationMatterCouplingResult FailMatterCoupling(
    OneGroupGrayRadiationMatterCouplingResult result,
    const std::string& reason) {
  const std::string nested = result.failure_diagnostics;
  result.success = false;
  result.updated_fields = 0u;
  result.failure_reason = reason;
  BuildMatterCouplingReport(
      result,
      false,
      "p3.radiation.one_group_gray_matter_coupling.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  if (!nested.empty()) {
    result.failure_diagnostics += "; nested_failure_diagnostics={" + nested + "}";
  }
  return result;
}

[[nodiscard]] bool CoefficientsAreValid(const OneGroupGrayRadiationCoefficients& coefficients) noexcept {
  return std::isfinite(coefficients.Dbar_cm2_per_s) &&
         std::isfinite(coefficients.kappaP_cm_inv) &&
         std::isfinite(coefficients.Bgray_erg_per_cm3) &&
         coefficients.Dbar_cm2_per_s >= 0.0 &&
         coefficients.kappaP_cm_inv >= 0.0 &&
         coefficients.Bgray_erg_per_cm3 >= 0.0;
}

[[nodiscard]] bool ScanRadiationGroup(
    const dec3d::core::Array3D<double>& group,
    double floor,
    double& min_value,
    double& max_value) noexcept {
  bool valid = true;
  min_value = std::numeric_limits<double>::infinity();
  max_value = -std::numeric_limits<double>::infinity();
  for (std::size_t r = 0u; r < group.extent_r(); ++r) {
    for (std::size_t theta = 0u; theta < group.extent_theta(); ++theta) {
      for (std::size_t phi = 0u; phi < group.extent_phi(); ++phi) {
        const double value = group(r, theta, phi);
        valid = valid && std::isfinite(value) && value >= floor;
        min_value = std::min(min_value, value);
        max_value = std::max(max_value, value);
      }
    }
  }
  return valid;
}

}  // namespace

bool OneGroupGrayRadiationResult::is_complete() const noexcept {
  return success && !report_line.empty() && cell_count > 0u &&
         Contains(report_line, "diagnostic_id=p3.radiation.one_group_gray_diffusion") &&
         Contains(report_line, "stage_id=R") &&
         Contains(report_line, "updated_fields=");
}

bool OneGroupGrayRadiationMatterCouplingResult::is_complete() const noexcept {
  return success && !report_line.empty() &&
         Contains(report_line, "diagnostic_id=p3.radiation.one_group_gray_matter_coupling") &&
         Contains(report_line, "stage_id=R") &&
         Contains(report_line, "radiation_electron_coupling_enabled=true") &&
         Contains(report_line, "updated_fields=");
}

dec3d::core::AuthoritativeFieldMask RadiationWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
}

dec3d::core::AuthoritativeFieldMask RadiationMatterWriteMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
         dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total);
}

OneGroupGrayRadiationResult ApplyOneGroupGrayRadiationDiffusion(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupGrayRadiationOptions& options) noexcept {
  OneGroupGrayRadiationResult result;
  result.dt_s = options.dt_s;
  result.Dbar_cm2_per_s = options.coefficients.Dbar_cm2_per_s;
  result.kappaP_cm_inv = options.coefficients.kappaP_cm_inv;
  result.Bgray_erg_per_cm3 = options.coefficients.Bgray_erg_per_cm3;
  result.radiation_energy_floor_erg_per_cm3 = options.radiation_energy_floor_erg_per_cm3;
  result.radiation_boundary_model = ToString(options.radiation_boundary_model);
  result.marshak_enabled =
      options.radiation_boundary_model ==
      OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;

  const auto layout_result = ValidateRadiationGroupStateStorage(options.group_layout, state);
  if (!layout_result.success) {
    return Fail(result, layout_result.failure_reason);
  }
  if (!geometry.is_valid()) {
    return Fail(result, "radiation diffusion requires valid spherical geometry");
  }
  if (!std::isfinite(options.dt_s) || options.dt_s < 0.0) {
    return Fail(result, "radiation diffusion dt_s must be finite and nonnegative");
  }
  if (!std::isfinite(options.radiation_energy_floor_erg_per_cm3)) {
    return Fail(result, "radiation energy floor must be finite");
  }
  if (!CoefficientsAreValid(options.coefficients)) {
    return Fail(result, "gray radiation coefficients must be finite and nonnegative");
  }

  bool input_valid = ScanRadiationGroup(
      state.radiation_groups[0],
      options.radiation_energy_floor_erg_per_cm3,
      result.min_Ug_before,
      result.max_Ug_before);
  result.cell_count = state.radiation_groups[0].size();
  if (!input_valid) {
    return Fail(result, "radiation group contains non-finite value or violates floor");
  }

  const bool no_change = options.dt_s == 0.0 ||
                         (options.coefficients.Dbar_cm2_per_s == 0.0 &&
                          options.coefficients.kappaP_cm_inv == 0.0 &&
                          options.coefficients.Bgray_erg_per_cm3 == 0.0);
  if (no_change) {
    result.success = true;
    result.updated_fields = 0u;
    result.min_Ug_after = result.min_Ug_before;
    result.max_Ug_after = result.max_Ug_before;
    result.assembly_report = "assembly_skipped=true; reason=no_change_contract_path";
    result.solve_report = "solve_skipped=true; reason=no_change_contract_path";
    BuildReport(result, false);
    return result;
  }

  dec3d::transport::GenericDiffusionProblem problem;
  problem.layout = dec3d::transport::DiffusionGridLayout{
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells};
  problem.geometry = geometry;
  problem.dt_s = options.dt_s;
  problem.boundary_policy = EffectiveBoundaryPolicy(options);
  problem.scalar_old = state.radiation_groups[0];
  problem.coefficient_A = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      1.0);
  problem.coefficient_D = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      options.coefficients.Dbar_cm2_per_s);

  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  problem.coefficient_C = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      -c * options.coefficients.kappaP_cm_inv);
  problem.coefficient_B = dec3d::core::Array3D<double>(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      c * options.coefficients.kappaP_cm_inv * options.coefficients.Bgray_erg_per_cm3);

  const auto assembly = dec3d::transport::AssembleGenericImplicitDiffusionSystem(problem);
  result.assembly_report = assembly.report_line;
  result.marshak_outer_face_count = assembly.marshak_outer_face_count;
  result.marshak_nonzero_diagonal_loss_count =
      assembly.marshak_nonzero_diagonal_loss_count;
  if (!assembly.success || !dec3d::transport::ValidateGenericDiffusionAssemblyDiagnostics(assembly)) {
    const std::string reason = assembly.failure_reason.empty()
                                   ? "radiation diffusion assembly diagnostics missing"
                                   : assembly.failure_reason;
    result.failure_diagnostics = assembly.failure_diagnostics;
    return Fail(result, reason);
  }

  const auto solve =
      dec3d::transport::SolveGenericDiffusionReference(assembly, options.serial_reference_options);
  result.solve_report = solve.report_line;
  result.backend_executed = solve.backend;
  if (!solve.success || !dec3d::transport::ValidateGenericDiffusionSolveDiagnostics(solve)) {
    const std::string reason = solve.failure_reason.empty()
                                   ? "radiation diffusion solve diagnostics missing"
                                   : solve.failure_reason;
    result.failure_diagnostics = solve.failure_diagnostics;
    return Fail(result, reason);
  }
  if (solve.scalar_new.size() != state.radiation_groups[0].size()) {
    return Fail(result, "radiation diffusion solve returned wrong scalar size");
  }

  dec3d::core::Array3D<double> staged(
      state.layout.radial_cells,
      state.layout.theta_cells,
      state.layout.phi_cells,
      0.0);
  result.min_Ug_after = std::numeric_limits<double>::infinity();
  result.max_Ug_after = -std::numeric_limits<double>::infinity();

  std::size_t offset = 0u;
  bool output_valid = true;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double value = solve.scalar_new[offset++];
        staged(r, theta, phi) = value;
        output_valid = output_valid &&
                       std::isfinite(value) &&
                       value >= options.radiation_energy_floor_erg_per_cm3;
        result.min_Ug_after = std::min(result.min_Ug_after, value);
        result.max_Ug_after = std::max(result.max_Ug_after, value);
      }
    }
  }
  if (!output_valid) {
    return Fail(result, "radiation solve produced non-finite value or violated floor");
  }

  state.radiation_groups[0] = staged;
  result.updated_fields = RadiationWriteMask();
  state.ApplyAuthoritativeWrite(result.updated_fields);
  result.success = true;
  BuildReport(result, true);
  return result;
}

OneGroupGrayRadiationMatterCouplingResult ApplyOneGroupGrayRadiationMatterCoupling(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupGrayRadiationMatterCouplingOptions& options) noexcept {
  OneGroupGrayRadiationMatterCouplingResult result;
  result.dt_s = options.radiation_options.dt_s;
  result.backend_requested = RequestedRadiationBackendName(options.radiation_options);
  result.marshak_enabled =
      options.radiation_options.radiation_boundary_model ==
      OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;

  if (!std::isfinite(options.electron_energy_floor_erg_per_cm3) ||
      !std::isfinite(options.ion_energy_floor_erg_per_cm3)) {
    return FailMatterCoupling(result, "radiation matter floors must be finite");
  }
  if (state.radiation_groups.empty()) {
    return FailMatterCoupling(result, "radiation matter coupling requires radiation group storage");
  }

  auto staged_state = state;
  const auto radiation_result =
      ApplyOneGroupGrayRadiationDiffusion(staged_state, geometry, options.radiation_options);
  result.radiation_report = radiation_result.report_line;
  result.backend_executed = radiation_result.backend_executed;
  result.cell_count = radiation_result.cell_count;
  result.min_Ug_after = radiation_result.min_Ug_after;

  if (!radiation_result.success || !ValidateOneGroupGrayRadiationDiagnostics(radiation_result)) {
    result.failure_diagnostics = radiation_result.failure_diagnostics;
    return FailMatterCoupling(result, "radiation solve failed before matter coupling");
  }

  if (options.radiation_options.dt_s == 0.0 || radiation_result.updated_fields == 0u) {
    result.success = true;
    result.updated_fields = 0u;
    result.thermodynamic_recovery_report =
        "recovery_skipped=true; reason=no_change_contract_path";
    BuildMatterCouplingReport(result, false);
    return result;
  }

  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double kappa = options.radiation_options.coefficients.kappaP_cm_inv;
  const double bgray = options.radiation_options.coefficients.Bgray_erg_per_cm3;

  bool valid = true;
  double min_electron = std::numeric_limits<double>::infinity();
  double min_ion = std::numeric_limits<double>::infinity();
  std::size_t linear = 0u;

  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double volume = CellVolume(geometry, linear++);
        const double u_old = state.radiation_groups[0](r, theta, phi);
        const double u_new = staged_state.radiation_groups[0](r, theta, phi);
        const double source_gain =
            options.radiation_options.dt_s * c * kappa * (bgray - u_new);
        const double delta_electron = -source_gain;

        staged_state.e_electron(r, theta, phi) =
            state.e_electron(r, theta, phi) + delta_electron;
        staged_state.e_fluid_total(r, theta, phi) =
            state.e_fluid_total(r, theta, phi) + delta_electron;

        double ion = -std::numeric_limits<double>::infinity();
        if (std::isfinite(staged_state.rho(r, theta, phi)) &&
            staged_state.rho(r, theta, phi) > 0.0) {
          const double kinetic =
              0.5 * ((staged_state.mom_r(r, theta, phi) *
                      staged_state.mom_r(r, theta, phi)) +
                     (staged_state.mom_theta(r, theta, phi) *
                      staged_state.mom_theta(r, theta, phi)) +
                     (staged_state.mom_phi(r, theta, phi) *
                      staged_state.mom_phi(r, theta, phi))) /
              staged_state.rho(r, theta, phi);
          ion = staged_state.e_fluid_total(r, theta, phi) - kinetic -
                staged_state.e_electron(r, theta, phi);
        }

        min_electron = std::min(min_electron, staged_state.e_electron(r, theta, phi));
        min_ion = std::min(min_ion, ion);
        valid = valid && std::isfinite(volume) && volume > 0.0 &&
                std::isfinite(u_old) && std::isfinite(u_new) &&
                std::isfinite(source_gain) &&
                std::isfinite(staged_state.e_electron(r, theta, phi)) &&
                std::isfinite(staged_state.e_fluid_total(r, theta, phi)) &&
                staged_state.e_electron(r, theta, phi) >=
                    options.electron_energy_floor_erg_per_cm3 &&
                std::isfinite(ion) &&
                ion >= options.ion_energy_floor_erg_per_cm3;

        result.delta_radiation_total += (u_new - u_old) * volume;
        result.source_gain_radiation_total += source_gain * volume;
        result.delta_electron_total += delta_electron * volume;
      }
    }
  }

  result.min_e_electron_after = min_electron;
  result.min_e_ion_after = min_ion;
  result.boundary_leak_total =
      result.source_gain_radiation_total - result.delta_radiation_total;
  result.radiation_electron_exchange_residual =
      std::abs(result.delta_electron_total + result.source_gain_radiation_total);
  result.global_radiation_plus_electron_residual =
      std::abs(result.delta_radiation_total + result.delta_electron_total +
               result.boundary_leak_total);

  if (!valid) {
    return FailMatterCoupling(result, "radiation matter coupling produced nonphysical staged state");
  }

  const auto recovery = dec3d::state::RecoverThermodynamicState(
      staged_state,
      dec3d::state::ThermodynamicRecoveryOptions{
          options.electron_energy_floor_erg_per_cm3,
          options.ion_energy_floor_erg_per_cm3});
  result.thermodynamic_recovery_report =
      recovery.success ? recovery.recovery_diagnostics : recovery.failure_diagnostics;
  if (!recovery.success || !dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovery)) {
    result.failure_diagnostics = recovery.failure_diagnostics;
    return FailMatterCoupling(result, "staged radiation matter thermodynamic recovery failed");
  }

  const double exchange_scale = BudgetScale(result.source_gain_radiation_total);
  const double total_scale =
      BudgetScale(result.delta_radiation_total + result.delta_electron_total);
  const double exchange_tolerance = 1.0e-12 * exchange_scale;
  const double total_tolerance = 1.0e-12 * total_scale;

  if (result.radiation_electron_exchange_residual > exchange_tolerance) {
    return FailMatterCoupling(result, "radiation electron exchange residual exceeds tolerance");
  }
  if (result.global_radiation_plus_electron_residual > total_tolerance) {
    return FailMatterCoupling(result, "radiation plus electron budget residual exceeds tolerance");
  }
  if (!result.marshak_enabled && std::abs(result.boundary_leak_total) > total_tolerance) {
    return FailMatterCoupling(result, "zero-flux radiation matter coupling reported boundary leak");
  }

  state.radiation_groups[0] = staged_state.radiation_groups[0];
  state.e_electron = staged_state.e_electron;
  state.e_fluid_total = staged_state.e_fluid_total;
  result.updated_fields = RadiationMatterWriteMask();
  state.ApplyAuthoritativeWrite(result.updated_fields);
  result.success = true;
  BuildMatterCouplingReport(result, true);
  return result;
}

bool ValidateOneGroupGrayRadiationDiagnostics(
    const OneGroupGrayRadiationResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }

  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.one_group_gray_diffusion") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "implementation_id=p3.radiation.one_group_gray_diffusion_v1") &&
         Contains(line, "radiation_group_count=1") &&
         Contains(line, "radiation_group_mode=gray_full_spectrum") &&
         Contains(line, "frequency_edges_used=false") &&
         Contains(line, "unknown=Ug") &&
         Contains(line, "Ug_unit=erg_per_cm3") &&
         Contains(line, "coefficient_time_level=old_time_lagged") &&
         Contains(line, "Dbar_cm2_per_s=") &&
         Contains(line, "min_Dbar_cm2_per_s=") &&
         Contains(line, "max_Dbar_cm2_per_s=") &&
         Contains(line, "kappaP_cm_inv=") &&
         Contains(line, "min_kappaP_cm_inv=") &&
         Contains(line, "max_kappaP_cm_inv=") &&
         Contains(line, "Bgray_erg_per_cm3=") &&
         Contains(line, "min_Bgray_erg_per_cm3=") &&
         Contains(line, "max_Bgray_erg_per_cm3=") &&
         Contains(line, "Bgray_source=fixed_user_supplied") &&
         Contains(line, "radiation_energy_floor_erg_per_cm3=") &&
         Contains(line, "matter_coupling=fixed_source_no_electron_writeback") &&
         Contains(line, "electron_writeback=false") &&
         Contains(line, "radiation_boundary_model=") &&
         Contains(line, "marshak_enabled=") &&
         Contains(line, "marshak_outer_face_count=") &&
         Contains(line, "marshak_nonzero_diagonal_loss_count=") &&
         Contains(line, "parity_claim_allowed=false") &&
         Contains(line, "fallback_used=false") &&
         Contains(line, "speed_of_light_cm_per_s=") &&
         Contains(line, "backend_requested=") &&
         Contains(line, "backend_executed=") &&
         Contains(line, "matrix_report_present=true") &&
         Contains(line, "solve_report_present=true") &&
         Contains(line, "updated_fields=") &&
         Contains(line, "metadata_written=") &&
         Contains(line, "canonical_state_mutated=");
}

bool ValidateOneGroupGrayRadiationMatterCouplingDiagnostics(
    const OneGroupGrayRadiationMatterCouplingResult& result) noexcept {
  if (!result.success || !result.is_complete()) {
    return false;
  }

  const std::string& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.one_group_gray_matter_coupling") &&
         Contains(line, "stage_id=R") &&
         Contains(line, "implementation_id=p3.radiation.one_group_gray_matter_coupling_v1") &&
         Contains(line, "radiation_electron_coupling_enabled=true") &&
         Contains(line, "radiation_group_mode=gray_full_spectrum") &&
         Contains(line, "group_count=1") &&
         Contains(line, "boundary_model=") &&
         Contains(line, "marshak_enabled=") &&
         Contains(line, "backend_requested=") &&
         Contains(line, "backend_executed=") &&
         Contains(line, "radiation_report_present=true") &&
         Contains(line, "thermodynamic_recovery_report_present=true") &&
         Contains(line, "electron_writeback=true") &&
         Contains(line, "updated_fields=") &&
         Contains(line, "opacity_provider=none") &&
         Contains(line, "Bgray_source=fixed_user_supplied") &&
         Contains(line, "Bgray_feedback=disabled") &&
         Contains(line, "advection_enabled=false") &&
         Contains(line, "radiation_pressure_work_enabled=false") &&
         Contains(line, "parity_claim_allowed=false") &&
         Contains(line, "delta_radiation_total=") &&
         Contains(line, "delta_electron_total=") &&
         Contains(line, "source_gain_radiation_total=") &&
         Contains(line, "boundary_leak_total=") &&
         Contains(line, "radiation_electron_exchange_residual=") &&
         Contains(line, "global_radiation_plus_electron_residual=") &&
         Contains(line, "min_Ug_after=") &&
         Contains(line, "min_e_electron_after=") &&
         Contains(line, "min_e_ion_after=") &&
         Contains(line, "canonical_state_mutated=");
}

}  // namespace dec3d::radiation
