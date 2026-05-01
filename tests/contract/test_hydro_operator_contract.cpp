#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] std::uint32_t HydroStageMask() noexcept {
  return 1u << static_cast<std::uint32_t>(dec3d::core::StageId::hydro);
}

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }

  return false;
}

[[nodiscard]] bool HasDiagnosticMessageFragment(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* fragment) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.message.find(fragment) != std::string::npos) {
      return true;
    }
  }

  return false;
}

}  // namespace

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::core::StageId;
    using dec3d::hydro::HydroOperator;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::IngestStageResult;
    using dec3d::runtime::IsStageRegistered;
    using dec3d::runtime::RegisterStage;
    using dec3d::runtime::ValidateStageRegistry;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::BuildHydroAuthorizedWriteMask;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto runtime = CreateRuntimeScaffold(PhaseId::p1, "p1-v0.1");
    DEC3D_CHECK(RegisterStage(runtime, StageId::hydro));
    DEC3D_CHECK(IsStageRegistered(runtime, StageId::hydro));
    DEC3D_CHECK(!IsStageRegistered(runtime, StageId::thermal));
    DEC3D_CHECK(!IsStageRegistered(runtime, StageId::equilibration));
    DEC3D_CHECK(!IsStageRegistered(runtime, StageId::radiation));
    DEC3D_CHECK(!IsStageRegistered(runtime, StageId::alpha));

    const auto validation = ValidateStageRegistry(runtime);
    DEC3D_CHECK(validation.success);
    DEC3D_CHECK_EQ(validation.required_stage_mask, HydroStageMask());
    DEC3D_CHECK_EQ(validation.registered_stage_mask, HydroStageMask());

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{2, 2, 2, 0.0, 1.0});
    DEC3D_CHECK(geometry.is_valid());

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 0});
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::rho));
    const dec3d::hydro::HydroPrimitiveState left{
        1.0,
        0.0,
        0.30,
        -0.20,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const dec3d::hydro::HydroPrimitiveState right{
        0.125,
        0.0,
        -0.10,
        0.25,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};
    const auto left_state = dec3d::hydro::MakeConservativeState(left);
    const auto right_state = dec3d::hydro::MakeConservativeState(right);

    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(0, theta, phi) = left_state.rho;
        state.mom_r(0, theta, phi) = left_state.mom_r;
        state.mom_theta(0, theta, phi) = left_state.mom_theta;
        state.mom_phi(0, theta, phi) = left_state.mom_phi;
        state.e_fluid_total(0, theta, phi) = left_state.e_fluid_total;
        state.e_electron(0, theta, phi) = ElectronEnergyDensityFromPressure(0.4);

        state.rho(1, theta, phi) = right_state.rho;
        state.mom_r(1, theta, phi) = right_state.mom_r;
        state.mom_theta(1, theta, phi) = right_state.mom_theta;
        state.mom_phi(1, theta, phi) = right_state.mom_phi;
        state.e_fluid_total(1, theta, phi) = right_state.e_fluid_total;
        state.e_electron(1, theta, phi) = ElectronEnergyDensityFromPressure(0.04);
      }
    }

    const auto initial_rho = state.rho(0, 0, 0);
    const auto initial_mom_r = state.mom_r(0, 0, 0);
    const auto initial_mom_theta = state.mom_theta(0, 0, 0);
    const auto initial_mom_phi = state.mom_phi(0, 0, 0);
    const auto initial_e_fluid_total = state.e_fluid_total(0, 0, 0);
    const auto initial_e_electron = state.e_electron(0, 0, 0);

    const StageContext context{
        0.0,
        1.0e-2,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1",
        "ownership:rank0",
        "diagnostics:p1.hydro"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, state));
    DEC3D_CHECK(hydro.is_bound());

    const auto dt_advice = hydro.estimate_dt();
    DEC3D_CHECK(dt_advice.is_complete());
    DEC3D_CHECK(std::isfinite(dt_advice.hard_cap_dt));
    DEC3D_CHECK(dt_advice.hard_cap_dt < std::numeric_limits<double>::infinity());
    DEC3D_CHECK(dt_advice.evidence != "p1.hydro.dt.scaffold.placeholder");
    DEC3D_CHECK(dt_advice.evidence.find("p1.hydro.dt.explicit_cfl") != std::string::npos);

    const auto result = hydro.advance();
    DEC3D_CHECK(result.is_semantically_complete());
    DEC3D_CHECK_EQ(result.updated_fields, BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask, BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::rho));
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::mom_r));
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::mom_theta));
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::mom_phi));
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::e_fluid_total));
    DEC3D_CHECK(dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::e_electron));
    DEC3D_CHECK(!dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::radiation_groups));
    DEC3D_CHECK(!dec3d::core::MaskContains(result.updated_fields, AuthoritativeField::alpha_state));
    DEC3D_CHECK(dec3d::core::MaskContains(state.last_invalidated_cached_mask, CachedField::electron_temperature));
    DEC3D_CHECK(dec3d::core::MaskContains(state.last_invalidated_cached_mask, CachedField::chi_e));
    DEC3D_CHECK(dec3d::core::MaskContains(state.last_invalidated_cached_mask, CachedField::velocity));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.stage.view_writeback"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.stage.hllc_static_grid"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.source.geometric_step"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.theta_pole_contract.executed"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.writeback.electron_channel"));
    DEC3D_CHECK(result.execution_evidence.has_value());
    DEC3D_CHECK_EQ(result.execution_evidence->touched_cell_count, state.rho.size());
    DEC3D_CHECK(std::abs(state.rho(0, 0, 0) - initial_rho) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.mom_r(0, 0, 0) - initial_mom_r) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.mom_theta(0, 0, 0) - initial_mom_theta) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.mom_phi(0, 0, 0) - initial_mom_phi) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.e_fluid_total(0, 0, 0) - initial_e_fluid_total) > 1.0e-8);
    DEC3D_CHECK(std::abs(state.e_electron(0, 0, 0) - initial_e_electron) > 1.0e-8);

    DEC3D_CHECK(IngestStageResult(runtime, StageId::hydro, result));

    const auto runtime_report = BuildRuntimeSubstrateReport(runtime);
    DEC3D_CHECK(runtime_report.is_complete());
    DEC3D_CHECK(runtime_report.required_stage_contract_satisfied);
    DEC3D_CHECK(runtime_report.stage_result_ingestion_complete);
    DEC3D_CHECK(runtime_report.execution_evidence_complete);
    DEC3D_CHECK_EQ(runtime_report.registered_stage_mask, HydroStageMask());
    DEC3D_CHECK_EQ(runtime_report.ingested_stage_mask, HydroStageMask());
    DEC3D_CHECK_EQ(runtime_report.missing_stage_report_mask, static_cast<std::uint32_t>(0));
    DEC3D_CHECK_EQ(runtime_report.authoritative_write_set_digest, BuildHydroAuthorizedWriteMask());

    {
      auto radiation_runtime = CreateRuntimeScaffold(PhaseId::p1, "p1-v0.1");
      DEC3D_CHECK(RegisterStage(radiation_runtime, StageId::hydro));

      auto radiation_state = CanonicalState::Create(CanonicalStateLayout{3, 2, 2, 1});
      const auto radiation_geometry =
          BuildSphericalGeometry(SphericalMeshDescriptor{3, 2, 2, 0.0, 1.0});
      DEC3D_CHECK(radiation_geometry.is_valid());

      for (std::size_t radial = 0; radial < radiation_state.rho.extent_r(); ++radial) {
        for (std::size_t theta = 0; theta < radiation_state.rho.extent_theta(); ++theta) {
          for (std::size_t phi = 0; phi < radiation_state.rho.extent_phi(); ++phi) {
            const double density = 1.0 + 0.1 * static_cast<double>(radial);
            const double velocity = 0.05 + 0.01 * static_cast<double>(radial);
            radiation_state.rho(radial, theta, phi) = density;
            radiation_state.mom_r(radial, theta, phi) = density * velocity;
            radiation_state.mom_theta(radial, theta, phi) = 0.0;
            radiation_state.mom_phi(radial, theta, phi) = 0.0;
            radiation_state.e_fluid_total(radial, theta, phi) = 2.5 + 0.1 * density;
            radiation_state.e_electron(radial, theta, phi) =
                ElectronEnergyDensityFromPressure(0.35 + 0.02 * density);
            radiation_state.radiation_groups[0](radial, theta, phi) =
                1.0 + 0.4 * static_cast<double>(radial);
          }
        }
      }

      const StageContext radiation_context{
          0.0,
          1.0e-3,
          1,
          PhaseId::p1,
          "p1-v0.1",
          "mesh:p1",
          "ownership:rank0",
          "diagnostics:p1.hydro"};
      DEC3D_CHECK(radiation_context.is_complete());

      HydroOperator radiation_hydro;
      DEC3D_CHECK(radiation_hydro.bind(radiation_context, radiation_geometry, radiation_state));
      dec3d::hydro::StaticGridHydroOptions radiation_options{};
      radiation_options.apply_theta_sweep = false;
      radiation_options.apply_phi_sweep = false;
      radiation_options.enable_radiation_hydro_terms = true;
      radiation_hydro.SetStaticGridOptions(radiation_options);

      const auto radiation_result = radiation_hydro.advance();
      DEC3D_CHECK(radiation_result.is_semantically_complete());
      DEC3D_CHECK(dec3d::core::MaskContains(
          radiation_result.updated_fields,
          AuthoritativeField::radiation_groups));
      DEC3D_CHECK(dec3d::core::MaskContains(
          radiation_state.last_authoritative_write_mask,
          AuthoritativeField::radiation_groups));
      DEC3D_CHECK(HasDiagnosticCode(
          radiation_result.diagnostics,
          "p3.radiation.hydro_terms.stage_H"));
      DEC3D_CHECK(HasDiagnosticCode(
          radiation_result.diagnostics,
          "p3.radiation.hydro_terms.writeback"));
      DEC3D_CHECK(HasDiagnosticMessageFragment(
          radiation_result.diagnostics,
          "advected_radiation_scalar=P_g_power_3_over_4"));
      DEC3D_CHECK(HasDiagnosticMessageFragment(
          radiation_result.diagnostics,
          "passive_Ug_advection=false"));

      DEC3D_CHECK(IngestStageResult(radiation_runtime, StageId::hydro, radiation_result));
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
