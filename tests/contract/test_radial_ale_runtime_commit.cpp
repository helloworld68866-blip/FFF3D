#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "runtime/runtime_scaffold.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace {

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

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::core::StageId;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::ApplyRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::runtime::BuildRuntimeSubstrateReport;
    using dec3d::runtime::CommitMeshUpdateProposal;
    using dec3d::runtime::CreateRuntimeScaffold;
    using dec3d::runtime::IngestStageResult;
    using dec3d::runtime::RegisterStage;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{4u, 2u, 2u, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());
    const auto initial_geometry = geometry;

    auto state = CanonicalState::Create(CanonicalStateLayout{4u, 2u, 2u, 0u});
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          const double rho = radial < 2u ? 1.0 : 0.5;
          const double v_r = radial < 2u ? -0.15 : 0.30;
          const double v_theta = 0.0;
          const double v_phi = 0.0;
          const double pressure = radial < 2u ? 1.0 : 0.3;
          state.rho(radial, theta, phi) = rho;
          state.mom_r(radial, theta, phi) = rho * v_r;
          state.mom_theta(radial, theta, phi) = rho * v_theta;
          state.mom_phi(radial, theta, phi) = rho * v_phi;
          state.e_fluid_total(radial, theta, phi) =
              pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) +
              0.5 * rho * (v_r * v_r + v_theta * v_theta + v_phi * v_phi);
          state.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(0.4 * pressure);
        }
      }
    }

    const StageContext context{
        0.0,
        1.0e-2,
        1u,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.ale",
        "ownership:rank0",
        "diagnostics:p1.ale"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    StaticGridHydroOptions options;
    options.request_radial_ale_proposal = true;
    options.apply_radial_ale_flux_correction = true;
    hydro.SetStaticGridOptions(options);
    DEC3D_CHECK(hydro.bind(context, geometry, state));

    const auto result = hydro.advance();
    DEC3D_CHECK(result.is_semantically_complete());
    DEC3D_CHECK(result.mesh_update_proposal.has_value());
    DEC3D_CHECK(result.mesh_update_proposal->is_complete(geometry.radial_faces.size()));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.ale.proposal.generated"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.ale.radial_flux_correction.executed"));
    for (std::size_t radial_face = 0; radial_face < geometry.radial_faces.size(); ++radial_face) {
      DEC3D_CHECK(std::abs(
          geometry.radial_faces[radial_face] - initial_geometry.radial_faces[radial_face]) < 1.0e-12);
    }

    auto runtime = CreateRuntimeScaffold(PhaseId::p1, "p1-v0.1");
    DEC3D_CHECK(RegisterStage(runtime, StageId::hydro));
    DEC3D_CHECK(IngestStageResult(runtime, StageId::hydro, result));

    const auto pending_report = BuildRuntimeSubstrateReport(runtime);
    DEC3D_CHECK(pending_report.is_complete());
    DEC3D_CHECK(pending_report.mesh_proposal_present);
    DEC3D_CHECK(pending_report.mesh_commit_required);
    DEC3D_CHECK(!pending_report.mesh_commit_performed);
    DEC3D_CHECK(!pending_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(!pending_report.mesh_commit_failure_reason.empty());
    DEC3D_CHECK(pending_report.report_line.find("mesh_commit_required=true") != std::string::npos);

    const auto commit_result = CommitMeshUpdateProposal(runtime, geometry, StageId::hydro, result);
    DEC3D_CHECK(commit_result.is_complete());
    DEC3D_CHECK(commit_result.success);
    DEC3D_CHECK(commit_result.commit_performed);

    bool any_face_changed = false;
    for (std::size_t radial_face = 0; radial_face < geometry.radial_faces.size(); ++radial_face) {
      if (std::abs(geometry.radial_faces[radial_face] - initial_geometry.radial_faces[radial_face]) > 1.0e-12) {
        any_face_changed = true;
        break;
      }
    }
    DEC3D_CHECK(any_face_changed);

    const auto committed_report = BuildRuntimeSubstrateReport(runtime);
    DEC3D_CHECK(committed_report.is_complete());
    DEC3D_CHECK(committed_report.mesh_proposal_present);
    DEC3D_CHECK(committed_report.mesh_commit_required);
    DEC3D_CHECK(committed_report.mesh_commit_performed);
    DEC3D_CHECK(committed_report.stage_report_aggregation_succeeded);
    DEC3D_CHECK(committed_report.mesh_commit_failure_reason.empty());
    DEC3D_CHECK(committed_report.report_line.find("mesh_commit_performed=true") != std::string::npos);

    auto standalone_geometry = initial_geometry;
    DEC3D_CHECK(ApplyRadialAleMeshUpdateProposal(*result.mesh_update_proposal, standalone_geometry).success);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
