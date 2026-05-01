#include "runtime/p0_acceptance.hpp"
#include "runtime/p1_acceptance.hpp"
#include "test_assert.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace {

[[nodiscard]] std::string ReadTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::string content;
  std::string line;
  while (std::getline(input, line)) {
    content += line;
    content.push_back('\n');
  }
  return content;
}

dec3d::runtime::BuildVerificationArtifact MakeBuildArtifact() {
  dec3d::runtime::BuildVerificationArtifact artifact;
  artifact.recorded = true;
  artifact.clean_rebuild_verified = true;
  artifact.report_line = "clean_rebuild=true; generator=vs2022";
  return artifact;
}

dec3d::runtime::FailingFirstEvidenceArtifact MakeFailingFirstArtifact() {
  dec3d::runtime::FailingFirstEvidenceArtifact artifact;
  artifact.recorded = true;
  artifact.failing_test_count = 7;
  artifact.report_line = "failing_tests=macro_zoning_case2_case3_case4; rerun=green";
  return artifact;
}

dec3d::runtime::HydroStageArtifact MakeHydroStageArtifact() {
  dec3d::runtime::HydroStageArtifact artifact;
  artifact.recorded = true;
  artifact.runtime_executed = true;
  artifact.diagnostics_complete = true;
  artifact.authorized_write_set = true;
  artifact.report_line = "implementation_id=p1.hydro.operator.hllc_static_grid; runtime_executed=true; diagnostics_complete=true";
  return artifact;
}

dec3d::runtime::AleRuntimeCommitArtifact MakeAleArtifact() {
  dec3d::runtime::AleRuntimeCommitArtifact artifact;
  artifact.recorded = true;
  artifact.contract_satisfied = true;
  artifact.diagnostics_complete = true;
  artifact.report_line = "proposal_commit_contract=true; radial_ale_flux_entry=true";
  return artifact;
}

dec3d::runtime::MacroZoningArtifact MakeMacroZoningArtifact() {
  dec3d::runtime::MacroZoningArtifact artifact;
  artifact.recorded = true;
  artifact.substrate_complete = true;
  artifact.single_rank_coupling_complete = true;
  artifact.mpi_case2_complete = true;
  artifact.mpi_case3_complete = true;
  artifact.mpi_case4_complete = true;
  artifact.diagnostics_complete = true;
  artifact.report_line = "workstream7a=true; workstream7b=true; case2_macro_on_mpi=true; case3_macro_on_mpi=true; case4_macro_on_mpi=true";
  return artifact;
}

dec3d::runtime::RadialMpiParityArtifact MakeRadialMpiArtifact() {
  dec3d::runtime::RadialMpiParityArtifact artifact;
  artifact.recorded = true;
  artifact.case2_complete = true;
  artifact.case3_complete = true;
  artifact.case4_complete = true;
  artifact.diagnostics_complete = true;
  artifact.report_line = "case2=true; case3=true; case4=true";
  return artifact;
}

dec3d::runtime::HydroNumericalCheckArtifact MakeNumericalArtifact() {
  dec3d::runtime::HydroNumericalCheckArtifact artifact;
  artifact.recorded = true;
  artifact.case2_complete = true;
  artifact.case3_complete = true;
  artifact.case4_complete = true;
  artifact.budget_complete = true;
  artifact.diagnostics_complete = true;
  artifact.report_line = "budget_complete=true; case2=true; case3=true; case4=true";
  return artifact;
}

dec3d::runtime::FieldUpdateSummaryArtifact MakeFieldUpdateArtifact() {
  dec3d::runtime::FieldUpdateSummaryArtifact artifact;
  artifact.recorded = true;
  artifact.authoritative_write_mask_nonzero = true;
  artifact.diagnostics_complete = true;
  artifact.report_line = "updated_fields=rho,mom_r,mom_theta,mom_phi,E_fluid_total,E_electron";
  return artifact;
}

dec3d::runtime::AssumptionLedgerDeltaArtifact MakeAssumptionDelta() {
  dec3d::runtime::AssumptionLedgerDeltaArtifact artifact;
  artifact.recorded = true;
  artifact.empty = false;
  artifact.entry_count = 3;
  artifact.report_line = "phase=p1; recorded=true; empty=false; entry_count=3";
  return artifact;
}

}  // namespace

int main() {
  try {
    using dec3d::runtime::BuildP1AcceptanceSummary;
    using dec3d::runtime::WriteP1AcceptanceArtifacts;

    const auto build_artifact = MakeBuildArtifact();
    const auto failing_first = MakeFailingFirstArtifact();
    const auto hydro_stage = MakeHydroStageArtifact();
    const auto ale_artifact = MakeAleArtifact();
    const auto macro_artifact = MakeMacroZoningArtifact();
    const auto radial_mpi_artifact = MakeRadialMpiArtifact();
    const auto numerical_artifact = MakeNumericalArtifact();
    const auto field_update_artifact = MakeFieldUpdateArtifact();
    const auto assumption_delta = MakeAssumptionDelta();

    const auto summary = BuildP1AcceptanceSummary(
        build_artifact,
        failing_first,
        hydro_stage,
        ale_artifact,
        macro_artifact,
        radial_mpi_artifact,
        numerical_artifact,
        field_update_artifact,
        assumption_delta);

    DEC3D_CHECK(summary.is_complete());
    DEC3D_CHECK(summary.success);
    DEC3D_CHECK(summary.workstream7_complete);
    DEC3D_CHECK(!summary.full_p1_closeout_ready);
    DEC3D_CHECK(summary.build_artifact_complete);
    DEC3D_CHECK(summary.macro_zoning_artifact_complete);
    DEC3D_CHECK(summary.radial_mpi_parity_complete);
    DEC3D_CHECK(summary.hydro_numerical_check_complete);
    DEC3D_CHECK(summary.field_update_summary_complete);
    DEC3D_CHECK(summary.assumption_delta_complete);

    const std::filesystem::path artifact_dir = "F:\\dec3d\\artifacts\\p1";
    std::filesystem::remove_all(artifact_dir);
    DEC3D_CHECK(WriteP1AcceptanceArtifacts(
        artifact_dir,
        build_artifact,
        failing_first,
        hydro_stage,
        ale_artifact,
        macro_artifact,
        radial_mpi_artifact,
        numerical_artifact,
        field_update_artifact,
        assumption_delta,
        summary));

    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "build-report.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "failing-first-evidence.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "hydro-stage-report.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "ale-runtime-commit-report.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "macro-zoning-report.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "radial-mpi-parity-report.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "hydro-numerical-check-summary.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "field-update-summary.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "assumption-ledger-delta.md"));
    DEC3D_CHECK(std::filesystem::exists(artifact_dir / "p1-closeout-summary.md"));

    const auto closeout_text = ReadTextFile(artifact_dir / "p1-closeout-summary.md");
    DEC3D_CHECK(closeout_text.find("workstream7_complete=true") != std::string::npos);
    DEC3D_CHECK(closeout_text.find("full_p1_closeout_ready=false") != std::string::npos);

    auto incomplete_macro = MakeMacroZoningArtifact();
    incomplete_macro.mpi_case4_complete = false;
    incomplete_macro.report_line = "workstream7a=true; workstream7b=true; case2_macro_on_mpi=true; case3_macro_on_mpi=true; case4_macro_on_mpi=false";
    const auto incomplete_summary = BuildP1AcceptanceSummary(
        build_artifact,
        failing_first,
        hydro_stage,
        ale_artifact,
        incomplete_macro,
        radial_mpi_artifact,
        numerical_artifact,
        field_update_artifact,
        assumption_delta);
    DEC3D_CHECK(incomplete_summary.is_complete());
    DEC3D_CHECK(!incomplete_summary.success);
    DEC3D_CHECK(!incomplete_summary.workstream7_complete);
    DEC3D_CHECK(!incomplete_summary.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
