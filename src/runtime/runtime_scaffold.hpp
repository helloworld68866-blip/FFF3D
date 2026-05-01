#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/ghost/ghost_topology.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <string>
#include <vector>

namespace dec3d::runtime {

struct StageRegistryEntry {
  dec3d::core::StageId stage_id{dec3d::core::StageId::hydro};
  std::string stage_name;
  std::string binding_id;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct IngestedStageResult {
  dec3d::core::StageId stage_id{dec3d::core::StageId::hydro};
  dec3d::core::StageResult stage_result;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RuntimeDiagnosticsSink {
  bool initialized{false};
  std::string sink_id;
  std::vector<IngestedStageResult> ingested_stage_results;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RuntimeMeshCommitState {
  bool proposal_present{false};
  bool commit_required{false};
  bool commit_performed{false};
  std::string proposal_summary;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RuntimeScaffold {
  bool initialized{false};
  dec3d::core::PhaseId phase_id{dec3d::core::PhaseId::p0};
  std::string contract_version{"p0-v0.1"};
  std::uint32_t registered_stage_mask{0};
  std::vector<StageRegistryEntry> stage_registry_entries;
  RuntimeDiagnosticsSink diagnostics_sink;
  RuntimeMeshCommitState mesh_commit_state;

  [[nodiscard]] bool has_registered_physics_stages() const noexcept {
    return registered_stage_mask != 0;
  }
};

struct RuntimeExecutionEvidence {
  bool stage_registry_initialized{false};
  bool required_stage_validation_performed{false};
  bool stage_machine_initialization_attempted{false};
  bool stage_machine_initialized{false};
  bool diagnostics_sink_initialized{false};
  std::size_t ingested_stage_result_count{0};
  std::size_t ingested_execution_evidence_count{0};
  std::string implementation_id;

  [[nodiscard]] bool is_present() const noexcept;
};

struct StageResultIngestionResult {
  bool attempt_performed{false};
  bool accepted{false};
  bool diagnostics_sink_initialized{false};
  bool stage_registered{false};
  bool stage_required_for_phase{false};
  bool stage_result_complete{false};
  bool duplicate_stage_result{false};
  std::uint32_t required_stage_mask{0};
  std::uint32_t registered_stage_mask{0};
  std::uint32_t ingested_stage_mask{0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct StageRegistryValidationResult {
  bool success{false};
  std::uint32_t required_stage_mask{0};
  std::uint32_t registered_stage_mask{0};
  std::uint32_t missing_required_stage_mask{0};
  std::uint32_t forbidden_registered_stage_mask{0};
  std::string failure_reason;
  std::string report_line;
};

struct StageMachineInitializationResult {
  bool attempt_performed{false};
  bool success{false};
  bool initialized{false};
  bool zero_physics_mode{false};
  std::size_t registry_entry_count{0};
  std::uint32_t required_stage_mask{0};
  std::uint32_t registered_stage_mask{0};
  std::uint32_t missing_required_stage_mask{0};
  std::uint32_t forbidden_registered_stage_mask{0};
  RuntimeExecutionEvidence execution_evidence;
  std::vector<StageRegistryEntry> registry_entries;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct StageReportAggregation {
  bool aggregation_performed{false};
  bool success{false};
  bool validation_succeeded{false};
  bool stage_machine_initialized{false};
  bool diagnostics_sink_initialized{false};
  bool stage_result_ingestion_complete{false};
  bool execution_evidence_complete{false};
  bool diagnostics_complete{false};
  bool mesh_proposal_present{false};
  bool mesh_commit_required{false};
  bool mesh_commit_performed{false};
  bool mesh_commit_complete{false};
  std::size_t registered_stage_count{0};
  std::size_t ingested_stage_result_count{0};
  std::size_t ingested_execution_evidence_count{0};
  std::size_t stage_report_count{0};
  std::uint32_t registered_stage_digest{0};
  std::uint32_t validated_required_stage_digest{0};
  std::uint32_t registered_but_not_ingested_mask{0};
  std::uint32_t ingested_stage_mask{0};
  std::uint32_t missing_stage_report_mask{0};
  dec3d::core::AuthoritativeFieldMask authoritative_write_set_digest{0};
  RuntimeExecutionEvidence execution_evidence_summary;
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string validation_failure_reason;
  std::string stage_machine_failure_reason;
  std::string ingestion_failure_reason;
  std::string mesh_commit_failure_reason;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RuntimeSubstrateReport {
  bool initialized{false};
  dec3d::core::PhaseId phase_id{dec3d::core::PhaseId::p0};
  std::size_t registered_stage_count{0};
  std::size_t stage_registry_entry_count{0};
  std::size_t stage_report_count{0};
  std::size_t ingested_stage_result_count{0};
  std::size_t ingested_execution_evidence_count{0};
  std::size_t diagnostics_entry_count{0};
  bool zero_physics_mode{false};
  bool phase_id_valid{false};
  bool validation_performed{false};
  bool stage_machine_initialization_attempted{false};
  bool stage_machine_initialized{false};
  bool stage_report_aggregation_performed{false};
  bool stage_report_aggregation_succeeded{false};
  bool diagnostics_sink_initialized{false};
  bool stage_result_ingestion_complete{false};
  bool execution_evidence_complete{false};
  bool diagnostics_complete{false};
  bool required_stage_contract_satisfied{false};
  bool mesh_proposal_present{false};
  bool mesh_commit_required{false};
  bool mesh_commit_performed{false};
  std::uint32_t required_stage_mask{0};
  std::uint32_t registered_stage_mask{0};
  std::uint32_t missing_required_stage_mask{0};
  std::uint32_t forbidden_registered_stage_mask{0};
  std::uint32_t registered_stage_digest{0};
  std::uint32_t validated_required_stage_digest{0};
  std::uint32_t registered_but_not_ingested_mask{0};
  std::uint32_t ingested_stage_mask{0};
  std::uint32_t missing_stage_report_mask{0};
  dec3d::core::AuthoritativeFieldMask authoritative_write_set_digest{0};
  RuntimeExecutionEvidence execution_evidence;
  std::string validation_failure_reason;
  std::string stage_machine_failure_reason;
  std::string ingestion_failure_reason;
  std::string mesh_commit_failure_reason;
  std::string stage_report_failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct RuntimeAcceptanceSnapshot {
  bool built{false};
  dec3d::core::PhaseId phase_id{dec3d::core::PhaseId::p0};
  bool runtime_report_complete{false};
  bool execution_evidence_present{false};
  bool diagnostics_complete{false};
  std::uint32_t registered_stage_digest{0};
  std::uint32_t validated_required_stage_digest{0};
  std::uint32_t ingested_stage_digest{0};
  std::uint32_t registered_but_not_ingested_mask{0};
  std::size_t stage_report_count{0};
  std::size_t diagnostics_entry_count{0};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MeshCommitResult {
  bool attempted{false};
  bool success{false};
  bool proposal_present{false};
  bool commit_performed{false};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MeshSubstrateRuntimeReport {
  bool mesh_geometry_built{false};
  bool ownership_built{false};
  bool ghost_topology_initialized{false};
  std::string implementation_id;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] MeshSubstrateRuntimeReport BuildMeshSubstrateRuntimeReport(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& ownership,
    const dec3d::mesh::GhostTopology& topology,
    const dec3d::mesh::MeshBuildEvidence& evidence) noexcept;

[[nodiscard]] RuntimeScaffold CreateRuntimeScaffold(
    dec3d::core::PhaseId phase_id,
    std::string contract_version) noexcept;

[[nodiscard]] bool RegisterStage(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept;

[[nodiscard]] bool IsStageRegistered(
    const RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept;

[[nodiscard]] bool IngestStageResult(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept;

[[nodiscard]] StageResultIngestionResult TryIngestStageResult(
    RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept;

[[nodiscard]] bool HasIngestedStageResult(
    const RuntimeScaffold& runtime,
    dec3d::core::StageId stage_id) noexcept;

[[nodiscard]] StageRegistryValidationResult ValidateStageRegistry(
    const RuntimeScaffold& runtime) noexcept;

[[nodiscard]] StageMachineInitializationResult InitializeStageMachine(
    const RuntimeScaffold& runtime) noexcept;

[[nodiscard]] StageReportAggregation BuildStageReportAggregation(
    const RuntimeScaffold& runtime) noexcept;

[[nodiscard]] RuntimeSubstrateReport BuildRuntimeSubstrateReport(
    const RuntimeScaffold& runtime) noexcept;

[[nodiscard]] MeshCommitResult CommitMeshUpdateProposal(
    RuntimeScaffold& runtime,
    dec3d::mesh::SphericalGeometryMetadata& geometry,
    dec3d::core::StageId stage_id,
    const dec3d::core::StageResult& stage_result) noexcept;

[[nodiscard]] RuntimeAcceptanceSnapshot BuildRuntimeAcceptanceSnapshot(
    const RuntimeScaffold& runtime) noexcept;

}  // namespace dec3d::runtime
