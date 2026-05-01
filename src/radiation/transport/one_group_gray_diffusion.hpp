#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cstddef>
#include <string>

namespace dec3d::radiation {

struct OneGroupGrayRadiationCoefficients {
  double Dbar_cm2_per_s{0.0};
  double kappaP_cm_inv{0.0};
  double Bgray_erg_per_cm3{0.0};
};

enum class OneGroupGrayRadiationBoundaryModel {
  contract_zero_flux_or_scalar_remap,
  thesis_marshak_vacuum
};

struct OneGroupGrayRadiationOptions {
  double dt_s{0.0};
  RadiationGroupLayout group_layout{MakeGrayFullSpectrumRadiationGroupLayout()};
  OneGroupGrayRadiationCoefficients coefficients;
  dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy;
  OneGroupGrayRadiationBoundaryModel radiation_boundary_model{
      OneGroupGrayRadiationBoundaryModel::contract_zero_flux_or_scalar_remap};
  double radiation_energy_floor_erg_per_cm3{0.0};
  dec3d::transport::GenericDiffusionReferenceSolveOptions serial_reference_options{};
};

struct OneGroupGrayRadiationResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string assembly_report;
  std::string solve_report;
  std::string backend_requested{"serial_dense_reference"};
  std::string backend_executed{"none"};
  std::size_t cell_count{0u};
  double dt_s{0.0};
  double Dbar_cm2_per_s{0.0};
  double kappaP_cm_inv{0.0};
  double Bgray_erg_per_cm3{0.0};
  double radiation_energy_floor_erg_per_cm3{0.0};
  double min_Ug_before{0.0};
  double max_Ug_before{0.0};
  double min_Ug_after{0.0};
  double max_Ug_after{0.0};
  std::string radiation_boundary_model{"contract_zero_flux_or_scalar_remap"};
  bool marshak_enabled{false};
  std::size_t marshak_outer_face_count{0};
  std::size_t marshak_nonzero_diagonal_loss_count{0};

  [[nodiscard]] bool is_complete() const noexcept;
};

struct OneGroupGrayRadiationMatterCouplingOptions {
  OneGroupGrayRadiationOptions radiation_options;
  double electron_energy_floor_erg_per_cm3{0.0};
  double ion_energy_floor_erg_per_cm3{0.0};
};

struct OneGroupGrayRadiationMatterCouplingResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
  std::string radiation_report;
  std::string thermodynamic_recovery_report;
  std::string backend_requested{"serial_dense_reference"};
  std::string backend_executed{"none"};
  std::size_t cell_count{0u};
  double dt_s{0.0};
  double delta_radiation_total{0.0};
  double delta_electron_total{0.0};
  double source_gain_radiation_total{0.0};
  double boundary_leak_total{0.0};
  double radiation_electron_exchange_residual{0.0};
  double global_radiation_plus_electron_residual{0.0};
  double min_Ug_after{0.0};
  double min_e_electron_after{0.0};
  double min_e_ion_after{0.0};
  bool marshak_enabled{false};

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] dec3d::core::AuthoritativeFieldMask RadiationWriteMask() noexcept;

[[nodiscard]] dec3d::core::AuthoritativeFieldMask RadiationMatterWriteMask() noexcept;

[[nodiscard]] OneGroupGrayRadiationResult ApplyOneGroupGrayRadiationDiffusion(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupGrayRadiationOptions& options) noexcept;

[[nodiscard]] OneGroupGrayRadiationMatterCouplingResult
ApplyOneGroupGrayRadiationMatterCoupling(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const OneGroupGrayRadiationMatterCouplingOptions& options) noexcept;

[[nodiscard]] bool ValidateOneGroupGrayRadiationDiagnostics(
    const OneGroupGrayRadiationResult& result) noexcept;

[[nodiscard]] bool ValidateOneGroupGrayRadiationMatterCouplingDiagnostics(
    const OneGroupGrayRadiationMatterCouplingResult& result) noexcept;

}  // namespace dec3d::radiation
