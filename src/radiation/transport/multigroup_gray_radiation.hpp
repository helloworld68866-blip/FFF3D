#pragma once

#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/transport/one_group_gray_diffusion.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace dec3d::radiation {

struct MultigroupGrayRadiationCoefficients {
  std::vector<dec3d::core::Array3D<double>> Dbar_cm2_per_s;
  std::vector<dec3d::core::Array3D<double>> kappaP_cm_inv;
  std::vector<dec3d::core::Array3D<double>> B_erg_per_cm3;
  std::string source{"fixed_user_supplied"};
};

struct MultigroupGrayRadiationMatterCouplingOptions {
  double dt_s{0.0};
  RadiationGroupLayout group_layout;
  MultigroupGrayRadiationCoefficients coefficients;
  dec3d::transport::GenericDiffusionBoundaryPolicy boundary_policy;
  OneGroupGrayRadiationBoundaryModel radiation_boundary_model{
      OneGroupGrayRadiationBoundaryModel::contract_zero_flux_or_scalar_remap};
  double radiation_energy_floor_erg_per_cm3{0.0};
  double electron_energy_floor_erg_per_cm3{0.0};
  double ion_energy_floor_erg_per_cm3{0.0};
  dec3d::transport::GenericDiffusionReferenceSolveOptions serial_reference_options{};
};

struct MultigroupGrayRadiationMatterCouplingResult {
  bool success{false};
  dec3d::core::AuthoritativeFieldMask updated_fields{0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
  std::size_t first_failing_group_index{std::numeric_limits<std::size_t>::max()};

  std::string group_layout_report;
  std::vector<std::string> per_group_radiation_reports;
  std::string thermodynamic_recovery_report;

  std::size_t group_count{0u};
  std::size_t per_group_solve_count{0u};
  bool all_groups_updated{false};
  std::string backend_requested{"serial_dense_reference"};
  std::string backend_executed{"none"};
  std::string per_group_coefficients_source{"fixed_user_supplied"};
  std::string opacity_provider{"none"};
  std::string Bg_source{"fixed_user_supplied"};
  bool marshak_enabled{false};

  double delta_radiation_total_all_groups{0.0};
  double delta_electron_total{0.0};
  double source_gain_radiation_total_all_groups{0.0};
  double boundary_leak_total_all_groups{0.0};
  double radiation_electron_exchange_residual{0.0};
  double global_radiation_plus_electron_residual{0.0};

  std::vector<double> delta_radiation_total_by_group;
  std::vector<double> source_gain_radiation_total_by_group;
  std::vector<double> boundary_leak_total_by_group;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] dec3d::core::AuthoritativeFieldMask MultigroupRadiationMatterWriteMask() noexcept;

[[nodiscard]] MultigroupGrayRadiationMatterCouplingResult
ApplyMultigroupGrayRadiationMatterCoupling(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const MultigroupGrayRadiationMatterCouplingOptions& options) noexcept;

[[nodiscard]] bool ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(
    const MultigroupGrayRadiationMatterCouplingResult& result) noexcept;

}  // namespace dec3d::radiation
