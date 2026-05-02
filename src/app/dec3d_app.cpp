#include "app/dec3d_app.hpp"

#include "alpha/alpha_operator.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "initialization/profile_initializer.hpp"
#include "io/input_deck.hpp"
#include "io/radial_profile.hpp"
#include "io/runtime_output.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "radiation/transport/radiation_flux_limiter.hpp"
#include "state/thermodynamics/electron_ion_equilibration.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "transport/thermal/electron_flux_limiter.hpp"
#include "transport/thermal/thermal_conduction.hpp"

#ifdef DEC3D_ENABLE_HYPRE
#include "alpha/distributed_alpha_operator.hpp"
#include "hydro/driver/macro_ale_hllc_mpi.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/ale/radial_ale_mpi.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "radiation/transport/distributed_multigroup_radiation.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"
#include "transport/thermal/distributed_thermal_conduction.hpp"

#include <mpi.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace dec3d::app {
namespace {

constexpr std::size_t kSerialImplicitRuntimeCellLimit = 64u;
constexpr double kRuntimeSerialReferenceResidualTolerance = 1.0e-2;

#ifdef DEC3D_ENABLE_HYPRE
struct MpiRuntimeContext {
  bool initialized{false};
  bool distributed{false};
  int rank{0};
  int rank_count{1};
};

[[nodiscard]] MpiRuntimeContext DetectMpiRuntime() noexcept {
  MpiRuntimeContext context;
  int initialized = 0;
  MPI_Initialized(&initialized);
  context.initialized = initialized != 0;
  if (!context.initialized) {
    return context;
  }
  MPI_Comm_rank(MPI_COMM_WORLD, &context.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &context.rank_count);
  context.distributed = context.rank_count > 1;
  return context;
}
#endif

Dec3DAppResult Fail(std::string reason, std::string diagnostics) {
  Dec3DAppResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = std::move(diagnostics);
  if (result.failure_diagnostics.empty()) {
    result.failure_diagnostics = "diagnostic_id=p5.io.runtime_entry.failure; failure_reason=" +
                                 result.failure_reason;
  }
  return result;
}

struct RuntimeStageResult {
  bool success{false};
  double hydro_view_wall_s{0.0};
  double hydro_halo_exchange_wall_s{0.0};
  double hydro_ale_proposal_wall_s{0.0};
  double hydro_stage_state_wall_s{0.0};
  double hydro_advance_wall_s{0.0};
  double hydro_snapshot_wall_s{0.0};
  double hydro_scratch_wall_s{0.0};
  double hydro_radial_sweep_wall_s{0.0};
  double hydro_macro_detect_wall_s{0.0};
  double hydro_macro_restrict_wall_s{0.0};
  double hydro_macro_update_wall_s{0.0};
  double hydro_macro_radial_update_wall_s{0.0};
  double hydro_macro_theta_update_wall_s{0.0};
  double hydro_macro_phi_update_wall_s{0.0};
  double hydro_macro_state_update_wall_s{0.0};
  double hydro_macro_prolong_wall_s{0.0};
  double hydro_theta_sweep_wall_s{0.0};
  double hydro_phi_sweep_wall_s{0.0};
  double hydro_commit_wall_s{0.0};
  double hydro_source_wall_s{0.0};
  double hydro_budget_wall_s{0.0};
  double hydro_diagnostics_wall_s{0.0};
  double hydro_writeback_wall_s{0.0};
  double hydro_global_gate_wall_s{0.0};
  double coefficient_provider_wall_s{0.0};
  double flux_limiter_wall_s{0.0};
  double assembly_wall_s{0.0};
  double hypre_setup_wall_s{0.0};
  double hypre_solve_wall_s{0.0};
  double writeback_wall_s{0.0};
  int solver_iterations{0};
  std::size_t lagged_amg_candidate_count{0u};
  std::size_t lagged_amg_reuse_attempted_count{0u};
  std::size_t lagged_amg_reuse_accepted_count{0u};
  std::size_t lagged_amg_rebuild_count{0u};
  std::size_t lagged_amg_fallback_rebuild_count{0u};
  double max_lagged_amg_global_matrix_rel_change{0.0};
  std::size_t global_phi_coupling_count{0u};
  std::size_t global_duplicate_column_row_count{0u};
  bool global_matrix_diagnostics_present{false};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

#ifdef DEC3D_ENABLE_HYPRE
struct RuntimeNohExactInflowState {
  bool enabled{false};
  bool local_outer_rank{false};
  bool complete{true};
  double rho0_g_cm3{0.0};
  double inward_speed_cm_s{0.0};
  double pressure0_dyn_cm2{0.0};
  double chi_e0{0.0};
  std::string report_line{
      "noh_exact_inflow_enabled=false; noh_exact_inflow_local_outer_rank=false"};
  std::string failure_reason;
};
#endif

[[nodiscard]] RuntimeStageResult FailStage(
    char stage,
    std::string reason,
    std::string nested = {}) {
  RuntimeStageResult result;
  result.failure_reason = std::move(reason);
  std::ostringstream out;
  out << "diagnostic_id=p5.runtime.stage.failure"
      << "; stage_id=" << stage
      << "; failure_reason=" << result.failure_reason;
  if (!nested.empty()) {
    out << "; nested_failure_diagnostics={" << nested << "}";
  }
  result.failure_diagnostics = out.str();
  return result;
}

#ifdef DEC3D_ENABLE_HYPRE
[[nodiscard]] bool HasToken(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] std::optional<std::size_t> MaxUnsignedDiagnosticValue(
    std::string_view text,
    std::string_view key) noexcept {
  std::optional<std::size_t> value;
  std::size_t pos = 0u;
  while ((pos = text.find(key, pos)) != std::string_view::npos) {
    pos += key.size();
    std::size_t parsed = 0u;
    bool has_digits = false;
    while (pos < text.size() && text[pos] >= '0' && text[pos] <= '9') {
      has_digits = true;
      parsed = parsed * 10u + static_cast<std::size_t>(text[pos] - '0');
      ++pos;
    }
    if (has_digits) {
      value = value.has_value() ? std::max(*value, parsed) : parsed;
    }
  }
  return value;
}

void AttachDistributedMatrixDiagnostics(
    RuntimeStageResult& result,
    const std::string& report_line) noexcept {
  const auto phi = MaxUnsignedDiagnosticValue(report_line, "global_phi_coupling_count=");
  const auto duplicates =
      MaxUnsignedDiagnosticValue(report_line, "global_duplicate_column_row_count=");
  if (phi.has_value() && duplicates.has_value()) {
    result.global_phi_coupling_count = *phi;
    result.global_duplicate_column_row_count = *duplicates;
    result.global_matrix_diagnostics_present = true;
  }
}
#endif

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start);

void AddHydroAdvanceTimingFromDiagnostics(
    RuntimeStageResult& result,
    const dec3d::core::DiagnosticsPayload& diagnostics);

void AddHydroAdvanceTimingFromStaticGrid(
    RuntimeStageResult& result,
    const dec3d::hydro::StaticGridHydroResult& hydro);

void AppendHydroAdvanceTimingFields(std::ostringstream& report,
                                    const RuntimeStageResult& result);

[[nodiscard]] std::string JoinStageOrder(const std::vector<char>& stages) {
  std::ostringstream out;
  for (std::size_t i = 0u; i < stages.size(); ++i) {
    if (i != 0u) {
      out << ',';
    }
    out << stages[i];
  }
  return out.str();
}

[[nodiscard]] bool StageEnabled(const dec3d::io::InputDeckConfig& config, char stage) noexcept {
  switch (stage) {
    case 'H':
      return config.physics.enable_hydro;
    case 'T':
      return config.physics.enable_thermal;
    case 'E':
      return config.physics.enable_equilibration;
    case 'R':
      return config.physics.enable_radiation;
    case 'A':
      return config.physics.enable_alpha;
    default:
      return false;
  }
}

[[nodiscard]] dec3d::transport::DiffusionBoundaryKind BoundaryKindFromDeck(
    const std::string& value) noexcept {
  if (value == "scalar_origin_remap_required") {
    return dec3d::transport::DiffusionBoundaryKind::scalar_origin_remap_required;
  }
  if (value == "scalar_pole_remap_required") {
    return dec3d::transport::DiffusionBoundaryKind::scalar_pole_remap_required;
  }
  if (value == "neumann_zero_flux") {
    return dec3d::transport::DiffusionBoundaryKind::neumann_zero_flux;
  }
  if (value == "periodic") {
    return dec3d::transport::DiffusionBoundaryKind::periodic;
  }
  if (value == "thesis_marshak_vacuum") {
    return dec3d::transport::DiffusionBoundaryKind::radiation_marshak_vacuum;
  }
  return dec3d::transport::DiffusionBoundaryKind::missing;
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy GenericBoundaryPolicy(
    const dec3d::io::InputDeckConfig& config) noexcept {
  dec3d::transport::GenericDiffusionBoundaryPolicy policy;
  policy.inner_radial = BoundaryKindFromDeck(config.boundaries.inner_radial);
  policy.outer_radial = BoundaryKindFromDeck(config.boundaries.outer_radial);
  policy.theta_lower = BoundaryKindFromDeck(config.boundaries.theta);
  policy.theta_upper = BoundaryKindFromDeck(config.boundaries.theta);
  policy.phi = BoundaryKindFromDeck(config.boundaries.phi);
  return policy;
}

[[nodiscard]] dec3d::radiation::OneGroupGrayRadiationBoundaryModel RadiationBoundaryModel(
    const dec3d::io::InputDeckConfig& config) noexcept {
  if (config.radiation.boundary_model == "thesis_marshak_vacuum") {
    return dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;
  }
  return dec3d::radiation::OneGroupGrayRadiationBoundaryModel::contract_zero_flux_or_scalar_remap;
}

[[nodiscard]] dec3d::transport::ThermalConductionKappaModel ThermalKappaModel(
    const std::string& value) noexcept {
  if (value == "spitzer_no_degeneracy") {
    return dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy;
  }
  if (value == "lee_more_with_degeneracy") {
    return dec3d::transport::ThermalConductionKappaModel::lee_more_with_degeneracy;
  }
  if (value == "constant_user_supplied") {
    return dec3d::transport::ThermalConductionKappaModel::constant_user_supplied;
  }
  return dec3d::transport::ThermalConductionKappaModel::unsupported;
}

[[nodiscard]] dec3d::transport::ElectronFluxLimiterModel ElectronFluxLimiterModelFromDeck(
    const std::string& value) noexcept {
  if (value.empty() || value == "disabled") {
    return dec3d::transport::ElectronFluxLimiterModel::disabled;
  }
  if (value == "minmax_old_time_face_effective_kappa") {
    return dec3d::transport::ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa;
  }
  return dec3d::transport::ElectronFluxLimiterModel::unsupported;
}

[[nodiscard]] dec3d::radiation::RadiationFluxLimiterModel RadiationFluxLimiterModelFromDeck(
    const std::string& value) noexcept {
  if (value.empty() || value == "disabled") {
    return dec3d::radiation::RadiationFluxLimiterModel::disabled;
  }
  if (value == "harmonic_eq_5_209") {
    return dec3d::radiation::RadiationFluxLimiterModel::harmonic_eq_5_209;
  }
  if (value == "larsen_eq_5_210_n2") {
    return dec3d::radiation::RadiationFluxLimiterModel::larsen_eq_5_210_n2;
  }
  if (value == "minmax_eq_5_211") {
    return dec3d::radiation::RadiationFluxLimiterModel::minmax_eq_5_211;
  }
  return dec3d::radiation::RadiationFluxLimiterModel::unsupported;
}

[[nodiscard]] std::filesystem::path FindDefaultTopsTableRoot(
    const std::filesystem::path& input_deck_path) {
  const std::filesystem::path relative("data/opacities/tops_dt_2026_04_27");
  std::vector<std::filesystem::path> candidates;
  candidates.push_back(std::filesystem::current_path() / relative);
  candidates.push_back(std::filesystem::current_path().parent_path() / relative);
  auto cursor = input_deck_path.parent_path();
  for (int i = 0; i < 5 && !cursor.empty(); ++i) {
    candidates.push_back(cursor / relative);
    cursor = cursor.parent_path();
  }
  for (const auto& candidate : candidates) {
    if (std::filesystem::exists(candidate / "metadata.json") &&
        std::filesystem::exists(candidate / "multigroup_opacities.csv")) {
      return candidate;
    }
  }
  return {};
}

[[nodiscard]] dec3d::radiation::TopsOpacityProviderOptions MakeRuntimeTopsProviderOptions(
    const std::filesystem::path& input_deck_path) {
  dec3d::radiation::TopsOpacityProviderOptions provider_options;
  provider_options.table_root = FindDefaultTopsTableRoot(input_deck_path).string();
  provider_options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  provider_options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  provider_options.density_clip_policy = dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return provider_options;
}

#ifdef DEC3D_ENABLE_HYPRE
struct RuntimeRadialGhostOverrideResult {
  bool success{false};
  dec3d::hydro::RadialGhostOverride override;
  std::string failure_reason;
};

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildLocalRuntimeGeometrySlice(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::RadialOwnershipSlice& slice,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  dec3d::mesh::SphericalGeometryMetadata local_geometry;
  if (!geometry.is_valid() ||
      !slice.initialized ||
      slice.local_cell_count() == 0u ||
      slice.end_index >= geometry.radial_faces.size()) {
    return local_geometry;
  }

  local_geometry.valid = true;
  local_geometry.radial_faces.assign(
      geometry.radial_faces.begin() + static_cast<std::ptrdiff_t>(slice.begin_index),
      geometry.radial_faces.begin() + static_cast<std::ptrdiff_t>(slice.end_index + 1u));
  local_geometry.theta_faces = geometry.theta_faces;
  local_geometry.phi_faces = geometry.phi_faces;
  local_geometry.cell_volumes.reserve(slice.local_cell_count() * theta_cells * phi_cells);
  for (std::size_t radial = slice.begin_index; radial < slice.end_index; ++radial) {
    for (std::size_t theta = 0u; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < phi_cells; ++phi) {
        const std::size_t global_index =
            ((radial * theta_cells) + theta) * phi_cells + phi;
        const double volume = geometry.cell_volumes[global_index];
        local_geometry.cell_volumes.push_back(volume);
        local_geometry.global_volume += volume;
      }
    }
  }
  return local_geometry;
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions RuntimeHydroOptions(
    const dec3d::io::InputDeckConfig& config) noexcept {
  dec3d::hydro::StaticGridHydroOptions options;
  options.use_ppm_reconstruction = config.mesh.macro_zoning;
  options.use_macro_zoning = config.mesh.macro_zoning;
  options.apply_radial_ale_flux_correction = config.mesh.moving_mesh;
  options.use_macro_ale_direct_moving_face_hllc = config.mesh.moving_mesh;
  options.enable_radiation_hydro_terms = config.physics.enable_radiation;
  options.enable_alpha_hydro_terms = config.physics.enable_alpha;
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality)) {
    options.apply_radial_sweep = true;
    options.apply_theta_sweep = true;
    options.apply_phi_sweep = false;
    options.macro_zoning_use_phi_ppm = false;
    options.apply_radial_ale_flux_correction = false;
    options.request_radial_ale_proposal = false;
    options.use_macro_ale_direct_moving_face_hllc = false;
  }
  return options;
}

[[nodiscard]] bool HaloExchangeComplete(
    const dec3d::mesh::RadialHaloFieldExchange& exchange,
    std::size_t theta_cells,
    std::size_t phi_cells,
    std::size_t ghost_layers) noexcept {
  return exchange.is_complete(theta_cells, phi_cells, ghost_layers);
}

void FillGhostStateFromExchanges(
    dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>& states,
    const dec3d::mesh::RadialHaloFieldExchange& rho,
    const dec3d::mesh::RadialHaloFieldExchange& mom_r,
    const dec3d::mesh::RadialHaloFieldExchange& mom_theta,
    const dec3d::mesh::RadialHaloFieldExchange& mom_phi,
    const dec3d::mesh::RadialHaloFieldExchange& e_fluid_total,
    const dec3d::mesh::RadialHaloFieldExchange& chi_e,
    const dec3d::mesh::RadialHaloFieldExchange* alpha_chi,
    const std::vector<dec3d::mesh::RadialHaloFieldExchange>& radiation_chi,
    bool inner) {
  for (std::size_t radial = 0u; radial < states.extent_r(); ++radial) {
    for (std::size_t theta = 0u; theta < states.extent_theta(); ++theta) {
      for (std::size_t phi = 0u; phi < states.extent_phi(); ++phi) {
        const auto& rho_values = inner ? rho.inner_ghost_values : rho.outer_ghost_values;
        const auto& mr_values = inner ? mom_r.inner_ghost_values : mom_r.outer_ghost_values;
        const auto& mt_values =
            inner ? mom_theta.inner_ghost_values : mom_theta.outer_ghost_values;
        const auto& mp_values = inner ? mom_phi.inner_ghost_values : mom_phi.outer_ghost_values;
        const auto& e_values =
            inner ? e_fluid_total.inner_ghost_values : e_fluid_total.outer_ghost_values;
        const auto& chi_values = inner ? chi_e.inner_ghost_values : chi_e.outer_ghost_values;
        dec3d::hydro::HydroConservativeState state;
        state.rho = rho_values(radial, theta, phi);
        state.mom_r = mr_values(radial, theta, phi);
        state.mom_theta = mt_values(radial, theta, phi);
        state.mom_phi = mp_values(radial, theta, phi);
        state.e_fluid_total = e_values(radial, theta, phi);
        state.chi_e = chi_values(radial, theta, phi);
        if (alpha_chi != nullptr) {
          const auto& alpha_values =
              inner ? alpha_chi->inner_ghost_values : alpha_chi->outer_ghost_values;
          state.alpha_chi = alpha_values(radial, theta, phi);
        }
        state.radiation_chi.resize(radiation_chi.size(), 0.0);
        for (std::size_t group = 0u; group < radiation_chi.size(); ++group) {
          const auto& group_values =
              inner ? radiation_chi[group].inner_ghost_values
                    : radiation_chi[group].outer_ghost_values;
          state.radiation_chi[group] = group_values(radial, theta, phi);
        }
        states(radial, theta, phi) = std::move(state);
      }
    }
  }
}

[[nodiscard]] RuntimeNohExactInflowState BuildRuntimeNohExactInflowState(
    const dec3d::io::InputDeckConfig& config,
    const dec3d::state::CanonicalState& local_state,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    int rank) {
  RuntimeNohExactInflowState state;
  if (config.boundaries.outer_radial != "noh_exact_inflow") {
    return state;
  }
  state.enabled = true;
  state.complete = false;
  if (rank < 0 || static_cast<std::size_t>(rank) >= decomposition.slices.size()) {
    state.failure_reason = "noh exact inflow requires a valid radial ownership rank";
    return state;
  }
  const auto& slice = decomposition.slices[static_cast<std::size_t>(rank)];
  state.local_outer_rank = slice.initialized &&
                           slice.end_index == decomposition.global_radial_cells;
  if (!state.local_outer_rank) {
    state.complete = true;
    state.report_line =
        "noh_exact_inflow_enabled=true"
        "; noh_exact_inflow_local_outer_rank=false"
        "; noh_exact_inflow_density_time_level=current_time"
        "; noh_exact_inflow_velocity_time_level=constant"
        "; noh_exact_inflow_pressure_time_level=current_time_adiabatic"
        "; noh_exact_inflow_chi_e_time_level=current_time_density_scaled";
    return state;
  }
  if (local_state.layout.radial_cells == 0u ||
      local_state.layout.theta_cells == 0u ||
      local_state.layout.phi_cells == 0u) {
    state.failure_reason = "noh exact inflow requires a non-empty local outer slab";
    return state;
  }

  const std::size_t r = local_state.layout.radial_cells - 1u;
  dec3d::hydro::HydroConservativeState conservative;
  conservative.rho = local_state.rho(r, 0u, 0u);
  conservative.mom_r = local_state.mom_r(r, 0u, 0u);
  conservative.mom_theta = local_state.mom_theta(r, 0u, 0u);
  conservative.mom_phi = local_state.mom_phi(r, 0u, 0u);
  conservative.e_fluid_total = local_state.e_fluid_total(r, 0u, 0u);
  conservative.chi_e = dec3d::state::ChiEFromElectronPressure(
      dec3d::state::ElectronPressureFromElectronEnergyDensity(
          local_state.e_electron(r, 0u, 0u)));
  const auto primitive = dec3d::hydro::RecoverPrimitiveState(conservative);
  if (!primitive.is_physical() || !(primitive.v_r < 0.0)) {
    state.failure_reason =
        "noh exact inflow requires physical outer profile with negative radial velocity";
    return state;
  }

  state.rho0_g_cm3 = primitive.rho;
  state.inward_speed_cm_s = -primitive.v_r;
  state.pressure0_dyn_cm2 = primitive.pressure;
  state.chi_e0 = primitive.chi_e;
  state.complete = true;
  std::ostringstream report;
  report << std::setprecision(17)
         << "noh_exact_inflow_enabled=true"
         << "; noh_exact_inflow_local_outer_rank=true"
         << "; noh_exact_inflow_density_time_level=current_time"
         << "; noh_exact_inflow_velocity_time_level=constant"
         << "; noh_exact_inflow_pressure_time_level=current_time_adiabatic"
         << "; noh_exact_inflow_chi_e_time_level=current_time_density_scaled"
         << "; noh_exact_inflow_rho0_g_cm3=" << state.rho0_g_cm3
         << "; noh_exact_inflow_u0_cm_s=" << state.inward_speed_cm_s
         << "; noh_exact_inflow_pressure0_dyn_cm2=" << state.pressure0_dyn_cm2;
  state.report_line = report.str();
  return state;
}

[[nodiscard]] dec3d::hydro::HydroConservativeState BuildNohExactInflowGhostState(
    const RuntimeNohExactInflowState& noh,
    double radius_cm,
    double time_s) noexcept {
  const double safe_radius = radius_cm > 0.0 ? radius_cm : 1.0;
  const double density_factor =
      (safe_radius + noh.inward_speed_cm_s * std::max(0.0, time_s)) / safe_radius;
  const double density_scale = density_factor * density_factor;
  const dec3d::hydro::HydroPrimitiveState primitive{
      noh.rho0_g_cm3 * density_scale,
      -noh.inward_speed_cm_s,
      0.0,
      0.0,
      noh.pressure0_dyn_cm2 * std::pow(density_scale, dec3d::state::HydroIdealGasGamma()),
      noh.chi_e0 * density_scale};
  return dec3d::hydro::MakeConservativeState(primitive);
}

[[nodiscard]] RuntimeRadialGhostOverrideResult BuildRuntimeRadialGhostOverride(
    const dec3d::state::HydroStateView& local_hydro_view,
    std::size_t ghost_layers,
    const RuntimeNohExactInflowState& noh_exact_inflow,
    const dec3d::mesh::SphericalGeometryMetadata& global_geometry,
    double time_s) {
  RuntimeRadialGhostOverrideResult result;
  if (!local_hydro_view.is_complete()) {
    result.failure_reason = local_hydro_view.failure_reason.empty()
                                ? "distributed hydro requires a complete local hydro view"
                                : local_hydro_view.failure_reason;
    return result;
  }
  const std::size_t theta_cells = local_hydro_view.rho->extent_theta();
  const std::size_t phi_cells = local_hydro_view.rho->extent_phi();
  const auto rho_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*local_hydro_view.rho, ghost_layers);
  const auto mom_r_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*local_hydro_view.mom_r, ghost_layers);
  const auto mom_theta_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*local_hydro_view.mom_theta, ghost_layers);
  const auto mom_phi_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*local_hydro_view.mom_phi, ghost_layers);
  const auto e_total_exchange =
      dec3d::mesh::ExchangeRadialHaloField(*local_hydro_view.e_fluid_total, ghost_layers);
  const auto chi_e_exchange =
      dec3d::mesh::ExchangeRadialHaloField(local_hydro_view.chi_e, ghost_layers);
  if (!HaloExchangeComplete(rho_exchange, theta_cells, phi_cells, ghost_layers) ||
      !HaloExchangeComplete(mom_r_exchange, theta_cells, phi_cells, ghost_layers) ||
      !HaloExchangeComplete(mom_theta_exchange, theta_cells, phi_cells, ghost_layers) ||
      !HaloExchangeComplete(mom_phi_exchange, theta_cells, phi_cells, ghost_layers) ||
      !HaloExchangeComplete(e_total_exchange, theta_cells, phi_cells, ghost_layers) ||
      !HaloExchangeComplete(chi_e_exchange, theta_cells, phi_cells, ghost_layers)) {
    result.failure_reason = "distributed hydro radial halo exchange failed";
    return result;
  }

  std::optional<dec3d::mesh::RadialHaloFieldExchange> alpha_exchange;
  if (local_hydro_view.operator_local_alpha_chi) {
    alpha_exchange = dec3d::mesh::ExchangeRadialHaloField(
        local_hydro_view.alpha_chi,
        ghost_layers);
    if (!HaloExchangeComplete(*alpha_exchange, theta_cells, phi_cells, ghost_layers)) {
      result.failure_reason = "distributed hydro alpha chi radial halo exchange failed";
      return result;
    }
  }

  std::vector<dec3d::mesh::RadialHaloFieldExchange> radiation_exchanges;
  radiation_exchanges.reserve(local_hydro_view.radiation_chi.size());
  for (const auto& group_chi : local_hydro_view.radiation_chi) {
    radiation_exchanges.push_back(dec3d::mesh::ExchangeRadialHaloField(group_chi, ghost_layers));
    if (!HaloExchangeComplete(radiation_exchanges.back(), theta_cells, phi_cells, ghost_layers)) {
      result.failure_reason = "distributed hydro radiation chi radial halo exchange failed";
      return result;
    }
  }

  result.override.has_inner_neighbor = rho_exchange.has_inner_neighbor;
  result.override.has_outer_neighbor =
      rho_exchange.has_outer_neighbor ||
      (noh_exact_inflow.enabled && noh_exact_inflow.local_outer_rank);
  result.override.ghost_layers = ghost_layers;
  result.override.report_line =
      "mpi_radial_halo_exchange=true; hydro_scalar_bundle_halo=true";
  if (result.override.has_inner_neighbor) {
    result.override.inner_ghost_states =
        dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
            ghost_layers,
            theta_cells,
            phi_cells);
    FillGhostStateFromExchanges(
        result.override.inner_ghost_states,
        rho_exchange,
        mom_r_exchange,
        mom_theta_exchange,
        mom_phi_exchange,
        e_total_exchange,
        chi_e_exchange,
        alpha_exchange.has_value() ? &*alpha_exchange : nullptr,
        radiation_exchanges,
        true);
  }
  if (result.override.has_outer_neighbor) {
    result.override.outer_ghost_states =
        dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
            ghost_layers,
            theta_cells,
            phi_cells);
    if (rho_exchange.has_outer_neighbor) {
      FillGhostStateFromExchanges(
          result.override.outer_ghost_states,
          rho_exchange,
          mom_r_exchange,
          mom_theta_exchange,
          mom_phi_exchange,
          e_total_exchange,
          chi_e_exchange,
          alpha_exchange.has_value() ? &*alpha_exchange : nullptr,
          radiation_exchanges,
          false);
    } else {
      const double outer_radius =
          global_geometry.radial_faces.empty() ? 0.0 : global_geometry.radial_faces.back();
      const double outer_dr = global_geometry.radial_faces.size() >= 2u
                                  ? global_geometry.radial_faces.back() -
                                        global_geometry.radial_faces[
                                            global_geometry.radial_faces.size() - 2u]
                                  : 0.0;
      for (std::size_t ghost = 0u; ghost < ghost_layers; ++ghost) {
        const double ghost_radius =
            outer_radius + (static_cast<double>(ghost) + 0.5) * outer_dr;
        const auto ghost_state =
            BuildNohExactInflowGhostState(noh_exact_inflow, ghost_radius, time_s);
        for (std::size_t theta = 0u; theta < theta_cells; ++theta) {
          for (std::size_t phi = 0u; phi < phi_cells; ++phi) {
            result.override.outer_ghost_states(ghost, theta, phi) = ghost_state;
          }
        }
      }
      result.override.outer_ghost_reuses_boundary_partition = true;
    }
  }
  if (noh_exact_inflow.enabled) {
    result.override.report_line += "; " + noh_exact_inflow.report_line;
  }
  result.success = result.override.is_complete(theta_cells, phi_cells, ghost_layers);
  if (!result.success) {
    result.failure_reason = "distributed hydro radial ghost override is incomplete";
  }
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteDistributedHydroStage(
    dec3d::state::CanonicalState& local_state,
    dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    dec3d::mesh::SphericalGeometryMetadata& global_geometry,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    const dec3d::io::InputDeckConfig& config,
    double dt_s,
    double time_s,
    const RuntimeNohExactInflowState& noh_exact_inflow) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  dec3d::state::HydroAlphaViewOptions alpha_options;
  alpha_options.enabled = config.physics.enable_alpha;
  RuntimeStageResult result;

  auto section_start = std::chrono::steady_clock::now();
  auto stage_entry_view = dec3d::state::BuildHydroStateView(local_state, alpha_options);
  result.hydro_view_wall_s += ElapsedSecondsSince(section_start);
  section_start = std::chrono::steady_clock::now();
  auto ghost = BuildRuntimeRadialGhostOverride(
      stage_entry_view,
      3u,
      noh_exact_inflow,
      global_geometry,
      time_s);
  result.hydro_halo_exchange_wall_s += ElapsedSecondsSince(section_start);
  if (!ghost.success) {
    return FailStage('H', ghost.failure_reason);
  }

  auto options = RuntimeHydroOptions(config);
  if (config.mesh.moving_mesh) {
    section_start = std::chrono::steady_clock::now();
    const auto proposal = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
        local_state.rho,
        local_state.mom_r,
        decomposition,
        global_geometry.radial_faces,
        dt_s,
        1.0,
        MPI_COMM_WORLD);
    result.hydro_ale_proposal_wall_s += ElapsedSecondsSince(section_start);
    if (!proposal.is_complete()) {
      return FailStage('H',
                       proposal.failure_reason.empty()
                           ? "distributed hydro ALE proposal failed"
                           : proposal.failure_reason,
                       proposal.report_line);
    }
    auto global_geometry_candidate = global_geometry;
    const auto global_geometry_commit =
        dec3d::mesh::ApplyRadialAleMeshUpdateProposal(proposal.proposal,
                                                      global_geometry_candidate);
    if (!global_geometry_commit.success) {
      return FailStage('H',
                       global_geometry_commit.failure_reason.empty()
                           ? "distributed hydro global geometry preview failed"
                           : global_geometry_commit.failure_reason,
                       global_geometry_commit.report_line);
    }
    section_start = std::chrono::steady_clock::now();
    auto hydro_view = dec3d::state::BuildHydroStateView(local_state, alpha_options);
    result.hydro_stage_state_wall_s += ElapsedSecondsSince(section_start);
    section_start = std::chrono::steady_clock::now();
    const auto hydro = dec3d::hydro::AdvanceMacroAleHllcMpiStep(
        hydro_view,
        local_geometry,
        decomposition,
        proposal.proposal,
        dt_s,
        ghost.override,
        options,
        MPI_COMM_WORLD);
    result.hydro_advance_wall_s += ElapsedSecondsSince(section_start);
    if (!hydro.success || !hydro.is_complete()) {
      return FailStage('H',
                       hydro.failure_reason.empty()
                           ? "distributed macro ALE HLLC hydro failed"
                           : hydro.failure_reason,
                       hydro.diagnostics.report_line);
    }
    global_geometry = std::move(global_geometry_candidate);
    result.success = true;
    std::ostringstream report;
    report << "diagnostic_id=p5.runtime.stage"
           << "; stage_id=H"
           << "; stage_backend=mpi_macro_ale_hllc"
           << "; hydro_runtime_mode=distributed_mpi_hydro"
           << "; hydro_mpi_global_stage_ok=" << (hydro.global_stage_ok ? "true" : "false")
           << "; hydro_mpi_global_publish_ok=" << (hydro.global_publish_ok ? "true" : "false")
           << "; hydro_halo_exchange_used=true"
           << "; hydro_scalar_bundle_halo=true"
           << "; moving_mesh=" << (config.mesh.moving_mesh ? "true" : "false")
           << "; macro_zoning=" << (config.mesh.macro_zoning ? "true" : "false")
           << "; ppm_reconstruction=true"
           << "; rank0_gather_solve_used=false"
           << "; h_view_wall_s=" << result.hydro_view_wall_s
           << "; h_halo_exchange_wall_s=" << result.hydro_halo_exchange_wall_s
           << "; h_ale_proposal_wall_s=" << result.hydro_ale_proposal_wall_s
           << "; h_stage_state_wall_s=" << result.hydro_stage_state_wall_s
           << "; h_hydro_advance_wall_s=" << result.hydro_advance_wall_s;
    AppendHydroAdvanceTimingFields(report, result);
    report << "; h_writeback_wall_s=" << result.hydro_writeback_wall_s
           << "; h_global_gate_wall_s=" << result.hydro_global_gate_wall_s
           << "; nested_ghost_report={" << ghost.override.report_line << "}"
           << "; alpha_hydro_terms_report_present=true"
           << "; alpha_hydro_terms_enabled=" << (config.physics.enable_alpha ? "true" : "false")
           << "; alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure"
           << "; radiation_hydro_terms_report_present=true"
           << "; radiation_hydro_terms_enabled="
           << (config.physics.enable_radiation ? "true" : "false")
           << "; radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"
           << "; nested_hydro_report={" << hydro.diagnostics.report_line << "}";
    result.report_line = report.str();
    return result;
  }

  section_start = std::chrono::steady_clock::now();
  auto staged_state = local_state;
  auto staged_view = dec3d::state::BuildHydroWorkView(staged_state, alpha_options);
  result.hydro_stage_state_wall_s += ElapsedSecondsSince(section_start);
  if (!staged_view.is_complete()) {
    return FailStage('H',
                     staged_view.failure_reason.empty()
                         ? "distributed static hydro work view failed"
                         : staged_view.failure_reason,
                     staged_view.report_line);
  }
  options.use_macro_ale_direct_moving_face_hllc = false;
  section_start = std::chrono::steady_clock::now();
    const auto hydro = dec3d::hydro::AdvanceStaticGridHydro(
        staged_view,
        local_geometry,
        dt_s,
        ghost.override,
        options);
    result.hydro_advance_wall_s += ElapsedSecondsSince(section_start);
    AddHydroAdvanceTimingFromStaticGrid(result, hydro);
    const int local_stage_ok = hydro.success ? 1 : 0;
  int global_stage_ok = 0;
  section_start = std::chrono::steady_clock::now();
  MPI_Allreduce(&local_stage_ok, &global_stage_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  result.hydro_global_gate_wall_s += ElapsedSecondsSince(section_start);
  if (global_stage_ok != 1) {
    return FailStage('H',
                     hydro.failure_reason.empty()
                         ? "distributed static hydro global stage failed"
                         : hydro.failure_reason,
                     hydro.report_line);
  }
  auto updated_fields = dec3d::state::BuildHydroAuthorizedWriteMask();
  if (config.physics.enable_radiation) {
    updated_fields |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
  }
  if (config.physics.enable_alpha) {
    updated_fields |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
  }
  section_start = std::chrono::steady_clock::now();
  const auto writeback =
      dec3d::state::CommitHydroWriteback(local_state, staged_view, updated_fields);
  result.hydro_writeback_wall_s += ElapsedSecondsSince(section_start);
  const int local_publish_ok = writeback.success ? 1 : 0;
  int global_publish_ok = 0;
  section_start = std::chrono::steady_clock::now();
  MPI_Allreduce(&local_publish_ok, &global_publish_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  result.hydro_global_gate_wall_s += ElapsedSecondsSince(section_start);
  if (global_publish_ok != 1) {
    return FailStage('H',
                     writeback.failure_reason.empty()
                         ? "distributed static hydro publish failed"
                         : writeback.failure_reason,
                     writeback.report_line);
  }

  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=H"
         << "; stage_backend=mpi_static_hllc"
         << "; hydro_runtime_mode=distributed_mpi_hydro"
         << "; hydro_mpi_global_stage_ok=true"
         << "; hydro_mpi_global_publish_ok=true"
         << "; hydro_halo_exchange_used=true"
         << "; hydro_scalar_bundle_halo=true"
         << "; moving_mesh=false"
         << "; macro_zoning=" << (config.mesh.macro_zoning ? "true" : "false")
         << "; ppm_reconstruction=" << (options.use_ppm_reconstruction ? "true" : "false")
         << "; rank0_gather_solve_used=false"
         << "; h_view_wall_s=" << result.hydro_view_wall_s
         << "; h_halo_exchange_wall_s=" << result.hydro_halo_exchange_wall_s
         << "; h_ale_proposal_wall_s=" << result.hydro_ale_proposal_wall_s
         << "; h_stage_state_wall_s=" << result.hydro_stage_state_wall_s
         << "; h_hydro_advance_wall_s=" << result.hydro_advance_wall_s;
  AppendHydroAdvanceTimingFields(report, result);
  report << "; h_writeback_wall_s=" << result.hydro_writeback_wall_s
         << "; h_global_gate_wall_s=" << result.hydro_global_gate_wall_s
         << "; nested_ghost_report={" << ghost.override.report_line << "}"
         << "; alpha_hydro_terms_report_present=true"
         << "; alpha_hydro_terms_enabled=" << (config.physics.enable_alpha ? "true" : "false")
         << "; alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure"
         << "; radiation_hydro_terms_report_present=true"
         << "; radiation_hydro_terms_enabled="
         << (config.physics.enable_radiation ? "true" : "false")
         << "; radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"
         << "; nested_hydro_report={" << hydro.report_line << "}"
         << "; rank=" << rank;
  result.report_line = report.str();
  return result;
}

template <typename Field>
void CopyFieldSlabFromGlobal(
    const Field& global,
    Field& local,
    std::size_t global_radial_begin) {
  for (std::size_t lr = 0u; lr < local.extent_r(); ++lr) {
    const std::size_t gr = global_radial_begin + lr;
    for (std::size_t theta = 0u; theta < local.extent_theta(); ++theta) {
      for (std::size_t phi = 0u; phi < local.extent_phi(); ++phi) {
        local(lr, theta, phi) = global(gr, theta, phi);
      }
    }
  }
}

[[nodiscard]] dec3d::state::CanonicalState MakeLocalSlabState(
    const dec3d::state::CanonicalState& global,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership) {
  const std::size_t local_radial =
      ownership.global_radial_end - ownership.global_radial_begin;
  auto local = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{
          local_radial,
          ownership.global_theta_cells,
          ownership.global_phi_cells,
          global.radiation_groups.size()});
  CopyFieldSlabFromGlobal(global.rho, local.rho, ownership.global_radial_begin);
  CopyFieldSlabFromGlobal(global.mom_r, local.mom_r, ownership.global_radial_begin);
  CopyFieldSlabFromGlobal(global.mom_theta, local.mom_theta, ownership.global_radial_begin);
  CopyFieldSlabFromGlobal(global.mom_phi, local.mom_phi, ownership.global_radial_begin);
  CopyFieldSlabFromGlobal(
      global.e_fluid_total, local.e_fluid_total, ownership.global_radial_begin);
  CopyFieldSlabFromGlobal(
      global.e_electron, local.e_electron, ownership.global_radial_begin);
  for (std::size_t group = 0u; group < global.radiation_groups.size(); ++group) {
    CopyFieldSlabFromGlobal(
        global.radiation_groups[group],
        local.radiation_groups[group],
        ownership.global_radial_begin);
  }
  CopyFieldSlabFromGlobal(
      global.alpha_state.storage,
      local.alpha_state.storage,
      ownership.global_radial_begin);
  return local;
}

template <typename Field>
void AllgatherFieldSlabToGlobal(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const Field& local,
    Field& global) {
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  const int local_count = static_cast<int>(local.size());
  const int local_displacement =
      static_cast<int>(ownership.global_radial_begin * plane_size);
  std::vector<int> counts(static_cast<std::size_t>(ownership.rank_count), 0);
  std::vector<int> displacements(static_cast<std::size_t>(ownership.rank_count), 0);
  MPI_Allgather(
      &local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, communicator);
  MPI_Allgather(
      &local_displacement, 1, MPI_INT, displacements.data(), 1, MPI_INT, communicator);
  MPI_Allgatherv(
      local.storage().data(),
      local_count,
      MPI_DOUBLE,
      global.storage().data(),
      counts.data(),
      displacements.data(),
      MPI_DOUBLE,
      communicator);
}

struct RuntimeSlabGatherPlan {
  bool success{false};
  int local_count{0};
  int local_displacement{0};
  std::vector<int> counts;
  std::vector<int> displacements;
};

[[nodiscard]] RuntimeSlabGatherPlan BuildRuntimeSlabGatherPlan(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership) {
  RuntimeSlabGatherPlan plan;
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  const std::size_t local_radial =
      ownership.global_radial_end - ownership.global_radial_begin;
  plan.local_count = static_cast<int>(local_radial * plane_size);
  plan.local_displacement = static_cast<int>(ownership.global_radial_begin * plane_size);
  plan.counts.assign(static_cast<std::size_t>(ownership.rank_count), 0);
  plan.displacements.assign(static_cast<std::size_t>(ownership.rank_count), 0);
  MPI_Allgather(
      &plan.local_count, 1, MPI_INT, plan.counts.data(), 1, MPI_INT, communicator);
  MPI_Allgather(
      &plan.local_displacement,
      1,
      MPI_INT,
      plan.displacements.data(),
      1,
      MPI_INT,
      communicator);
  plan.success = true;
  return plan;
}

void AllgatherLocalStateToGlobal(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::state::CanonicalState& local,
    dec3d::state::CanonicalState& global) {
  AllgatherFieldSlabToGlobal(communicator, ownership, local.rho, global.rho);
  AllgatherFieldSlabToGlobal(communicator, ownership, local.mom_r, global.mom_r);
  AllgatherFieldSlabToGlobal(communicator, ownership, local.mom_theta, global.mom_theta);
  AllgatherFieldSlabToGlobal(communicator, ownership, local.mom_phi, global.mom_phi);
  AllgatherFieldSlabToGlobal(
      communicator, ownership, local.e_fluid_total, global.e_fluid_total);
  AllgatherFieldSlabToGlobal(
      communicator, ownership, local.e_electron, global.e_electron);
  for (std::size_t group = 0u; group < global.radiation_groups.size(); ++group) {
    AllgatherFieldSlabToGlobal(
        communicator, ownership, local.radiation_groups[group], global.radiation_groups[group]);
  }
  AllgatherFieldSlabToGlobal(
      communicator, ownership, local.alpha_state.storage, global.alpha_state.storage);
}

template <typename Field>
void GatherFieldSlabToRoot(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const RuntimeSlabGatherPlan& plan,
    const Field& local,
    Field& global,
    int root_rank) {
  const bool is_root = ownership.rank == root_rank;
  MPI_Gatherv(
      local.storage().data(),
      plan.local_count,
      MPI_DOUBLE,
      is_root ? global.storage().data() : nullptr,
      is_root ? plan.counts.data() : nullptr,
      is_root ? plan.displacements.data() : nullptr,
      MPI_DOUBLE,
      root_rank,
      communicator);
}

void GatherLocalStateToRoot(
    MPI_Comm communicator,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const RuntimeSlabGatherPlan& gather_plan,
    const dec3d::state::CanonicalState& local,
    dec3d::state::CanonicalState& global,
    int root_rank) {
  GatherFieldSlabToRoot(communicator, ownership, gather_plan, local.rho, global.rho, root_rank);
  GatherFieldSlabToRoot(
      communicator, ownership, gather_plan, local.mom_r, global.mom_r, root_rank);
  GatherFieldSlabToRoot(
      communicator, ownership, gather_plan, local.mom_theta, global.mom_theta, root_rank);
  GatherFieldSlabToRoot(
      communicator, ownership, gather_plan, local.mom_phi, global.mom_phi, root_rank);
  GatherFieldSlabToRoot(
      communicator, ownership, gather_plan, local.e_fluid_total, global.e_fluid_total, root_rank);
  GatherFieldSlabToRoot(
      communicator, ownership, gather_plan, local.e_electron, global.e_electron, root_rank);
  for (std::size_t group = 0u; group < global.radiation_groups.size(); ++group) {
    GatherFieldSlabToRoot(
        communicator,
        ownership,
        gather_plan,
        local.radiation_groups[group],
        global.radiation_groups[group],
        root_rank);
  }
  GatherFieldSlabToRoot(
      communicator,
      ownership,
      gather_plan,
      local.alpha_state.storage,
      global.alpha_state.storage,
      root_rank);
}

[[nodiscard]] dec3d::io::RuntimeOutputExtrema ReduceRuntimeExtremaToRoot(
    MPI_Comm communicator,
    const dec3d::io::RuntimeOutputExtrema& local,
    int root_rank) {
  double mins[3] = {local.rho_min, local.te_min_keV, local.ti_min_keV};
  double maxs[3] = {local.rho_max, local.te_max_keV, local.ti_max_keV};
  double global_mins[3] = {0.0, 0.0, 0.0};
  double global_maxs[3] = {0.0, 0.0, 0.0};
  MPI_Reduce(mins, global_mins, 3, MPI_DOUBLE, MPI_MIN, root_rank, communicator);
  MPI_Reduce(maxs, global_maxs, 3, MPI_DOUBLE, MPI_MAX, root_rank, communicator);
  dec3d::io::RuntimeOutputExtrema global;
  global.rho_min = global_mins[0];
  global.te_min_keV = global_mins[1];
  global.ti_min_keV = global_mins[2];
  global.rho_max = global_maxs[0];
  global.te_max_keV = global_maxs[1];
  global.ti_max_keV = global_maxs[2];
  return global;
}

[[nodiscard]] dec3d::radiation::TopsOpacityProviderOptions RuntimeTopsProviderOptions(
    const std::filesystem::path& input_deck_path) {
  dec3d::radiation::TopsOpacityProviderOptions provider_options;
  provider_options.table_root = FindDefaultTopsTableRoot(input_deck_path).string();
  provider_options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  provider_options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  provider_options.density_clip_policy = dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return provider_options;
}
#endif

[[nodiscard]] RuntimeStageResult ExecuteHydroStage(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::io::InputDeckConfig& config,
    double dt_s,
    std::uint64_t step,
    double time_s) {
  dec3d::hydro::HydroOperator hydro;
  auto options = RuntimeHydroOptions(config);
  hydro.SetStaticGridOptions(options);

  dec3d::core::StageContext context;
  context.time_s = time_s;
  context.dt_s = dt_s;
  context.step = step;
  context.phase_id = dec3d::core::PhaseId::p1;
  context.contract_version = "p5.runtime.serial_callable.v1";
  context.mesh_snapshot_handle = "p5.runtime.mesh.current";
  context.ownership_handle = "single_rank_full_domain";
  context.diagnostics_sink_handle = "p5.runtime.dec3d_out";
  if (!hydro.bind(context, geometry, state)) {
    return FailStage('H', "hydro bind failed");
  }
  RuntimeStageResult result;
  auto section_start = std::chrono::steady_clock::now();
  const auto hydro_result = hydro.advance();
  result.hydro_advance_wall_s += ElapsedSecondsSince(section_start);
  if (!hydro_result.success) {
    std::ostringstream nested;
    for (const auto& entry : hydro_result.diagnostics.entries) {
      nested << '[' << entry.code << ':' << entry.message << ']';
    }
    return FailStage('H',
                     hydro_result.failure_reason.empty() ? "hydro advance failed"
                                                          : hydro_result.failure_reason,
                     nested.str());
  }
  AddHydroAdvanceTimingFromDiagnostics(result, hydro_result.diagnostics);
  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=H"
         << "; stage_backend=hydro_operator"
         << "; updated_fields_mask=" << hydro_result.updated_fields
         << "; alpha_hydro_terms_report_present=true"
         << "; alpha_hydro_terms_enabled="
         << (config.physics.enable_alpha ? "true" : "false")
         << "; alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure"
         << "; radiation_hydro_terms_report_present=true"
         << "; radiation_hydro_terms_enabled="
         << (config.physics.enable_radiation ? "true" : "false")
         << "; radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure"
         << "; h_hydro_advance_wall_s=" << result.hydro_advance_wall_s;
  AppendHydroAdvanceTimingFields(report, result);
  result.report_line = report.str();
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteThermalStage(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::io::InputDeckConfig& config,
    double dt_s) {
  if (state.rho.size() > kSerialImplicitRuntimeCellLimit) {
    return FailStage('T', "serial T runtime exceeds implicit cell limit; distributed P5 runtime required");
  }
  dec3d::transport::ThermalConductionOptions options;
  options.dt_s = dt_s;
  options.kappa_model = ThermalKappaModel(config.thermal.kappa_model);
  options.backend = dec3d::transport::ThermalConductionBackend::serial_dense_reference;
  options.boundary_policy = GenericBoundaryPolicy(config);
  options.serial_reference_options.row_limit = kSerialImplicitRuntimeCellLimit;
  options.serial_reference_options.residual_tolerance = kRuntimeSerialReferenceResidualTolerance;
  const auto thermal = dec3d::transport::ApplyThermalConduction(state, geometry, options);
  if (!thermal.success || !dec3d::transport::ValidateThermalConductionDiagnostics(thermal)) {
    return FailStage('T',
                     thermal.failure_reason.empty() ? "thermal conduction failed"
                                                    : thermal.failure_reason,
                     thermal.failure_diagnostics.empty() ? thermal.report_line
                                                         : thermal.failure_diagnostics);
  }
  RuntimeStageResult result;
  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=T"
         << "; stage_backend=" << thermal.backend_executed
         << "; thermal_stage_report_present=true"
         << "; runtime_serial_reference_residual_tolerance="
         << kRuntimeSerialReferenceResidualTolerance
         << "; updated_fields_mask=" << thermal.updated_fields;
  result.report_line = report.str();
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteEquilibrationStage(
    dec3d::state::CanonicalState& state,
    double dt_s) {
  dec3d::state::ElectronIonEquilibrationOptions options;
  options.dt_s = dt_s;
  options.tau_model = dec3d::state::ElectronIonTauModel::thesis_spitzer_eq_5_241;
  const auto equilibration = dec3d::state::ApplyLocalElectronIonEquilibration(state, options);
  if (!equilibration.success ||
      !dec3d::state::ValidateElectronIonEquilibrationDiagnostics(equilibration)) {
    return FailStage('E',
                     equilibration.failure_reason.empty() ? "electron-ion equilibration failed"
                                                          : equilibration.failure_reason,
                     equilibration.failure_diagnostics.empty() ? equilibration.report_line
                                                               : equilibration.failure_diagnostics);
  }
  RuntimeStageResult result;
  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=E"
         << "; stage_backend=local_exact_exponential"
         << "; equilibration_stage_report_present=true"
         << "; updated_fields_mask=" << equilibration.updated_fields;
  result.report_line = report.str();
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteRadiationStage(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const dec3d::io::InputDeckConfig& config,
    const std::filesystem::path& input_deck_path,
    const dec3d::radiation::TopsOpacityTable* runtime_opacity_table,
    const dec3d::radiation::TopsOpacityProviderOptions* runtime_provider_options,
    double dt_s) {
  if (state.rho.size() > kSerialImplicitRuntimeCellLimit) {
    return FailStage('R', "serial R runtime exceeds implicit cell limit; distributed P5 runtime required");
  }
  if (config.radiation.opacity_provider != "tops_dt_tabulated") {
    return FailStage('R', "unsupported runtime opacity_provider");
  }
  dec3d::radiation::TopsOpacityProviderOptions provider_options =
      runtime_provider_options != nullptr
          ? *runtime_provider_options
          : MakeRuntimeTopsProviderOptions(input_deck_path);
  dec3d::radiation::TopsOpacityTable loaded_table;
  const dec3d::radiation::TopsOpacityTable* table = runtime_opacity_table;
  if (provider_options.table_root.empty()) {
    return FailStage('R', "TOPS opacity table root not found");
  }
  if (table == nullptr) {
    const auto table_result = dec3d::radiation::LoadTopsOpacityTable(provider_options);
    if (!table_result.success) {
      return FailStage('R',
                       table_result.failure_reason.empty() ? "TOPS opacity table load failed"
                                                           : table_result.failure_reason,
                       table_result.failure_diagnostics);
    }
    loaded_table = table_result.table;
    table = &loaded_table;
  }
  const auto provider =
      dec3d::radiation::BuildRadiationCoefficientArrays(
          state, geometry, group_layout, *table, provider_options);
  if (!provider.success ||
      !dec3d::radiation::ValidateRadiationCoefficientProviderDiagnostics(provider)) {
    return FailStage('R',
                     provider.failure_reason.empty() ? "radiation coefficient provider failed"
                                                     : provider.failure_reason,
                     provider.failure_diagnostics.empty() ? provider.report_line
                                                          : provider.failure_diagnostics);
  }
  dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions options;
  options.dt_s = dt_s;
  options.group_layout = group_layout;
  options.coefficients = provider.coefficients;
  options.boundary_policy = GenericBoundaryPolicy(config);
  options.radiation_boundary_model = RadiationBoundaryModel(config);
  options.serial_reference_options.row_limit = kSerialImplicitRuntimeCellLimit;
  const auto radiation =
      dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(state, geometry, options);
  if (!radiation.success ||
      !dec3d::radiation::ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(radiation)) {
    return FailStage('R',
                     radiation.failure_reason.empty() ? "radiation coupling failed"
                                                      : radiation.failure_reason,
                     radiation.failure_diagnostics.empty() ? radiation.report_line
                                                           : radiation.failure_diagnostics);
  }
  RuntimeStageResult result;
  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=R"
         << "; stage_backend=" << radiation.backend_executed
         << "; radiation_stage_report_present=true"
         << "; radiation_provider_report_present=true"
         << "; opacity_provider=tops_dt_tabulated"
         << "; radiation_groups_epoch=post_H_committed"
         << "; electron_thermal_state_epoch=post_E_committed"
         << "; radiation_boundary_model=" << config.radiation.boundary_model
         << "; radiation_flux_limiter_model=" << config.radiation.flux_limiter
         << "; updated_fields_mask=" << radiation.updated_fields;
  result.report_line = report.str();
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteAlphaStage(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::io::InputDeckConfig& config,
    double dt_s) {
  if (state.rho.size() > kSerialImplicitRuntimeCellLimit) {
    return FailStage('A', "serial A runtime exceeds implicit cell limit; distributed P5 runtime required");
  }
  dec3d::alpha::OneGroupAlphaTransportOptions options;
  options.dt_s = dt_s;
  options.provider_options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.provider_options.tau_model = dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.provider_options.reactivity_model =
      config.alpha.reactivity_model == "bosch_hale_dt"
          ? dec3d::alpha::AlphaReactivityModel::bosch_hale_dt
          : dec3d::alpha::AlphaReactivityModel::missing;
  options.boundary_policy = GenericBoundaryPolicy(config);
  options.serial_reference_options.row_limit = kSerialImplicitRuntimeCellLimit;
  const auto alpha = dec3d::alpha::ApplyOneGroupAlphaTransport(state, geometry, options);
  if (!alpha.success || !dec3d::alpha::ValidateAlphaOneGroupOperatorDiagnostics(alpha)) {
    return FailStage('A',
                     alpha.failure_reason.empty() ? "alpha transport failed"
                                                  : alpha.failure_reason,
                     alpha.failure_diagnostics.empty() ? alpha.report_line
                                                       : alpha.failure_diagnostics);
  }
  RuntimeStageResult result;
  result.success = true;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=A"
         << "; stage_backend=" << alpha.backend_executed
         << "; alpha_stage_report_present=true"
         << "; alpha_transport_model=atzeni_one_group"
         << "; reactivity_model_executed=" << alpha.reactivity_model_executed
         << "; fuel_depletion_enabled=false"
         << "; updated_fields_mask=" << alpha.updated_fields;
  result.report_line = report.str();
  return result;
}

#ifdef DEC3D_ENABLE_HYPRE
[[nodiscard]] RuntimeStageResult ExecuteDistributedThermalStage(
    dec3d::state::CanonicalState& local_state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::io::InputDeckConfig& config,
    double dt_s) {
  dec3d::transport::DistributedThermalConductionProblem problem;
  problem.ownership = ownership;
  problem.local_state = local_state;
  problem.global_geometry = geometry;
  problem.boundary_policy = GenericBoundaryPolicy(config);

  dec3d::transport::DistributedThermalConductionOptions options;
  options.dt_s = dt_s;
  options.kappa_model = ThermalKappaModel(config.thermal.kappa_model);
  options.electron_flux_limiter_model =
      ElectronFluxLimiterModelFromDeck(config.thermal.electron_flux_limiter);
  options.electron_flux_limiter_alpha_e =
      options.electron_flux_limiter_model ==
              dec3d::transport::ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa
          ? 1.0e-8
          : 0.0;
  const auto thermal =
      dec3d::transport::ApplyDistributedVariableKappaThermalConduction(problem, options);
  if (!thermal.success ||
      !dec3d::transport::ValidateDistributedThermalConductionDiagnostics(thermal)) {
    return FailStage('T',
                     thermal.failure_reason.empty() ? "distributed thermal conduction failed"
                                                    : thermal.failure_reason,
                     thermal.failure_diagnostics.empty() ? thermal.report_line
                                                         : thermal.failure_diagnostics);
  }
  local_state = problem.local_state;
  RuntimeStageResult result;
  result.success = true;
  result.coefficient_provider_wall_s = thermal.coefficient_provider_wall_s;
  result.assembly_wall_s = thermal.assembly_wall_s;
  result.hypre_setup_wall_s = thermal.hypre_setup_wall_s;
  result.hypre_solve_wall_s = thermal.hypre_solve_wall_s;
  result.writeback_wall_s = thermal.writeback_wall_s;
  result.solver_iterations = thermal.solver_iterations;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=T"
         << "; stage_backend=hypre_parcsr_gmres_boomeramg"
         << "; thermal_stage_report_present=true"
         << "; distributed_writeback=" << (thermal.distributed_writeback ? "true" : "false")
         << "; gathered_writeback_used=" << (thermal.gathered_writeback_used ? "true" : "false")
         << "; updated_fields_mask=" << thermal.updated_fields
         << "; nested_thermal_report={" << thermal.report_line << "}";
  result.report_line = report.str();
  AttachDistributedMatrixDiagnostics(result,
                                     thermal.electron_matrix_report + ";" +
                                         thermal.ion_matrix_report);
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) &&
      !result.global_matrix_diagnostics_present) {
    return FailStage('T', "distributed implicit matrix diagnostics missing", result.report_line);
  }
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteDistributedRadiationStage(
    dec3d::state::CanonicalState& local_state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const dec3d::radiation::TopsOpacityTable& opacity_table,
    const dec3d::io::InputDeckConfig& config,
    const dec3d::radiation::TopsOpacityProviderOptions& provider_options,
    dec3d::transport::DistributedLaggedBoomerAmgCache* radiation_amg_cache,
    double dt_s) {
  dec3d::radiation::DistributedMultigroupRadiationProblem problem;
  problem.ownership = ownership;
  problem.local_state = local_state;
  problem.global_geometry = geometry;
  problem.boundary_policy = GenericBoundaryPolicy(config);
  problem.group_layout = group_layout;
  problem.opacity_table = &opacity_table;
  problem.dt_s = dt_s;

  dec3d::radiation::DistributedMultigroupRadiationOptions options;
  options.provider_options = provider_options;
  options.boundary_model =
      config.radiation.boundary_model == "thesis_marshak_vacuum"
          ? dec3d::radiation::DistributedRadiationBoundaryModel::thesis_marshak_vacuum
          : dec3d::radiation::DistributedRadiationBoundaryModel::
                contract_zero_flux_or_scalar_remap;
  options.radiation_flux_limiter.model =
      RadiationFluxLimiterModelFromDeck(config.radiation.flux_limiter);
  options.solve_options.relative_tolerance = 1.0e-8;
  options.solve_options.max_iterations = 100;
  options.lagged_amg_cache = radiation_amg_cache;
  options.lagged_amg_enabled = radiation_amg_cache != nullptr;
  options.lagged_amg_rebuild_every = 4;
  options.lagged_amg_max_matrix_relative_change = 0.1;
  options.lagged_amg_max_iteration_growth = 1.5;
  const auto radiation =
      dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(problem, options);
  if (!radiation.success ||
      !dec3d::radiation::ValidateDistributedMultigroupRadiationDiagnostics(radiation)) {
    return FailStage('R',
                     radiation.failure_reason.empty() ? "distributed radiation failed"
                                                      : radiation.failure_reason,
                     radiation.failure_diagnostics.empty() ? radiation.report_line
                                                           : radiation.failure_diagnostics);
  }
  local_state = problem.local_state;
  RuntimeStageResult result;
  result.success = true;
  result.coefficient_provider_wall_s = radiation.coefficient_provider_wall_s;
  result.flux_limiter_wall_s = radiation.flux_limiter_wall_s;
  result.assembly_wall_s = radiation.assembly_wall_s;
  result.hypre_setup_wall_s = radiation.hypre_setup_wall_s;
  result.hypre_solve_wall_s = radiation.hypre_solve_wall_s;
  result.writeback_wall_s = radiation.writeback_wall_s;
  result.solver_iterations = radiation.solver_iterations;
  result.lagged_amg_candidate_count = radiation.lagged_amg_candidate_count;
  result.lagged_amg_reuse_attempted_count =
      radiation.lagged_amg_reuse_attempted_count;
  result.lagged_amg_reuse_accepted_count =
      radiation.lagged_amg_reuse_accepted_count;
  result.lagged_amg_rebuild_count = radiation.lagged_amg_rebuild_count;
  result.lagged_amg_fallback_rebuild_count =
      radiation.lagged_amg_fallback_rebuild_count;
  result.max_lagged_amg_global_matrix_rel_change =
      radiation.max_lagged_amg_global_matrix_rel_change;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=R"
         << "; stage_backend=hypre_parcsr_gmres_boomeramg"
         << "; radiation_stage_report_present=true"
         << "; distributed_per_group=true"
         << "; rank0_gather_solve_used=false"
         << "; radiation_groups_epoch=post_H_committed"
         << "; electron_thermal_state_epoch=post_E_committed"
         << "; radiation_boundary_model=" << config.radiation.boundary_model
         << "; radiation_flux_limiter_model=" << config.radiation.flux_limiter
         << "; updated_fields_mask=" << radiation.updated_fields
         << "; nested_radiation_report={" << radiation.report_line << "}";
  result.report_line = report.str();
  std::string radiation_matrix_reports;
  for (const auto& matrix_report : radiation.per_group_matrix_reports) {
    radiation_matrix_reports += matrix_report;
    radiation_matrix_reports += ';';
  }
  AttachDistributedMatrixDiagnostics(result, radiation_matrix_reports);
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) &&
      !result.global_matrix_diagnostics_present) {
    return FailStage('R', "distributed implicit matrix diagnostics missing", result.report_line);
  }
  return result;
}

[[nodiscard]] RuntimeStageResult ExecuteDistributedAlphaStage(
    dec3d::state::CanonicalState& local_state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::io::InputDeckConfig& config,
    double dt_s) {
  dec3d::alpha::DistributedOneGroupAlphaTransportProblem problem;
  problem.communicator = ownership.communicator;
  problem.ownership = ownership;
  problem.local_state = &local_state;
  problem.global_geometry = geometry;
  problem.dt_s = dt_s;

  dec3d::alpha::DistributedOneGroupAlphaTransportOptions options;
  options.provider_options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.provider_options.tau_model = dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.provider_options.reactivity_model =
      config.alpha.reactivity_model == "bosch_hale_dt"
          ? dec3d::alpha::AlphaReactivityModel::bosch_hale_dt
          : dec3d::alpha::AlphaReactivityModel::missing;
  options.boundary_policy = GenericBoundaryPolicy(config);
  options.solve_options.relative_tolerance = 1.0e-8;
  options.solve_options.max_iterations = 100;
  const auto alpha = dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(problem, options);
  if (!alpha.success || !dec3d::alpha::ValidateDistributedOneGroupAlphaDiagnostics(alpha)) {
    return FailStage('A',
                     alpha.failure_reason.empty() ? "distributed alpha transport failed"
                                                  : alpha.failure_reason,
                     alpha.failure_diagnostics.empty() ? alpha.report_line
                                                       : alpha.failure_diagnostics);
  }
  RuntimeStageResult result;
  result.success = true;
  result.coefficient_provider_wall_s = alpha.coefficient_provider_wall_s;
  result.assembly_wall_s = alpha.assembly_wall_s;
  result.hypre_setup_wall_s = alpha.hypre_setup_wall_s;
  result.hypre_solve_wall_s = alpha.hypre_solve_wall_s;
  result.writeback_wall_s = alpha.writeback_wall_s;
  result.solver_iterations = alpha.solver_iterations;
  std::ostringstream report;
  report << "diagnostic_id=p5.runtime.stage"
         << "; stage_id=A"
         << "; stage_backend=hypre_parcsr_gmres_boomeramg"
         << "; alpha_stage_report_present=true"
         << "; alpha_transport_model=atzeni_one_group"
         << "; fuel_depletion_enabled=false"
         << "; rank0_gather_solve_used=false"
         << "; serial_dense_fallback_used=false"
         << "; updated_fields_mask=" << alpha.updated_fields
         << "; nested_alpha_report={" << alpha.report_line << "}";
  result.report_line = report.str();
  AttachDistributedMatrixDiagnostics(result, alpha.assembly_report);
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) &&
      !result.global_matrix_diagnostics_present) {
    return FailStage('A', "distributed implicit matrix diagnostics missing", result.report_line);
  }
  return result;
}
#endif

struct RuntimeLoopResult {
  bool success{false};
  std::size_t steps_executed{0u};
  double final_time_s{0.0};
  double dt_estimator_wall_s{0.0};
  double h_stage_wall_s{0.0};
  double h_view_wall_s{0.0};
  double h_halo_exchange_wall_s{0.0};
  double h_ale_proposal_wall_s{0.0};
  double h_stage_state_wall_s{0.0};
  double h_hydro_advance_wall_s{0.0};
  double h_hydro_snapshot_wall_s{0.0};
  double h_hydro_scratch_wall_s{0.0};
  double h_hydro_radial_sweep_wall_s{0.0};
  double h_hydro_macro_detect_wall_s{0.0};
  double h_hydro_macro_restrict_wall_s{0.0};
  double h_hydro_macro_update_wall_s{0.0};
  double h_hydro_macro_radial_update_wall_s{0.0};
  double h_hydro_macro_theta_update_wall_s{0.0};
  double h_hydro_macro_phi_update_wall_s{0.0};
  double h_hydro_macro_state_update_wall_s{0.0};
  double h_hydro_macro_prolong_wall_s{0.0};
  double h_hydro_theta_sweep_wall_s{0.0};
  double h_hydro_phi_sweep_wall_s{0.0};
  double h_hydro_commit_wall_s{0.0};
  double h_hydro_source_wall_s{0.0};
  double h_hydro_budget_wall_s{0.0};
  double h_hydro_diagnostics_wall_s{0.0};
  double h_writeback_wall_s{0.0};
  double h_global_gate_wall_s{0.0};
  double t_stage_wall_s{0.0};
  double e_stage_wall_s{0.0};
  double r_stage_wall_s{0.0};
  double a_stage_wall_s{0.0};
  double runtime_output_wall_s{0.0};
  double t_coefficient_provider_wall_s{0.0};
  double t_assembly_wall_s{0.0};
  double t_hypre_setup_wall_s{0.0};
  double t_hypre_solve_wall_s{0.0};
  double t_writeback_wall_s{0.0};
  int t_solver_iterations{0};
  double r_coefficient_provider_wall_s{0.0};
  double r_flux_limiter_wall_s{0.0};
  double r_assembly_wall_s{0.0};
  double r_hypre_setup_wall_s{0.0};
  double r_hypre_solve_wall_s{0.0};
  double r_writeback_wall_s{0.0};
  int r_solver_iterations{0};
  std::size_t r_lagged_amg_candidate_count{0u};
  std::size_t r_lagged_amg_reuse_attempted_count{0u};
  std::size_t r_lagged_amg_reuse_accepted_count{0u};
  std::size_t r_lagged_amg_rebuild_count{0u};
  std::size_t r_lagged_amg_fallback_rebuild_count{0u};
  double r_max_lagged_amg_global_matrix_rel_change{0.0};
  double a_coefficient_provider_wall_s{0.0};
  double a_assembly_wall_s{0.0};
  double a_hypre_setup_wall_s{0.0};
  double a_hypre_solve_wall_s{0.0};
  double a_writeback_wall_s{0.0};
  int a_solver_iterations{0};
  std::size_t global_phi_coupling_count{0u};
  std::size_t global_duplicate_column_row_count{0u};
  bool global_matrix_diagnostics_present{false};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct AxisymmetricInvariantSummary {
  double max_abs_mom_phi{0.0};
  double max_abs_v_phi{0.0};
  double max_abs_mom_total{0.0};
  double max_abs_velocity{0.0};
  double mom_phi_tol{1.0e-30};
  double v_phi_tol{1.0e-30};
  bool ok{true};
};

void FinalizeAxisymmetricInvariantSummary(AxisymmetricInvariantSummary& summary) noexcept {
  summary.mom_phi_tol = std::max(1.0e-30, 1.0e-14 * summary.max_abs_mom_total);
  summary.v_phi_tol = std::max(1.0e-30, 1.0e-14 * summary.max_abs_velocity);
  summary.ok = summary.max_abs_mom_phi <= summary.mom_phi_tol &&
               summary.max_abs_v_phi <= summary.v_phi_tol;
}

[[nodiscard]] AxisymmetricInvariantSummary EvaluateAxisymmetricInvariants(
    const dec3d::state::CanonicalState& state) noexcept {
  AxisymmetricInvariantSummary summary;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double rho = state.rho(r, t, p);
        const double mom_r = state.mom_r(r, t, p);
        const double mom_theta = state.mom_theta(r, t, p);
        const double mom_phi = state.mom_phi(r, t, p);
        const double mom_total =
            std::sqrt(mom_r * mom_r + mom_theta * mom_theta + mom_phi * mom_phi);
        summary.max_abs_mom_phi = std::max(summary.max_abs_mom_phi, std::abs(mom_phi));
        summary.max_abs_mom_total = std::max(summary.max_abs_mom_total, mom_total);
        if (rho > 0.0 && std::isfinite(rho)) {
          const double v_r = mom_r / rho;
          const double v_theta = mom_theta / rho;
          const double v_phi = mom_phi / rho;
          const double velocity =
              std::sqrt(v_r * v_r + v_theta * v_theta + v_phi * v_phi);
          summary.max_abs_v_phi = std::max(summary.max_abs_v_phi, std::abs(v_phi));
          summary.max_abs_velocity = std::max(summary.max_abs_velocity, velocity);
        }
      }
    }
  }
  FinalizeAxisymmetricInvariantSummary(summary);
  return summary;
}

#ifdef DEC3D_ENABLE_HYPRE
[[nodiscard]] AxisymmetricInvariantSummary ReduceAxisymmetricInvariants(
    MPI_Comm comm,
    const AxisymmetricInvariantSummary& local) noexcept {
  double local_values[4] = {
      local.max_abs_mom_phi,
      local.max_abs_v_phi,
      local.max_abs_mom_total,
      local.max_abs_velocity,
  };
  double global_values[4] = {};
  MPI_Allreduce(local_values, global_values, 4, MPI_DOUBLE, MPI_MAX, comm);
  AxisymmetricInvariantSummary summary;
  summary.max_abs_mom_phi = global_values[0];
  summary.max_abs_v_phi = global_values[1];
  summary.max_abs_mom_total = global_values[2];
  summary.max_abs_velocity = global_values[3];
  FinalizeAxisymmetricInvariantSummary(summary);
  return summary;
}
#endif

void AppendAxisymmetricInvariantDiagnostics(
    std::ostringstream& report,
    const dec3d::io::InputDeckConfig& config,
    const AxisymmetricInvariantSummary& axisym) {
  report << "; max_abs_mom_phi=" << axisym.max_abs_mom_phi
         << "; max_abs_v_phi=" << axisym.max_abs_v_phi
         << "; axisymmetric_mom_phi_tol=" << axisym.mom_phi_tol
         << "; axisymmetric_v_phi_tol=" << axisym.v_phi_tol
         << "; axisymmetric_invariant_ok="
         << (!dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) || axisym.ok ? "true"
                                                                                   : "false");
}

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] bool TargetTimeEnabled(const dec3d::io::InputDeckConfig& config) noexcept {
  return config.run.target_time_s > 0.0;
}

[[nodiscard]] bool TargetTimeReached(const dec3d::io::InputDeckConfig& config,
                                     const double time_s) noexcept {
  if (!TargetTimeEnabled(config)) {
    return false;
  }
  const double tolerance = std::max(1.0e-30, std::abs(config.run.target_time_s) * 1.0e-12);
  return time_s + tolerance >= config.run.target_time_s;
}

[[nodiscard]] double ClampDtToTargetTime(const dec3d::io::InputDeckConfig& config,
                                         const double time_s,
                                         const double dt_s) noexcept {
  if (!TargetTimeEnabled(config)) {
    return dt_s;
  }
  const double remaining_s = config.run.target_time_s - time_s;
  if (!(remaining_s > 0.0)) {
    return 0.0;
  }
  return std::min(dt_s, remaining_s);
}

[[nodiscard]] bool ParseRuntimeDoubleField(
    std::string_view report_line,
    std::string_view key,
    double* value_out) {
  if (value_out == nullptr) {
    return false;
  }
  const std::string needle = std::string(key) + "=";
  const auto key_position = report_line.find(needle);
  if (key_position == std::string_view::npos) {
    return false;
  }
  const auto value_begin = key_position + needle.size();
  auto value_end = report_line.find(';', value_begin);
  if (value_end == std::string_view::npos) {
    value_end = report_line.size();
  }
  std::istringstream in(std::string(report_line.substr(value_begin, value_end - value_begin)));
  double parsed = 0.0;
  in >> parsed;
  if (!in) {
    return false;
  }
  *value_out = parsed;
  return true;
}

void AddHydroAdvanceTimingFromReport(RuntimeStageResult& result,
                                     std::string_view report_line) {
  const auto add_field = [&](std::string_view key, double& target) {
    double value = 0.0;
    if (ParseRuntimeDoubleField(report_line, key, &value)) {
      target += value;
    }
  };
  add_field("h_hydro_snapshot_wall_s", result.hydro_snapshot_wall_s);
  add_field("h_hydro_scratch_wall_s", result.hydro_scratch_wall_s);
  add_field("h_hydro_radial_sweep_wall_s", result.hydro_radial_sweep_wall_s);
  add_field("h_hydro_macro_detect_wall_s", result.hydro_macro_detect_wall_s);
  add_field("h_hydro_macro_restrict_wall_s", result.hydro_macro_restrict_wall_s);
  add_field("h_hydro_macro_update_wall_s", result.hydro_macro_update_wall_s);
  add_field("h_hydro_macro_radial_update_wall_s",
            result.hydro_macro_radial_update_wall_s);
  add_field("h_hydro_macro_theta_update_wall_s",
            result.hydro_macro_theta_update_wall_s);
  add_field("h_hydro_macro_phi_update_wall_s",
            result.hydro_macro_phi_update_wall_s);
  add_field("h_hydro_macro_state_update_wall_s",
            result.hydro_macro_state_update_wall_s);
  add_field("h_hydro_macro_prolong_wall_s", result.hydro_macro_prolong_wall_s);
  add_field("h_hydro_theta_sweep_wall_s", result.hydro_theta_sweep_wall_s);
  add_field("h_hydro_phi_sweep_wall_s", result.hydro_phi_sweep_wall_s);
  add_field("h_hydro_commit_wall_s", result.hydro_commit_wall_s);
  add_field("h_hydro_source_wall_s", result.hydro_source_wall_s);
  add_field("h_hydro_budget_wall_s", result.hydro_budget_wall_s);
  add_field("h_hydro_diagnostics_wall_s", result.hydro_diagnostics_wall_s);
}

void AddHydroAdvanceTimingFromDiagnostics(
    RuntimeStageResult& result,
    const dec3d::core::DiagnosticsPayload& diagnostics) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == "p1.hydro.static_grid.timing") {
      AddHydroAdvanceTimingFromReport(result, entry.message);
    }
  }
}

void AddHydroAdvanceTimingFromStaticGrid(
    RuntimeStageResult& result,
    const dec3d::hydro::StaticGridHydroResult& hydro) {
  result.hydro_snapshot_wall_s += hydro.timing_snapshot_wall_s;
  result.hydro_scratch_wall_s += hydro.timing_scratch_wall_s;
  result.hydro_radial_sweep_wall_s += hydro.timing_radial_sweep_wall_s;
  result.hydro_macro_detect_wall_s += hydro.timing_macro_detect_wall_s;
  result.hydro_macro_restrict_wall_s += hydro.timing_macro_restrict_wall_s;
  result.hydro_macro_update_wall_s += hydro.timing_macro_update_wall_s;
  result.hydro_macro_radial_update_wall_s += hydro.timing_macro_radial_update_wall_s;
  result.hydro_macro_theta_update_wall_s += hydro.timing_macro_theta_update_wall_s;
  result.hydro_macro_phi_update_wall_s += hydro.timing_macro_phi_update_wall_s;
  result.hydro_macro_state_update_wall_s += hydro.timing_macro_state_update_wall_s;
  result.hydro_macro_prolong_wall_s += hydro.timing_macro_prolong_wall_s;
  result.hydro_theta_sweep_wall_s += hydro.timing_theta_sweep_wall_s;
  result.hydro_phi_sweep_wall_s += hydro.timing_phi_sweep_wall_s;
  result.hydro_commit_wall_s += hydro.timing_commit_wall_s;
  result.hydro_source_wall_s += hydro.timing_source_wall_s;
  result.hydro_budget_wall_s += hydro.timing_budget_wall_s;
  result.hydro_diagnostics_wall_s += hydro.timing_diagnostics_wall_s;
}

void AppendHydroAdvanceTimingFields(std::ostringstream& report,
                                    const RuntimeStageResult& result) {
  report << "; h_hydro_snapshot_wall_s=" << result.hydro_snapshot_wall_s
         << "; h_hydro_scratch_wall_s=" << result.hydro_scratch_wall_s
         << "; h_hydro_radial_sweep_wall_s=" << result.hydro_radial_sweep_wall_s
         << "; h_hydro_macro_detect_wall_s=" << result.hydro_macro_detect_wall_s
         << "; h_hydro_macro_restrict_wall_s=" << result.hydro_macro_restrict_wall_s
         << "; h_hydro_macro_update_wall_s=" << result.hydro_macro_update_wall_s
         << "; h_hydro_macro_radial_update_wall_s="
         << result.hydro_macro_radial_update_wall_s
         << "; h_hydro_macro_theta_update_wall_s="
         << result.hydro_macro_theta_update_wall_s
         << "; h_hydro_macro_phi_update_wall_s="
         << result.hydro_macro_phi_update_wall_s
         << "; h_hydro_macro_state_update_wall_s="
         << result.hydro_macro_state_update_wall_s
         << "; h_hydro_macro_prolong_wall_s=" << result.hydro_macro_prolong_wall_s
         << "; h_hydro_theta_sweep_wall_s=" << result.hydro_theta_sweep_wall_s
         << "; h_hydro_phi_sweep_wall_s=" << result.hydro_phi_sweep_wall_s
         << "; h_hydro_commit_wall_s=" << result.hydro_commit_wall_s
         << "; h_hydro_source_wall_s=" << result.hydro_source_wall_s
         << "; h_hydro_budget_wall_s=" << result.hydro_budget_wall_s
         << "; h_hydro_diagnostics_wall_s=" << result.hydro_diagnostics_wall_s;
}

void AddRuntimeStageWallTime(RuntimeLoopResult& loop, const char stage, const double seconds) {
  switch (stage) {
    case 'H':
      loop.h_stage_wall_s += seconds;
      break;
    case 'T':
      loop.t_stage_wall_s += seconds;
      break;
    case 'E':
      loop.e_stage_wall_s += seconds;
      break;
    case 'R':
      loop.r_stage_wall_s += seconds;
      break;
    case 'A':
      loop.a_stage_wall_s += seconds;
      break;
    default:
      break;
  }
}

void AddRuntimeStagePerformance(RuntimeLoopResult& loop,
                                const char stage,
                                const RuntimeStageResult& result) {
  switch (stage) {
    case 'H':
      loop.h_view_wall_s += result.hydro_view_wall_s;
      loop.h_halo_exchange_wall_s += result.hydro_halo_exchange_wall_s;
      loop.h_ale_proposal_wall_s += result.hydro_ale_proposal_wall_s;
      loop.h_stage_state_wall_s += result.hydro_stage_state_wall_s;
      loop.h_hydro_advance_wall_s += result.hydro_advance_wall_s;
      loop.h_hydro_snapshot_wall_s += result.hydro_snapshot_wall_s;
      loop.h_hydro_scratch_wall_s += result.hydro_scratch_wall_s;
      loop.h_hydro_radial_sweep_wall_s += result.hydro_radial_sweep_wall_s;
      loop.h_hydro_macro_detect_wall_s += result.hydro_macro_detect_wall_s;
      loop.h_hydro_macro_restrict_wall_s += result.hydro_macro_restrict_wall_s;
      loop.h_hydro_macro_update_wall_s += result.hydro_macro_update_wall_s;
      loop.h_hydro_macro_radial_update_wall_s += result.hydro_macro_radial_update_wall_s;
      loop.h_hydro_macro_theta_update_wall_s += result.hydro_macro_theta_update_wall_s;
      loop.h_hydro_macro_phi_update_wall_s += result.hydro_macro_phi_update_wall_s;
      loop.h_hydro_macro_state_update_wall_s += result.hydro_macro_state_update_wall_s;
      loop.h_hydro_macro_prolong_wall_s += result.hydro_macro_prolong_wall_s;
      loop.h_hydro_theta_sweep_wall_s += result.hydro_theta_sweep_wall_s;
      loop.h_hydro_phi_sweep_wall_s += result.hydro_phi_sweep_wall_s;
      loop.h_hydro_commit_wall_s += result.hydro_commit_wall_s;
      loop.h_hydro_source_wall_s += result.hydro_source_wall_s;
      loop.h_hydro_budget_wall_s += result.hydro_budget_wall_s;
      loop.h_hydro_diagnostics_wall_s += result.hydro_diagnostics_wall_s;
      loop.h_writeback_wall_s += result.hydro_writeback_wall_s;
      loop.h_global_gate_wall_s += result.hydro_global_gate_wall_s;
      break;
    case 'T':
      loop.t_coefficient_provider_wall_s += result.coefficient_provider_wall_s;
      loop.t_assembly_wall_s += result.assembly_wall_s;
      loop.t_hypre_setup_wall_s += result.hypre_setup_wall_s;
      loop.t_hypre_solve_wall_s += result.hypre_solve_wall_s;
      loop.t_writeback_wall_s += result.writeback_wall_s;
      loop.t_solver_iterations += result.solver_iterations;
      break;
    case 'R':
      loop.r_coefficient_provider_wall_s += result.coefficient_provider_wall_s;
      loop.r_flux_limiter_wall_s += result.flux_limiter_wall_s;
      loop.r_assembly_wall_s += result.assembly_wall_s;
      loop.r_hypre_setup_wall_s += result.hypre_setup_wall_s;
      loop.r_hypre_solve_wall_s += result.hypre_solve_wall_s;
      loop.r_writeback_wall_s += result.writeback_wall_s;
      loop.r_solver_iterations += result.solver_iterations;
      loop.r_lagged_amg_candidate_count += result.lagged_amg_candidate_count;
      loop.r_lagged_amg_reuse_attempted_count +=
          result.lagged_amg_reuse_attempted_count;
      loop.r_lagged_amg_reuse_accepted_count +=
          result.lagged_amg_reuse_accepted_count;
      loop.r_lagged_amg_rebuild_count += result.lagged_amg_rebuild_count;
      loop.r_lagged_amg_fallback_rebuild_count +=
          result.lagged_amg_fallback_rebuild_count;
      loop.r_max_lagged_amg_global_matrix_rel_change =
          std::max(loop.r_max_lagged_amg_global_matrix_rel_change,
                   result.max_lagged_amg_global_matrix_rel_change);
      break;
    case 'A':
      loop.a_coefficient_provider_wall_s += result.coefficient_provider_wall_s;
      loop.a_assembly_wall_s += result.assembly_wall_s;
      loop.a_hypre_setup_wall_s += result.hypre_setup_wall_s;
      loop.a_hypre_solve_wall_s += result.hypre_solve_wall_s;
      loop.a_writeback_wall_s += result.writeback_wall_s;
      loop.a_solver_iterations += result.solver_iterations;
      break;
    default:
      break;
  }
  if (result.global_matrix_diagnostics_present) {
    loop.global_phi_coupling_count =
        std::max(loop.global_phi_coupling_count, result.global_phi_coupling_count);
    loop.global_duplicate_column_row_count =
        std::max(loop.global_duplicate_column_row_count,
                 result.global_duplicate_column_row_count);
    loop.global_matrix_diagnostics_present = true;
  }
}

void AppendRuntimeTimingDiagnostics(std::ostringstream& report,
                                    const RuntimeLoopResult& loop) {
  report << "; dt_estimator_wall_s=" << loop.dt_estimator_wall_s
         << "; h_stage_wall_s=" << loop.h_stage_wall_s
         << "; h_view_wall_s=" << loop.h_view_wall_s
         << "; h_halo_exchange_wall_s=" << loop.h_halo_exchange_wall_s
         << "; h_ale_proposal_wall_s=" << loop.h_ale_proposal_wall_s
         << "; h_stage_state_wall_s=" << loop.h_stage_state_wall_s
         << "; h_hydro_advance_wall_s=" << loop.h_hydro_advance_wall_s
         << "; h_hydro_snapshot_wall_s=" << loop.h_hydro_snapshot_wall_s
         << "; h_hydro_scratch_wall_s=" << loop.h_hydro_scratch_wall_s
         << "; h_hydro_radial_sweep_wall_s=" << loop.h_hydro_radial_sweep_wall_s
         << "; h_hydro_macro_detect_wall_s=" << loop.h_hydro_macro_detect_wall_s
         << "; h_hydro_macro_restrict_wall_s=" << loop.h_hydro_macro_restrict_wall_s
         << "; h_hydro_macro_update_wall_s=" << loop.h_hydro_macro_update_wall_s
         << "; h_hydro_macro_radial_update_wall_s="
         << loop.h_hydro_macro_radial_update_wall_s
         << "; h_hydro_macro_theta_update_wall_s="
         << loop.h_hydro_macro_theta_update_wall_s
         << "; h_hydro_macro_phi_update_wall_s="
         << loop.h_hydro_macro_phi_update_wall_s
         << "; h_hydro_macro_state_update_wall_s="
         << loop.h_hydro_macro_state_update_wall_s
         << "; h_hydro_macro_prolong_wall_s=" << loop.h_hydro_macro_prolong_wall_s
         << "; h_hydro_theta_sweep_wall_s=" << loop.h_hydro_theta_sweep_wall_s
         << "; h_hydro_phi_sweep_wall_s=" << loop.h_hydro_phi_sweep_wall_s
         << "; h_hydro_commit_wall_s=" << loop.h_hydro_commit_wall_s
         << "; h_hydro_source_wall_s=" << loop.h_hydro_source_wall_s
         << "; h_hydro_budget_wall_s=" << loop.h_hydro_budget_wall_s
         << "; h_hydro_diagnostics_wall_s=" << loop.h_hydro_diagnostics_wall_s
         << "; h_writeback_wall_s=" << loop.h_writeback_wall_s
         << "; h_global_gate_wall_s=" << loop.h_global_gate_wall_s
         << "; t_stage_wall_s=" << loop.t_stage_wall_s
         << "; e_stage_wall_s=" << loop.e_stage_wall_s
         << "; r_stage_wall_s=" << loop.r_stage_wall_s
         << "; a_stage_wall_s=" << loop.a_stage_wall_s
         << "; runtime_output_wall_s=" << loop.runtime_output_wall_s
         << "; t_coefficient_provider_wall_s=" << loop.t_coefficient_provider_wall_s
         << "; t_assembly_wall_s=" << loop.t_assembly_wall_s
         << "; t_hypre_setup_wall_s=" << loop.t_hypre_setup_wall_s
         << "; t_hypre_solve_wall_s=" << loop.t_hypre_solve_wall_s
         << "; t_writeback_wall_s=" << loop.t_writeback_wall_s
         << "; t_solver_iterations=" << loop.t_solver_iterations
         << "; r_coefficient_provider_wall_s=" << loop.r_coefficient_provider_wall_s
         << "; r_flux_limiter_wall_s=" << loop.r_flux_limiter_wall_s
         << "; r_assembly_wall_s=" << loop.r_assembly_wall_s
         << "; r_hypre_setup_wall_s=" << loop.r_hypre_setup_wall_s
         << "; r_hypre_solve_wall_s=" << loop.r_hypre_solve_wall_s
         << "; r_writeback_wall_s=" << loop.r_writeback_wall_s
         << "; r_solver_iterations=" << loop.r_solver_iterations
         << "; r_lagged_amg_candidate_count=" << loop.r_lagged_amg_candidate_count
         << "; r_lagged_amg_reuse_attempted_count="
         << loop.r_lagged_amg_reuse_attempted_count
         << "; r_lagged_amg_reuse_accepted_count="
         << loop.r_lagged_amg_reuse_accepted_count
         << "; r_lagged_amg_rebuild_count=" << loop.r_lagged_amg_rebuild_count
         << "; r_lagged_amg_fallback_rebuild_count="
         << loop.r_lagged_amg_fallback_rebuild_count
         << "; r_max_lagged_amg_global_matrix_rel_change="
         << loop.r_max_lagged_amg_global_matrix_rel_change
         << "; a_coefficient_provider_wall_s=" << loop.a_coefficient_provider_wall_s
         << "; a_assembly_wall_s=" << loop.a_assembly_wall_s
         << "; a_hypre_setup_wall_s=" << loop.a_hypre_setup_wall_s
         << "; a_hypre_solve_wall_s=" << loop.a_hypre_solve_wall_s
         << "; a_writeback_wall_s=" << loop.a_writeback_wall_s
         << "; a_solver_iterations=" << loop.a_solver_iterations;
}

struct RestartInitialState {
  bool success{false};
  dec3d::state::CanonicalState state;
  std::size_t step{0u};
  double time_s{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct RuntimeStartPoint {
  std::size_t step{0u};
  double time_s{0.0};
};

[[nodiscard]] std::vector<std::string> SplitCsvLine(const std::string& line) {
  std::vector<std::string> parts;
  std::string current;
  std::istringstream in(line);
  while (std::getline(in, current, ',')) {
    parts.push_back(current);
  }
  return parts;
}

[[nodiscard]] bool ParseRestartBool(const std::string& value) noexcept {
  return value == "true" || value == "1";
}

[[nodiscard]] RestartInitialState FailRestartLoad(
    std::string reason,
    const std::filesystem::path& restart_path) {
  RestartInitialState result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics =
      "diagnostic_id=p5.io.restart_load.failure; failure_reason=" +
      result.failure_reason + "; restart_file=" + restart_path.string();
  return result;
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask RestartAuthoritativeMask() noexcept {
  using dec3d::core::AuthoritativeField;
  return dec3d::core::ToMask(AuthoritativeField::rho) |
         dec3d::core::ToMask(AuthoritativeField::mom_r) |
         dec3d::core::ToMask(AuthoritativeField::mom_theta) |
         dec3d::core::ToMask(AuthoritativeField::mom_phi) |
         dec3d::core::ToMask(AuthoritativeField::e_fluid_total) |
         dec3d::core::ToMask(AuthoritativeField::e_electron) |
         dec3d::core::ToMask(AuthoritativeField::radiation_groups) |
         dec3d::core::ToMask(AuthoritativeField::alpha_state);
}

[[nodiscard]] RestartInitialState LoadRestartInitialState(
    const std::filesystem::path& restart_path,
    const dec3d::io::InputDeckConfig& config,
    const dec3d::radiation::RadiationGroupLayout& group_layout) {
  std::ifstream in(restart_path);
  if (!in) {
    return FailRestartLoad("missing restart file", restart_path);
  }

  bool saw_header = false;
  bool restart_compatible = false;
  std::size_t step = 0u;
  double time_s = 0.0;
  dec3d::state::CanonicalStateLayout layout;
  std::string line;
  while (std::getline(in, line)) {
    if (line.rfind("#", 0) != 0) {
      saw_header = line == "field,group,radial,theta,phi,value";
      break;
    }
    const auto eq = line.find('=');
    if (eq == std::string::npos) {
      continue;
    }
    const auto key_begin = line.find_first_not_of("# \t");
    if (key_begin == std::string::npos || key_begin >= eq) {
      continue;
    }
    const auto key = line.substr(key_begin, eq - key_begin);
    const auto value = line.substr(eq + 1u);
    try {
      if (key == "restart_compatible") {
        restart_compatible = ParseRestartBool(value);
      } else if (key == "step") {
        step = static_cast<std::size_t>(std::stoull(value));
      } else if (key == "time_s") {
        time_s = std::stod(value);
      } else if (key == "radial_cells") {
        layout.radial_cells = static_cast<std::size_t>(std::stoull(value));
      } else if (key == "theta_cells") {
        layout.theta_cells = static_cast<std::size_t>(std::stoull(value));
      } else if (key == "phi_cells") {
        layout.phi_cells = static_cast<std::size_t>(std::stoull(value));
      } else if (key == "radiation_group_count") {
        layout.radiation_group_count = static_cast<std::size_t>(std::stoull(value));
      }
    } catch (const std::exception&) {
      return FailRestartLoad("invalid restart metadata value", restart_path);
    }
  }

  if (!restart_compatible) {
    return FailRestartLoad("restart file is not restart compatible", restart_path);
  }
  if (!saw_header) {
    return FailRestartLoad("restart field header is missing", restart_path);
  }
  if (layout.radial_cells != config.mesh.radial_cells ||
      layout.theta_cells != config.mesh.theta_cells ||
      layout.phi_cells != config.mesh.phi_cells ||
      layout.radiation_group_count != group_layout.group_count) {
    return FailRestartLoad("restart layout does not match input deck", restart_path);
  }
  if (!(time_s >= 0.0) || !std::isfinite(time_s)) {
    return FailRestartLoad("restart time is not finite", restart_path);
  }

  auto state = dec3d::state::CanonicalState::Create(layout);
  const std::size_t cell_count =
      layout.radial_cells * layout.theta_cells * layout.phi_cells;
  std::size_t rho_count = 0u;
  std::size_t mom_r_count = 0u;
  std::size_t mom_theta_count = 0u;
  std::size_t mom_phi_count = 0u;
  std::size_t e_fluid_total_count = 0u;
  std::size_t e_electron_count = 0u;
  std::size_t radiation_count = 0u;
  std::size_t alpha_count = 0u;

  while (std::getline(in, line)) {
    if (line.empty() || line.rfind("#", 0) == 0) {
      continue;
    }
    const auto parts = SplitCsvLine(line);
    if (parts.size() != 6u) {
      return FailRestartLoad("invalid restart field row", restart_path);
    }
    try {
      const auto& field = parts[0];
      const int group = std::stoi(parts[1]);
      const std::size_t r = static_cast<std::size_t>(std::stoull(parts[2]));
      const std::size_t t = static_cast<std::size_t>(std::stoull(parts[3]));
      const std::size_t p = static_cast<std::size_t>(std::stoull(parts[4]));
      const double value = std::stod(parts[5]);
      if (r >= layout.radial_cells || t >= layout.theta_cells ||
          p >= layout.phi_cells || !std::isfinite(value)) {
        return FailRestartLoad("restart field row is out of range", restart_path);
      }
      if (field == "rho" && group == -1) {
        state.rho(r, t, p) = value;
        ++rho_count;
      } else if (field == "mom_r" && group == -1) {
        state.mom_r(r, t, p) = value;
        ++mom_r_count;
      } else if (field == "mom_theta" && group == -1) {
        state.mom_theta(r, t, p) = value;
        ++mom_theta_count;
      } else if (field == "mom_phi" && group == -1) {
        state.mom_phi(r, t, p) = value;
        ++mom_phi_count;
      } else if (field == "e_electron" && group == -1) {
        state.e_electron(r, t, p) = value;
        ++e_electron_count;
      } else if (field == "e_fluid_total" && group == -1) {
        state.e_fluid_total(r, t, p) = value;
        ++e_fluid_total_count;
      } else if (field == "radiation_groups" && group >= 0 &&
                 static_cast<std::size_t>(group) < state.radiation_groups.size()) {
        state.radiation_groups[static_cast<std::size_t>(group)](r, t, p) = value;
        ++radiation_count;
      } else if (field == "alpha_state" && group == -1) {
        state.alpha_state.storage(r, t, p) = value;
        ++alpha_count;
      }
    } catch (const std::exception&) {
      return FailRestartLoad("invalid restart field value", restart_path);
    }
  }

  if (rho_count != cell_count || mom_r_count != cell_count ||
      mom_theta_count != cell_count || mom_phi_count != cell_count ||
      e_electron_count != cell_count || e_fluid_total_count != cell_count ||
      alpha_count != cell_count ||
      radiation_count != cell_count * group_layout.group_count) {
    return FailRestartLoad("restart payload is incomplete", restart_path);
  }

  state.ApplyAuthoritativeWrite(RestartAuthoritativeMask());
  RestartInitialState result;
  result.success = true;
  result.state = std::move(state);
  result.step = step;
  result.time_s = time_s;
  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p5.io.restart_load"
         << "; restart_file=" << restart_path.string()
         << "; restart_step=" << step
         << "; restart_time_s=" << time_s
         << "; restart_layout_matches_deck=true"
         << "; restart_payload_complete=true"
         << "; restart_authoritative_state_restored=true";
  result.report_line = report.str();
  return result;
}

[[nodiscard]] RuntimeLoopResult ExecuteP5RuntimeLoop(
    const dec3d::io::InputDeckConfig& config,
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& input_deck_path,
    const std::filesystem::path& profile_path,
    const RuntimeStartPoint& start) {
  RuntimeLoopResult loop;
  if (config.run.step_count <= 0) {
    loop.success = true;
    loop.report_line = "diagnostic_id=p5.runtime.loop; runtime_time_loop_executed=false";
    return loop;
  }
  for (char stage : config.run.stage_order) {
    if (!StageEnabled(config, stage)) {
      loop.failure_reason = "stage_order includes disabled physics stage";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason +
          "; disabled_stage=" + std::string(1, stage);
      return loop;
    }
  }

  const auto stage_order = JoinStageOrder(config.run.stage_order);
  dec3d::radiation::TopsOpacityTable serial_opacity_table;
  dec3d::radiation::TopsOpacityProviderOptions serial_provider_options;
  bool serial_opacity_table_loaded = false;
  if (config.physics.enable_radiation) {
    if (config.radiation.opacity_provider != "tops_dt_tabulated") {
      loop.failure_reason = "unsupported runtime opacity_provider";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
      return loop;
    }
    serial_provider_options = MakeRuntimeTopsProviderOptions(input_deck_path);
    if (serial_provider_options.table_root.empty()) {
      loop.failure_reason = "TOPS opacity table root not found";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
      return loop;
    }
    const auto table_result = dec3d::radiation::LoadTopsOpacityTable(serial_provider_options);
    if (!table_result.success || table_result.report_line.empty()) {
      loop.failure_reason =
          table_result.failure_reason.empty() ? "TOPS opacity table load failed"
                                              : table_result.failure_reason;
      loop.failure_diagnostics = table_result.failure_diagnostics.empty()
                                     ? table_result.report_line
                                     : table_result.failure_diagnostics;
      return loop;
    }
    serial_opacity_table = table_result.table;
    serial_opacity_table_loaded = true;
  }
  double time_s = start.time_s;
  dec3d::io::RuntimeOutputStepState output_state;
  output_state.last_field_checkpoint_time_s = start.time_s;
  output_state.last_restart_checkpoint_time_s = start.time_s;
  output_state.last_history_profile_time_s = start.time_s;
  const auto wall_start = std::chrono::steady_clock::now();
  for (int local_step_index = 1; local_step_index <= config.run.step_count; ++local_step_index) {
    const std::size_t step_index = start.step + static_cast<std::size_t>(local_step_index);
    if (TargetTimeReached(config, time_s)) {
      break;
    }
    const auto dt_timer = std::chrono::steady_clock::now();
    dec3d::hydro::HydroOperator hydro_for_dt;
    auto hydro_options = RuntimeHydroOptions(config);
    hydro_for_dt.SetStaticGridOptions(hydro_options);
    dec3d::core::StageContext dt_context;
    dt_context.time_s = time_s;
    dt_context.dt_s = 0.0;
    dt_context.step = static_cast<std::uint64_t>(step_index);
    dt_context.phase_id = dec3d::core::PhaseId::p1;
    dt_context.contract_version = "p5.runtime.serial_callable.v1";
    dt_context.mesh_snapshot_handle = "p5.runtime.mesh.current";
    dt_context.ownership_handle = "single_rank_full_domain";
    dt_context.diagnostics_sink_handle = "p5.runtime.dec3d_out";
    if (!hydro_for_dt.bind(dt_context, geometry, state)) {
      loop.failure_reason = "hydro dt estimator bind failed";
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason;
      return loop;
    }
    const auto dt_advice = hydro_for_dt.estimate_dt();
    if (!dt_advice.is_complete()) {
      loop.failure_reason =
          dt_advice.reason.empty() ? "hydro dt estimation failed" : dt_advice.reason;
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason + "; dt_evidence=" + dt_advice.evidence;
      return loop;
    }
    double dt_s = std::max(0.0, config.run.cfl) * dt_advice.hard_cap_dt;
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      loop.failure_reason = "runtime dt is not positive finite";
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason;
      return loop;
    }
    dt_s = ClampDtToTargetTime(config, time_s, dt_s);
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      break;
    }
    loop.dt_estimator_wall_s += ElapsedSecondsSince(dt_timer);

    for (char stage : config.run.stage_order) {
      const auto stage_timer = std::chrono::steady_clock::now();
      RuntimeStageResult stage_result;
      switch (stage) {
        case 'H':
          stage_result = ExecuteHydroStage(state,
                                           geometry,
                                           config,
                                           dt_s,
                                           static_cast<std::uint64_t>(step_index),
                                           time_s);
          break;
        case 'T':
          stage_result = ExecuteThermalStage(state, geometry, config, dt_s);
          break;
        case 'E':
          stage_result = ExecuteEquilibrationStage(state, dt_s);
          break;
        case 'R':
          stage_result =
              ExecuteRadiationStage(state,
                                    geometry,
                                    group_layout,
                                    config,
                                    input_deck_path,
                                    serial_opacity_table_loaded ? &serial_opacity_table : nullptr,
                                    serial_opacity_table_loaded ? &serial_provider_options : nullptr,
                                    dt_s);
          break;
        case 'A':
          stage_result = ExecuteAlphaStage(state, geometry, config, dt_s);
          break;
        default:
          stage_result = FailStage(stage, "unsupported runtime stage");
          break;
      }
      AddRuntimeStageWallTime(loop, stage, ElapsedSecondsSince(stage_timer));
      AddRuntimeStagePerformance(loop, stage, stage_result);
      if (!stage_result.success) {
        loop.failure_reason = stage_result.failure_reason;
        loop.failure_diagnostics = stage_result.failure_diagnostics;
        return loop;
      }
    }

    time_s += dt_s;
    loop.steps_executed = static_cast<std::size_t>(local_step_index);
    loop.final_time_s = time_s;
    const auto wall_elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    const auto output_timer = std::chrono::steady_clock::now();
    const auto output = dec3d::io::WriteRuntimeStepOutputs(
        config, state, geometry, group_layout, profile_path, static_cast<std::size_t>(step_index),
        time_s, dt_s, wall_elapsed, "advanced", output_state);
    loop.runtime_output_wall_s += ElapsedSecondsSince(output_timer);
    if (!output.success) {
      loop.failure_reason = output.failure_reason.empty() ? "runtime step output failed"
                                                          : output.failure_reason;
      loop.failure_diagnostics = output.failure_diagnostics;
      return loop;
    }
  }

  const auto axisym = EvaluateAxisymmetricInvariants(state);
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) && !axisym.ok) {
    loop.failure_reason = "axisymmetric invariant violated";
    std::ostringstream failure;
    failure << "diagnostic_id=p5.runtime.loop.failure"
            << "; failure_reason=" << loop.failure_reason;
    AppendAxisymmetricInvariantDiagnostics(failure, config, axisym);
    loop.failure_diagnostics = failure.str();
    return loop;
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p5.runtime.loop"
         << "; runtime_time_loop_executed=true"
         << "; runtime_stage_order_executed=" << stage_order
         << "; runtime_stage_backend=serial_callable_runtime"
         << "; runtime_steps_executed=" << loop.steps_executed
         << "; final_time_s=" << loop.final_time_s
         << "; mesh_dimensionality=" << dec3d::io::ToString(config.mesh.dimensionality)
         << "; active_hydro_directions="
         << (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) ? "r,theta"
                                                                     : "r,theta,phi")
         << "; phi_sweep_executed="
         << (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality)
                 ? "false"
                 : (loop.h_hydro_phi_sweep_wall_s > 0.0 ? "true" : "false"))
         << "; target_time_trigger_enabled="
         << (TargetTimeEnabled(config) ? "true" : "false")
         << "; target_time_s=" << config.run.target_time_s
         << "; target_time_reached="
         << (TargetTimeReached(config, loop.final_time_s) ? "true" : "false")
         << "; h_stage_executed=" << (stage_order.find('H') != std::string::npos ? "true" : "false")
         << "; t_stage_executed=" << (stage_order.find('T') != std::string::npos ? "true" : "false")
         << "; e_stage_executed=" << (stage_order.find('E') != std::string::npos ? "true" : "false")
         << "; r_stage_executed=" << (stage_order.find('R') != std::string::npos ? "true" : "false")
         << "; a_stage_executed=" << (stage_order.find('A') != std::string::npos ? "true" : "false")
         << "; distributed_mpi_runtime=false"
         << "; serial_radiation_table_loaded_once="
         << (serial_opacity_table_loaded ? "true" : "false")
         << "; serial_implicit_cell_limit=" << kSerialImplicitRuntimeCellLimit;
  AppendAxisymmetricInvariantDiagnostics(report, config, axisym);
  report << "; global_matrix_diagnostics_present="
         << (loop.global_matrix_diagnostics_present ? "true" : "false")
         << "; global_phi_coupling_count=" << loop.global_phi_coupling_count
         << "; global_duplicate_column_row_count="
         << loop.global_duplicate_column_row_count;
  AppendRuntimeTimingDiagnostics(report, loop);
  loop.success = true;
  loop.report_line = report.str();
  return loop;
}

#ifdef DEC3D_ENABLE_HYPRE
[[nodiscard]] RuntimeLoopResult ExecuteP5DistributedRuntimeLoop(
    const dec3d::io::InputDeckConfig& config,
    dec3d::state::CanonicalState& global_state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& input_deck_path,
    const std::filesystem::path& profile_path,
    const MpiRuntimeContext& mpi_context,
    const RuntimeStartPoint& start) {
  RuntimeLoopResult loop;
  if (config.run.step_count <= 0) {
    loop.success = true;
    loop.report_line =
        "diagnostic_id=p5.runtime.loop; runtime_time_loop_executed=false; "
        "distributed_mpi_runtime=true";
    return loop;
  }
  for (char stage : config.run.stage_order) {
    if (!StageEnabled(config, stage)) {
      loop.failure_reason = "stage_order includes disabled physics stage";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason +
          "; disabled_stage=" + std::string(1, stage);
      return loop;
    }
  }

  const auto ownership = dec3d::transport::BuildDistributedDiffusionRowOwnership(
      MPI_COMM_WORLD,
      global_state.layout.radial_cells,
      global_state.layout.theta_cells,
      global_state.layout.phi_cells);
  if (!ownership.success) {
    loop.failure_reason = "distributed runtime ownership build failed";
    loop.failure_diagnostics = ownership.failure_diagnostics;
    return loop;
  }
  const auto hydro_decomposition = dec3d::mesh::BuildRadialOwnership(
      global_state.layout.radial_cells,
      static_cast<std::size_t>(mpi_context.rank_count));
  if (!hydro_decomposition.is_valid() ||
      mpi_context.rank < 0 ||
      static_cast<std::size_t>(mpi_context.rank) >= hydro_decomposition.slices.size()) {
    loop.failure_reason = hydro_decomposition.failure_reason.empty()
                              ? "distributed runtime hydro ownership build failed"
                              : hydro_decomposition.failure_reason;
    loop.failure_diagnostics =
        "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
    return loop;
  }
  const auto& hydro_slice =
      hydro_decomposition.slices[static_cast<std::size_t>(mpi_context.rank)];
  if (ownership.global_radial_begin != hydro_slice.begin_index ||
      ownership.global_radial_end != hydro_slice.end_index) {
    loop.failure_reason =
        "distributed runtime hydro and diffusion ownership layouts differ";
    loop.failure_diagnostics =
        "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason +
        "; hydro_begin=" + std::to_string(hydro_slice.begin_index) +
        "; hydro_end=" + std::to_string(hydro_slice.end_index) +
        "; diffusion_begin=" + std::to_string(ownership.global_radial_begin) +
        "; diffusion_end=" + std::to_string(ownership.global_radial_end);
    return loop;
  }
  auto runtime_global_geometry = geometry;
  auto local_geometry = BuildLocalRuntimeGeometrySlice(runtime_global_geometry,
                                                       hydro_slice,
                                                       global_state.layout.theta_cells,
                                                       global_state.layout.phi_cells);
  if (!local_geometry.is_valid()) {
    loop.failure_reason = "distributed runtime local hydro geometry slice failed";
    loop.failure_diagnostics =
        "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
    return loop;
  }
  auto local_state = MakeLocalSlabState(global_state, ownership);
  const auto gather_plan = BuildRuntimeSlabGatherPlan(MPI_COMM_WORLD, ownership);
  if (!gather_plan.success) {
    loop.failure_reason = "distributed runtime gather metadata precompute failed";
    loop.failure_diagnostics =
        "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
    return loop;
  }

  dec3d::radiation::TopsOpacityTable opacity_table;
  dec3d::radiation::TopsOpacityProviderOptions provider_options;
  if (config.physics.enable_radiation) {
    if (config.radiation.opacity_provider != "tops_dt_tabulated") {
      loop.failure_reason = "unsupported runtime opacity_provider";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
      return loop;
    }
    provider_options = RuntimeTopsProviderOptions(input_deck_path);
    if (provider_options.table_root.empty()) {
      loop.failure_reason = "TOPS opacity table root not found";
      loop.failure_diagnostics =
          "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
      return loop;
    }
    const auto table_result = dec3d::radiation::LoadTopsOpacityTable(provider_options);
    if (!table_result.success || table_result.report_line.empty()) {
      loop.failure_reason =
          table_result.failure_reason.empty() ? "TOPS opacity table load failed"
                                              : table_result.failure_reason;
      loop.failure_diagnostics = table_result.failure_diagnostics.empty()
                                     ? table_result.report_line
                                     : table_result.failure_diagnostics;
      return loop;
    }
    opacity_table = table_result.table;
  }

  const auto stage_order = JoinStageOrder(config.run.stage_order);
  std::string h_stage_backend = "none";
  std::string t_stage_backend = "none";
  std::string r_stage_backend = "none";
  std::string a_stage_backend = "none";
  bool allgather_after_distributed_stages = false;
  bool root_gather_for_checkpoint_outputs = false;
  bool runtime_scalar_output_uses_mpi_reduce = false;
  dec3d::transport::DistributedLaggedBoomerAmgCache radiation_amg_cache;
  double time_s = start.time_s;
  dec3d::io::RuntimeOutputStepState output_state;
  output_state.last_field_checkpoint_time_s = start.time_s;
  output_state.last_restart_checkpoint_time_s = start.time_s;
  output_state.last_history_profile_time_s = start.time_s;
  const auto wall_start = std::chrono::steady_clock::now();
  const auto noh_exact_inflow =
      BuildRuntimeNohExactInflowState(config, local_state, hydro_decomposition, mpi_context.rank);
  if (!noh_exact_inflow.complete) {
    loop.failure_reason = noh_exact_inflow.failure_reason.empty()
                              ? "noh exact inflow initialization failed"
                              : noh_exact_inflow.failure_reason;
    loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                               loop.failure_reason + "; " + noh_exact_inflow.report_line;
    return loop;
  }

  for (int local_step_index = 1; local_step_index <= config.run.step_count; ++local_step_index) {
    const std::size_t step_index = start.step + static_cast<std::size_t>(local_step_index);
    if (TargetTimeReached(config, time_s)) {
      break;
    }
    const auto dt_timer = std::chrono::steady_clock::now();
    dec3d::hydro::HydroOperator hydro_for_dt;
    auto hydro_options = RuntimeHydroOptions(config);
    hydro_for_dt.SetStaticGridOptions(hydro_options);
    dec3d::core::StageContext dt_context;
    dt_context.time_s = time_s;
    dt_context.dt_s = 0.0;
    dt_context.step = static_cast<std::uint64_t>(step_index);
    dt_context.phase_id = dec3d::core::PhaseId::p1;
    dt_context.contract_version = "p5.runtime.mpi_hypre_callable.v1";
    dt_context.mesh_snapshot_handle = "p5.runtime.mesh.current";
    dt_context.ownership_handle = "distributed_mpi_hydro";
    dt_context.diagnostics_sink_handle = "p5.runtime.dec3d_out";
    if (!hydro_for_dt.bind(dt_context, local_geometry, local_state)) {
      loop.failure_reason = "hydro dt estimator bind failed";
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason;
      return loop;
    }
    const auto dt_advice = hydro_for_dt.estimate_dt();
    if (!dt_advice.is_complete()) {
      loop.failure_reason =
          dt_advice.reason.empty() ? "hydro dt estimation failed" : dt_advice.reason;
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason + "; dt_evidence=" + dt_advice.evidence;
      return loop;
    }
    double global_hard_cap_dt = dt_advice.hard_cap_dt;
    MPI_Allreduce(
        MPI_IN_PLACE,
        &global_hard_cap_dt,
        1,
        MPI_DOUBLE,
        MPI_MIN,
        MPI_COMM_WORLD);
    double dt_s = std::max(0.0, config.run.cfl) * global_hard_cap_dt;
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      loop.failure_reason = "runtime dt is not positive finite";
      loop.failure_diagnostics = "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                 loop.failure_reason;
      return loop;
    }
    dt_s = ClampDtToTargetTime(config, time_s, dt_s);
    if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
      break;
    }
    loop.dt_estimator_wall_s += ElapsedSecondsSince(dt_timer);

    for (char stage : config.run.stage_order) {
      const auto stage_timer = std::chrono::steady_clock::now();
      RuntimeStageResult stage_result;
      switch (stage) {
        case 'H':
          stage_result = ExecuteDistributedHydroStage(local_state,
                                                      local_geometry,
                                                      runtime_global_geometry,
                                                      hydro_decomposition,
                                                      config,
                                                      dt_s,
                                                      time_s,
                                                      noh_exact_inflow);
          h_stage_backend =
              stage_result.success
                  ? (config.mesh.moving_mesh ? "mpi_macro_ale_hllc" : "mpi_static_hllc")
                  : h_stage_backend;
          break;
        case 'T':
          stage_result =
              ExecuteDistributedThermalStage(
                  local_state, ownership, runtime_global_geometry, config, dt_s);
          t_stage_backend =
              stage_result.success ? "hypre_parcsr_gmres_boomeramg" : t_stage_backend;
          break;
        case 'E':
          stage_result = ExecuteEquilibrationStage(local_state, dt_s);
          break;
        case 'R':
          stage_result = ExecuteDistributedRadiationStage(local_state,
                                                          ownership,
                                                          runtime_global_geometry,
                                                          group_layout,
                                                          opacity_table,
                                                          config,
                                                          provider_options,
                                                          &radiation_amg_cache,
                                                          dt_s);
          r_stage_backend =
              stage_result.success ? "hypre_parcsr_gmres_boomeramg" : r_stage_backend;
          break;
        case 'A':
          stage_result =
              ExecuteDistributedAlphaStage(
                  local_state, ownership, runtime_global_geometry, config, dt_s);
          a_stage_backend =
              stage_result.success ? "hypre_parcsr_gmres_boomeramg" : a_stage_backend;
          break;
        default:
          stage_result = FailStage(stage, "unsupported runtime stage");
          break;
      }
      AddRuntimeStageWallTime(loop, stage, ElapsedSecondsSince(stage_timer));
      AddRuntimeStagePerformance(loop, stage, stage_result);
      if (!stage_result.success) {
        loop.failure_reason = stage_result.failure_reason;
        loop.failure_diagnostics = stage_result.failure_diagnostics;
        return loop;
      }
    }

    time_s += dt_s;
    loop.steps_executed = static_cast<std::size_t>(local_step_index);
    loop.final_time_s = time_s;
    const auto wall_elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
    const auto output_timer = std::chrono::steady_clock::now();
    const bool field_checkpoint_due =
        dec3d::io::ShouldWriteCheckpoint(
            dec3d::io::FieldCheckpointTrigger(config.output),
            static_cast<std::size_t>(step_index),
            time_s,
            output_state.last_field_checkpoint_time_s);
    const bool restart_checkpoint_due =
        dec3d::io::ShouldWriteCheckpoint(
            dec3d::io::RestartCheckpointTrigger(config.output),
            static_cast<std::size_t>(step_index),
            time_s,
            output_state.last_restart_checkpoint_time_s);
    const bool history_profile_due =
        dec3d::io::ShouldWriteCheckpoint(
            dec3d::io::HistoryProfileTrigger(config.output),
            static_cast<std::size_t>(step_index),
            time_s,
            output_state.last_history_profile_time_s);
    const bool full_state_output_due =
        field_checkpoint_due || restart_checkpoint_due || history_profile_due;
    const bool dec3d_out_due =
        config.output.write_dec3d_out_every_steps > 0 &&
        static_cast<std::size_t>(step_index) %
            static_cast<std::size_t>(config.output.write_dec3d_out_every_steps) == 0u;

    dec3d::io::RuntimeOutputWriteResult output;
    bool output_attempted = false;
    if (full_state_output_due) {
      GatherLocalStateToRoot(
          MPI_COMM_WORLD, ownership, gather_plan, local_state, global_state, 0);
      root_gather_for_checkpoint_outputs = true;
      output_attempted = true;
      if (mpi_context.rank == 0) {
        output = dec3d::io::WriteRuntimeStepOutputs(
            config, global_state, runtime_global_geometry, group_layout, profile_path,
            static_cast<std::size_t>(step_index), time_s, dt_s, wall_elapsed, "advanced",
            output_state);
      } else {
        output.success = true;
      }
    } else if (dec3d_out_due) {
      const auto local_extrema = dec3d::io::ComputeRuntimeOutputExtrema(local_state);
      int local_extrema_ok = local_extrema.success ? 1 : 0;
      int global_extrema_ok = 0;
      MPI_Allreduce(
          &local_extrema_ok, &global_extrema_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if (global_extrema_ok == 0) {
        loop.failure_reason = local_extrema.failure_reason.empty()
                                  ? "runtime output extrema failed"
                                  : local_extrema.failure_reason;
        loop.failure_diagnostics = local_extrema.failure_diagnostics.empty()
                                       ? "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                             loop.failure_reason
                                       : local_extrema.failure_diagnostics;
        return loop;
      }
      const auto global_extrema =
          ReduceRuntimeExtremaToRoot(MPI_COMM_WORLD, local_extrema.extrema, 0);
      runtime_scalar_output_uses_mpi_reduce = true;
      output_attempted = true;
      if (mpi_context.rank == 0) {
        output = dec3d::io::WriteRuntimeStepDec3DOut(
            config,
            static_cast<std::size_t>(step_index),
            time_s,
            dt_s,
            wall_elapsed,
            "advanced",
            global_extrema);
      } else {
        output.success = true;
      }
    }

    if (output_attempted) {
      int local_output_ok = output.success ? 1 : 0;
      int global_output_ok = 0;
      MPI_Allreduce(
          &local_output_ok, &global_output_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if (global_output_ok == 0) {
        loop.failure_reason = output.failure_reason.empty() ? "runtime step output failed"
                                                            : output.failure_reason;
        loop.failure_diagnostics = output.failure_diagnostics.empty()
                                       ? "diagnostic_id=p5.runtime.loop.failure; failure_reason=" +
                                             loop.failure_reason
                                       : output.failure_diagnostics;
        return loop;
      }
      if (field_checkpoint_due) {
        output_state.last_field_checkpoint_time_s = time_s;
      }
      if (restart_checkpoint_due) {
        output_state.last_restart_checkpoint_time_s = time_s;
      }
      if (history_profile_due) {
        output_state.last_history_profile_time_s = time_s;
      }
    }
    loop.runtime_output_wall_s += ElapsedSecondsSince(output_timer);
  }

  const auto axisym =
      ReduceAxisymmetricInvariants(MPI_COMM_WORLD, EvaluateAxisymmetricInvariants(local_state));
  if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) && !axisym.ok) {
    loop.failure_reason = "axisymmetric invariant violated";
    std::ostringstream failure;
    failure << "diagnostic_id=p5.runtime.loop.failure"
            << "; failure_reason=" << loop.failure_reason;
    AppendAxisymmetricInvariantDiagnostics(failure, config, axisym);
    loop.failure_diagnostics = failure.str();
    return loop;
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p5.runtime.loop"
         << "; runtime_time_loop_executed=true"
         << "; runtime_stage_order_executed=" << stage_order
         << "; runtime_stage_backend=mpi_hypre_callable_runtime"
         << "; runtime_steps_executed=" << loop.steps_executed
         << "; final_time_s=" << loop.final_time_s
         << "; mesh_dimensionality=" << dec3d::io::ToString(config.mesh.dimensionality)
         << "; active_hydro_directions="
         << (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) ? "r,theta"
                                                                     : "r,theta,phi")
         << "; phi_sweep_executed="
         << (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality)
                 ? "false"
                 : (loop.h_hydro_phi_sweep_wall_s > 0.0 ? "true" : "false"))
         << "; target_time_trigger_enabled="
         << (TargetTimeEnabled(config) ? "true" : "false")
         << "; target_time_s=" << config.run.target_time_s
         << "; target_time_reached="
         << (TargetTimeReached(config, loop.final_time_s) ? "true" : "false")
         << "; h_stage_executed=" << (stage_order.find('H') != std::string::npos ? "true" : "false")
         << "; t_stage_executed=" << (stage_order.find('T') != std::string::npos ? "true" : "false")
         << "; e_stage_executed=" << (stage_order.find('E') != std::string::npos ? "true" : "false")
         << "; r_stage_executed=" << (stage_order.find('R') != std::string::npos ? "true" : "false")
         << "; a_stage_executed=" << (stage_order.find('A') != std::string::npos ? "true" : "false")
         << "; distributed_mpi_runtime=true"
         << "; mpi_rank_count=" << mpi_context.rank_count
         << "; hydro_runtime_mode=distributed_mpi_hydro"
         << "; h_stage_backend=" << h_stage_backend
         << "; hydro_halo_exchange_used=true"
         << "; hydro_scalar_bundle_halo=true"
         << "; " << noh_exact_inflow.report_line
         << "; distributed_implicit_runtime=true"
         << "; t_stage_backend=" << t_stage_backend
         << "; r_stage_backend=" << r_stage_backend
         << "; a_stage_backend=" << a_stage_backend
         << "; rank0_gather_solve_used=false"
         << "; local_slab_allgather_after_distributed_stages="
         << (allgather_after_distributed_stages ? "true" : "false")
         << "; root_gather_for_checkpoint_outputs="
         << (root_gather_for_checkpoint_outputs ? "true" : "false")
         << "; runtime_scalar_output_uses_mpi_reduce="
         << (runtime_scalar_output_uses_mpi_reduce ? "true" : "false")
         << "; gather_metadata_precomputed=true"
         << "; serial_implicit_cell_limit=" << kSerialImplicitRuntimeCellLimit;
  AppendAxisymmetricInvariantDiagnostics(report, config, axisym);
  report << "; global_matrix_diagnostics_present="
         << (loop.global_matrix_diagnostics_present ? "true" : "false")
         << "; global_phi_coupling_count=" << loop.global_phi_coupling_count
         << "; global_duplicate_column_row_count="
         << loop.global_duplicate_column_row_count;
  AppendRuntimeTimingDiagnostics(report, loop);
  loop.success = true;
  loop.report_line = report.str();
  return loop;
}
#endif

}  // namespace

Dec3DAppResult RunDec3DCommandLine(const std::vector<std::string>& argv) noexcept {
#ifdef DEC3D_ENABLE_HYPRE
  const auto mpi_context = DetectMpiRuntime();
#endif
  const auto args = dec3d::io::ParseP5IORuntimeArguments(argv);
  if (!args.success) {
    return Fail(args.failure_reason, args.failure_diagnostics);
  }
  const auto deck = dec3d::io::LoadInputDeck(args.input_deck_path);
  if (!deck.success || !dec3d::io::ValidateInputDeckDiagnostics(deck)) {
    return Fail(deck.failure_reason.empty() ? "input deck validation failed" : deck.failure_reason,
                deck.failure_diagnostics);
  }
  const auto profile =
      dec3d::io::LoadRadialProfile(args.profile_path,
                                   deck.config.initial_condition.profile_radius_unit);
  if (!profile.success || !dec3d::io::ValidateRadialProfileDiagnostics(profile)) {
    return Fail(profile.failure_reason.empty() ? "profile validation failed" : profile.failure_reason,
                profile.failure_diagnostics);
  }
  auto init =
      dec3d::initialization::InitializeFromRadialProfile(deck.config, profile.profile,
                                                         args.profile_path);
  if (!init.success ||
      !dec3d::initialization::ValidateProfileInitializationDiagnostics(init)) {
    return Fail(init.failure_reason.empty() ? "profile initialization failed" : init.failure_reason,
                init.failure_diagnostics);
  }

  RuntimeStartPoint runtime_start;
  bool restart_initialization_used = false;
  std::string restart_report_line =
      "diagnostic_id=p5.io.restart_load; restart_file_present=false";
  if (!args.restart_path.empty()) {
    const auto restart =
        LoadRestartInitialState(args.restart_path, deck.config, init.group_layout);
    if (!restart.success) {
      return Fail(restart.failure_reason.empty() ? "restart load failed"
                                                 : restart.failure_reason,
                  restart.failure_diagnostics);
    }
    init.state = std::move(restart.state);
    runtime_start.step = restart.step;
    runtime_start.time_s = restart.time_s;
    restart_initialization_used = true;
    restart_report_line = restart.report_line;
  }

#ifdef DEC3D_ENABLE_HYPRE
  dec3d::io::RuntimeOutputWriteResult outputs;
  if (restart_initialization_used) {
    outputs.success = true;
    outputs.report_line =
        "diagnostic_id=p5.io.runtime_output; initial_outputs_skipped_for_restart=true";
  } else if (!mpi_context.distributed || mpi_context.rank == 0) {
    outputs = dec3d::io::WriteInitialRuntimeOutputs(
        deck.config,
        init.state,
        init.geometry,
        init.group_layout,
        args.profile_path);
  } else {
    outputs.success = true;
    outputs.report_line =
        "diagnostic_id=p5.io.runtime_output; output_skipped_on_nonzero_mpi_rank=true";
  }
#else
  const auto outputs =
      restart_initialization_used
          ? [] {
              dec3d::io::RuntimeOutputWriteResult skipped;
              skipped.success = true;
              skipped.report_line =
                  "diagnostic_id=p5.io.runtime_output; initial_outputs_skipped_for_restart=true";
              return skipped;
            }()
          : dec3d::io::WriteInitialRuntimeOutputs(
                deck.config,
                init.state,
                init.geometry,
                init.group_layout,
                args.profile_path);
#endif
  if (!outputs.success) {
    return Fail(outputs.failure_reason.empty() ? "runtime output failed"
                                               : outputs.failure_reason,
                outputs.failure_diagnostics);
  }

#ifdef DEC3D_ENABLE_HYPRE
  const auto loop =
      mpi_context.distributed
          ? ExecuteP5DistributedRuntimeLoop(
                deck.config,
                init.state,
                init.geometry,
                init.group_layout,
                args.input_deck_path,
                args.profile_path,
                mpi_context,
                runtime_start)
          : ExecuteP5RuntimeLoop(
                deck.config,
                init.state,
                init.geometry,
                init.group_layout,
                args.input_deck_path,
                args.profile_path,
                runtime_start);
#else
  const auto loop = ExecuteP5RuntimeLoop(
      deck.config,
      init.state,
      init.geometry,
      init.group_layout,
      args.input_deck_path,
      args.profile_path,
      runtime_start);
#endif
  if (!loop.success) {
    return Fail(loop.failure_reason.empty() ? "runtime loop failed" : loop.failure_reason,
                loop.failure_diagnostics);
  }

  std::ostringstream report;
  report << std::setprecision(17);
  report << "diagnostic_id=p5.io.runtime_entry"
         << "; input_deck_report_present=true"
         << "; radial_profile_report_present=true"
         << "; initialization_report_present=true"
         << "; profile_file_source=command_line.profile"
         << "; implicit_profile_used=false"
         << "; canonical_state_initialized=true"
         << "; restart_initialization_used="
         << (restart_initialization_used ? "true" : "false")
         << "; " << restart_report_line
         << "; runtime_output_report_present=true"
         << "; dec3d_out_written=" << (outputs.dec3d_out_written ? "true" : "false")
         << "; field_checkpoint_written="
         << (outputs.field_checkpoint_written ? "true" : "false")
         << "; restart_checkpoint_written="
         << (outputs.restart_checkpoint_written ? "true" : "false")
         << "; case_name=" << deck.config.run.case_name
         << "; " << loop.report_line;

  Dec3DAppResult result;
  result.exit_code = 0;
  result.report_line = report.str();
  result.output_report_line = outputs.report_line;
  return result;
}

}  // namespace dec3d::app
