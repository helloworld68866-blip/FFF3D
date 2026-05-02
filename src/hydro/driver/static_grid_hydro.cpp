#include "hydro/driver/static_grid_hydro.hpp"

#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/driver/radial_overlap_remap.hpp"
#include "hydro/driver/radial_boundary_contract.hpp"
#include "hydro/driver/theta_pole_boundary.hpp"
#include "hydro/reconstruction/ppm_reconstruction.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "hydro/source/geometric_source.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/macro_zoning/macro_zoning.hpp"

#include <algorithm>
#include <charconv>
#include <cctype>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace dec3d::hydro {

namespace {

enum class SweepDirection {
  radial,
  theta,
  phi,
};

[[nodiscard]] HydroDirection ToHydroDirection(SweepDirection direction) noexcept;
[[nodiscard]] const char* SweepDirectionName(SweepDirection direction) noexcept;

[[nodiscard]] std::vector<double> ScaleRadiationBundle(
    const std::vector<double>& values,
    double factor) {
  std::vector<double> scaled;
  scaled.reserve(values.size());
  for (const double value : values) {
    scaled.push_back(value * factor);
  }
  return scaled;
}

[[nodiscard]] std::vector<double> AddRadiationBundles(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) {
  if (lhs.size() != rhs.size()) {
    return {};
  }

  std::vector<double> sum(lhs.size(), 0.0);
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    sum[index] = lhs[index] + rhs[index];
  }
  return sum;
}

[[nodiscard]] std::vector<double> SubtractRadiationBundles(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) {
  if (lhs.size() != rhs.size()) {
    return {};
  }

  std::vector<double> delta(lhs.size(), 0.0);
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    delta[index] = lhs[index] - rhs[index];
  }
  return delta;
}

struct GridShape {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
};

struct HydroStateSnapshot {
  GridShape shape;
  std::vector<HydroConservativeState> cells;

  [[nodiscard]] bool is_complete() const noexcept {
    return shape.radial_cells > 0u &&
           shape.theta_cells > 0u &&
           shape.phi_cells > 0u &&
           cells.size() == shape.radial_cells * shape.theta_cells * shape.phi_cells;
  }
};

struct DirectionalPpmSummary {
  bool executed{false};
  std::size_t ghost_layers{0};
  std::size_t line_count{0};
  std::size_t downgraded_interface_count{0};
  std::string report_line;
};

struct RadialAleFluxSummary {
  bool executed{false};
  double max_face_speed{0.0};
  std::string report_line;
  std::vector<std::string> debug_report_lines;
};

struct MacroZonedHydroSummary {
  bool executed{false};
  bool radial_executed{false};
  bool theta_executed{false};
  bool phi_executed{false};
  std::size_t coarse_cells{0};
  std::size_t theta_bands{0};
  bool coarse_ghost_bootstrap_executed{false};
  bool coarse_source_budget_executed{false};
  DirectionalPpmSummary radial_ppm_summary;
  DirectionalPpmSummary theta_ppm_summary;
  DirectionalPpmSummary phi_ppm_summary;
  std::string coarse_ghost_bootstrap_report_line;
  std::string coarse_source_budget_report_line;
  std::string detect_report_line;
  std::string restrict_report_line;
  std::string coarse_update_report_line;
  std::string prolong_report_line;
  std::vector<std::string> angular_stage_report_lines;
  double detect_wall_s{0.0};
  double restrict_wall_s{0.0};
  double update_wall_s{0.0};
  double radial_update_wall_s{0.0};
  double theta_update_wall_s{0.0};
  double phi_update_wall_s{0.0};
  double state_update_wall_s{0.0};
  double prolong_wall_s{0.0};
};

[[nodiscard]] double ElapsedHydroSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now() - start)
      .count();
}

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] std::string_view Trim(std::string_view text) noexcept {
  while (!text.empty() && std::isspace(static_cast<unsigned char>(text.front())) != 0) {
    text.remove_prefix(1u);
  }
  while (!text.empty() && std::isspace(static_cast<unsigned char>(text.back())) != 0) {
    text.remove_suffix(1u);
  }
  return text;
}

[[nodiscard]] bool ParseDoubleField(
    std::string_view report_line,
    std::string_view key,
    double* value_out) noexcept {
  if (value_out == nullptr) {
    return false;
  }

  const std::string needle = std::string(key) + "=";
  const std::size_t key_position = report_line.find(needle);
  if (key_position == std::string_view::npos) {
    return false;
  }

  const std::size_t value_begin = key_position + needle.size();
  const std::size_t value_end = report_line.find(';', value_begin);
  const auto raw_value = Trim(report_line.substr(
      value_begin,
      value_end == std::string_view::npos ? std::string_view::npos : value_end - value_begin));
  if (raw_value.empty()) {
    return false;
  }

  const auto* begin = raw_value.data();
  const auto* end = raw_value.data() + raw_value.size();
  const auto parsed = std::from_chars(begin, end, *value_out);
  return parsed.ec == std::errc{} && parsed.ptr == end;
}

[[nodiscard]] std::size_t LinearIndex(
    const GridShape& shape,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((radial * shape.theta_cells) + theta) * shape.phi_cells + phi;
}

[[nodiscard]] HydroConservativeState ToDirectionalState(
    const HydroConservativeState& state,
    SweepDirection direction) noexcept {
  switch (direction) {
    case SweepDirection::radial:
      return state;
    case SweepDirection::theta:
      return {
          state.rho,
          state.mom_theta,
          state.mom_r,
          state.mom_phi,
          state.e_fluid_total,
          state.chi_e,
          state.alpha_chi,
          state.radiation_chi};
    case SweepDirection::phi:
      return {
          state.rho,
          state.mom_phi,
          state.mom_r,
          state.mom_theta,
          state.e_fluid_total,
          state.chi_e,
          state.alpha_chi,
          state.radiation_chi};
  }

  return {};
}

[[nodiscard]] HydroConservativeState FromDirectionalFlux(
    const HydroConservativeState& flux,
    SweepDirection direction) noexcept {
  switch (direction) {
    case SweepDirection::radial:
      return flux;
    case SweepDirection::theta:
      return {
          flux.rho,
          flux.mom_theta,
          flux.mom_r,
          flux.mom_phi,
          flux.e_fluid_total,
          flux.chi_e,
          flux.alpha_chi,
          flux.radiation_chi};
    case SweepDirection::phi:
      return {
          flux.rho,
          flux.mom_theta,
          flux.mom_phi,
          flux.mom_r,
          flux.e_fluid_total,
          flux.chi_e,
          flux.alpha_chi,
          flux.radiation_chi};
  }

  return {};
}

[[nodiscard]] HydroConservativeState LoadCellState(
    const dec3d::state::HydroStateView& hydro_view,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return {
      (*hydro_view.rho)(radial, theta, phi),
      (*hydro_view.mom_r)(radial, theta, phi),
      (*hydro_view.mom_theta)(radial, theta, phi),
      (*hydro_view.mom_phi)(radial, theta, phi),
      (*hydro_view.e_fluid_total)(radial, theta, phi),
      hydro_view.chi_e(radial, theta, phi),
      hydro_view.operator_local_alpha_chi ? hydro_view.alpha_chi(radial, theta, phi) : 0.0,
      [&]() {
        std::vector<double> bundle;
        bundle.reserve(hydro_view.radiation_chi.size());
        for (const auto& group_chi : hydro_view.radiation_chi) {
          bundle.push_back(group_chi(radial, theta, phi));
        }
        return bundle;
      }()};
}

void StoreCellState(
    const HydroConservativeState& cell_state,
    dec3d::state::HydroStateView& hydro_view,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  (*hydro_view.rho)(radial, theta, phi) = cell_state.rho;
  (*hydro_view.mom_r)(radial, theta, phi) = cell_state.mom_r;
  (*hydro_view.mom_theta)(radial, theta, phi) = cell_state.mom_theta;
  (*hydro_view.mom_phi)(radial, theta, phi) = cell_state.mom_phi;
  (*hydro_view.e_fluid_total)(radial, theta, phi) = cell_state.e_fluid_total;
  hydro_view.chi_e(radial, theta, phi) = cell_state.chi_e;
  if (hydro_view.operator_local_alpha_chi &&
      !hydro_view.alpha_chi.empty()) {
    hydro_view.alpha_chi(radial, theta, phi) = cell_state.alpha_chi;
  }
  for (std::size_t group = 0; group < hydro_view.radiation_chi.size() &&
                              group < cell_state.radiation_chi.size(); ++group) {
    hydro_view.radiation_chi[group](radial, theta, phi) =
        cell_state.radiation_chi[group];
  }
}

[[nodiscard]] dec3d::core::Array3D<HydroConservativeState> CaptureCellArray(
    const dec3d::state::HydroStateView& hydro_view) {
  dec3d::core::Array3D<HydroConservativeState> cells(
      hydro_view.rho->extent_r(),
      hydro_view.rho->extent_theta(),
      hydro_view.rho->extent_phi(),
      HydroConservativeState{});
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < cells.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < cells.extent_phi(); ++phi) {
        cells(radial, theta, phi) = LoadCellState(hydro_view, radial, theta, phi);
      }
    }
  }
  return cells;
}

void StoreCellArray(
    const dec3d::core::Array3D<HydroConservativeState>& cells,
    dec3d::state::HydroStateView& hydro_view) {
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < cells.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < cells.extent_phi(); ++phi) {
        StoreCellState(cells(radial, theta, phi), hydro_view, radial, theta, phi);
      }
    }
  }
}

[[nodiscard]] dec3d::state::HydroStateView BuildOwnedHydroWorkViewFromSnapshot(
    const HydroStateSnapshot& snapshot) {
  dec3d::state::HydroStateView view;
  view.owned_rho = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.owned_mom_r = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.owned_mom_theta = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.owned_mom_phi = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.owned_e_fluid_total = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.owned_e_electron = std::make_unique<dec3d::core::Array3D<double>>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);

  view.rho = view.owned_rho.get();
  view.mom_r = view.owned_mom_r.get();
  view.mom_theta = view.owned_mom_theta.get();
  view.mom_phi = view.owned_mom_phi.get();
  view.e_fluid_total = view.owned_e_fluid_total.get();
  view.e_electron = view.owned_e_electron.get();
  view.chi_e = dec3d::core::Array3D<double>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.alpha_chi = dec3d::core::Array3D<double>(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);
  view.operator_local_alpha_chi = true;
  std::size_t radiation_group_count = 0u;
  for (const auto& cell : snapshot.cells) {
    radiation_group_count = std::max(radiation_group_count, cell.radiation_chi.size());
  }
  view.radiation_chi.clear();
  view.radiation_chi.reserve(radiation_group_count);
  for (std::size_t group = 0; group < radiation_group_count; ++group) {
    view.radiation_chi.emplace_back(
        snapshot.shape.radial_cells,
        snapshot.shape.theta_cells,
        snapshot.shape.phi_cells,
        0.0);
  }
  view.transactional_work_copy = true;

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(snapshot.shape, radial, theta, phi);
        const auto& cell = snapshot.cells[index];
        (*view.rho)(radial, theta, phi) = cell.rho;
        (*view.mom_r)(radial, theta, phi) = cell.mom_r;
        (*view.mom_theta)(radial, theta, phi) = cell.mom_theta;
        (*view.mom_phi)(radial, theta, phi) = cell.mom_phi;
        (*view.e_fluid_total)(radial, theta, phi) = cell.e_fluid_total;
        view.chi_e(radial, theta, phi) = cell.chi_e;
        view.alpha_chi(radial, theta, phi) = cell.alpha_chi;
        for (std::size_t group = 0; group < view.radiation_chi.size() &&
                                    group < cell.radiation_chi.size(); ++group) {
          view.radiation_chi[group](radial, theta, phi) = cell.radiation_chi[group];
        }
        (*view.e_electron)(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(
                dec3d::state::ElectronPressureFromChiE(cell.chi_e));
      }
    }
  }

  std::ostringstream report;
  report << "hydro_view_complete=true"
         << "; operator_local_chi_e=true"
         << "; operator_local_alpha_chi=true"
         << "; operator_local_radiation_chi=true"
         << "; radiation_group_count=" << view.radiation_chi.size()
         << "; transactional_work_copy=true"
         << "; cell_count=" << view.rho->size();
  view.report_line = report.str();
  return view;
}

[[nodiscard]] HydroConservativeState Subtract(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  return {
      lhs.rho - rhs.rho,
      lhs.mom_r - rhs.mom_r,
      lhs.mom_theta - rhs.mom_theta,
      lhs.mom_phi - rhs.mom_phi,
      lhs.e_fluid_total - rhs.e_fluid_total,
      lhs.chi_e - rhs.chi_e,
      lhs.alpha_chi - rhs.alpha_chi,
      SubtractRadiationBundles(lhs.radiation_chi, rhs.radiation_chi)};
}

[[nodiscard]] HydroConservativeState Add(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  return {
      lhs.rho + rhs.rho,
      lhs.mom_r + rhs.mom_r,
      lhs.mom_theta + rhs.mom_theta,
      lhs.mom_phi + rhs.mom_phi,
      lhs.e_fluid_total + rhs.e_fluid_total,
      lhs.chi_e + rhs.chi_e,
      lhs.alpha_chi + rhs.alpha_chi,
      AddRadiationBundles(lhs.radiation_chi, rhs.radiation_chi)};
}

[[nodiscard]] double MovingInterfaceNormalPressure(
    const HydroConservativeState& moving_flux,
    const HydroConservativeState& selected_state,
    double face_speed) noexcept {
  if (!moving_flux.is_finite() ||
      !selected_state.is_finite() ||
      !(selected_state.rho > 0.0) ||
      !std::isfinite(face_speed)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double selected_velocity = selected_state.mom_r / selected_state.rho;
  const double advective_momentum_flux =
      selected_state.mom_r * (selected_velocity - face_speed);
  return moving_flux.mom_r - advective_momentum_flux;
}

[[nodiscard]] double FaceConsistentPressureBalance(
    double lower_face_pressure,
    double upper_face_pressure) noexcept {
  if (!std::isfinite(lower_face_pressure) ||
      !std::isfinite(upper_face_pressure)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 0.5 * (lower_face_pressure + upper_face_pressure);
}

[[nodiscard]] HydroConservativeState Scale(
    const HydroConservativeState& state,
    double factor) noexcept {
  return {
      state.rho * factor,
      state.mom_r * factor,
      state.mom_theta * factor,
      state.mom_phi * factor,
      state.e_fluid_total * factor,
      state.chi_e * factor,
      state.alpha_chi * factor,
      ScaleRadiationBundle(state.radiation_chi, factor)};
}

[[nodiscard]] HydroConservativeState ToExtensiveState(
    const HydroConservativeState& state,
    double volume) noexcept {
  return Scale(state, volume);
}

[[nodiscard]] HydroConservativeState FromExtensiveState(
    const HydroConservativeState& extensive_state,
    double volume) noexcept {
  if (!(volume > 0.0) || !std::isfinite(volume)) {
    return {};
  }

  return Scale(extensive_state, 1.0 / volume);
}

[[nodiscard]] HydroConservativeState SubtractScaled(
    const HydroConservativeState& state,
    const HydroConservativeState& flux_delta,
    double factor) noexcept {
  const auto scaled = Scale(flux_delta, factor);
  return {
      state.rho - scaled.rho,
      state.mom_r - scaled.mom_r,
      state.mom_theta - scaled.mom_theta,
      state.mom_phi - scaled.mom_phi,
      state.e_fluid_total - scaled.e_fluid_total,
      state.chi_e - scaled.chi_e,
      state.alpha_chi - scaled.alpha_chi,
      SubtractRadiationBundles(state.radiation_chi, scaled.radiation_chi)};
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const GridShape& shape,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return geometry.cell_volumes[LinearIndex(shape, radial, theta, phi)];
}

[[nodiscard]] double CellSolidAngleFactor(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta,
    std::size_t phi) noexcept {
  return (std::cos(geometry.theta_faces[theta]) -
          std::cos(geometry.theta_faces[theta + 1u])) *
         (geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi]);
}

[[nodiscard]] double RadialShellVolumeFactorFromFaces(
    const std::vector<double>& radial_faces,
    std::size_t radial) noexcept {
  const double radial_inner = radial_faces[radial];
  const double radial_outer = radial_faces[radial + 1u];
  return ((radial_outer * radial_outer * radial_outer) -
          (radial_inner * radial_inner * radial_inner)) / 3.0;
}

[[nodiscard]] HydroConservativeState BuildAleFactorizedUpdatedState(
    const HydroConservativeState& old_state,
    const HydroConservativeState& flux_delta_extensive,
    const HydroConservativeState& source_delta_extensive,
    const HydroConservativeState& radial_delta_per_solid_angle,
    const HydroConservativeState& radial_delta_extensive,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double angular_factor = CellSolidAngleFactor(old_geometry, theta, phi);
  const double old_radial_volume =
      RadialShellVolumeFactorFromFaces(old_geometry.radial_faces, radial);
  const double new_radial_volume =
      RadialShellVolumeFactorFromFaces(new_geometry.radial_faces, radial);
  if (!(angular_factor > 0.0) ||
      !(old_radial_volume > 0.0) ||
      !(new_radial_volume > 0.0) ||
      !std::isfinite(angular_factor) ||
      !std::isfinite(old_radial_volume) ||
      !std::isfinite(new_radial_volume)) {
    return {};
  }

  const auto non_radial_extensive_delta =
      Add(Subtract(flux_delta_extensive, radial_delta_extensive),
          source_delta_extensive);
  const auto per_solid_angle_state = Add(
      Add(Scale(old_state, old_radial_volume),
          radial_delta_per_solid_angle),
      Scale(non_radial_extensive_delta, 1.0 / angular_factor));
  return Scale(per_solid_angle_state, 1.0 / new_radial_volume);
}

[[nodiscard]] HydroConservativeState BuildMacroRadialFactorizedCoarseState(
    const HydroConservativeState& old_state,
    const HydroConservativeState& coarse_extensive_delta,
    const HydroConservativeState& radial_delta_per_solid_angle_reference,
    const HydroConservativeState& radial_delta_per_solid_angle_weighted_delta_sum,
    const HydroConservativeState& radial_delta_extensive,
    double radial_delta_solid_angle_weight,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry,
    const dec3d::mesh::MacroZoneCell& cell) noexcept {
  const double old_radial_volume =
      RadialShellVolumeFactorFromFaces(old_geometry.radial_faces, cell.radial_index);
  const double new_radial_volume =
      RadialShellVolumeFactorFromFaces(new_geometry.radial_faces, cell.radial_index);
  const double angular_factor = cell.total_volume / old_radial_volume;
  if (!(old_radial_volume > 0.0) ||
      !(new_radial_volume > 0.0) ||
      !(angular_factor > 0.0) ||
      !std::isfinite(old_radial_volume) ||
      !std::isfinite(new_radial_volume) ||
      !std::isfinite(angular_factor)) {
    return {};
  }

  HydroConservativeState radial_delta_per_solid_angle{};
  if ((radial_delta_solid_angle_weight > 0.0) &&
      std::isfinite(radial_delta_solid_angle_weight)) {
    radial_delta_per_solid_angle = Add(
        radial_delta_per_solid_angle_reference,
        Scale(
            radial_delta_per_solid_angle_weighted_delta_sum,
            1.0 / radial_delta_solid_angle_weight));
  }
  const auto non_radial_extensive_delta =
      Subtract(coarse_extensive_delta, radial_delta_extensive);
  const auto new_per_solid_angle_extensive = Add(
      Add(Scale(old_state, old_radial_volume), radial_delta_per_solid_angle),
      Scale(non_radial_extensive_delta, 1.0 / angular_factor));
  return Scale(new_per_solid_angle_extensive, 1.0 / new_radial_volume);
}

[[nodiscard]] double CellMass(
    const HydroConservativeState& state,
    double volume) noexcept {
  return state.rho * volume;
}

[[nodiscard]] double CellRadialMomentum(
    const HydroConservativeState& state,
    double volume) noexcept {
  return state.mom_r * volume;
}

[[nodiscard]] double CellFluidEnergy(
    const HydroConservativeState& state,
    double volume) noexcept {
  return state.e_fluid_total * volume;
}

struct AngularStageMetrics {
  bool valid{false};
  double max_rho_angular_relative_spread{0.0};
  std::size_t max_spread_radial_index{0};
  double max_rho_absolute_spread{0.0};
  std::size_t max_rho_absolute_spread_radial_index{0};
  double max_mom_r_absolute_spread{0.0};
  std::size_t max_mom_r_absolute_spread_radial_index{0};
  double max_e_total_absolute_spread{0.0};
  std::size_t max_e_total_absolute_spread_radial_index{0};
  double max_chi_e_absolute_spread{0.0};
  std::size_t max_chi_e_absolute_spread_radial_index{0};
  double max_pressure_angular_relative_spread{0.0};
  std::size_t max_pressure_spread_radial_index{0};
  double max_pressure_absolute_spread{0.0};
  std::size_t max_pressure_absolute_spread_radial_index{0};
  double max_abs_mom_theta{0.0};
  double max_abs_mom_phi{0.0};
  std::size_t max_abs_mom_theta_radial{0};
  std::size_t max_abs_mom_theta_theta{0};
  std::size_t max_abs_mom_theta_phi{0};
  std::size_t max_abs_mom_phi_radial{0};
  std::size_t max_abs_mom_phi_theta{0};
  std::size_t max_abs_mom_phi_phi{0};
  std::size_t max_unweighted_tangential_radial_index{0};
  double max_unweighted_tangential_radial_sum{0.0};
  double tangential_abs_momentum{0.0};
  double radial_abs_momentum{0.0};
  double tangential_to_radial_momentum_ratio{0.0};
  double unweighted_tangential_abs_momentum{0.0};
  double unweighted_radial_abs_momentum{0.0};
  double unweighted_tangential_to_radial_momentum_ratio{0.0};
};

[[nodiscard]] std::string FormatAngularStageMetrics(
    std::string_view stage,
    std::string_view representation,
    const AngularStageMetrics& metrics) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "stage=" << stage
         << "; representation=" << representation
         << "; valid=" << (metrics.valid ? "true" : "false")
         << "; max_rho_angular_relative_spread=" << metrics.max_rho_angular_relative_spread
         << "; max_spread_radial_index=" << metrics.max_spread_radial_index
         << "; max_rho_absolute_spread=" << metrics.max_rho_absolute_spread
         << "; max_rho_absolute_spread_radial_index="
         << metrics.max_rho_absolute_spread_radial_index
         << "; max_mom_r_absolute_spread=" << metrics.max_mom_r_absolute_spread
         << "; max_mom_r_absolute_spread_radial_index="
         << metrics.max_mom_r_absolute_spread_radial_index
         << "; max_e_total_absolute_spread=" << metrics.max_e_total_absolute_spread
         << "; max_e_total_absolute_spread_radial_index="
         << metrics.max_e_total_absolute_spread_radial_index
         << "; max_chi_e_absolute_spread=" << metrics.max_chi_e_absolute_spread
         << "; max_chi_e_absolute_spread_radial_index="
         << metrics.max_chi_e_absolute_spread_radial_index
         << "; max_pressure_angular_relative_spread="
         << metrics.max_pressure_angular_relative_spread
         << "; max_pressure_spread_radial_index="
         << metrics.max_pressure_spread_radial_index
         << "; max_pressure_absolute_spread=" << metrics.max_pressure_absolute_spread
         << "; max_pressure_absolute_spread_radial_index="
         << metrics.max_pressure_absolute_spread_radial_index
         << "; max_abs_mom_theta=" << metrics.max_abs_mom_theta
         << "; max_abs_mom_theta_index="
         << metrics.max_abs_mom_theta_radial << ','
         << metrics.max_abs_mom_theta_theta << ','
         << metrics.max_abs_mom_theta_phi
         << "; max_abs_mom_phi=" << metrics.max_abs_mom_phi
         << "; max_abs_mom_phi_index="
         << metrics.max_abs_mom_phi_radial << ','
         << metrics.max_abs_mom_phi_theta << ','
         << metrics.max_abs_mom_phi_phi
         << "; max_unweighted_tangential_radial_index="
         << metrics.max_unweighted_tangential_radial_index
         << "; max_unweighted_tangential_radial_sum="
         << metrics.max_unweighted_tangential_radial_sum
         << "; tangential_abs_momentum=" << metrics.tangential_abs_momentum
         << "; radial_abs_momentum=" << metrics.radial_abs_momentum
         << "; tangential_to_radial_momentum_ratio="
         << metrics.tangential_to_radial_momentum_ratio
         << "; unweighted_tangential_abs_momentum="
         << metrics.unweighted_tangential_abs_momentum
         << "; unweighted_radial_abs_momentum="
         << metrics.unweighted_radial_abs_momentum
         << "; unweighted_tangential_to_radial_momentum_ratio="
         << metrics.unweighted_tangential_to_radial_momentum_ratio;
  return report.str();
}

void AppendPrimitiveReport(
    std::ostringstream& report,
    std::string_view prefix,
    const DirectionalPrimitiveState& primitive) {
  report << "; " << prefix << "_rho=" << primitive.rho
         << "; " << prefix << "_v_n=" << primitive.v_n
         << "; " << prefix << "_v_t1=" << primitive.v_t1
         << "; " << prefix << "_v_t2=" << primitive.v_t2
         << "; " << prefix << "_pressure=" << primitive.pressure
         << "; " << prefix << "_chi_e=" << primitive.chi_e;
}

void AppendMacroCellReport(
    std::ostringstream& report,
    std::string_view prefix,
    const dec3d::mesh::MacroZoneCell& cell) {
  report << "; " << prefix << "_radial=" << cell.radial_index
         << "; " << prefix << "_theta_begin=" << cell.theta_begin
         << "; " << prefix << "_theta_end=" << cell.theta_end
         << "; " << prefix << "_phi_begin=" << cell.phi_begin
         << "; " << prefix << "_phi_end=" << cell.phi_end
         << "; " << prefix << "_theta_factor=" << cell.theta_factor
         << "; " << prefix << "_phi_factor=" << cell.phi_factor
         << "; " << prefix << "_volume=" << cell.total_volume;
}

[[nodiscard]] AngularStageMetrics BuildFineAngularStageMetrics(
    const HydroStateSnapshot& snapshot,
    const std::vector<HydroConservativeState>& accumulated_delta,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  AngularStageMetrics metrics;
  if (!snapshot.is_complete() ||
      accumulated_delta.size() != snapshot.cells.size() ||
      !geometry.is_valid()) {
    return metrics;
  }

  std::vector<double> rho_min(
      snapshot.shape.radial_cells,
      std::numeric_limits<double>::infinity());
  std::vector<double> rho_max(
      snapshot.shape.radial_cells,
      -std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_min(
      snapshot.shape.radial_cells,
      std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_max(
      snapshot.shape.radial_cells,
      -std::numeric_limits<double>::infinity());
  std::vector<double> e_total_min(
      snapshot.shape.radial_cells,
      std::numeric_limits<double>::infinity());
  std::vector<double> e_total_max(
      snapshot.shape.radial_cells,
      -std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_min(
      snapshot.shape.radial_cells,
      std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_max(
      snapshot.shape.radial_cells,
      -std::numeric_limits<double>::infinity());
  std::vector<double> rho_volume_sum(snapshot.shape.radial_cells, 0.0);
  std::vector<double> pressure_min(
      snapshot.shape.radial_cells,
      std::numeric_limits<double>::infinity());
  std::vector<double> pressure_max(
      snapshot.shape.radial_cells,
      -std::numeric_limits<double>::infinity());
  std::vector<double> pressure_volume_sum(snapshot.shape.radial_cells, 0.0);
  std::vector<double> volume_sum(snapshot.shape.radial_cells, 0.0);
  std::vector<double> unweighted_tangential_radial_sum(snapshot.shape.radial_cells, 0.0);

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(snapshot.shape, radial, theta, phi);
        const auto state = Add(snapshot.cells[index], accumulated_delta[index]);
        const auto primitive = RecoverPrimitiveState(state);
        const double volume = CellVolume(geometry, snapshot.shape, radial, theta, phi);
        const double abs_mom_theta = std::abs(state.mom_theta);
        const double abs_mom_phi = std::abs(state.mom_phi);
        rho_min[radial] = std::min(rho_min[radial], state.rho);
        rho_max[radial] = std::max(rho_max[radial], state.rho);
        mom_r_min[radial] = std::min(mom_r_min[radial], state.mom_r);
        mom_r_max[radial] = std::max(mom_r_max[radial], state.mom_r);
        e_total_min[radial] = std::min(e_total_min[radial], state.e_fluid_total);
        e_total_max[radial] = std::max(e_total_max[radial], state.e_fluid_total);
        chi_e_min[radial] = std::min(chi_e_min[radial], state.chi_e);
        chi_e_max[radial] = std::max(chi_e_max[radial], state.chi_e);
        rho_volume_sum[radial] += state.rho * volume;
        if (primitive.is_physical()) {
          pressure_min[radial] = std::min(pressure_min[radial], primitive.pressure);
          pressure_max[radial] = std::max(pressure_max[radial], primitive.pressure);
          pressure_volume_sum[radial] += primitive.pressure * volume;
        }
        volume_sum[radial] += volume;
        if (abs_mom_theta > metrics.max_abs_mom_theta) {
          metrics.max_abs_mom_theta = abs_mom_theta;
          metrics.max_abs_mom_theta_radial = radial;
          metrics.max_abs_mom_theta_theta = theta;
          metrics.max_abs_mom_theta_phi = phi;
        }
        if (abs_mom_phi > metrics.max_abs_mom_phi) {
          metrics.max_abs_mom_phi = abs_mom_phi;
          metrics.max_abs_mom_phi_radial = radial;
          metrics.max_abs_mom_phi_theta = theta;
          metrics.max_abs_mom_phi_phi = phi;
        }
        metrics.tangential_abs_momentum +=
            (abs_mom_theta + abs_mom_phi) * volume;
        metrics.radial_abs_momentum += std::abs(state.mom_r) * volume;
        metrics.unweighted_tangential_abs_momentum += abs_mom_theta + abs_mom_phi;
        metrics.unweighted_radial_abs_momentum += std::abs(state.mom_r);
        unweighted_tangential_radial_sum[radial] += abs_mom_theta + abs_mom_phi;
      }
    }
  }

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    if (unweighted_tangential_radial_sum[radial] >
        metrics.max_unweighted_tangential_radial_sum) {
      metrics.max_unweighted_tangential_radial_sum =
          unweighted_tangential_radial_sum[radial];
      metrics.max_unweighted_tangential_radial_index = radial;
    }
    if (!(volume_sum[radial] > 0.0)) {
      continue;
    }
    const auto update_absolute_spread =
        [&](double spread, double& max_spread, std::size_t& max_radial_index) {
          if (spread > max_spread) {
            max_spread = spread;
            max_radial_index = radial;
          }
        };
    update_absolute_spread(
        rho_max[radial] - rho_min[radial],
        metrics.max_rho_absolute_spread,
        metrics.max_rho_absolute_spread_radial_index);
    update_absolute_spread(
        mom_r_max[radial] - mom_r_min[radial],
        metrics.max_mom_r_absolute_spread,
        metrics.max_mom_r_absolute_spread_radial_index);
    update_absolute_spread(
        e_total_max[radial] - e_total_min[radial],
        metrics.max_e_total_absolute_spread,
        metrics.max_e_total_absolute_spread_radial_index);
    update_absolute_spread(
        chi_e_max[radial] - chi_e_min[radial],
        metrics.max_chi_e_absolute_spread,
        metrics.max_chi_e_absolute_spread_radial_index);
    const double mean_rho = rho_volume_sum[radial] / volume_sum[radial];
    const double spread = std::abs(mean_rho) > 0.0
                              ? (rho_max[radial] - rho_min[radial]) / std::abs(mean_rho)
                              : 0.0;
    if (spread > metrics.max_rho_angular_relative_spread) {
      metrics.max_rho_angular_relative_spread = spread;
      metrics.max_spread_radial_index = radial;
    }
    const double mean_pressure = pressure_volume_sum[radial] / volume_sum[radial];
    const double pressure_spread =
        std::abs(mean_pressure) > 0.0
            ? (pressure_max[radial] - pressure_min[radial]) / std::abs(mean_pressure)
            : 0.0;
    update_absolute_spread(
        pressure_max[radial] - pressure_min[radial],
        metrics.max_pressure_absolute_spread,
        metrics.max_pressure_absolute_spread_radial_index);
    if (pressure_spread > metrics.max_pressure_angular_relative_spread) {
      metrics.max_pressure_angular_relative_spread = pressure_spread;
      metrics.max_pressure_spread_radial_index = radial;
    }
  }

  metrics.tangential_to_radial_momentum_ratio =
      metrics.radial_abs_momentum > 0.0
          ? metrics.tangential_abs_momentum / metrics.radial_abs_momentum
          : 0.0;
  metrics.unweighted_tangential_to_radial_momentum_ratio =
      metrics.unweighted_radial_abs_momentum > 0.0
          ? metrics.unweighted_tangential_abs_momentum /
                metrics.unweighted_radial_abs_momentum
          : 0.0;
  metrics.valid = true;
  return metrics;
}

[[nodiscard]] AngularStageMetrics BuildFineAngularStageMetricsFromView(
    const dec3d::state::HydroStateView& hydro_view,
    const GridShape& shape,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  AngularStageMetrics metrics;
  if (!hydro_view.is_complete() || !geometry.is_valid()) {
    return metrics;
  }

  std::vector<double> rho_min(shape.radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> rho_max(shape.radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_min(shape.radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_max(shape.radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> e_total_min(shape.radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> e_total_max(shape.radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_min(shape.radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_max(shape.radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> rho_volume_sum(shape.radial_cells, 0.0);
  std::vector<double> pressure_min(shape.radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> pressure_max(shape.radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> pressure_volume_sum(shape.radial_cells, 0.0);
  std::vector<double> volume_sum(shape.radial_cells, 0.0);
  std::vector<double> unweighted_tangential_radial_sum(shape.radial_cells, 0.0);

  for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
        const auto state = LoadCellState(hydro_view, radial, theta, phi);
        const auto primitive = RecoverPrimitiveState(state);
        const double volume = CellVolume(geometry, shape, radial, theta, phi);
        const double abs_mom_theta = std::abs(state.mom_theta);
        const double abs_mom_phi = std::abs(state.mom_phi);
        rho_min[radial] = std::min(rho_min[radial], state.rho);
        rho_max[radial] = std::max(rho_max[radial], state.rho);
        mom_r_min[radial] = std::min(mom_r_min[radial], state.mom_r);
        mom_r_max[radial] = std::max(mom_r_max[radial], state.mom_r);
        e_total_min[radial] = std::min(e_total_min[radial], state.e_fluid_total);
        e_total_max[radial] = std::max(e_total_max[radial], state.e_fluid_total);
        chi_e_min[radial] = std::min(chi_e_min[radial], state.chi_e);
        chi_e_max[radial] = std::max(chi_e_max[radial], state.chi_e);
        rho_volume_sum[radial] += state.rho * volume;
        if (primitive.is_physical()) {
          pressure_min[radial] = std::min(pressure_min[radial], primitive.pressure);
          pressure_max[radial] = std::max(pressure_max[radial], primitive.pressure);
          pressure_volume_sum[radial] += primitive.pressure * volume;
        }
        volume_sum[radial] += volume;
        if (abs_mom_theta > metrics.max_abs_mom_theta) {
          metrics.max_abs_mom_theta = abs_mom_theta;
          metrics.max_abs_mom_theta_radial = radial;
          metrics.max_abs_mom_theta_theta = theta;
          metrics.max_abs_mom_theta_phi = phi;
        }
        if (abs_mom_phi > metrics.max_abs_mom_phi) {
          metrics.max_abs_mom_phi = abs_mom_phi;
          metrics.max_abs_mom_phi_radial = radial;
          metrics.max_abs_mom_phi_theta = theta;
          metrics.max_abs_mom_phi_phi = phi;
        }
        metrics.tangential_abs_momentum +=
            (abs_mom_theta + abs_mom_phi) * volume;
        metrics.radial_abs_momentum += std::abs(state.mom_r) * volume;
        metrics.unweighted_tangential_abs_momentum += abs_mom_theta + abs_mom_phi;
        metrics.unweighted_radial_abs_momentum += std::abs(state.mom_r);
        unweighted_tangential_radial_sum[radial] += abs_mom_theta + abs_mom_phi;
      }
    }
  }

  for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
    if (unweighted_tangential_radial_sum[radial] >
        metrics.max_unweighted_tangential_radial_sum) {
      metrics.max_unweighted_tangential_radial_sum =
          unweighted_tangential_radial_sum[radial];
      metrics.max_unweighted_tangential_radial_index = radial;
    }
    if (!(volume_sum[radial] > 0.0)) {
      continue;
    }
    const auto update_absolute_spread =
        [&](double spread, double& max_spread, std::size_t& max_radial_index) {
          if (spread > max_spread) {
            max_spread = spread;
            max_radial_index = radial;
          }
        };
    update_absolute_spread(
        rho_max[radial] - rho_min[radial],
        metrics.max_rho_absolute_spread,
        metrics.max_rho_absolute_spread_radial_index);
    update_absolute_spread(
        mom_r_max[radial] - mom_r_min[radial],
        metrics.max_mom_r_absolute_spread,
        metrics.max_mom_r_absolute_spread_radial_index);
    update_absolute_spread(
        e_total_max[radial] - e_total_min[radial],
        metrics.max_e_total_absolute_spread,
        metrics.max_e_total_absolute_spread_radial_index);
    update_absolute_spread(
        chi_e_max[radial] - chi_e_min[radial],
        metrics.max_chi_e_absolute_spread,
        metrics.max_chi_e_absolute_spread_radial_index);
    const double mean_rho = rho_volume_sum[radial] / volume_sum[radial];
    const double spread = std::abs(mean_rho) > 0.0
                              ? (rho_max[radial] - rho_min[radial]) / std::abs(mean_rho)
                              : 0.0;
    if (spread > metrics.max_rho_angular_relative_spread) {
      metrics.max_rho_angular_relative_spread = spread;
      metrics.max_spread_radial_index = radial;
    }
    const double mean_pressure = pressure_volume_sum[radial] / volume_sum[radial];
    const double pressure_spread =
        std::abs(mean_pressure) > 0.0
            ? (pressure_max[radial] - pressure_min[radial]) / std::abs(mean_pressure)
            : 0.0;
    update_absolute_spread(
        pressure_max[radial] - pressure_min[radial],
        metrics.max_pressure_absolute_spread,
        metrics.max_pressure_absolute_spread_radial_index);
    if (pressure_spread > metrics.max_pressure_angular_relative_spread) {
      metrics.max_pressure_angular_relative_spread = pressure_spread;
      metrics.max_pressure_spread_radial_index = radial;
    }
  }

  metrics.tangential_to_radial_momentum_ratio =
      metrics.radial_abs_momentum > 0.0
          ? metrics.tangential_abs_momentum / metrics.radial_abs_momentum
          : 0.0;
  metrics.unweighted_tangential_to_radial_momentum_ratio =
      metrics.unweighted_radial_abs_momentum > 0.0
          ? metrics.unweighted_tangential_abs_momentum /
                metrics.unweighted_radial_abs_momentum
          : 0.0;
  metrics.valid = true;
  return metrics;
}

[[nodiscard]] AngularStageMetrics BuildCoarseAngularStageMetrics(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const std::vector<HydroConservativeState>& coarse_extensive_delta) noexcept {
  AngularStageMetrics metrics;
  if (!map.is_complete() ||
      coarse_states.size() != map.coarse_cells.size() ||
      coarse_extensive_delta.size() != map.coarse_cells.size()) {
    return metrics;
  }

  std::vector<double> rho_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> rho_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> e_total_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> e_total_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> rho_volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> pressure_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> pressure_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> pressure_volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> unweighted_tangential_radial_sum(map.fine_radial_cells, 0.0);

  for (std::size_t coarse_index = 0; coarse_index < map.coarse_cells.size(); ++coarse_index) {
    const auto& cell = map.coarse_cells[coarse_index];
    const double volume = cell.total_volume;
    if (!(volume > 0.0)) {
      continue;
    }
    const auto state = Add(coarse_states[coarse_index],
                           Scale(coarse_extensive_delta[coarse_index], 1.0 / volume));
    const auto primitive = RecoverPrimitiveState(state);
    const std::size_t radial = cell.radial_index;
    const double abs_mom_theta = std::abs(state.mom_theta);
    const double abs_mom_phi = std::abs(state.mom_phi);
    rho_min[radial] = std::min(rho_min[radial], state.rho);
    rho_max[radial] = std::max(rho_max[radial], state.rho);
    mom_r_min[radial] = std::min(mom_r_min[radial], state.mom_r);
    mom_r_max[radial] = std::max(mom_r_max[radial], state.mom_r);
    e_total_min[radial] = std::min(e_total_min[radial], state.e_fluid_total);
    e_total_max[radial] = std::max(e_total_max[radial], state.e_fluid_total);
    chi_e_min[radial] = std::min(chi_e_min[radial], state.chi_e);
    chi_e_max[radial] = std::max(chi_e_max[radial], state.chi_e);
    rho_volume_sum[radial] += state.rho * volume;
    if (primitive.is_physical()) {
      pressure_min[radial] = std::min(pressure_min[radial], primitive.pressure);
      pressure_max[radial] = std::max(pressure_max[radial], primitive.pressure);
      pressure_volume_sum[radial] += primitive.pressure * volume;
    }
    volume_sum[radial] += volume;
    if (abs_mom_theta > metrics.max_abs_mom_theta) {
      metrics.max_abs_mom_theta = abs_mom_theta;
      metrics.max_abs_mom_theta_radial = radial;
      metrics.max_abs_mom_theta_theta = cell.theta_begin;
      metrics.max_abs_mom_theta_phi = cell.phi_begin;
    }
    if (abs_mom_phi > metrics.max_abs_mom_phi) {
      metrics.max_abs_mom_phi = abs_mom_phi;
      metrics.max_abs_mom_phi_radial = radial;
      metrics.max_abs_mom_phi_theta = cell.theta_begin;
      metrics.max_abs_mom_phi_phi = cell.phi_begin;
    }
    metrics.tangential_abs_momentum +=
        (abs_mom_theta + abs_mom_phi) * volume;
    metrics.radial_abs_momentum += std::abs(state.mom_r) * volume;
    metrics.unweighted_tangential_abs_momentum += abs_mom_theta + abs_mom_phi;
    metrics.unweighted_radial_abs_momentum += std::abs(state.mom_r);
    unweighted_tangential_radial_sum[radial] += abs_mom_theta + abs_mom_phi;
  }

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    if (unweighted_tangential_radial_sum[radial] >
        metrics.max_unweighted_tangential_radial_sum) {
      metrics.max_unweighted_tangential_radial_sum =
          unweighted_tangential_radial_sum[radial];
      metrics.max_unweighted_tangential_radial_index = radial;
    }
    if (!(volume_sum[radial] > 0.0)) {
      continue;
    }
    const auto update_absolute_spread =
        [&](double spread, double& max_spread, std::size_t& max_radial_index) {
          if (spread > max_spread) {
            max_spread = spread;
            max_radial_index = radial;
          }
        };
    update_absolute_spread(
        rho_max[radial] - rho_min[radial],
        metrics.max_rho_absolute_spread,
        metrics.max_rho_absolute_spread_radial_index);
    update_absolute_spread(
        mom_r_max[radial] - mom_r_min[radial],
        metrics.max_mom_r_absolute_spread,
        metrics.max_mom_r_absolute_spread_radial_index);
    update_absolute_spread(
        e_total_max[radial] - e_total_min[radial],
        metrics.max_e_total_absolute_spread,
        metrics.max_e_total_absolute_spread_radial_index);
    update_absolute_spread(
        chi_e_max[radial] - chi_e_min[radial],
        metrics.max_chi_e_absolute_spread,
        metrics.max_chi_e_absolute_spread_radial_index);
    const double mean_rho = rho_volume_sum[radial] / volume_sum[radial];
    const double spread = std::abs(mean_rho) > 0.0
                              ? (rho_max[radial] - rho_min[radial]) / std::abs(mean_rho)
                              : 0.0;
    if (spread > metrics.max_rho_angular_relative_spread) {
      metrics.max_rho_angular_relative_spread = spread;
      metrics.max_spread_radial_index = radial;
    }
    const double mean_pressure = pressure_volume_sum[radial] / volume_sum[radial];
    const double pressure_spread =
        std::abs(mean_pressure) > 0.0
            ? (pressure_max[radial] - pressure_min[radial]) / std::abs(mean_pressure)
            : 0.0;
    update_absolute_spread(
        pressure_max[radial] - pressure_min[radial],
        metrics.max_pressure_absolute_spread,
        metrics.max_pressure_absolute_spread_radial_index);
    if (pressure_spread > metrics.max_pressure_angular_relative_spread) {
      metrics.max_pressure_angular_relative_spread = pressure_spread;
      metrics.max_pressure_spread_radial_index = radial;
    }
  }

  metrics.tangential_to_radial_momentum_ratio =
      metrics.radial_abs_momentum > 0.0
          ? metrics.tangential_abs_momentum / metrics.radial_abs_momentum
          : 0.0;
  metrics.unweighted_tangential_to_radial_momentum_ratio =
      metrics.unweighted_radial_abs_momentum > 0.0
          ? metrics.unweighted_tangential_abs_momentum /
                metrics.unweighted_radial_abs_momentum
          : 0.0;
  metrics.valid = true;
  return metrics;
}

void AppendFineAngularStageDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string_view stage,
    const HydroStateSnapshot& snapshot,
    const std::vector<HydroConservativeState>& accumulated_delta,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.debug.angular.stage",
      FormatAngularStageMetrics(
          stage,
          "fine_snapshot_plus_accumulated_delta",
          BuildFineAngularStageMetrics(snapshot, accumulated_delta, geometry)));
}

void AppendFineAngularStageDiagnosticFromView(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string_view stage,
    const dec3d::state::HydroStateView& hydro_view,
    const GridShape& shape,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.debug.angular.stage",
      FormatAngularStageMetrics(
          stage,
          "fine_hydro_view",
          BuildFineAngularStageMetricsFromView(hydro_view, shape, geometry)));
}

void AppendAleCommitDensityBalanceDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    const HydroStateSnapshot& snapshot,
    const std::vector<HydroConservativeState>& flux_delta_extensive,
    const std::vector<HydroConservativeState>& source_delta_extensive,
    const std::vector<HydroConservativeState>& radial_delta_per_solid_angle,
    const std::vector<HydroConservativeState>& radial_delta_extensive,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) {
  if (!snapshot.is_complete() ||
      flux_delta_extensive.size() != snapshot.cells.size() ||
      source_delta_extensive.size() != snapshot.cells.size() ||
      radial_delta_per_solid_angle.size() != snapshot.cells.size() ||
      radial_delta_extensive.size() != snapshot.cells.size() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.debug.ale.commit_density_balance",
        "stage=ale_commit_density_balance; valid=false");
    return;
  }

  struct ShellBalance {
    double old_rho_min{std::numeric_limits<double>::infinity()};
    double old_rho_max{-std::numeric_limits<double>::infinity()};
    double flux_density_min{std::numeric_limits<double>::infinity()};
    double flux_density_max{-std::numeric_limits<double>::infinity()};
    double source_density_min{std::numeric_limits<double>::infinity()};
    double source_density_max{-std::numeric_limits<double>::infinity()};
    double new_volume_min{std::numeric_limits<double>::infinity()};
    double new_volume_max{-std::numeric_limits<double>::infinity()};
    double updated_rho_min{std::numeric_limits<double>::infinity()};
    double updated_rho_max{-std::numeric_limits<double>::infinity()};
    double updated_rho_sum{0.0};
    std::size_t updated_rho_min_theta{0u};
    std::size_t updated_rho_min_phi{0u};
    std::size_t updated_rho_max_theta{0u};
    std::size_t updated_rho_max_phi{0u};
    std::size_t count{0u};
  };

  std::vector<ShellBalance> shells(snapshot.shape.radial_cells);
  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(snapshot.shape, radial, theta, phi);
        const double old_volume = CellVolume(old_geometry, snapshot.shape, radial, theta, phi);
        const double new_volume = CellVolume(new_geometry, snapshot.shape, radial, theta, phi);
        if (!(old_volume > 0.0) || !(new_volume > 0.0)) {
          continue;
        }
        auto& shell = shells[radial];
        const double old_rho = snapshot.cells[index].rho;
        const double flux_density = flux_delta_extensive[index].rho / new_volume;
        const double source_density = source_delta_extensive[index].rho / new_volume;
        const auto updated_state = BuildAleFactorizedUpdatedState(
            snapshot.cells[index],
            flux_delta_extensive[index],
            source_delta_extensive[index],
            radial_delta_per_solid_angle[index],
            radial_delta_extensive[index],
            old_geometry,
            new_geometry,
            radial,
            theta,
            phi);
        const double updated_rho = updated_state.rho;

        shell.old_rho_min = std::min(shell.old_rho_min, old_rho);
        shell.old_rho_max = std::max(shell.old_rho_max, old_rho);
        shell.flux_density_min = std::min(shell.flux_density_min, flux_density);
        shell.flux_density_max = std::max(shell.flux_density_max, flux_density);
        shell.source_density_min = std::min(shell.source_density_min, source_density);
        shell.source_density_max = std::max(shell.source_density_max, source_density);
        shell.new_volume_min = std::min(shell.new_volume_min, new_volume);
        shell.new_volume_max = std::max(shell.new_volume_max, new_volume);
        if (updated_rho < shell.updated_rho_min) {
          shell.updated_rho_min = updated_rho;
          shell.updated_rho_min_theta = theta;
          shell.updated_rho_min_phi = phi;
        }
        if (updated_rho > shell.updated_rho_max) {
          shell.updated_rho_max = updated_rho;
          shell.updated_rho_max_theta = theta;
          shell.updated_rho_max_phi = phi;
        }
        shell.updated_rho_sum += updated_rho;
        ++shell.count;
      }
    }
  }

  std::size_t max_radial = 0u;
  double max_spread = 0.0;
  for (std::size_t radial = 0; radial < shells.size(); ++radial) {
    const auto& shell = shells[radial];
    if (shell.count == 0u) {
      continue;
    }
    const double mean = shell.updated_rho_sum / static_cast<double>(shell.count);
    const double spread =
        std::abs(mean) > 0.0 ? (shell.updated_rho_max - shell.updated_rho_min) / std::abs(mean) : 0.0;
    if (spread > max_spread) {
      max_spread = spread;
      max_radial = radial;
    }
  }

  const auto& shell = shells[max_radial];
  const double mean = shell.count > 0u
                          ? shell.updated_rho_sum / static_cast<double>(shell.count)
                          : 0.0;
  const auto relative_spread = [](double min_value, double max_value) {
    const double mean_value = 0.5 * (min_value + max_value);
    return std::abs(mean_value) > 0.0
               ? (max_value - min_value) / std::abs(mean_value)
               : 0.0;
  };

  std::ostringstream report;
  report << std::setprecision(17)
         << "stage=ale_commit_density_balance"
         << "; valid=" << (shell.count > 0u ? "true" : "false")
         << "; max_radial_index=" << max_radial
         << "; updated_rho_spread=" << max_spread
         << "; updated_rho_mean=" << mean
         << "; updated_rho_min=" << shell.updated_rho_min
         << "; updated_rho_min_theta_phi="
         << shell.updated_rho_min_theta << ',' << shell.updated_rho_min_phi
         << "; updated_rho_max=" << shell.updated_rho_max
         << "; updated_rho_max_theta_phi="
         << shell.updated_rho_max_theta << ',' << shell.updated_rho_max_phi
         << "; old_rho_min=" << shell.old_rho_min
         << "; old_rho_max=" << shell.old_rho_max
         << "; old_rho_relative_spread="
         << relative_spread(shell.old_rho_min, shell.old_rho_max)
         << "; flux_density_min=" << shell.flux_density_min
         << "; flux_density_max=" << shell.flux_density_max
         << "; flux_density_relative_spread="
         << relative_spread(shell.flux_density_min, shell.flux_density_max)
         << "; source_density_min=" << shell.source_density_min
         << "; source_density_max=" << shell.source_density_max
         << "; new_volume_min=" << shell.new_volume_min
         << "; new_volume_max=" << shell.new_volume_max
         << "; new_volume_relative_spread="
         << relative_spread(shell.new_volume_min, shell.new_volume_max);
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.debug.ale.commit_density_balance",
      report.str());
}

[[nodiscard]] const char* HllcActiveRegionName(
    HllcActiveRegion region) noexcept {
  switch (region) {
    case HllcActiveRegion::left_flux:
      return "left_flux";
    case HllcActiveRegion::left_star:
      return "left_star";
    case HllcActiveRegion::right_star:
      return "right_star";
    case HllcActiveRegion::right_flux:
      return "right_flux";
    case HllcActiveRegion::invalid:
    default:
      return "invalid";
  }
}

struct RadialAleFaceFluxDetail {
  std::size_t theta{0u};
  std::size_t phi{0u};
  double mass_flux{0.0};
  double face_speed{0.0};
  bool downgraded_to_first_order{false};
  bool left_profile_troubled{false};
  bool right_profile_troubled{false};
  bool left_trace_failed{false};
  bool right_trace_failed{false};
  int left_troubled_component{-1};
  int right_troubled_component{-1};
  double left_rho{0.0};
  double left_vn{0.0};
  double left_pressure{0.0};
  double right_rho{0.0};
  double right_vn{0.0};
  double right_pressure{0.0};
  double s_left{0.0};
  double s_star{0.0};
  double s_right{0.0};
  HllcActiveRegion moving_branch{HllcActiveRegion::invalid};
  PpmSmoothnessProbe left_component1_smoothness;
  PpmSmoothnessProbe right_component1_smoothness;
};

struct RadialAleFaceFluxBucket {
  bool valid{false};
  double mass_flux_min{std::numeric_limits<double>::infinity()};
  double mass_flux_max{-std::numeric_limits<double>::infinity()};
  double mass_flux_sum{0.0};
  RadialAleFaceFluxDetail min_detail;
  RadialAleFaceFluxDetail max_detail;
  std::size_t count{0u};
  std::size_t downgraded_count{0u};
  std::size_t left_troubled_count{0u};
  std::size_t right_troubled_count{0u};
  std::size_t trace_failed_count{0u};
  std::size_t hllc_failed_count{0u};
  std::size_t branch_left_flux_count{0u};
  std::size_t branch_left_star_count{0u};
  std::size_t branch_right_star_count{0u};
  std::size_t branch_right_flux_count{0u};
  std::size_t branch_invalid_count{0u};
  std::vector<std::size_t> line_count_by_theta;
  std::vector<std::size_t> downgraded_count_by_theta;
  std::vector<std::size_t> left_troubled_count_by_theta;
  std::vector<std::size_t> right_troubled_count_by_theta;
};

void CountMovingBranch(
    RadialAleFaceFluxBucket& bucket,
    HllcActiveRegion branch) noexcept {
  switch (branch) {
    case HllcActiveRegion::left_flux:
      ++bucket.branch_left_flux_count;
      break;
    case HllcActiveRegion::left_star:
      ++bucket.branch_left_star_count;
      break;
    case HllcActiveRegion::right_star:
      ++bucket.branch_right_star_count;
      break;
    case HllcActiveRegion::right_flux:
      ++bucket.branch_right_flux_count;
      break;
    case HllcActiveRegion::invalid:
    default:
      ++bucket.branch_invalid_count;
      break;
  }
}

[[nodiscard]] double RelativeSpreadFromMinMaxMean(
    double min_value,
    double max_value,
    double mean_value) noexcept {
  return std::abs(mean_value) > 0.0
             ? (max_value - min_value) / std::abs(mean_value)
             : 0.0;
}

void AccumulateRadialAleFaceFluxDebug(
    std::vector<RadialAleFaceFluxBucket>& buckets,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi,
    double face_speed,
    const PpmInterfaceState& interface_state,
    bool used_first_order,
    const HydroConservativeState& interface_flux) noexcept {
  if (radial_face >= buckets.size()) {
    return;
  }

  auto& bucket = buckets[radial_face];
  RadialAleFaceFluxDetail detail;
  detail.theta = theta;
  detail.phi = phi;
  detail.mass_flux = interface_flux.rho;
  detail.face_speed = face_speed;
  detail.downgraded_to_first_order = used_first_order;
  detail.left_profile_troubled = interface_state.left_profile_troubled;
  detail.right_profile_troubled = interface_state.right_profile_troubled;
  detail.left_trace_failed = interface_state.left_trace_failed;
  detail.right_trace_failed = interface_state.right_trace_failed;
  detail.left_troubled_component = interface_state.left_troubled_component;
  detail.right_troubled_component = interface_state.right_troubled_component;
  detail.left_rho = interface_state.left.rho;
  detail.left_vn = interface_state.left.v_n;
  detail.left_pressure = interface_state.left.pressure;
  detail.right_rho = interface_state.right.rho;
  detail.right_vn = interface_state.right.v_n;
  detail.right_pressure = interface_state.right.pressure;

  const auto left_state = MakeDirectionalConservativeState(interface_state.left);
  const auto right_state = MakeDirectionalConservativeState(interface_state.right);
  const auto hllc = SolveHllcRiemann(left_state, right_state);
  if (hllc.is_complete()) {
    detail.s_left = hllc.waves.s_left;
    detail.s_star = hllc.waves.s_star;
    detail.s_right = hllc.waves.s_right;
    const auto branch =
        SelectMovingInterfaceBranch(hllc, left_state, right_state, face_speed);
    detail.moving_branch = branch.active_region;
  } else {
    ++bucket.hllc_failed_count;
  }
  detail.left_component1_smoothness =
      interface_state.left_component1_smoothness;
  detail.right_component1_smoothness =
      interface_state.right_component1_smoothness;

  bucket.valid = true;
  bucket.mass_flux_sum += detail.mass_flux;
  ++bucket.count;
  if (theta >= bucket.line_count_by_theta.size()) {
    bucket.line_count_by_theta.resize(theta + 1u, 0u);
    bucket.downgraded_count_by_theta.resize(theta + 1u, 0u);
    bucket.left_troubled_count_by_theta.resize(theta + 1u, 0u);
    bucket.right_troubled_count_by_theta.resize(theta + 1u, 0u);
  }
  ++bucket.line_count_by_theta[theta];
  if (detail.downgraded_to_first_order) {
    ++bucket.downgraded_count;
    ++bucket.downgraded_count_by_theta[theta];
  }
  if (detail.left_profile_troubled) {
    ++bucket.left_troubled_count;
    ++bucket.left_troubled_count_by_theta[theta];
  }
  if (detail.right_profile_troubled) {
    ++bucket.right_troubled_count;
    ++bucket.right_troubled_count_by_theta[theta];
  }
  if (detail.left_trace_failed || detail.right_trace_failed) {
    ++bucket.trace_failed_count;
  }
  CountMovingBranch(bucket, detail.moving_branch);

  if (detail.mass_flux < bucket.mass_flux_min) {
    bucket.mass_flux_min = detail.mass_flux;
    bucket.min_detail = detail;
  }
  if (detail.mass_flux > bucket.mass_flux_max) {
    bucket.mass_flux_max = detail.mass_flux;
    bucket.max_detail = detail;
  }
}

[[nodiscard]] std::string FormatRadialAleFaceDetail(
    std::string_view prefix,
    const RadialAleFaceFluxDetail& detail) {
  std::ostringstream report;
  report << std::setprecision(17)
         << prefix << "_theta_phi=" << detail.theta << ',' << detail.phi
         << "; " << prefix << "_mass_flux=" << detail.mass_flux
         << "; " << prefix << "_branch="
         << HllcActiveRegionName(detail.moving_branch)
         << "; " << prefix << "_downgraded="
         << (detail.downgraded_to_first_order ? "true" : "false")
         << "; " << prefix << "_left_troubled="
         << (detail.left_profile_troubled ? "true" : "false")
         << "; " << prefix << "_right_troubled="
         << (detail.right_profile_troubled ? "true" : "false")
         << "; " << prefix << "_left_component="
         << detail.left_troubled_component
         << "; " << prefix << "_right_component="
         << detail.right_troubled_component
         << "; " << prefix << "_left_rho=" << detail.left_rho
         << "; " << prefix << "_left_vn=" << detail.left_vn
         << "; " << prefix << "_left_pressure=" << detail.left_pressure
         << "; " << prefix << "_right_rho=" << detail.right_rho
         << "; " << prefix << "_right_vn=" << detail.right_vn
         << "; " << prefix << "_right_pressure=" << detail.right_pressure
         << "; " << prefix << "_s_left=" << detail.s_left
         << "; " << prefix << "_s_star=" << detail.s_star
         << "; " << prefix << "_s_right=" << detail.s_right
         << "; " << prefix << "_left_c1_reversal="
         << (detail.left_component1_smoothness.sign_reversal ? "true" : "false")
         << "; " << prefix << "_left_c1_reversed_slope="
         << detail.left_component1_smoothness.reversed_slope
         << "; " << prefix << "_left_c1_tolerance="
         << detail.left_component1_smoothness.tolerance
         << "; " << prefix << "_left_c1_primary_slope="
         << detail.left_component1_smoothness.primary_slope
         << "; " << prefix << "_left_c1_q_im1="
         << detail.left_component1_smoothness.q_im1
         << "; " << prefix << "_left_c1_q_i="
         << detail.left_component1_smoothness.q_i
         << "; " << prefix << "_left_c1_q_ip1="
         << detail.left_component1_smoothness.q_ip1
         << "; " << prefix << "_right_c1_reversal="
         << (detail.right_component1_smoothness.sign_reversal ? "true" : "false")
         << "; " << prefix << "_right_c1_reversed_slope="
         << detail.right_component1_smoothness.reversed_slope
         << "; " << prefix << "_right_c1_tolerance="
         << detail.right_component1_smoothness.tolerance
         << "; " << prefix << "_right_c1_primary_slope="
         << detail.right_component1_smoothness.primary_slope
         << "; " << prefix << "_right_c1_q_im1="
         << detail.right_component1_smoothness.q_im1
         << "; " << prefix << "_right_c1_q_i="
         << detail.right_component1_smoothness.q_i
         << "; " << prefix << "_right_c1_q_ip1="
         << detail.right_component1_smoothness.q_ip1;
  return report.str();
}

[[nodiscard]] std::string FormatRadialAleThetaCounts(
    std::string_view prefix,
    const std::vector<std::size_t>& numerators,
    const std::vector<std::size_t>& denominators) {
  std::ostringstream report;
  report << prefix << '=';
  for (std::size_t theta = 0; theta < denominators.size(); ++theta) {
    if (theta != 0u) {
      report << ',';
    }
    const std::size_t numerator =
        theta < numerators.size() ? numerators[theta] : 0u;
    report << theta << ':' << numerator << '/' << denominators[theta];
  }
  return report.str();
}

[[nodiscard]] std::string BuildRadialAleFaceFluxDebugReport(
    const std::vector<RadialAleFaceFluxBucket>& buckets,
    const std::vector<double>& face_speeds) {
  std::size_t max_face = 0u;
  double max_spread = 0.0;
  for (std::size_t face = 0; face < buckets.size(); ++face) {
    const auto& bucket = buckets[face];
    if (!bucket.valid || bucket.count == 0u) {
      continue;
    }
    const double mean =
        bucket.mass_flux_sum / static_cast<double>(bucket.count);
    const double spread = RelativeSpreadFromMinMaxMean(
        bucket.mass_flux_min,
        bucket.mass_flux_max,
        mean);
    if (spread > max_spread) {
      max_spread = spread;
      max_face = face;
    }
  }

  const auto& bucket = buckets[max_face];
  const double mean = bucket.count > 0u
                          ? bucket.mass_flux_sum / static_cast<double>(bucket.count)
                          : 0.0;
  std::ostringstream report;
  report << std::setprecision(17)
         << "stage=radial_ale_face_flux"
         << "; valid=" << (bucket.valid ? "true" : "false")
         << "; max_radial_face=" << max_face
         << "; face_speed="
         << (max_face < face_speeds.size() ? face_speeds[max_face] : 0.0)
         << "; mass_flux_relative_spread=" << max_spread
         << "; mass_flux_mean=" << mean
         << "; mass_flux_min=" << bucket.mass_flux_min
         << "; mass_flux_max=" << bucket.mass_flux_max
         << "; count=" << bucket.count
         << "; downgraded_count=" << bucket.downgraded_count
         << "; left_troubled_count=" << bucket.left_troubled_count
         << "; right_troubled_count=" << bucket.right_troubled_count
         << "; trace_failed_count=" << bucket.trace_failed_count
         << "; hllc_failed_count=" << bucket.hllc_failed_count
         << "; branch_left_flux_count=" << bucket.branch_left_flux_count
         << "; branch_left_star_count=" << bucket.branch_left_star_count
         << "; branch_right_star_count=" << bucket.branch_right_star_count
         << "; branch_right_flux_count=" << bucket.branch_right_flux_count
         << "; branch_invalid_count=" << bucket.branch_invalid_count
         << "; " << FormatRadialAleThetaCounts(
                "theta_downgraded_counts",
                bucket.downgraded_count_by_theta,
                bucket.line_count_by_theta)
         << "; " << FormatRadialAleThetaCounts(
                "theta_left_troubled_counts",
                bucket.left_troubled_count_by_theta,
                bucket.line_count_by_theta)
         << "; " << FormatRadialAleThetaCounts(
                "theta_right_troubled_counts",
                bucket.right_troubled_count_by_theta,
                bucket.line_count_by_theta)
         << "; " << FormatRadialAleFaceDetail("min", bucket.min_detail)
         << "; " << FormatRadialAleFaceDetail("max", bucket.max_detail);
  return report.str();
}

[[nodiscard]] double ExtensiveMass(
    const HydroConservativeState& state) noexcept {
  return state.rho;
}

[[nodiscard]] double ExtensiveRadialMomentum(
    const HydroConservativeState& state) noexcept {
  return state.mom_r;
}

[[nodiscard]] double ExtensiveFluidEnergy(
    const HydroConservativeState& state) noexcept {
  return state.e_fluid_total;
}

[[nodiscard]] double RadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  const double theta_lower = geometry.theta_faces[theta];
  const double theta_upper = geometry.theta_faces[theta + 1u];
  const double phi_lower = geometry.phi_faces[phi];
  const double phi_upper = geometry.phi_faces[phi + 1u];
  return radius * radius *
         (std::cos(theta_lower) - std::cos(theta_upper)) *
         (phi_upper - phi_lower);
}

[[nodiscard]] double RadialShellVolumeFactor(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return RadialShellVolumeFactorFromFaces(geometry.radial_faces, radial);
}

[[nodiscard]] double RadialFaceAreaPerSolidAngle(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  return radius * radius;
}

[[nodiscard]] double ThetaFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  const double phi_lower = geometry.phi_faces[phi];
  const double phi_upper = geometry.phi_faces[phi + 1u];
  return 0.5 *
         (radial_outer * radial_outer - radial_inner * radial_inner) *
         std::sin(geometry.theta_faces[theta_face]) *
         (phi_upper - phi_lower);
}

[[nodiscard]] double PhiFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  const double theta_lower = geometry.theta_faces[theta];
  const double theta_upper = geometry.theta_faces[theta + 1u];
  return 0.5 *
         (radial_outer * radial_outer - radial_inner * radial_inner) *
         (theta_upper - theta_lower);
}

[[nodiscard]] double CellCenteredRadius(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

[[nodiscard]] double CellCenteredTheta(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] HydroConservativeState ComputeDirectionalInterfaceFlux(
    const HydroConservativeState& left_state,
    const HydroConservativeState& right_state,
    SweepDirection direction) noexcept {
  const auto hllc = SolveHllcRiemann(
      ToDirectionalState(left_state, direction),
      ToDirectionalState(right_state, direction));
  if (!hllc.is_complete()) {
    return {};
  }

  return FromDirectionalFlux(hllc.interface_flux, direction);
}

[[nodiscard]] std::vector<HydroConservativeState> ExtractGhostedLine(
    const dec3d::core::Array3D<HydroConservativeState>& ghosted_states,
    SweepDirection direction,
    std::size_t fixed_a,
    std::size_t fixed_b) {
  std::vector<HydroConservativeState> line;
  if (direction == SweepDirection::radial) {
    line.reserve(ghosted_states.extent_r());
    for (std::size_t radial = 0; radial < ghosted_states.extent_r(); ++radial) {
      line.push_back(ghosted_states(radial, fixed_a, fixed_b));
    }
    return line;
  }

  if (direction == SweepDirection::theta) {
    line.reserve(ghosted_states.extent_theta());
    for (std::size_t theta = 0; theta < ghosted_states.extent_theta(); ++theta) {
      line.push_back(ghosted_states(fixed_a, theta, fixed_b));
    }
    return line;
  }

  line.reserve(ghosted_states.extent_phi());
  for (std::size_t phi = 0; phi < ghosted_states.extent_phi(); ++phi) {
    line.push_back(ghosted_states(fixed_a, fixed_b, phi));
  }
  return line;
}

[[nodiscard]] std::vector<double> ExtractLineEffectiveWidths(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const GridShape& shape,
    SweepDirection direction,
    std::size_t fixed_a,
    std::size_t fixed_b,
    std::size_t ghost_layers) {
  std::vector<double> widths;
  if (direction == SweepDirection::radial) {
    widths.resize(shape.radial_cells + (2u * ghost_layers), 0.0);
    for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
      widths[ghost_layers + radial] =
          geometry.radial_faces[radial + 1u] - geometry.radial_faces[radial];
    }
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      widths[ghost] = widths[ghost_layers];
      widths[ghost_layers + shape.radial_cells + ghost] =
          widths[ghost_layers + shape.radial_cells - 1u];
    }
    return widths;
  }

  if (direction == SweepDirection::theta) {
    widths.resize(shape.theta_cells + (2u * ghost_layers), 0.0);
    const double radius = CellCenteredRadius(geometry, fixed_a);
    for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
      widths[ghost_layers + theta] =
          radius * (geometry.theta_faces[theta + 1u] - geometry.theta_faces[theta]);
    }
    for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
      widths[ghost_layers - ghost] = widths[ghost_layers + ghost - 1u];
      widths[ghost_layers + shape.theta_cells + ghost - 1u] =
          widths[ghost_layers + shape.theta_cells - ghost];
    }
    return widths;
  }

  widths.resize(shape.phi_cells + (2u * ghost_layers), 0.0);
  const double radius = CellCenteredRadius(geometry, fixed_a);
  const double theta_center = CellCenteredTheta(geometry, fixed_b);
  for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
    widths[ghost_layers + phi] =
        radius * std::sin(theta_center) *
        (geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi]);
  }
  for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
    const std::size_t lower_source =
        (shape.phi_cells - (ghost % shape.phi_cells)) % shape.phi_cells;
    const std::size_t upper_source = (ghost - 1u) % shape.phi_cells;
    widths[ghost_layers - ghost] = widths[ghost_layers + lower_source];
    widths[ghost_layers + shape.phi_cells + ghost - 1u] =
        widths[ghost_layers + upper_source];
  }
  return widths;
}

[[nodiscard]] bool BuildPpmInterfaceFluxes(
    const std::vector<HydroConservativeState>& ghosted_line,
    SweepDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers,
    double dt_s,
    const std::vector<double>& effective_cell_widths,
    std::vector<HydroConservativeState>& interface_fluxes,
    std::size_t& downgraded_interface_count,
    std::string& failure_reason,
    const std::vector<double>* radial_face_speeds = nullptr,
    const std::vector<bool>* zero_flux_interfaces = nullptr,
    PpmLineReconstructionResult* reconstruction_out = nullptr) {
  const auto reconstruction = ReconstructPpmLine(
      ghosted_line,
      ToHydroDirection(direction),
      interior_cells,
      ghost_layers,
      dt_s,
      effective_cell_widths,
      direction == SweepDirection::radial ? radial_face_speeds : nullptr);
  if (!reconstruction.success) {
    failure_reason = reconstruction.failure_reason.empty()
                         ? "PPM line reconstruction failed"
                         : reconstruction.failure_reason;
    return false;
  }

  downgraded_interface_count += reconstruction.downgraded_interface_count;
  if (reconstruction_out != nullptr) {
    *reconstruction_out = reconstruction;
  }
  for (std::size_t interface_index = 0; interface_index < reconstruction.interfaces.size(); ++interface_index) {
    if (zero_flux_interfaces != nullptr &&
        interface_index < zero_flux_interfaces->size() &&
        (*zero_flux_interfaces)[interface_index]) {
      interface_fluxes[interface_index] = HydroConservativeState{};
      continue;
    }

    const auto& interface_state = reconstruction.interfaces[interface_index];
    const auto left_state = MakeDirectionalConservativeState(interface_state.left);
    const auto right_state = MakeDirectionalConservativeState(interface_state.right);
    const auto hllc = SolveHllcRiemann(left_state, right_state);
    if (!hllc.is_complete()) {
      std::ostringstream reason;
      reason << "PPM-reconstructed HLLC solve failed"
             << "; direction=" << SweepDirectionName(direction)
             << "; interface_index=" << interface_index;
      if (direction == SweepDirection::radial && radial_face_speeds != nullptr) {
        reason << "; radial_face_speed=" << (*radial_face_speeds)[interface_index];
      }
      if (!hllc.failure_reason.empty()) {
        reason << "; " << hllc.failure_reason;
      }
      failure_reason = reason.str();
      return false;
    }
    auto interface_flux = hllc.interface_flux;
    if (direction == SweepDirection::radial && radial_face_speeds != nullptr) {
      interface_flux = ComputeAleCorrectedFlux(
          hllc,
          left_state,
          right_state,
          (*radial_face_speeds)[interface_index]);
      if (!interface_flux.is_finite()) {
        failure_reason = "ALE-corrected radial PPM interface flux is non-finite";
        return false;
      }
    }
    interface_fluxes[interface_index] = FromDirectionalFlux(interface_flux, direction);
  }

  return true;
}

[[nodiscard]] bool BuildPpmInterfaceFluxesFromReconstruction(
    const PpmLineReconstructionResult& reconstruction,
    const std::vector<HydroConservativeState>& ghosted_line,
    SweepDirection direction,
    std::size_t ghost_layers,
    const std::vector<bool>& shared_fallback_mask,
    std::vector<HydroConservativeState>& interface_fluxes,
    std::size_t& downgraded_interface_count,
    std::string& failure_reason,
    const std::vector<double>* radial_face_speeds = nullptr,
    const std::vector<bool>* zero_flux_interfaces = nullptr,
    MacroAleHllcDiagnostics* macro_ale_hllc_diagnostics = nullptr,
    std::vector<double>* moving_interface_pressure_terms = nullptr) {
  if (!reconstruction.success) {
    failure_reason = reconstruction.failure_reason.empty()
                         ? "PPM line reconstruction failed"
                         : reconstruction.failure_reason;
    return false;
  }
  if (shared_fallback_mask.size() != reconstruction.interfaces.size()) {
    failure_reason = "shell-wide PPM fallback mask size does not match interface count";
    return false;
  }
  if (interface_fluxes.size() != reconstruction.interfaces.size()) {
    failure_reason = "PPM interface flux buffer size does not match reconstruction";
    return false;
  }
  if (moving_interface_pressure_terms != nullptr &&
      moving_interface_pressure_terms->size() != reconstruction.interfaces.size()) {
    failure_reason = "moving-interface pressure buffer size does not match reconstruction";
    return false;
  }
  if (moving_interface_pressure_terms != nullptr &&
      !(direction == SweepDirection::radial &&
        radial_face_speeds != nullptr &&
        macro_ale_hllc_diagnostics != nullptr)) {
    failure_reason =
        "moving-interface pressure terms require direct radial moving-face HLLC diagnostics";
    return false;
  }

  for (std::size_t interface_index = 0;
       interface_index < reconstruction.interfaces.size();
       ++interface_index) {
    if (zero_flux_interfaces != nullptr &&
        interface_index < zero_flux_interfaces->size() &&
        (*zero_flux_interfaces)[interface_index]) {
      interface_fluxes[interface_index] = HydroConservativeState{};
      if (moving_interface_pressure_terms != nullptr) {
        (*moving_interface_pressure_terms)[interface_index] = 0.0;
      }
      continue;
    }

    const auto& interface_state = reconstruction.interfaces[interface_index];
    const bool use_first_order =
        interface_state.downgraded_to_first_order ||
        shared_fallback_mask[interface_index];

    DirectionalPrimitiveState left_primitive = interface_state.left;
    DirectionalPrimitiveState right_primitive = interface_state.right;
    if (use_first_order) {
      const std::size_t left_cell = ghost_layers + interface_index - 1u;
      const std::size_t right_cell = ghost_layers + interface_index;
      if (right_cell >= ghosted_line.size()) {
        failure_reason = "PPM first-order fallback cell index is outside ghosted line";
        return false;
      }
      left_primitive = ToDirectionalPrimitive(
          ghosted_line[left_cell],
          ToHydroDirection(direction));
      right_primitive = ToDirectionalPrimitive(
          ghosted_line[right_cell],
          ToHydroDirection(direction));
      ++downgraded_interface_count;
    }

    const auto left_state = MakeDirectionalConservativeState(left_primitive);
    const auto right_state = MakeDirectionalConservativeState(right_primitive);
    const auto hllc = SolveHllcRiemann(left_state, right_state);
    if (!hllc.is_complete()) {
      std::ostringstream reason;
      reason << "PPM-reconstructed HLLC solve failed"
             << "; direction=" << SweepDirectionName(direction)
             << "; interface_index=" << interface_index
             << "; fallback_mask="
             << (use_first_order ? "true" : "false");
      if (direction == SweepDirection::radial && radial_face_speeds != nullptr) {
        reason << "; radial_face_speed=" << (*radial_face_speeds)[interface_index];
      }
      if (!hllc.failure_reason.empty()) {
        reason << "; " << hllc.failure_reason;
      }
      failure_reason = reason.str();
      return false;
    }

    auto interface_flux = hllc.interface_flux;
    if (direction == SweepDirection::radial && radial_face_speeds != nullptr) {
      const double face_speed = (*radial_face_speeds)[interface_index];
      if (macro_ale_hllc_diagnostics != nullptr) {
        const auto moving_flux = SolveMovingInterfaceHllcFlux(
            left_state,
            right_state,
            face_speed);
        if (!moving_flux.is_complete()) {
          std::ostringstream reason;
          reason << "direct moving-face HLLC flux assembly failed"
                 << "; direction=" << SweepDirectionName(direction)
                 << "; interface_index=" << interface_index
                 << "; radial_face_speed=" << face_speed;
          if (!moving_flux.failure_reason.empty()) {
            reason << "; " << moving_flux.failure_reason;
          }
          failure_reason = reason.str();
          macro_ale_hllc_diagnostics->failure_class =
              MacroAleHllcFailureClass::moving_hllc_branch_unclassified;
          ++macro_ale_hllc_diagnostics->branch_invalid_count;
          return false;
        }
        switch (moving_flux.active_region) {
          case HllcActiveRegion::left_flux:
            ++macro_ale_hllc_diagnostics->branch_left_flux_count;
            break;
          case HllcActiveRegion::left_star:
            ++macro_ale_hllc_diagnostics->branch_left_star_count;
            break;
          case HllcActiveRegion::right_star:
            ++macro_ale_hllc_diagnostics->branch_right_star_count;
            break;
          case HllcActiveRegion::right_flux:
            ++macro_ale_hllc_diagnostics->branch_right_flux_count;
            break;
          case HllcActiveRegion::invalid:
            ++macro_ale_hllc_diagnostics->branch_invalid_count;
            break;
        }
        if (moving_flux.active_region != moving_flux.static_zero_region) {
          ++macro_ale_hllc_diagnostics->branch_static_zero_mismatch_count;
        }
        macro_ale_hllc_diagnostics->max_abs_w_face =
            std::max(macro_ale_hllc_diagnostics->max_abs_w_face, std::abs(face_speed));
        interface_flux = moving_flux.flux;
        if (moving_interface_pressure_terms != nullptr) {
          const double pressure = MovingInterfaceNormalPressure(
              moving_flux.flux,
              moving_flux.selected_state,
              face_speed);
          if (!(pressure > 0.0) || !std::isfinite(pressure)) {
            failure_reason =
                "direct moving-face HLLC pressure split produced non-physical interface pressure";
            macro_ale_hllc_diagnostics->failure_class =
                MacroAleHllcFailureClass::moving_hllc_nonphysical;
            return false;
          }
          (*moving_interface_pressure_terms)[interface_index] = pressure;
        }
      } else {
        interface_flux = ComputeAleCorrectedFlux(
            hllc,
            left_state,
            right_state,
            face_speed);
      }
      if (!interface_flux.is_finite()) {
        failure_reason = "ALE-corrected radial PPM interface flux is non-finite";
        if (macro_ale_hllc_diagnostics != nullptr) {
          macro_ale_hllc_diagnostics->failure_class =
              MacroAleHllcFailureClass::moving_hllc_nonphysical;
        }
        return false;
      }
    }
    interface_fluxes[interface_index] = FromDirectionalFlux(interface_flux, direction);
  }

  return true;
}

[[nodiscard]] HydroConservativeState ComputeDirectionalPhysicalFlux(
    const HydroConservativeState& canonical_state,
    SweepDirection direction) noexcept {
  return FromDirectionalFlux(
      ComputePhysicalFlux(ToDirectionalState(canonical_state, direction)),
      direction);
}

[[nodiscard]] bool IsPhysicalCellState(const HydroConservativeState& state) noexcept {
  const auto primitive = RecoverPrimitiveState(state);
  return state.is_finite() && primitive.is_physical();
}

[[nodiscard]] HydroDirection ToHydroDirection(SweepDirection direction) noexcept {
  switch (direction) {
    case SweepDirection::radial:
      return HydroDirection::radial;
    case SweepDirection::theta:
      return HydroDirection::theta;
    case SweepDirection::phi:
      return HydroDirection::phi;
  }

  return HydroDirection::radial;
}

[[nodiscard]] const char* SweepDirectionName(SweepDirection direction) noexcept {
  switch (direction) {
    case SweepDirection::radial:
      return "radial";
    case SweepDirection::theta:
      return "theta";
    case SweepDirection::phi:
      return "phi";
  }

  return "unknown";
}

[[nodiscard]] RadialLowerBoundaryMode InferRadialLowerBoundaryMode(
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  if (!geometry.radial_faces.empty() &&
      std::abs(geometry.radial_faces.front()) <= 1.0e-14) {
    return RadialLowerBoundaryMode::origin_remap;
  }
  return RadialLowerBoundaryMode::boundary_contract;
}

[[nodiscard]] HydroStateSnapshot CaptureSnapshot(
    const dec3d::state::HydroStateView& hydro_view) {
  HydroStateSnapshot snapshot;
  snapshot.shape = {
      hydro_view.rho->extent_r(),
      hydro_view.rho->extent_theta(),
      hydro_view.rho->extent_phi()};
  snapshot.cells.resize(snapshot.shape.radial_cells * snapshot.shape.theta_cells * snapshot.shape.phi_cells);

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        snapshot.cells[LinearIndex(snapshot.shape, radial, theta, phi)] =
            LoadCellState(hydro_view, radial, theta, phi);
      }
    }
  }

  return snapshot;
}

[[nodiscard]] dec3d::mesh::FineHydroPackageView BuildFineHydroPackageView(
    const dec3d::state::HydroStateView& hydro_view) noexcept {
  return {
      hydro_view.rho,
      hydro_view.mom_r,
      hydro_view.mom_theta,
      hydro_view.mom_phi,
      hydro_view.e_fluid_total,
      &hydro_view.chi_e};
}

[[nodiscard]] std::vector<HydroConservativeState> BuildMacroZoneStates(
    const dec3d::mesh::MacroZonedHydroPackage& coarse_package) {
  std::vector<HydroConservativeState> states(coarse_package.map.coarse_cells.size());
  for (std::size_t index = 0; index < states.size(); ++index) {
    states[index] = {
        coarse_package.rho[index],
        coarse_package.mom_r[index],
        coarse_package.mom_theta[index],
        coarse_package.mom_phi[index],
        coarse_package.e_fluid_total[index],
        coarse_package.chi_e[index]};
  }
  return states;
}

void StoreMacroZoneStates(
    const std::vector<HydroConservativeState>& states,
    dec3d::mesh::MacroZonedHydroPackage& coarse_package) {
  for (std::size_t index = 0; index < states.size(); ++index) {
    coarse_package.rho[index] = states[index].rho;
    coarse_package.mom_r[index] = states[index].mom_r;
    coarse_package.mom_theta[index] = states[index].mom_theta;
    coarse_package.mom_phi[index] = states[index].mom_phi;
    coarse_package.e_fluid_total[index] = states[index].e_fluid_total;
    coarse_package.chi_e[index] = states[index].chi_e;
  }
}

[[nodiscard]] std::vector<std::size_t> CoarseCellIndicesForBand(
    const dec3d::mesh::MacroZoneMap& map,
    std::size_t radial,
    std::size_t theta_begin,
    std::size_t theta_end) {
  std::vector<std::size_t> indices;
  if (radial >= map.fine_radial_cells ||
      theta_begin >= theta_end ||
      theta_end > map.fine_theta_cells) {
    return indices;
  }

  std::size_t phi = 0u;
  while (phi < map.fine_phi_cells) {
    const std::size_t coarse_index = map.fine_to_coarse(radial, theta_begin, phi);
    if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex() ||
        coarse_index >= map.coarse_cells.size()) {
      indices.clear();
      return indices;
    }
    const auto& cell = map.coarse_cells[coarse_index];
    if (cell.radial_index != radial ||
        cell.theta_begin != theta_begin ||
        cell.theta_end != theta_end ||
        cell.phi_begin != phi ||
        cell.phi_end <= phi ||
        cell.phi_end > map.fine_phi_cells) {
      indices.clear();
      return indices;
    }
    indices.push_back(coarse_index);
    phi = cell.phi_end;
  }
  return indices;
}

[[nodiscard]] double CoarsePhiFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta_begin,
    std::size_t theta_end) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  return 0.5 *
         (radial_outer * radial_outer - radial_inner * radial_inner) *
         (geometry.theta_faces[theta_end] - geometry.theta_faces[theta_begin]);
}

[[nodiscard]] double CoarseThetaFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi_begin,
    std::size_t phi_end) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  return 0.5 *
         (radial_outer * radial_outer - radial_inner * radial_inner) *
         std::sin(geometry.theta_faces[theta_face]) *
         (geometry.phi_faces[phi_end] - geometry.phi_faces[phi_begin]);
}

[[nodiscard]] double CoarseRadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face,
    std::size_t theta_begin,
    std::size_t theta_end,
    std::size_t phi_begin,
    std::size_t phi_end) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  return radius * radius *
         (std::cos(geometry.theta_faces[theta_begin]) - std::cos(geometry.theta_faces[theta_end])) *
         (geometry.phi_faces[phi_end] - geometry.phi_faces[phi_begin]);
}

[[nodiscard]] double CoarseCellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::MacroZoneCell& cell) noexcept {
  const double radial_inner = geometry.radial_faces[cell.radial_index];
  const double radial_outer = geometry.radial_faces[cell.radial_index + 1u];
  const double radial_factor =
      (radial_outer * radial_outer * radial_outer -
       radial_inner * radial_inner * radial_inner) /
      3.0;
  const double polar_factor =
      std::cos(geometry.theta_faces[cell.theta_begin]) -
      std::cos(geometry.theta_faces[cell.theta_end]);
  const double azimuthal_factor =
      geometry.phi_faces[cell.phi_end] - geometry.phi_faces[cell.phi_begin];
  return radial_factor * polar_factor * azimuthal_factor;
}

[[nodiscard]] AngularStageMetrics BuildCoarseStateAngularStageMetrics(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  AngularStageMetrics metrics;
  if (!map.is_complete() ||
      coarse_states.size() != map.coarse_cells.size() ||
      !geometry.is_valid()) {
    return metrics;
  }

  std::vector<double> rho_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> rho_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> mom_r_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> e_total_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> e_total_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> chi_e_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> rho_volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> pressure_min(map.fine_radial_cells, std::numeric_limits<double>::infinity());
  std::vector<double> pressure_max(map.fine_radial_cells, -std::numeric_limits<double>::infinity());
  std::vector<double> pressure_volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> volume_sum(map.fine_radial_cells, 0.0);
  std::vector<double> unweighted_tangential_radial_sum(map.fine_radial_cells, 0.0);

  for (std::size_t coarse_index = 0; coarse_index < map.coarse_cells.size(); ++coarse_index) {
    const auto& cell = map.coarse_cells[coarse_index];
    const double volume = CoarseCellVolume(geometry, cell);
    if (!(volume > 0.0) || !std::isfinite(volume)) {
      continue;
    }
    const auto& state = coarse_states[coarse_index];
    const auto primitive = RecoverPrimitiveState(state);
    const std::size_t radial = cell.radial_index;
    const double abs_mom_theta = std::abs(state.mom_theta);
    const double abs_mom_phi = std::abs(state.mom_phi);
    rho_min[radial] = std::min(rho_min[radial], state.rho);
    rho_max[radial] = std::max(rho_max[radial], state.rho);
    mom_r_min[radial] = std::min(mom_r_min[radial], state.mom_r);
    mom_r_max[radial] = std::max(mom_r_max[radial], state.mom_r);
    e_total_min[radial] = std::min(e_total_min[radial], state.e_fluid_total);
    e_total_max[radial] = std::max(e_total_max[radial], state.e_fluid_total);
    chi_e_min[radial] = std::min(chi_e_min[radial], state.chi_e);
    chi_e_max[radial] = std::max(chi_e_max[radial], state.chi_e);
    rho_volume_sum[radial] += state.rho * volume;
    if (primitive.is_physical()) {
      pressure_min[radial] = std::min(pressure_min[radial], primitive.pressure);
      pressure_max[radial] = std::max(pressure_max[radial], primitive.pressure);
      pressure_volume_sum[radial] += primitive.pressure * volume;
    }
    volume_sum[radial] += volume;
    if (abs_mom_theta > metrics.max_abs_mom_theta) {
      metrics.max_abs_mom_theta = abs_mom_theta;
      metrics.max_abs_mom_theta_radial = radial;
      metrics.max_abs_mom_theta_theta = cell.theta_begin;
      metrics.max_abs_mom_theta_phi = cell.phi_begin;
    }
    if (abs_mom_phi > metrics.max_abs_mom_phi) {
      metrics.max_abs_mom_phi = abs_mom_phi;
      metrics.max_abs_mom_phi_radial = radial;
      metrics.max_abs_mom_phi_theta = cell.theta_begin;
      metrics.max_abs_mom_phi_phi = cell.phi_begin;
    }
    metrics.tangential_abs_momentum +=
        (abs_mom_theta + abs_mom_phi) * volume;
    metrics.radial_abs_momentum += std::abs(state.mom_r) * volume;
    metrics.unweighted_tangential_abs_momentum += abs_mom_theta + abs_mom_phi;
    metrics.unweighted_radial_abs_momentum += std::abs(state.mom_r);
    unweighted_tangential_radial_sum[radial] += abs_mom_theta + abs_mom_phi;
  }

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    if (unweighted_tangential_radial_sum[radial] >
        metrics.max_unweighted_tangential_radial_sum) {
      metrics.max_unweighted_tangential_radial_sum =
          unweighted_tangential_radial_sum[radial];
      metrics.max_unweighted_tangential_radial_index = radial;
    }
    if (!(volume_sum[radial] > 0.0)) {
      continue;
    }
    const auto update_absolute_spread =
        [&](double spread, double& max_spread, std::size_t& max_radial_index) {
          if (spread > max_spread) {
            max_spread = spread;
            max_radial_index = radial;
          }
        };
    update_absolute_spread(
        rho_max[radial] - rho_min[radial],
        metrics.max_rho_absolute_spread,
        metrics.max_rho_absolute_spread_radial_index);
    update_absolute_spread(
        mom_r_max[radial] - mom_r_min[radial],
        metrics.max_mom_r_absolute_spread,
        metrics.max_mom_r_absolute_spread_radial_index);
    update_absolute_spread(
        e_total_max[radial] - e_total_min[radial],
        metrics.max_e_total_absolute_spread,
        metrics.max_e_total_absolute_spread_radial_index);
    update_absolute_spread(
        chi_e_max[radial] - chi_e_min[radial],
        metrics.max_chi_e_absolute_spread,
        metrics.max_chi_e_absolute_spread_radial_index);
    const double mean_rho = rho_volume_sum[radial] / volume_sum[radial];
    const double spread = std::abs(mean_rho) > 0.0
                              ? (rho_max[radial] - rho_min[radial]) / std::abs(mean_rho)
                              : 0.0;
    if (spread > metrics.max_rho_angular_relative_spread) {
      metrics.max_rho_angular_relative_spread = spread;
      metrics.max_spread_radial_index = radial;
    }
    const double mean_pressure = pressure_volume_sum[radial] / volume_sum[radial];
    const double pressure_spread =
        std::abs(mean_pressure) > 0.0
            ? (pressure_max[radial] - pressure_min[radial]) / std::abs(mean_pressure)
            : 0.0;
    update_absolute_spread(
        pressure_max[radial] - pressure_min[radial],
        metrics.max_pressure_absolute_spread,
        metrics.max_pressure_absolute_spread_radial_index);
    if (pressure_spread > metrics.max_pressure_angular_relative_spread) {
      metrics.max_pressure_angular_relative_spread = pressure_spread;
      metrics.max_pressure_spread_radial_index = radial;
    }
  }

  metrics.tangential_to_radial_momentum_ratio =
      metrics.radial_abs_momentum > 0.0
          ? metrics.tangential_abs_momentum / metrics.radial_abs_momentum
          : 0.0;
  metrics.unweighted_tangential_to_radial_momentum_ratio =
      metrics.unweighted_radial_abs_momentum > 0.0
          ? metrics.unweighted_tangential_abs_momentum /
                metrics.unweighted_radial_abs_momentum
          : 0.0;
  metrics.valid = true;
  return metrics;
}

[[nodiscard]] AngularStageMetrics BuildCoarseA4ExtensiveStageMetrics(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const std::vector<HydroConservativeState>& coarse_extensive_delta,
    const std::vector<HydroConservativeState>& radial_delta_per_solid_angle_reference,
    const std::vector<HydroConservativeState>& radial_delta_per_solid_angle_weighted_delta_sum,
    const std::vector<HydroConservativeState>& radial_delta_extensive,
    const std::vector<double>& radial_delta_solid_angle_weight,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) noexcept {
  std::vector<HydroConservativeState> predicted_states;
  if (!map.is_complete() ||
      coarse_states.size() != map.coarse_cells.size() ||
      coarse_extensive_delta.size() != map.coarse_cells.size() ||
      radial_delta_per_solid_angle_reference.size() != map.coarse_cells.size() ||
      radial_delta_per_solid_angle_weighted_delta_sum.size() != map.coarse_cells.size() ||
      radial_delta_extensive.size() != map.coarse_cells.size() ||
      radial_delta_solid_angle_weight.size() != map.coarse_cells.size() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return {};
  }

  predicted_states.resize(coarse_states.size());
  for (std::size_t coarse_index = 0; coarse_index < coarse_states.size(); ++coarse_index) {
    predicted_states[coarse_index] = BuildMacroRadialFactorizedCoarseState(
        coarse_states[coarse_index],
        coarse_extensive_delta[coarse_index],
        radial_delta_per_solid_angle_reference[coarse_index],
        radial_delta_per_solid_angle_weighted_delta_sum[coarse_index],
        radial_delta_extensive[coarse_index],
        radial_delta_solid_angle_weight[coarse_index],
        old_geometry,
        new_geometry,
        map.coarse_cells[coarse_index]);
    if (!IsPhysicalCellState(predicted_states[coarse_index])) {
      return {};
    }
  }

  return BuildCoarseStateAngularStageMetrics(map, predicted_states, new_geometry);
}

[[nodiscard]] double CoarsePhiCellWidth(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::MacroZoneCell& cell) noexcept {
  const double radial_center =
      0.5 * (geometry.radial_faces[cell.radial_index] + geometry.radial_faces[cell.radial_index + 1u]);
  const double theta_center =
      0.5 * (geometry.theta_faces[cell.theta_begin] + geometry.theta_faces[cell.theta_end]);
  return radial_center *
         std::sin(theta_center) *
         (geometry.phi_faces[cell.phi_end] - geometry.phi_faces[cell.phi_begin]);
}

[[nodiscard]] double CoarseThetaBandWidth(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::MacroZoneThetaBand& band) noexcept {
  const double radial_center =
      0.5 * (geometry.radial_faces[band.radial_index] + geometry.radial_faces[band.radial_index + 1u]);
  return radial_center *
         (geometry.theta_faces[band.theta_end] - geometry.theta_faces[band.theta_begin]);
}

[[nodiscard]] HydroConservativeState ComputeCoarseGeometricSource(
    const HydroConservativeState& coarse_state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::mesh::MacroZoneCell& cell) noexcept {
  HydroConservativeState source;
  const auto primitive = RecoverPrimitiveState(coarse_state);
  if (!primitive.is_physical() || !(cell.total_volume > 0.0)) {
    return source;
  }

  const double radial_area_difference =
      CoarseRadialFaceArea(
          geometry,
          cell.radial_index + 1u,
          cell.theta_begin,
          cell.theta_end,
          cell.phi_begin,
          cell.phi_end) -
      CoarseRadialFaceArea(
          geometry,
          cell.radial_index,
          cell.theta_begin,
          cell.theta_end,
          cell.phi_begin,
          cell.phi_end);
  const double theta_area_difference =
      CoarseThetaFaceArea(
          geometry,
          cell.radial_index,
          cell.theta_end,
          cell.phi_begin,
          cell.phi_end) -
      CoarseThetaFaceArea(
          geometry,
          cell.radial_index,
          cell.theta_begin,
          cell.phi_begin,
          cell.phi_end);
  const double inverse_radius_factor =
      0.5 * radial_area_difference / cell.total_volume;
  const double theta_curvature_factor =
      theta_area_difference / cell.total_volume;
  if (!(inverse_radius_factor > 0.0) ||
      !std::isfinite(inverse_radius_factor) ||
      !std::isfinite(theta_curvature_factor)) {
    return source;
  }

  const double radial_pressure_source =
      primitive.pressure * radial_area_difference / cell.total_volume;
  const double theta_pressure_source =
      primitive.pressure * theta_curvature_factor;

  source.mom_r =
      radial_pressure_source +
      ((coarse_state.mom_theta * coarse_state.mom_theta +
        coarse_state.mom_phi * coarse_state.mom_phi) /
       coarse_state.rho) *
          inverse_radius_factor;
  source.mom_theta =
      (-(coarse_state.mom_theta * coarse_state.mom_r) /
       coarse_state.rho) *
          inverse_radius_factor +
      (theta_curvature_factor *
       ((coarse_state.mom_phi * coarse_state.mom_phi) / coarse_state.rho)) +
      theta_pressure_source;
  source.mom_phi =
      (-(coarse_state.mom_r * coarse_state.mom_phi) /
       coarse_state.rho) *
          inverse_radius_factor -
      (theta_curvature_factor *
       ((coarse_state.mom_theta * coarse_state.mom_phi) / coarse_state.rho));
  return source;
}

[[nodiscard]] std::vector<HydroConservativeState> BuildCoarseGeometricSourceExtensiveDelta(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    bool exclude_balanced_pressure_terms = false) {
  std::vector<HydroConservativeState> source_extensive_delta(
      coarse_states.size(),
      HydroConservativeState{});
  if (!map.is_complete() || coarse_states.size() != map.coarse_cells.size()) {
    return source_extensive_delta;
  }

  for (std::size_t coarse_index = 0; coarse_index < coarse_states.size(); ++coarse_index) {
    const auto& cell = map.coarse_cells[coarse_index];
    const auto primitive = RecoverPrimitiveState(coarse_states[coarse_index]);
    if (!primitive.is_physical() || !(dt_s > 0.0) || !std::isfinite(dt_s)) {
      continue;
    }
    const double radial_area_difference =
        CoarseRadialFaceArea(
            geometry,
            cell.radial_index + 1u,
            cell.theta_begin,
            cell.theta_end,
            cell.phi_begin,
            cell.phi_end) -
        CoarseRadialFaceArea(
            geometry,
            cell.radial_index,
            cell.theta_begin,
            cell.theta_end,
            cell.phi_begin,
            cell.phi_end);
    const double theta_area_difference =
        CoarseThetaFaceArea(
            geometry,
            cell.radial_index,
            cell.theta_end,
            cell.phi_begin,
            cell.phi_end) -
        CoarseThetaFaceArea(
            geometry,
            cell.radial_index,
            cell.theta_begin,
            cell.phi_begin,
            cell.phi_end);
    if (!(radial_area_difference > 0.0) ||
        !std::isfinite(radial_area_difference) ||
        !std::isfinite(theta_area_difference)) {
      continue;
    }
    const double half_radial_area_difference = 0.5 * radial_area_difference;
    const double radial_pressure_source =
        exclude_balanced_pressure_terms
            ? 0.0
            : primitive.pressure * radial_area_difference;
    const double theta_pressure_source =
        exclude_balanced_pressure_terms
            ? 0.0
            : primitive.pressure * theta_area_difference;
    const double tangential_pressure =
        (coarse_states[coarse_index].mom_theta *
         coarse_states[coarse_index].mom_theta +
         coarse_states[coarse_index].mom_phi *
         coarse_states[coarse_index].mom_phi) /
        coarse_states[coarse_index].rho;
    const double phi_pressure =
        (coarse_states[coarse_index].mom_phi *
         coarse_states[coarse_index].mom_phi) /
        coarse_states[coarse_index].rho;
    HydroConservativeState source_extensive;
    source_extensive.mom_r =
        dt_s *
        (radial_pressure_source +
         tangential_pressure * half_radial_area_difference);
    source_extensive.mom_theta =
        dt_s *
        ((-(coarse_states[coarse_index].mom_theta *
            coarse_states[coarse_index].mom_r) /
          coarse_states[coarse_index].rho) *
             half_radial_area_difference +
         phi_pressure * theta_area_difference +
         theta_pressure_source);
    source_extensive.mom_phi =
        dt_s *
        ((-(coarse_states[coarse_index].mom_r *
            coarse_states[coarse_index].mom_phi) /
          coarse_states[coarse_index].rho) *
             half_radial_area_difference -
         ((coarse_states[coarse_index].mom_theta *
           coarse_states[coarse_index].mom_phi) /
          coarse_states[coarse_index].rho) *
             theta_area_difference);
    source_extensive_delta[coarse_index] = source_extensive;
  }
  return source_extensive_delta;
}

[[nodiscard]] std::string FormatCoarseGeometricSourceBalance(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const std::vector<HydroConservativeState>& flux_extensive_delta,
    const std::vector<HydroConservativeState>& source_extensive_delta) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "stage=macro_geometric_source_coarse_balance"
         << "; representation=macro_coarse_flux_source_delta"
         << "; valid=false";
  if (!map.is_complete() ||
      coarse_states.size() != map.coarse_cells.size() ||
      flux_extensive_delta.size() != coarse_states.size() ||
      source_extensive_delta.size() != coarse_states.size()) {
    return report.str();
  }

  bool valid = false;
  double max_abs_before_source_mom_theta = 0.0;
  double max_abs_after_source_mom_theta = 0.0;
  double max_abs_flux_update_mom_theta = 0.0;
  double max_abs_source_update_mom_theta = 0.0;
  double max_relative_after_source_mom_theta = 0.0;
  double max_after_old_mom_theta = 0.0;
  double max_after_flux_update_mom_theta = 0.0;
  double max_after_source_update_mom_theta = 0.0;
  double max_after_final_mom_theta = 0.0;
  std::size_t max_before_index = 0u;
  std::size_t max_after_index = 0u;
  std::size_t max_relative_index = 0u;

  for (std::size_t coarse_index = 0; coarse_index < coarse_states.size(); ++coarse_index) {
    const auto& cell = map.coarse_cells[coarse_index];
    const double volume = cell.total_volume;
    if (!(volume > 0.0) || !std::isfinite(volume)) {
      continue;
    }

    const auto flux_update = Scale(flux_extensive_delta[coarse_index], 1.0 / volume);
    const auto source_update = Scale(source_extensive_delta[coarse_index], 1.0 / volume);
    const auto before_source = Add(coarse_states[coarse_index], flux_update);
    const auto after_source = Add(before_source, source_update);
    const double abs_before = std::abs(before_source.mom_theta);
    const double abs_after = std::abs(after_source.mom_theta);
    const double abs_flux_update = std::abs(flux_update.mom_theta);
    const double abs_source_update = std::abs(source_update.mom_theta);
    const double cancellation_scale =
        std::abs(coarse_states[coarse_index].mom_theta) +
        abs_flux_update +
        abs_source_update;
    const double relative_after =
        cancellation_scale > 0.0 ? abs_after / cancellation_scale : 0.0;

    valid = true;
    if (abs_before > max_abs_before_source_mom_theta) {
      max_abs_before_source_mom_theta = abs_before;
      max_before_index = coarse_index;
    }
    if (abs_after > max_abs_after_source_mom_theta) {
      max_abs_after_source_mom_theta = abs_after;
      max_after_index = coarse_index;
      max_after_old_mom_theta = coarse_states[coarse_index].mom_theta;
      max_after_flux_update_mom_theta = flux_update.mom_theta;
      max_after_source_update_mom_theta = source_update.mom_theta;
      max_after_final_mom_theta = after_source.mom_theta;
    }
    if (abs_flux_update > max_abs_flux_update_mom_theta) {
      max_abs_flux_update_mom_theta = abs_flux_update;
    }
    if (abs_source_update > max_abs_source_update_mom_theta) {
      max_abs_source_update_mom_theta = abs_source_update;
    }
    if (relative_after > max_relative_after_source_mom_theta) {
      max_relative_after_source_mom_theta = relative_after;
      max_relative_index = coarse_index;
    }
  }

  const auto append_cell = [&](std::string_view prefix, std::size_t coarse_index) {
    const auto& cell = map.coarse_cells[coarse_index];
    report << "; " << prefix << "_coarse_index=" << coarse_index
           << "; " << prefix << "_radial=" << cell.radial_index
           << "; " << prefix << "_theta_begin=" << cell.theta_begin
           << "; " << prefix << "_theta_end=" << cell.theta_end
           << "; " << prefix << "_phi_begin=" << cell.phi_begin
           << "; " << prefix << "_phi_end=" << cell.phi_end
           << "; " << prefix << "_volume=" << cell.total_volume;
  };

  report.str(std::string{});
  report.clear();
  report << std::setprecision(17)
         << "stage=macro_geometric_source_coarse_balance"
         << "; representation=macro_coarse_flux_source_delta"
         << "; valid=" << (valid ? "true" : "false")
         << "; max_abs_before_source_mom_theta=" << max_abs_before_source_mom_theta
         << "; max_abs_after_source_mom_theta=" << max_abs_after_source_mom_theta
         << "; max_abs_flux_update_mom_theta=" << max_abs_flux_update_mom_theta
         << "; max_abs_source_update_mom_theta=" << max_abs_source_update_mom_theta
         << "; max_relative_after_source_mom_theta="
         << max_relative_after_source_mom_theta
         << "; max_after_old_mom_theta=" << max_after_old_mom_theta
         << "; max_after_flux_update_mom_theta=" << max_after_flux_update_mom_theta
         << "; max_after_source_update_mom_theta=" << max_after_source_update_mom_theta
         << "; max_after_final_mom_theta=" << max_after_final_mom_theta;
  if (valid) {
    append_cell("max_before", max_before_index);
    append_cell("max_after", max_after_index);
    append_cell("max_relative", max_relative_index);
  }
  return report.str();
}

[[nodiscard]] GeometricSourceStepResult ApplyMacroZonedGeometricSourceStep(
    dec3d::state::HydroStateView& hydro_view,
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double coarse_factor,
    double dt_s) noexcept {
  GeometricSourceStepResult result;
  if (!hydro_view.is_complete()) {
    result.failure_reason = "hydro view is incomplete for macro-zoned geometric source step";
  } else if (!source_snapshot.is_complete()) {
    result.failure_reason = "macro-zoned geometric source snapshot is incomplete";
  } else if (!geometry.is_valid()) {
    result.failure_reason = "spherical geometry metadata is invalid for macro-zoned geometric source step";
  } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
    result.failure_reason = "time step must be finite and positive for macro-zoned geometric source step";
  } else {
    const GridShape shape{
        source_snapshot.radial_cells,
        source_snapshot.theta_cells,
        source_snapshot.phi_cells};
    dec3d::core::Array3D<double> rho(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);
    dec3d::core::Array3D<double> mom_r(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);
    dec3d::core::Array3D<double> mom_theta(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);
    dec3d::core::Array3D<double> mom_phi(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);
    dec3d::core::Array3D<double> e_fluid_total(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);
    dec3d::core::Array3D<double> chi_e(shape.radial_cells, shape.theta_cells, shape.phi_cells, 0.0);

    for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
      for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
          const auto state = source_snapshot.states[LinearIndex(shape, radial, theta, phi)];
          rho(radial, theta, phi) = state.rho;
          mom_r(radial, theta, phi) = state.mom_r;
          mom_theta(radial, theta, phi) = state.mom_theta;
          mom_phi(radial, theta, phi) = state.mom_phi;
          e_fluid_total(radial, theta, phi) = state.e_fluid_total;
          chi_e(radial, theta, phi) = state.chi_e;
        }
      }
    }

    const auto map = dec3d::mesh::DetectMacroZones(
        geometry,
        dec3d::mesh::MacroZoningDetectionOptions{coarse_factor});
    if (!map.is_complete()) {
      result.failure_reason = map.failure_reason.empty()
                                  ? "macro-zoned geometric source detection failed"
                                  : map.failure_reason;
    } else {
      const dec3d::mesh::FineHydroPackageView fine_view{
          &rho,
          &mom_r,
          &mom_theta,
          &mom_phi,
          &e_fluid_total,
          &chi_e};
      const auto coarse_package = dec3d::mesh::RestrictFineHydroPackage(
          fine_view,
          geometry,
          map);
      if (!coarse_package.is_complete()) {
        result.failure_reason = coarse_package.failure_reason.empty()
                                    ? "macro-zoned geometric source restrict failed"
                                    : coarse_package.failure_reason;
      } else {
        const auto coarse_states = BuildMacroZoneStates(coarse_package);
        std::vector<HydroConservativeState> coarse_source_terms(
            coarse_states.size(),
            HydroConservativeState{});
        for (std::size_t coarse_index = 0; coarse_index < coarse_states.size(); ++coarse_index) {
          coarse_source_terms[coarse_index] = ComputeCoarseGeometricSource(
              coarse_states[coarse_index],
              geometry,
              map.coarse_cells[coarse_index]);
          if (!coarse_source_terms[coarse_index].is_finite()) {
            result.failure_reason = "macro-zoned coarse geometric source produced a non-finite source";
            break;
          }
        }

        for (std::size_t radial = 0; result.failure_reason.empty() && radial < shape.radial_cells; ++radial) {
          for (std::size_t theta = 0; result.failure_reason.empty() && theta < shape.theta_cells; ++theta) {
            for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
              const std::size_t coarse_index = map.fine_to_coarse(radial, theta, phi);
              const auto source = coarse_source_terms[coarse_index];
              HydroConservativeState updated = LoadCellState(hydro_view, radial, theta, phi);
              updated.mom_r += dt_s * source.mom_r;
              updated.mom_theta += dt_s * source.mom_theta;
              updated.mom_phi += dt_s * source.mom_phi;
              if (!updated.is_finite() || !RecoverPrimitiveState(updated).is_physical()) {
                result.failure_reason =
                    "macro-zoned geometric source step produced a non-physical hydro state";
                break;
              }
              (*hydro_view.mom_r)(radial, theta, phi) = updated.mom_r;
              (*hydro_view.mom_theta)(radial, theta, phi) = updated.mom_theta;
              (*hydro_view.mom_phi)(radial, theta, phi) = updated.mom_phi;
              ++result.updated_cell_count;
            }
          }
        }

        if (result.failure_reason.empty()) {
          result.success = result.updated_cell_count == source_snapshot.states.size();
        }
      }
    }
  }

  if (result.success) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.source.geometric_step",
        "hydro geometric source step used macro-zoned coarse-consistent source terms and conservative fine-cell distribution");
  } else {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.source.geometric_step.failed",
        result.failure_reason.empty()
            ? "macro-zoned geometric source step failed"
            : result.failure_reason);
  }

  std::ostringstream report;
  report << "geometric_source_step_success=" << (result.success ? "true" : "false")
         << "; macro_zoned_coarse_consistent=true"
         << "; updated_cell_count=" << result.updated_cell_count;
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();
  return result;
}

struct MacroGhostBootstrapSummary {
  bool executed{false};
  bool inner_consumed{false};
  bool outer_consumed{false};
  std::size_t ghost_layers{0};
  std::size_t inner_coarse_cells{0};
  std::size_t outer_coarse_cells{0};
  std::string inner_topology_source{"none"};
  std::string outer_topology_source{"none"};
  std::string report_line;
};

struct MacroCoarseGhostSide {
  bool consumed{false};
  dec3d::mesh::SphericalGeometryMetadata geometry;
  dec3d::mesh::MacroZoneMap map;
  std::vector<HydroConservativeState> states;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept {
    return consumed &&
           geometry.is_valid() &&
           map.is_complete() &&
           states.size() == map.coarse_cells.size();
  }
};

struct MacroGhostBootstrap {
  MacroGhostBootstrapSummary summary;
  MacroCoarseGhostSide inner;
  MacroCoarseGhostSide outer;
};

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildExtrapolatedGhostGeometry(
    const dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    std::size_t ghost_layers,
    bool inner_side) noexcept {
  dec3d::mesh::SphericalGeometryMetadata ghost_geometry;
  if (!local_geometry.is_valid() ||
      local_geometry.radial_faces.size() < 2u ||
      ghost_layers == 0u) {
    return ghost_geometry;
  }

  const double boundary_width = inner_side
                                    ? local_geometry.radial_faces[1u] - local_geometry.radial_faces[0u]
                                    : local_geometry.radial_faces.back() -
                                          local_geometry.radial_faces[local_geometry.radial_faces.size() - 2u];
  if (!(boundary_width > 0.0) || !std::isfinite(boundary_width)) {
    return ghost_geometry;
  }

  ghost_geometry.valid = true;
  ghost_geometry.radial_faces.resize(ghost_layers + 1u, 0.0);
  if (inner_side) {
    const double boundary_face = local_geometry.radial_faces.front();
    for (std::size_t face = 0; face <= ghost_layers; ++face) {
      ghost_geometry.radial_faces[face] =
          boundary_face - boundary_width * static_cast<double>(ghost_layers - face);
    }
  } else {
    const double boundary_face = local_geometry.radial_faces.back();
    for (std::size_t face = 0; face <= ghost_layers; ++face) {
      ghost_geometry.radial_faces[face] =
          boundary_face + boundary_width * static_cast<double>(face);
    }
  }

  for (double face : ghost_geometry.radial_faces) {
    if (!(face >= 0.0) || !std::isfinite(face)) {
      return {};
    }
  }

  ghost_geometry.theta_faces = local_geometry.theta_faces;
  ghost_geometry.phi_faces = local_geometry.phi_faces;
  ghost_geometry.cell_volumes.reserve(
      ghost_layers *
      (ghost_geometry.theta_faces.size() - 1u) *
      (ghost_geometry.phi_faces.size() - 1u));

  for (std::size_t radial = 0; radial < ghost_layers; ++radial) {
    const double radial_factor =
        (std::pow(ghost_geometry.radial_faces[radial + 1u], 3.0) -
         std::pow(ghost_geometry.radial_faces[radial], 3.0)) / 3.0;
    for (std::size_t theta = 0; theta + 1u < ghost_geometry.theta_faces.size(); ++theta) {
      const double theta_factor =
          std::cos(ghost_geometry.theta_faces[theta]) -
          std::cos(ghost_geometry.theta_faces[theta + 1u]);
      for (std::size_t phi = 0; phi + 1u < ghost_geometry.phi_faces.size(); ++phi) {
        const double phi_factor =
            ghost_geometry.phi_faces[phi + 1u] - ghost_geometry.phi_faces[phi];
        const double volume = radial_factor * theta_factor * phi_factor;
        ghost_geometry.cell_volumes.push_back(volume);
        ghost_geometry.global_volume += volume;
      }
    }
  }

  return ghost_geometry;
}

[[nodiscard]] dec3d::mesh::MacroZoneMap BuildBoundaryLockedGhostMacroZoneMap(
    const dec3d::mesh::SphericalGeometryMetadata& ghost_geometry,
    const dec3d::mesh::MacroZoneMap& interior_map,
    bool inner_side) noexcept {
  dec3d::mesh::MacroZoneMap map;
  if (!ghost_geometry.is_valid() ||
      !interior_map.is_complete() ||
      ghost_geometry.radial_faces.size() < 2u ||
      ghost_geometry.theta_faces.size() != interior_map.fine_theta_cells + 1u ||
      ghost_geometry.phi_faces.size() != interior_map.fine_phi_cells + 1u) {
    map.failure_reason =
        "boundary-locked macro ghost topology requires valid ghost geometry and interior map";
    return map;
  }

  const std::size_t ghost_radial_cells = ghost_geometry.radial_faces.size() - 1u;
  const std::size_t theta_cells = interior_map.fine_theta_cells;
  const std::size_t phi_cells = interior_map.fine_phi_cells;
  if (ghost_geometry.cell_volumes.size() != ghost_radial_cells * theta_cells * phi_cells) {
    map.failure_reason =
        "boundary-locked macro ghost topology requires ghost volumes matching ghost extents";
    return map;
  }

  const std::size_t reference_radial =
      inner_side ? 0u : interior_map.fine_radial_cells - 1u;
  map.fine_radial_cells = ghost_radial_cells;
  map.fine_theta_cells = theta_cells;
  map.fine_phi_cells = phi_cells;
  map.fine_to_coarse = dec3d::core::Array3D<std::size_t>(
      ghost_radial_cells,
      theta_cells,
      phi_cells,
      dec3d::mesh::InvalidMacroZoneIndex());

  bool found_reference_band = false;
  for (std::size_t ghost_radial = 0; ghost_radial < ghost_radial_cells; ++ghost_radial) {
    for (const auto& band : interior_map.theta_bands) {
      if (band.radial_index != reference_radial) {
        continue;
      }
      found_reference_band = true;
      auto ghost_band = band;
      ghost_band.radial_index = ghost_radial;
      map.theta_bands.push_back(ghost_band);
    }

    bool found_reference_cell = false;
    for (const auto& cell : interior_map.coarse_cells) {
      if (cell.radial_index != reference_radial) {
        continue;
      }
      found_reference_cell = true;
      auto ghost_cell = cell;
      ghost_cell.radial_index = ghost_radial;
      ghost_cell.total_volume = 0.0;
      const std::size_t coarse_index = map.coarse_cells.size();
      for (std::size_t theta = ghost_cell.theta_begin; theta < ghost_cell.theta_end; ++theta) {
        for (std::size_t phi = ghost_cell.phi_begin; phi < ghost_cell.phi_end; ++phi) {
          ghost_cell.total_volume += ghost_geometry.cell_volumes[LinearIndex(
              GridShape{ghost_radial_cells, theta_cells, phi_cells},
              ghost_radial,
              theta,
              phi)];
          map.fine_to_coarse(ghost_radial, theta, phi) = coarse_index;
        }
      }
      map.coarse_cells.push_back(ghost_cell);
    }

    if (!found_reference_cell) {
      map.failure_reason =
          "boundary-locked macro ghost topology could not find reference boundary cells";
      return map;
    }
  }

  if (!found_reference_band) {
    map.failure_reason =
        "boundary-locked macro ghost topology could not find reference boundary bands";
    return map;
  }

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    bool have_band = false;
    std::size_t last_phi_factor = 0u;
    for (const auto& band : map.theta_bands) {
      if (band.radial_index != radial) {
        continue;
      }
      if (have_band && band.phi_factor != last_phi_factor) {
        map.phi_factor_varies_by_theta_band = true;
      }
      have_band = true;
      last_phi_factor = band.phi_factor;
    }
  }

  std::ostringstream report;
  report << "macro_zoning_detected=true"
         << " radial_cells=" << map.fine_radial_cells
         << " theta_cells=" << map.fine_theta_cells
         << " phi_cells=" << map.fine_phi_cells
         << " coarse_cells=" << map.coarse_cells.size()
         << " theta_bands=" << map.theta_bands.size()
         << " phi_factor_varies_by_theta_band="
         << (map.phi_factor_varies_by_theta_band ? "true" : "false")
         << " topology_source=boundary_interior_layer";
  map.report_line = report.str();
  map.valid = true;
  if (!map.is_complete()) {
    map.valid = false;
    map.failure_reason =
        "boundary-locked macro ghost topology produced an incomplete coarse map";
  }
  return map;
}

[[nodiscard]] MacroCoarseGhostSide RestrictGhostStatesToCoarseSide(
    const dec3d::core::Array3D<HydroConservativeState>& ghost_states,
    const dec3d::mesh::SphericalGeometryMetadata& ghost_geometry,
    double coarse_factor,
    const dec3d::mesh::MacroZoneMap* locked_map = nullptr) noexcept {
  MacroCoarseGhostSide side;
  side.geometry = ghost_geometry;
  if (ghost_states.empty() ||
      !ghost_geometry.is_valid() ||
      ghost_states.extent_r() + 1u != ghost_geometry.radial_faces.size() ||
      ghost_states.extent_theta() + 1u != ghost_geometry.theta_faces.size() ||
      ghost_states.extent_phi() + 1u != ghost_geometry.phi_faces.size()) {
    side.failure_reason = "macro-zoning ghost restrict requires complete ghost states and geometry";
    return side;
  }

  side.map = locked_map == nullptr
                 ? dec3d::mesh::DetectMacroZones(
                       ghost_geometry,
                       dec3d::mesh::MacroZoningDetectionOptions{coarse_factor})
                 : *locked_map;
  if (!side.map.is_complete()) {
    side.failure_reason = side.map.failure_reason.empty()
                              ? "macro-zoning ghost detect failed"
                              : side.map.failure_reason;
    return side;
  }

  std::vector<double> rho(side.map.coarse_cells.size(), 0.0);
  std::vector<double> mom_r(side.map.coarse_cells.size(), 0.0);
  std::vector<double> mom_theta(side.map.coarse_cells.size(), 0.0);
  std::vector<double> mom_phi(side.map.coarse_cells.size(), 0.0);
  std::vector<double> e_fluid_total(side.map.coarse_cells.size(), 0.0);
  std::vector<double> chi_e(side.map.coarse_cells.size(), 0.0);
  std::vector<bool> has_reference(side.map.coarse_cells.size(), false);
  std::vector<HydroConservativeState> reference(side.map.coarse_cells.size());

  const std::size_t theta_cells = ghost_states.extent_theta();
  const std::size_t phi_cells = ghost_states.extent_phi();
  for (std::size_t radial = 0; radial < ghost_states.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const std::size_t coarse_index = side.map.fine_to_coarse(radial, theta, phi);
        const std::size_t linear_index = ((radial * theta_cells) + theta) * phi_cells + phi;
        const double cell_volume = ghost_geometry.cell_volumes[linear_index];
        const auto& state = ghost_states(radial, theta, phi);
        if (!has_reference[coarse_index]) {
          has_reference[coarse_index] = true;
          reference[coarse_index] = state;
        }
        rho[coarse_index] +=
            (state.rho - reference[coarse_index].rho) * cell_volume;
        mom_r[coarse_index] +=
            (state.mom_r - reference[coarse_index].mom_r) * cell_volume;
        mom_theta[coarse_index] +=
            (state.mom_theta - reference[coarse_index].mom_theta) * cell_volume;
        mom_phi[coarse_index] +=
            (state.mom_phi - reference[coarse_index].mom_phi) * cell_volume;
        e_fluid_total[coarse_index] +=
            (state.e_fluid_total - reference[coarse_index].e_fluid_total) *
            cell_volume;
        chi_e[coarse_index] +=
            (state.chi_e - reference[coarse_index].chi_e) * cell_volume;
      }
    }
  }

  side.states.assign(side.map.coarse_cells.size(), HydroConservativeState{});
  for (std::size_t coarse_index = 0; coarse_index < side.map.coarse_cells.size(); ++coarse_index) {
    const double coarse_volume = side.map.coarse_cells[coarse_index].total_volume;
    if (!(coarse_volume > 0.0) || !std::isfinite(coarse_volume)) {
      side.failure_reason = "macro-zoning ghost restrict produced a non-positive coarse volume";
      return side;
    }
    side.states[coarse_index] = {
        reference[coarse_index].rho + (rho[coarse_index] / coarse_volume),
        reference[coarse_index].mom_r + (mom_r[coarse_index] / coarse_volume),
        reference[coarse_index].mom_theta + (mom_theta[coarse_index] / coarse_volume),
        reference[coarse_index].mom_phi + (mom_phi[coarse_index] / coarse_volume),
        reference[coarse_index].e_fluid_total +
            (e_fluid_total[coarse_index] / coarse_volume),
        reference[coarse_index].chi_e + (chi_e[coarse_index] / coarse_volume)};
    if (!IsPhysicalCellState(side.states[coarse_index])) {
      side.failure_reason = "macro-zoning ghost restrict produced a non-physical coarse state";
      return side;
    }
  }

  side.consumed = true;
  return side;
}

[[nodiscard]] MacroGhostBootstrap BuildMacroGhostBootstrap(
    const RadialGhostOverride& radial_ghost_override,
    const dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    const dec3d::mesh::MacroZoneMap& interior_map,
    double coarse_factor) noexcept {
  MacroGhostBootstrap bootstrap;
  if (radial_ghost_override.ghost_layers == 0u) {
    return bootstrap;
  }

  bootstrap.summary.executed = true;
  bootstrap.summary.ghost_layers = radial_ghost_override.ghost_layers;

  if (radial_ghost_override.has_inner_neighbor) {
    const auto inner_geometry = BuildExtrapolatedGhostGeometry(
        local_geometry,
        radial_ghost_override.ghost_layers,
        true);
    const auto inner_locked_map =
        radial_ghost_override.inner_ghost_reuses_boundary_partition
            ? BuildBoundaryLockedGhostMacroZoneMap(inner_geometry, interior_map, true)
            : dec3d::mesh::MacroZoneMap{};
    bootstrap.inner = RestrictGhostStatesToCoarseSide(
        radial_ghost_override.inner_ghost_states,
        inner_geometry,
        coarse_factor,
        radial_ghost_override.inner_ghost_reuses_boundary_partition ? &inner_locked_map : nullptr);
    bootstrap.summary.inner_consumed = bootstrap.inner.is_complete();
    bootstrap.summary.inner_coarse_cells = bootstrap.inner.states.size();
    bootstrap.summary.inner_topology_source =
        radial_ghost_override.inner_ghost_reuses_boundary_partition
            ? "boundary_interior_layer"
            : "ghost_extrapolated_geometry";
  }

  if (radial_ghost_override.has_outer_neighbor) {
    const auto outer_geometry = BuildExtrapolatedGhostGeometry(
        local_geometry,
        radial_ghost_override.ghost_layers,
        false);
    const auto outer_locked_map =
        radial_ghost_override.outer_ghost_reuses_boundary_partition
            ? BuildBoundaryLockedGhostMacroZoneMap(outer_geometry, interior_map, false)
            : dec3d::mesh::MacroZoneMap{};
    bootstrap.outer = RestrictGhostStatesToCoarseSide(
        radial_ghost_override.outer_ghost_states,
        outer_geometry,
        coarse_factor,
        radial_ghost_override.outer_ghost_reuses_boundary_partition ? &outer_locked_map : nullptr);
    bootstrap.summary.outer_consumed = bootstrap.outer.is_complete();
    bootstrap.summary.outer_coarse_cells = bootstrap.outer.states.size();
    bootstrap.summary.outer_topology_source =
        radial_ghost_override.outer_ghost_reuses_boundary_partition
            ? "boundary_interior_layer"
            : "ghost_extrapolated_geometry";
  }

  std::ostringstream report;
  report << "ghost_layers=" << bootstrap.summary.ghost_layers
         << "; inner_neighbor=" << (radial_ghost_override.has_inner_neighbor ? "true" : "false")
         << "; outer_neighbor=" << (radial_ghost_override.has_outer_neighbor ? "true" : "false")
         << "; inner_consumed=" << (bootstrap.summary.inner_consumed ? "true" : "false")
         << "; outer_consumed=" << (bootstrap.summary.outer_consumed ? "true" : "false")
         << "; inner_coarse_cells=" << bootstrap.summary.inner_coarse_cells
         << "; outer_coarse_cells=" << bootstrap.summary.outer_coarse_cells
         << "; inner_topology_source=" << bootstrap.summary.inner_topology_source
         << "; outer_topology_source=" << bootstrap.summary.outer_topology_source
         << "; source=mpi_radial_halo_override"
         << "; target=macro_zoned_coarse_ghost_work_view";
  bootstrap.summary.report_line = report.str();
  return bootstrap;
}

[[nodiscard]] std::size_t FindBandCellCoveringFinePhi(
    const std::vector<std::size_t>& band_cells,
    const dec3d::mesh::MacroZoneMap& map,
    std::size_t fine_phi_index) noexcept {
  if (band_cells.empty() || fine_phi_index >= map.fine_phi_cells) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  const std::size_t reference_index = band_cells.front();
  if (reference_index >= map.coarse_cells.size()) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  const auto& reference_cell = map.coarse_cells[reference_index];
  const std::size_t coarse_index =
      map.fine_to_coarse(reference_cell.radial_index,
                         reference_cell.theta_begin,
                         fine_phi_index);
  if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex() ||
      coarse_index >= map.coarse_cells.size()) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  const auto& cell = map.coarse_cells[coarse_index];
  if (cell.radial_index != reference_cell.radial_index ||
      cell.theta_begin != reference_cell.theta_begin ||
      cell.theta_end != reference_cell.theta_end ||
      !(cell.phi_begin <= fine_phi_index && fine_phi_index < cell.phi_end)) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  return coarse_index;
}

[[nodiscard]] std::size_t FindCoarseCellCoveringFineThetaPhi(
    const dec3d::mesh::MacroZoneMap& map,
    std::size_t radial,
    std::size_t fine_theta_index,
    std::size_t fine_phi_index) noexcept {
  if (radial >= map.fine_radial_cells ||
      fine_theta_index >= map.fine_theta_cells ||
      fine_phi_index >= map.fine_phi_cells) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  const std::size_t coarse_index =
      map.fine_to_coarse(radial, fine_theta_index, fine_phi_index);
  if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex() ||
      coarse_index >= map.coarse_cells.size()) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  const auto& cell = map.coarse_cells[coarse_index];
  if (cell.radial_index != radial ||
      !(cell.theta_begin <= fine_theta_index && fine_theta_index < cell.theta_end) ||
      !(cell.phi_begin <= fine_phi_index && fine_phi_index < cell.phi_end)) {
    return dec3d::mesh::InvalidMacroZoneIndex();
  }
  return coarse_index;
}

[[nodiscard]] bool BuildMacroRadialPpmLine(
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const MacroGhostBootstrap& ghost_bootstrap,
    const RadialGhostOverride* radial_ghost_override,
    std::size_t fine_theta_index,
    std::size_t fine_phi_index,
    std::size_t ghost_layers,
    std::vector<std::size_t>& line_cell_indices,
    std::vector<HydroConservativeState>& ghosted_line,
    std::vector<double>& effective_widths,
    std::string& failure_reason) {
  const std::size_t interior_cells = map.fine_radial_cells;
  if (interior_cells == 0u) {
    failure_reason = "macro-zoning radial PPM line requires at least one coarse radial cell";
    return false;
  }

  line_cell_indices.assign(interior_cells, dec3d::mesh::InvalidMacroZoneIndex());
  ghosted_line.assign(interior_cells + (2u * ghost_layers), HydroConservativeState{});
  effective_widths.assign(interior_cells + (2u * ghost_layers), 0.0);

  for (std::size_t radial = 0; radial < interior_cells; ++radial) {
    const std::size_t coarse_index = FindCoarseCellCoveringFineThetaPhi(
        map,
        radial,
        fine_theta_index,
        fine_phi_index);
    if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
      failure_reason = "macro-zoning radial PPM line could not find an interior coarse cell";
      return false;
    }
    line_cell_indices[radial] = coarse_index;
    ghosted_line[ghost_layers + radial] = coarse_states[coarse_index];
    effective_widths[ghost_layers + radial] =
        geometry.radial_faces[radial + 1u] - geometry.radial_faces[radial];
  }

  const bool has_inner_neighbor =
      radial_ghost_override != nullptr && radial_ghost_override->has_inner_neighbor;
  const bool has_outer_neighbor =
      radial_ghost_override != nullptr && radial_ghost_override->has_outer_neighbor;

  if (has_inner_neighbor) {
    if (!ghost_bootstrap.inner.is_complete()) {
      failure_reason = ghost_bootstrap.inner.failure_reason.empty()
                           ? "macro-zoning radial PPM requires complete inner coarse ghost states"
                           : ghost_bootstrap.inner.failure_reason;
      return false;
    }
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      const std::size_t coarse_index = FindCoarseCellCoveringFineThetaPhi(
          ghost_bootstrap.inner.map,
          ghost,
          fine_theta_index,
          fine_phi_index);
      if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
        failure_reason = "macro-zoning radial PPM line could not find an inner coarse ghost cell";
        return false;
      }
      ghosted_line[ghost] = ghost_bootstrap.inner.states[coarse_index];
      effective_widths[ghost] =
          ghost_bootstrap.inner.geometry.radial_faces[ghost + 1u] -
          ghost_bootstrap.inner.geometry.radial_faces[ghost];
    }
  } else {
    if (interior_cells < ghost_layers) {
      failure_reason =
          "macro-zoning radial origin ghost remap requires at least ghost_layers interior radial cells";
      return false;
    }
    if ((map.fine_phi_cells % 2u) != 0u) {
      failure_reason =
          "macro-zoning radial origin ghost remap requires an even phi cell count";
      return false;
    }
    const std::size_t mapped_theta =
        MapThetaAcrossOrigin(fine_theta_index, map.fine_theta_cells);
    const std::size_t mapped_phi =
        MapPhiAcrossOrigin(fine_phi_index, map.fine_phi_cells);
    if (mapped_theta >= map.fine_theta_cells || mapped_phi >= map.fine_phi_cells) {
      failure_reason =
          "macro-zoning radial origin ghost remap could not map antipodal angular cell";
      return false;
    }
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      const std::size_t requested_radial = ghost_layers - 1u - ghost;
      const std::size_t coarse_index = FindCoarseCellCoveringFineThetaPhi(
          map,
          requested_radial,
          mapped_theta,
          mapped_phi);
      if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
        failure_reason =
            "macro-zoning radial origin ghost remap could not find an interior coarse cell";
        return false;
      }
      ghosted_line[ghost] = BuildOriginRadialGhostState(coarse_states[coarse_index]);
      effective_widths[ghost] = effective_widths[ghost_layers];
    }
  }

  if (has_outer_neighbor) {
    if (!ghost_bootstrap.outer.is_complete()) {
      failure_reason = ghost_bootstrap.outer.failure_reason.empty()
                           ? "macro-zoning radial PPM requires complete outer coarse ghost states"
                           : ghost_bootstrap.outer.failure_reason;
      return false;
    }
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      const std::size_t coarse_index = FindCoarseCellCoveringFineThetaPhi(
          ghost_bootstrap.outer.map,
          ghost,
          fine_theta_index,
          fine_phi_index);
      if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
        failure_reason = "macro-zoning radial PPM line could not find an outer coarse ghost cell";
        return false;
      }
      ghosted_line[ghost_layers + interior_cells + ghost] =
          ghost_bootstrap.outer.states[coarse_index];
      effective_widths[ghost_layers + interior_cells + ghost] =
          ghost_bootstrap.outer.geometry.radial_faces[ghost + 1u] -
          ghost_bootstrap.outer.geometry.radial_faces[ghost];
    }
  } else {
    const auto prepared = PrepareRadialBoundaryStates(
        BoundaryFace::upper,
        ghosted_line[ghost_layers + interior_cells - 1u]);
    if (!prepared.is_complete()) {
      failure_reason = prepared.failure_reason.empty()
                           ? "macro-zoning radial upper boundary ghost fill failed"
                           : prepared.failure_reason;
      return false;
    }
    for (std::size_t ghost = 0; ghost < ghost_layers; ++ghost) {
      ghosted_line[ghost_layers + interior_cells + ghost] = prepared.right_state;
      effective_widths[ghost_layers + interior_cells + ghost] =
          effective_widths[ghost_layers + interior_cells - 1u];
    }
  }

  return true;
}

[[nodiscard]] bool BuildMacroPhiPpmLine(
    const std::vector<std::size_t>& band_cells,
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t ghost_layers,
    std::vector<HydroConservativeState>& ghosted_line,
    std::vector<double>& effective_widths,
    std::string& failure_reason) {
  const std::size_t interior_cells = band_cells.size();
  if (interior_cells == 0u) {
    failure_reason = "macro-zoning phi PPM line requires at least one coarse cell";
    return false;
  }

  ghosted_line.assign(interior_cells + (2u * ghost_layers), HydroConservativeState{});
  effective_widths.assign(interior_cells + (2u * ghost_layers), 0.0);

  for (std::size_t local = 0; local < interior_cells; ++local) {
    const std::size_t coarse_index = band_cells[local];
    ghosted_line[ghost_layers + local] = coarse_states[coarse_index];
    effective_widths[ghost_layers + local] =
        CoarsePhiCellWidth(geometry, map.coarse_cells[coarse_index]);
  }

  for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
    const std::size_t lower_source =
        (interior_cells - (ghost % interior_cells)) % interior_cells;
    const std::size_t upper_source = (ghost - 1u) % interior_cells;
    ghosted_line[ghost_layers - ghost] = coarse_states[band_cells[lower_source]];
    ghosted_line[ghost_layers + interior_cells + ghost - 1u] =
        coarse_states[band_cells[upper_source]];
    effective_widths[ghost_layers - ghost] = effective_widths[ghost_layers + lower_source];
    effective_widths[ghost_layers + interior_cells + ghost - 1u] =
        effective_widths[ghost_layers + upper_source];
  }

  return true;
}

[[nodiscard]] bool BuildMacroThetaPpmLine(
    const std::vector<dec3d::mesh::MacroZoneThetaBand>& radial_bands,
    const std::vector<std::vector<std::size_t>>& radial_band_cells,
    const dec3d::mesh::MacroZoneMap& map,
    const std::vector<HydroConservativeState>& coarse_states,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t fine_phi_index,
    std::size_t ghost_layers,
    std::vector<std::size_t>& line_cell_indices,
    std::vector<HydroConservativeState>& ghosted_line,
    std::vector<double>& effective_widths,
    std::string& failure_reason) {
  const std::size_t interior_cells = radial_bands.size();
  if (interior_cells == 0u) {
    failure_reason = "macro-zoning theta PPM line requires at least one coarse theta band";
    return false;
  }

  line_cell_indices.assign(interior_cells, dec3d::mesh::InvalidMacroZoneIndex());
  ghosted_line.assign(interior_cells + (2u * ghost_layers), HydroConservativeState{});
  effective_widths.assign(interior_cells + (2u * ghost_layers), 0.0);

  for (std::size_t band_index = 0; band_index < interior_cells; ++band_index) {
    const std::size_t coarse_index = FindBandCellCoveringFinePhi(
        radial_band_cells[band_index],
        map,
        fine_phi_index);
    if (coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
      failure_reason = "macro-zoning theta PPM line could not find a coarse cell for the fine-phi track";
      return false;
    }

    line_cell_indices[band_index] = coarse_index;
    ghosted_line[ghost_layers + band_index] = coarse_states[coarse_index];
    effective_widths[ghost_layers + band_index] =
        CoarseThetaBandWidth(geometry, radial_bands[band_index]);
  }

  const std::size_t mapped_phi = MapPhiAcrossPole(fine_phi_index, map.fine_phi_cells);
  if (mapped_phi >= map.fine_phi_cells) {
    failure_reason = "macro-zoning theta PPM line produced an invalid mapped phi index";
    return false;
  }

  for (std::size_t ghost = 1u; ghost <= ghost_layers; ++ghost) {
    const std::size_t lower_band = ghost - 1u;
    const std::size_t upper_band = interior_cells - ghost;
    const std::size_t lower_coarse_index = FindBandCellCoveringFinePhi(
        radial_band_cells[lower_band],
        map,
        mapped_phi);
    const std::size_t upper_coarse_index = FindBandCellCoveringFinePhi(
        radial_band_cells[upper_band],
        map,
        mapped_phi);
    if (lower_coarse_index == dec3d::mesh::InvalidMacroZoneIndex() ||
        upper_coarse_index == dec3d::mesh::InvalidMacroZoneIndex()) {
      failure_reason = "macro-zoning theta PPM line could not build pole-remapped ghost states";
      return false;
    }

    ghosted_line[ghost_layers - ghost] =
        BuildThetaPoleGhostState(coarse_states[lower_coarse_index]);
    ghosted_line[ghost_layers + interior_cells + ghost - 1u] =
        BuildThetaPoleGhostState(coarse_states[upper_coarse_index]);
    effective_widths[ghost_layers - ghost] = effective_widths[ghost_layers + lower_band];
    effective_widths[ghost_layers + interior_cells + ghost - 1u] =
        effective_widths[ghost_layers + upper_band];
  }

  return true;
}

[[nodiscard]] bool AccumulateMacroZonedHydroDelta(
    const HydroStateSnapshot& snapshot,
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const StaticGridHydroOptions& options,
    const RadialGhostOverride* radial_ghost_override,
    std::vector<HydroConservativeState>& accumulated_delta,
    std::string& failure_reason,
    bool allow_post_macro_ale_remap,
    MacroZonedHydroSummary* summary,
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal = nullptr,
    const dec3d::mesh::SphericalGeometryMetadata* ale_updated_geometry = nullptr,
    MacroAleHllcDiagnostics* macro_ale_hllc_diagnostics = nullptr,
    MacroAleHllcLocalFaceWindow* macro_ale_hllc_local_face_window = nullptr) {
  if (!snapshot.is_complete() ||
      !hydro_view.is_complete() ||
      accumulated_delta.size() != snapshot.cells.size()) {
    failure_reason = "macro-zoning angular update requires complete fine hydro state";
    return false;
  }

  const bool macro_ale_direct_hllc_mode =
      options.apply_radial_ale_flux_correction &&
      options.use_macro_ale_direct_moving_face_hllc;
  if (options.apply_radial_ale_flux_correction &&
      !allow_post_macro_ale_remap &&
      !macro_ale_direct_hllc_mode) {
    failure_reason =
        "macro-zoning coupling currently does not support ALE-aware angular updates";
    return false;
  }
  if (macro_ale_direct_hllc_mode) {
    if (radial_ale_proposal == nullptr ||
        ale_updated_geometry == nullptr ||
        macro_ale_hllc_diagnostics == nullptr ||
        macro_ale_hllc_local_face_window == nullptr) {
      failure_reason =
          "direct moving-face HLLC macro ALE requires proposal, preview geometry, diagnostics, and face window";
      if (macro_ale_hllc_diagnostics != nullptr) {
        macro_ale_hllc_diagnostics->failure_class =
            MacroAleHllcFailureClass::proposal_window_missing;
      }
      return false;
    }
    if (!radial_ale_proposal->radial_face_indexing_is_global ||
        !radial_ale_proposal->has_radial_face_window(
            options.radial_ale_global_face_begin_index,
            geometry.radial_faces.size())) {
      failure_reason =
          "direct moving-face HLLC macro ALE requires a complete local global-face proposal window";
      macro_ale_hllc_diagnostics->failure_class =
          MacroAleHllcFailureClass::proposal_window_missing;
      return false;
    }
  }

  auto timing_start = std::chrono::steady_clock::now();
  const auto map = dec3d::mesh::DetectMacroZones(
      geometry,
      dec3d::mesh::MacroZoningDetectionOptions{options.macro_zoning_coarse_factor});
  if (summary != nullptr) {
    summary->detect_wall_s += ElapsedHydroSecondsSince(timing_start);
  }
  if (!map.is_complete()) {
    failure_reason = map.failure_reason.empty()
                         ? "macro-zoning detection failed"
                         : map.failure_reason;
    return false;
  }

  const auto fine_view = BuildFineHydroPackageView(hydro_view);
  timing_start = std::chrono::steady_clock::now();
  const auto coarse_package = dec3d::mesh::RestrictFineHydroPackage(
      fine_view,
      geometry,
      map);
  if (summary != nullptr) {
    summary->restrict_wall_s += ElapsedHydroSecondsSince(timing_start);
  }
  if (!coarse_package.is_complete()) {
    failure_reason = coarse_package.failure_reason.empty()
                         ? "macro-zoning restrict failed"
                         : coarse_package.failure_reason;
    return false;
  }

  const auto macro_update_start = std::chrono::steady_clock::now();
  MacroGhostBootstrap ghost_bootstrap;
  if (radial_ghost_override != nullptr &&
      radial_ghost_override->ghost_layers != 0u) {
    ghost_bootstrap = BuildMacroGhostBootstrap(
        *radial_ghost_override,
        geometry,
        map,
        options.macro_zoning_coarse_factor);
    if (summary != nullptr) {
      summary->coarse_ghost_bootstrap_executed = ghost_bootstrap.summary.executed;
      summary->coarse_ghost_bootstrap_report_line = ghost_bootstrap.summary.report_line;
    }
  }

  auto updated_coarse = coarse_package;
  auto coarse_states = BuildMacroZoneStates(coarse_package);
  std::vector<HydroConservativeState> coarse_extensive_delta(
      coarse_states.size(),
      HydroConservativeState{});
  std::vector<HydroConservativeState> coarse_radial_delta_per_solid_angle_reference(
      coarse_states.size(),
      HydroConservativeState{});
  std::vector<HydroConservativeState> coarse_radial_delta_per_solid_angle_weighted_delta_sum(
      coarse_states.size(),
      HydroConservativeState{});
  std::vector<HydroConservativeState> coarse_radial_delta_extensive(
      coarse_states.size(),
      HydroConservativeState{});
  std::vector<double> coarse_radial_delta_solid_angle_weight(
      coarse_states.size(),
      0.0);
  std::vector<bool> coarse_radial_delta_has_reference(coarse_states.size(), false);
  std::vector<HydroConservativeState> coarse_source_extensive_delta(
      coarse_states.size(),
      HydroConservativeState{});
  bool coarse_source_budget_integrated = false;
  const std::size_t ppm_ghost_layers =
      options.use_ppm_reconstruction ? options.reconstruction_ghost_layers : 1u;
  dec3d::mesh::RadialAleLocalProposalWindow macro_radial_ale_window;
  if (macro_ale_direct_hllc_mode) {
    macro_radial_ale_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
        *radial_ale_proposal,
        options.radial_ale_global_face_begin_index,
        geometry.radial_faces.size());
    if (!macro_radial_ale_window.is_complete()) {
      failure_reason =
          macro_radial_ale_window.failure_reason.empty()
              ? "direct moving-face HLLC macro ALE could not build local proposal window"
              : macro_radial_ale_window.failure_reason;
      macro_ale_hllc_diagnostics->failure_class =
          MacroAleHllcFailureClass::proposal_window_missing;
      return false;
    }
    macro_ale_hllc_local_face_window->available = true;
    macro_ale_hllc_local_face_window->global_face_begin =
        options.radial_ale_global_face_begin_index;
    macro_ale_hllc_local_face_window->local_face_count =
        geometry.radial_faces.size();
    macro_ale_hllc_local_face_window->face_fluxes.assign(
        geometry.radial_faces.size(),
        HydroConservativeState{});
    macro_ale_hllc_local_face_window->face_speeds =
        macro_radial_ale_window.radial_face_velocities;
    macro_ale_hllc_local_face_window->proposed_face_radii =
        macro_radial_ale_window.proposed_radial_faces;
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "a4_after_restrict_raw_coarse_state",
        "macro_coarse_state_raw_average",
        BuildCoarseStateAngularStageMetrics(map, coarse_states, geometry)));
    if (macro_ale_direct_hllc_mode) {
      summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
          "a4_after_restrict_zero_delta_vq_prediction",
          "macro_coarse_vq_prediction_new_geometry",
          BuildCoarseA4ExtensiveStageMetrics(
              map,
              coarse_states,
              coarse_extensive_delta,
              coarse_radial_delta_per_solid_angle_reference,
              coarse_radial_delta_per_solid_angle_weighted_delta_sum,
              coarse_radial_delta_extensive,
              coarse_radial_delta_solid_angle_weight,
              geometry,
              *ale_updated_geometry)));
    }
  }

  if (macro_ale_direct_hllc_mode && options.apply_geometric_source) {
    coarse_source_extensive_delta =
        BuildCoarseGeometricSourceExtensiveDelta(
            map,
            coarse_states,
            geometry,
            dt_s,
            true);
    if (coarse_source_extensive_delta.size() != coarse_extensive_delta.size()) {
      failure_reason =
          "direct moving-face HLLC macro ALE could not build coarse source extensive delta";
      macro_ale_hllc_diagnostics->failure_class =
          MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
      return false;
    }
    double max_abs_source_delta = 0.0;
    for (std::size_t coarse_index = 0; coarse_index < coarse_extensive_delta.size(); ++coarse_index) {
      if (!coarse_source_extensive_delta[coarse_index].is_finite()) {
        failure_reason =
            "direct moving-face HLLC macro ALE coarse source extensive delta is non-finite";
        macro_ale_hllc_diagnostics->failure_class =
            MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
        return false;
      }
      max_abs_source_delta = std::max(
          max_abs_source_delta,
          std::max({
              std::abs(coarse_source_extensive_delta[coarse_index].mom_r),
              std::abs(coarse_source_extensive_delta[coarse_index].mom_theta),
              std::abs(coarse_source_extensive_delta[coarse_index].mom_phi)}));
      coarse_extensive_delta[coarse_index] = Add(
          coarse_extensive_delta[coarse_index],
          coarse_source_extensive_delta[coarse_index]);
    }
    coarse_source_budget_integrated = true;
    if (summary != nullptr) {
      summary->coarse_source_budget_executed = true;
      std::ostringstream source_report;
      source_report << std::setprecision(17)
                    << "mode=direct_moving_face_hllc"
                    << "; source_budget=coarse_extensive_accumulator"
                    << "; source_order=old_time_coarse_state"
                    << "; source_geometry=geometry_n"
                    << "; pressure_source_balance=moving_interface_face_pressure_quadrature"
                    << "; source_components=mom_r,mom_theta,mom_phi"
                    << "; max_abs_source_extensive_delta="
                    << max_abs_source_delta
                    << "; coarse_cells=" << coarse_extensive_delta.size();
      summary->coarse_source_budget_report_line = source_report.str();
    }
  }

  if (summary != nullptr &&
      options.debug_angular_stage_diagnostics &&
      macro_ale_direct_hllc_mode) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "a4_after_source_budget_vq_prediction",
        "macro_coarse_vq_prediction_new_geometry",
        BuildCoarseA4ExtensiveStageMetrics(
            map,
            coarse_states,
            coarse_extensive_delta,
            coarse_radial_delta_per_solid_angle_reference,
            coarse_radial_delta_per_solid_angle_weighted_delta_sum,
            coarse_radial_delta_extensive,
            coarse_radial_delta_solid_angle_weight,
            geometry,
            *ale_updated_geometry)));
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "post_macro_restrict_coarse_state",
        "macro_coarse_state",
        BuildCoarseAngularStageMetrics(map, coarse_states, coarse_extensive_delta)));
  }

  if (options.apply_radial_sweep &&
      options.use_ppm_reconstruction &&
      options.macro_zoning_use_radial_ppm) {
    const auto macro_radial_start = std::chrono::steady_clock::now();
    struct MacroRadialFluxLineDetail {
      bool valid{false};
      std::size_t theta{0};
      std::size_t phi{0};
      std::size_t left_coarse_index{dec3d::mesh::InvalidMacroZoneIndex()};
      std::size_t right_coarse_index{dec3d::mesh::InvalidMacroZoneIndex()};
      double left_width{0.0};
      double right_width{0.0};
      bool downgraded_to_first_order{false};
      bool shared_shell_fallback{false};
      bool left_profile_troubled{false};
      bool right_profile_troubled{false};
      bool left_trace_failed{false};
      bool right_trace_failed{false};
      int left_troubled_component{-1};
      int right_troubled_component{-1};
      double left_trouble_q_im1{0.0};
      double left_trouble_q_i{0.0};
      double left_trouble_q_ip1{0.0};
      double left_trouble_total_variation{0.0};
      double left_trouble_end_to_end{0.0};
      double right_trouble_q_im1{0.0};
      double right_trouble_q_i{0.0};
      double right_trouble_q_ip1{0.0};
      double right_trouble_total_variation{0.0};
      double right_trouble_end_to_end{0.0};
      double mass_flux{0.0};
      double mom_r_flux{0.0};
      double e_total_flux{0.0};
      double face_speed{0.0};
      double s_left{0.0};
      double s_star{0.0};
      double s_right{0.0};
      HllcActiveRegion moving_branch{HllcActiveRegion::invalid};
      DirectionalPrimitiveState left_traced;
      DirectionalPrimitiveState right_traced;
      DirectionalPrimitiveState left_cell_center;
      DirectionalPrimitiveState right_cell_center;
      PpmSmoothnessProbe left_component1_smoothness;
      PpmSmoothnessProbe right_component1_smoothness;
      std::vector<int> stencil_logical_radial;
      std::vector<double> stencil_rho;
      std::vector<double> stencil_v_n;
      std::vector<double> stencil_pressure;
    };

    struct MacroRadialPpmLineWork {
      std::size_t theta{0};
      std::size_t phi{0};
      std::vector<std::size_t> line_cell_indices;
      std::vector<HydroConservativeState> ghosted_line;
      std::vector<double> effective_widths;
      PpmLineReconstructionResult line_reconstruction;
      std::size_t downgraded_interface_count{0u};
    };

    const bool collect_radial_flux_diagnostics =
        summary != nullptr && options.debug_angular_stage_diagnostics;
    const std::size_t radial_face_diagnostic_count =
        collect_radial_flux_diagnostics ? map.fine_radial_cells + 1u : 0u;
    const std::size_t radial_cell_diagnostic_count =
        collect_radial_flux_diagnostics ? map.fine_radial_cells : 0u;
    std::vector<double> radial_mass_flux_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_mass_flux_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> radial_mass_flux_sum(radial_face_diagnostic_count, 0.0);
    std::vector<std::size_t> radial_mass_flux_count(radial_face_diagnostic_count, 0u);
    std::vector<MacroRadialFluxLineDetail> radial_mass_flux_min_detail(
        radial_face_diagnostic_count);
    std::vector<MacroRadialFluxLineDetail> radial_mass_flux_max_detail(
        radial_face_diagnostic_count);
    std::vector<double> radial_mom_r_flux_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_mom_r_flux_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> radial_mom_r_flux_sum(radial_face_diagnostic_count, 0.0);
    std::vector<std::size_t> radial_mom_r_flux_count(radial_face_diagnostic_count, 0u);
    std::vector<MacroRadialFluxLineDetail> radial_mom_r_flux_min_detail(
        radial_face_diagnostic_count);
    std::vector<MacroRadialFluxLineDetail> radial_mom_r_flux_max_detail(
        radial_face_diagnostic_count);
    std::vector<double> radial_e_total_flux_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_e_total_flux_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> radial_e_total_flux_sum(radial_face_diagnostic_count, 0.0);
    std::vector<std::size_t> radial_e_total_flux_count(radial_face_diagnostic_count, 0u);
    std::vector<MacroRadialFluxLineDetail> radial_e_total_flux_min_detail(
        radial_face_diagnostic_count);
    std::vector<MacroRadialFluxLineDetail> radial_e_total_flux_max_detail(
        radial_face_diagnostic_count);
    const auto record_radial_flux_spread =
        [](std::size_t radial_face,
           double flux,
           const MacroRadialFluxLineDetail& detail,
           std::vector<double>& min_flux,
           std::vector<double>& max_flux,
           std::vector<double>& sum_flux,
           std::vector<std::size_t>& count_flux,
           std::vector<MacroRadialFluxLineDetail>& min_detail,
           std::vector<MacroRadialFluxLineDetail>& max_detail) {
          if (flux < min_flux[radial_face]) {
            min_flux[radial_face] = flux;
            min_detail[radial_face] = detail;
          }
          if (flux > max_flux[radial_face]) {
            max_flux[radial_face] = flux;
            max_detail[radial_face] = detail;
          }
          sum_flux[radial_face] += flux;
          count_flux[radial_face] += 1u;
    };
    std::vector<MacroRadialFluxLineDetail> debug_face_smoothness_details;
    std::vector<MacroRadialPpmLineWork> radial_ppm_lines;
    radial_ppm_lines.reserve(map.fine_theta_cells * map.fine_phi_cells);
    std::vector<bool> shell_fallback_mask(map.fine_radial_cells + 1u, false);
    std::vector<double> radial_input_rho_min(
        radial_cell_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_input_rho_max(
        radial_cell_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> radial_input_vn_min(
        radial_cell_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_input_vn_max(
        radial_cell_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> radial_input_pressure_min(
        radial_cell_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> radial_input_pressure_max(
        radial_cell_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<std::size_t> radial_input_count(radial_cell_diagnostic_count, 0u);
    std::vector<double> traced_left_rho_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_left_rho_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_rho_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_rho_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> traced_left_vn_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_left_vn_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_vn_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_vn_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_rho_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_rho_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_rho_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_rho_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_vn_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_vn_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_vn_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_vn_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_pressure_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_left_pressure_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_pressure_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> edge_right_pressure_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<std::size_t> edge_interface_count(radial_face_diagnostic_count, 0u);
    std::vector<double> traced_left_pressure_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_left_pressure_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_pressure_min(
        radial_face_diagnostic_count,
        std::numeric_limits<double>::infinity());
    std::vector<double> traced_right_pressure_max(
        radial_face_diagnostic_count,
        -std::numeric_limits<double>::infinity());
    std::vector<std::size_t> traced_interface_count(radial_face_diagnostic_count, 0u);
    const auto append_primitive_spread_report =
        [&](std::string_view stage,
            std::string_view representation,
            const std::vector<double>& rho_min,
            const std::vector<double>& rho_max,
            const std::vector<double>& vn_min,
            const std::vector<double>& vn_max,
            const std::vector<double>& pressure_min,
            const std::vector<double>& pressure_max,
            const std::vector<std::size_t>& counts) {
          std::size_t max_rho_index = 0u;
          double max_rho_abs = 0.0;
          double max_rho_relative = 0.0;
          std::size_t max_vn_index = 0u;
          double max_vn_abs = 0.0;
          double max_vn_relative = 0.0;
          std::size_t max_pressure_index = 0u;
          double max_pressure_abs = 0.0;
          double max_pressure_relative = 0.0;
          std::size_t total_count = 0u;
          for (std::size_t index = 0; index < counts.size(); ++index) {
            if (counts[index] == 0u) {
              continue;
            }
            total_count += counts[index];
            const double rho_abs = rho_max[index] - rho_min[index];
            const double rho_mean = 0.5 * (rho_max[index] + rho_min[index]);
            const double rho_relative =
                std::abs(rho_mean) > 0.0 ? rho_abs / std::abs(rho_mean) : 0.0;
            if (rho_relative > max_rho_relative) {
              max_rho_relative = rho_relative;
              max_rho_abs = rho_abs;
              max_rho_index = index;
            }
            const double vn_abs = vn_max[index] - vn_min[index];
            const double vn_mean = 0.5 * (vn_max[index] + vn_min[index]);
            const double vn_relative =
                std::abs(vn_mean) > 0.0 ? vn_abs / std::abs(vn_mean) : 0.0;
            if (vn_relative > max_vn_relative) {
              max_vn_relative = vn_relative;
              max_vn_abs = vn_abs;
              max_vn_index = index;
            }
            const double pressure_abs = pressure_max[index] - pressure_min[index];
            const double pressure_mean =
                0.5 * (pressure_max[index] + pressure_min[index]);
            const double pressure_relative =
                std::abs(pressure_mean) > 0.0
                    ? pressure_abs / std::abs(pressure_mean)
                    : 0.0;
            if (pressure_relative > max_pressure_relative) {
              max_pressure_relative = pressure_relative;
              max_pressure_abs = pressure_abs;
              max_pressure_index = index;
            }
          }
          std::ostringstream report;
          report << std::setprecision(17)
                 << "stage=" << stage
                 << "; representation=" << representation
                 << "; valid=" << (total_count > 0u ? "true" : "false")
                 << "; sample_count=" << total_count
                 << "; max_rho_relative_spread=" << max_rho_relative
                 << "; max_rho_spread_index=" << max_rho_index
                 << "; max_rho_absolute_spread=" << max_rho_abs
                 << "; max_vn_relative_spread=" << max_vn_relative
                 << "; max_vn_spread_index=" << max_vn_index
                 << "; max_vn_absolute_spread=" << max_vn_abs
                 << "; max_pressure_relative_spread=" << max_pressure_relative
                 << "; max_pressure_spread_index=" << max_pressure_index
                 << "; max_pressure_absolute_spread=" << max_pressure_abs;
          summary->angular_stage_report_lines.push_back(report.str());
        };

    for (std::size_t theta = 0; theta < map.fine_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < map.fine_phi_cells; ++phi) {
        MacroRadialPpmLineWork line_work;
        line_work.theta = theta;
        line_work.phi = phi;
        if (!BuildMacroRadialPpmLine(
                map,
                coarse_states,
                geometry,
                ghost_bootstrap,
                radial_ghost_override,
                theta,
                phi,
                ppm_ghost_layers,
                line_work.line_cell_indices,
                line_work.ghosted_line,
                line_work.effective_widths,
                failure_reason)) {
          return false;
        }

        if (collect_radial_flux_diagnostics) {
          for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
            const auto primitive = ToDirectionalPrimitive(
                line_work.ghosted_line[ppm_ghost_layers + radial],
                HydroDirection::radial);
            radial_input_rho_min[radial] =
                std::min(radial_input_rho_min[radial], primitive.rho);
            radial_input_rho_max[radial] =
                std::max(radial_input_rho_max[radial], primitive.rho);
            radial_input_vn_min[radial] =
                std::min(radial_input_vn_min[radial], primitive.v_n);
            radial_input_vn_max[radial] =
                std::max(radial_input_vn_max[radial], primitive.v_n);
            radial_input_pressure_min[radial] =
                std::min(radial_input_pressure_min[radial], primitive.pressure);
            radial_input_pressure_max[radial] =
                std::max(radial_input_pressure_max[radial], primitive.pressure);
            ++radial_input_count[radial];
          }
        }

        line_work.line_reconstruction = ReconstructPpmLine(
            line_work.ghosted_line,
            HydroDirection::radial,
            map.fine_radial_cells,
            ppm_ghost_layers,
            dt_s,
            line_work.effective_widths,
            macro_ale_direct_hllc_mode
                ? &macro_radial_ale_window.radial_face_velocities
                : nullptr);
        if (!line_work.line_reconstruction.success) {
          failure_reason = line_work.line_reconstruction.failure_reason.empty()
                               ? "macro radial PPM line reconstruction failed"
                               : line_work.line_reconstruction.failure_reason;
          return false;
        }
        AccumulateShellWidePpmFallbackMask(
            line_work.line_reconstruction,
            shell_fallback_mask);
        radial_ppm_lines.push_back(std::move(line_work));
      }
    }

    std::vector<HydroConservativeState> radial_interface_fluxes(
        map.fine_radial_cells + 1u,
        HydroConservativeState{});
    std::vector<double> radial_moving_interface_pressure_terms(
        macro_ale_direct_hllc_mode ? map.fine_radial_cells + 1u : 0u,
        std::numeric_limits<double>::quiet_NaN());
    for (auto& line_work : radial_ppm_lines) {
      if (macro_ale_direct_hllc_mode) {
        std::fill(
            radial_moving_interface_pressure_terms.begin(),
            radial_moving_interface_pressure_terms.end(),
            std::numeric_limits<double>::quiet_NaN());
      }
      if (!BuildPpmInterfaceFluxesFromReconstruction(
              line_work.line_reconstruction,
              line_work.ghosted_line,
              SweepDirection::radial,
              ppm_ghost_layers,
              shell_fallback_mask,
              radial_interface_fluxes,
              line_work.downgraded_interface_count,
              failure_reason,
              macro_ale_direct_hllc_mode
                  ? &macro_radial_ale_window.radial_face_velocities
                  : nullptr,
              nullptr,
              macro_ale_direct_hllc_mode
                  ? macro_ale_hllc_diagnostics
                  : nullptr,
              macro_ale_direct_hllc_mode
                  ? &radial_moving_interface_pressure_terms
                  : nullptr)) {
        return false;
      }

      for (std::size_t radial_face = 0; radial_face <= map.fine_radial_cells; ++radial_face) {
        const double mass_flux = radial_interface_fluxes[radial_face].rho;
        const double mom_r_flux = radial_interface_fluxes[radial_face].mom_r;
        const double e_total_flux =
            radial_interface_fluxes[radial_face].e_fluid_total;
        if (macro_ale_direct_hllc_mode) {
          const double angular_face_area =
              RadialFaceAreaPerSolidAngle(geometry, radial_face) *
              CellSolidAngleFactor(geometry, line_work.theta, line_work.phi);
          if (!(angular_face_area >= 0.0) || !std::isfinite(angular_face_area)) {
            failure_reason =
                "direct moving-face HLLC macro ALE encountered invalid radial face area";
            macro_ale_hllc_diagnostics->failure_class =
                MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
            return false;
          }
          macro_ale_hllc_local_face_window->face_fluxes[radial_face] =
              Add(
                  macro_ale_hllc_local_face_window->face_fluxes[radial_face],
                  Scale(radial_interface_fluxes[radial_face], angular_face_area));
        }
        if (collect_radial_flux_diagnostics) {
          MacroRadialFluxLineDetail detail;
          if (line_work.line_reconstruction.success &&
              radial_face < line_work.line_reconstruction.interfaces.size()) {
          const auto& interface_state =
              line_work.line_reconstruction.interfaces[radial_face];
          const std::size_t left_line_cell = ppm_ghost_layers + radial_face - 1u;
          const std::size_t right_line_cell = ppm_ghost_layers + radial_face;
          detail.valid = true;
          detail.theta = line_work.theta;
          detail.phi = line_work.phi;
          detail.mass_flux = mass_flux;
          detail.mom_r_flux = mom_r_flux;
          detail.e_total_flux = e_total_flux;
          if (macro_ale_direct_hllc_mode &&
              radial_face < macro_radial_ale_window.radial_face_velocities.size()) {
            detail.face_speed =
                macro_radial_ale_window.radial_face_velocities[radial_face];
          }
          detail.shared_shell_fallback = shell_fallback_mask[radial_face];
          detail.downgraded_to_first_order =
              interface_state.downgraded_to_first_order ||
              detail.shared_shell_fallback;
          detail.left_profile_troubled =
              interface_state.left_profile_troubled;
          detail.right_profile_troubled =
              interface_state.right_profile_troubled;
          detail.left_trace_failed =
              interface_state.left_trace_failed;
          detail.right_trace_failed =
              interface_state.right_trace_failed;
          detail.left_troubled_component =
              interface_state.left_troubled_component;
          detail.right_troubled_component =
              interface_state.right_troubled_component;
          detail.left_trouble_q_im1 =
              interface_state.left_trouble_q_im1;
          detail.left_trouble_q_i =
              interface_state.left_trouble_q_i;
          detail.left_trouble_q_ip1 =
              interface_state.left_trouble_q_ip1;
          detail.left_trouble_total_variation =
              interface_state.left_trouble_total_variation;
          detail.left_trouble_end_to_end =
              interface_state.left_trouble_end_to_end;
          detail.right_trouble_q_im1 =
              interface_state.right_trouble_q_im1;
          detail.right_trouble_q_i =
              interface_state.right_trouble_q_i;
          detail.right_trouble_q_ip1 =
              interface_state.right_trouble_q_ip1;
          detail.right_trouble_total_variation =
              interface_state.right_trouble_total_variation;
          detail.right_trouble_end_to_end =
              interface_state.right_trouble_end_to_end;
          detail.left_cell_center = ToDirectionalPrimitive(
              line_work.ghosted_line[left_line_cell],
              HydroDirection::radial);
          detail.right_cell_center = ToDirectionalPrimitive(
              line_work.ghosted_line[right_line_cell],
              HydroDirection::radial);
          detail.left_traced = detail.downgraded_to_first_order
                                   ? detail.left_cell_center
                                   : interface_state.left;
          detail.right_traced = detail.downgraded_to_first_order
                                     ? detail.right_cell_center
                                     : interface_state.right;
          if (interface_state.edge_state_diagnostics_available) {
            edge_left_rho_min[radial_face] =
                std::min(edge_left_rho_min[radial_face],
                         interface_state.left_edge_before_trace.rho);
            edge_left_rho_max[radial_face] =
                std::max(edge_left_rho_max[radial_face],
                         interface_state.left_edge_before_trace.rho);
            edge_right_rho_min[radial_face] =
                std::min(edge_right_rho_min[radial_face],
                         interface_state.right_edge_before_trace.rho);
            edge_right_rho_max[radial_face] =
                std::max(edge_right_rho_max[radial_face],
                         interface_state.right_edge_before_trace.rho);
            edge_left_vn_min[radial_face] =
                std::min(edge_left_vn_min[radial_face],
                         interface_state.left_edge_before_trace.v_n);
            edge_left_vn_max[radial_face] =
                std::max(edge_left_vn_max[radial_face],
                         interface_state.left_edge_before_trace.v_n);
            edge_right_vn_min[radial_face] =
                std::min(edge_right_vn_min[radial_face],
                         interface_state.right_edge_before_trace.v_n);
            edge_right_vn_max[radial_face] =
                std::max(edge_right_vn_max[radial_face],
                         interface_state.right_edge_before_trace.v_n);
            edge_left_pressure_min[radial_face] =
                std::min(edge_left_pressure_min[radial_face],
                         interface_state.left_edge_before_trace.pressure);
            edge_left_pressure_max[radial_face] =
                std::max(edge_left_pressure_max[radial_face],
                         interface_state.left_edge_before_trace.pressure);
            edge_right_pressure_min[radial_face] =
                std::min(edge_right_pressure_min[radial_face],
                         interface_state.right_edge_before_trace.pressure);
            edge_right_pressure_max[radial_face] =
                std::max(edge_right_pressure_max[radial_face],
                         interface_state.right_edge_before_trace.pressure);
            ++edge_interface_count[radial_face];
          }
          traced_left_rho_min[radial_face] =
              std::min(traced_left_rho_min[radial_face], detail.left_traced.rho);
          traced_left_rho_max[radial_face] =
              std::max(traced_left_rho_max[radial_face], detail.left_traced.rho);
          traced_right_rho_min[radial_face] =
              std::min(traced_right_rho_min[radial_face], detail.right_traced.rho);
          traced_right_rho_max[radial_face] =
              std::max(traced_right_rho_max[radial_face], detail.right_traced.rho);
          traced_left_vn_min[radial_face] =
              std::min(traced_left_vn_min[radial_face], detail.left_traced.v_n);
          traced_left_vn_max[radial_face] =
              std::max(traced_left_vn_max[radial_face], detail.left_traced.v_n);
          traced_right_vn_min[radial_face] =
              std::min(traced_right_vn_min[radial_face], detail.right_traced.v_n);
          traced_right_vn_max[radial_face] =
              std::max(traced_right_vn_max[radial_face], detail.right_traced.v_n);
          traced_left_pressure_min[radial_face] =
              std::min(traced_left_pressure_min[radial_face], detail.left_traced.pressure);
          traced_left_pressure_max[radial_face] =
              std::max(traced_left_pressure_max[radial_face], detail.left_traced.pressure);
          traced_right_pressure_min[radial_face] =
              std::min(traced_right_pressure_min[radial_face], detail.right_traced.pressure);
          traced_right_pressure_max[radial_face] =
              std::max(traced_right_pressure_max[radial_face], detail.right_traced.pressure);
          ++traced_interface_count[radial_face];
          if (macro_ale_direct_hllc_mode) {
            const auto left_directional_state =
                MakeDirectionalConservativeState(detail.left_traced);
            const auto right_directional_state =
                MakeDirectionalConservativeState(detail.right_traced);
            const auto hllc =
                SolveHllcRiemann(left_directional_state, right_directional_state);
            if (hllc.is_complete()) {
              detail.s_left = hllc.waves.s_left;
              detail.s_star = hllc.waves.s_star;
              detail.s_right = hllc.waves.s_right;
              detail.moving_branch =
                  SelectMovingInterfaceBranch(
                      hllc,
                      left_directional_state,
                      right_directional_state,
                      detail.face_speed)
                      .active_region;
            }
          }
          detail.left_component1_smoothness =
              interface_state.left_component1_smoothness;
          detail.right_component1_smoothness =
              interface_state.right_component1_smoothness;
          detail.left_width = line_work.effective_widths[left_line_cell];
          detail.right_width = line_work.effective_widths[right_line_cell];
          const std::size_t stencil_begin = left_line_cell >= 2u ? left_line_cell - 2u : 0u;
          const std::size_t stencil_end = std::min(
              right_line_cell + 2u,
              line_work.ghosted_line.size() - 1u);
          for (std::size_t stencil = stencil_begin; stencil <= stencil_end; ++stencil) {
            const auto primitive = ToDirectionalPrimitive(
                line_work.ghosted_line[stencil],
                HydroDirection::radial);
            detail.stencil_logical_radial.push_back(
                static_cast<int>(stencil) - static_cast<int>(ppm_ghost_layers));
            detail.stencil_rho.push_back(primitive.rho);
            detail.stencil_v_n.push_back(primitive.v_n);
            detail.stencil_pressure.push_back(primitive.pressure);
          }
          if (radial_face > 0u) {
            detail.left_coarse_index = line_work.line_cell_indices[radial_face - 1u];
          }
          if (radial_face < line_work.line_cell_indices.size()) {
            detail.right_coarse_index = line_work.line_cell_indices[radial_face];
          }
          }
          record_radial_flux_spread(
              radial_face,
              mass_flux,
              detail,
              radial_mass_flux_min,
              radial_mass_flux_max,
              radial_mass_flux_sum,
              radial_mass_flux_count,
              radial_mass_flux_min_detail,
              radial_mass_flux_max_detail);
          record_radial_flux_spread(
              radial_face,
              mom_r_flux,
              detail,
              radial_mom_r_flux_min,
              radial_mom_r_flux_max,
              radial_mom_r_flux_sum,
              radial_mom_r_flux_count,
              radial_mom_r_flux_min_detail,
              radial_mom_r_flux_max_detail);
          record_radial_flux_spread(
              radial_face,
              e_total_flux,
              detail,
              radial_e_total_flux_min,
              radial_e_total_flux_max,
              radial_e_total_flux_sum,
              radial_e_total_flux_count,
              radial_e_total_flux_min_detail,
              radial_e_total_flux_max_detail);
          if (options.debug_macro_radial_ppm_face_smoothness &&
              radial_face == options.debug_macro_radial_ppm_face &&
              detail.valid) {
            debug_face_smoothness_details.push_back(detail);
          }
        }
      }

      if (summary != nullptr) {
        summary->radial_ppm_summary.executed = true;
        summary->radial_ppm_summary.ghost_layers = ppm_ghost_layers;
        summary->radial_ppm_summary.line_count += 1u;
        summary->radial_ppm_summary.downgraded_interface_count +=
            line_work.downgraded_interface_count;
      }

      for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
        const std::size_t coarse_index = line_work.line_cell_indices[radial];
        auto face_minus_per_solid_angle = Scale(
            radial_interface_fluxes[radial],
            RadialFaceAreaPerSolidAngle(geometry, radial));
        auto face_plus_per_solid_angle = Scale(
            radial_interface_fluxes[radial + 1u],
            RadialFaceAreaPerSolidAngle(geometry, radial + 1u));
        if (macro_ale_direct_hllc_mode && options.apply_geometric_source) {
          const auto primitive = RecoverPrimitiveState(coarse_states[coarse_index]);
          if (!primitive.is_physical()) {
            failure_reason =
                "direct moving-face HLLC pressure-balanced radial flux encountered non-physical coarse state";
            macro_ale_hllc_diagnostics->failure_class =
                MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
            return false;
          }
          if (radial_moving_interface_pressure_terms.size() !=
              radial_interface_fluxes.size()) {
            failure_reason =
                "direct moving-face HLLC pressure-balanced radial flux is missing interface pressure terms";
            macro_ale_hllc_diagnostics->failure_class =
                MacroAleHllcFailureClass::required_diagnostic_missing;
            return false;
          }
          const bool cell_has_moving_radial_face =
              macro_radial_ale_window.radial_face_velocities[radial] != 0.0 ||
              macro_radial_ale_window.radial_face_velocities[radial + 1u] != 0.0;
          const double pressure_balance =
              cell_has_moving_radial_face
                  ? FaceConsistentPressureBalance(
                        radial_moving_interface_pressure_terms[radial],
                        radial_moving_interface_pressure_terms[radial + 1u])
                  : primitive.pressure;
          if (!(pressure_balance > 0.0) || !std::isfinite(pressure_balance)) {
            failure_reason =
                "direct moving-face HLLC pressure-balanced radial flux encountered non-physical face-consistent pressure";
            macro_ale_hllc_diagnostics->failure_class =
                MacroAleHllcFailureClass::moving_hllc_nonphysical;
            return false;
          }
          face_minus_per_solid_angle.mom_r -=
              pressure_balance * RadialFaceAreaPerSolidAngle(geometry, radial);
          face_plus_per_solid_angle.mom_r -=
              pressure_balance * RadialFaceAreaPerSolidAngle(geometry, radial + 1u);
        }
        const auto per_solid_angle_delta = Scale(
            Subtract(face_minus_per_solid_angle, face_plus_per_solid_angle),
            dt_s);
        const double angular_factor =
            CellSolidAngleFactor(geometry, line_work.theta, line_work.phi);
        if (!(angular_factor > 0.0) || !std::isfinite(angular_factor)) {
          failure_reason = "macro radial PPM encountered invalid angular factor";
          return false;
        }
        const auto extensive_delta =
            Scale(per_solid_angle_delta, angular_factor);
        coarse_extensive_delta[coarse_index] =
            Add(coarse_extensive_delta[coarse_index], extensive_delta);
        if (!coarse_radial_delta_has_reference[coarse_index]) {
          coarse_radial_delta_per_solid_angle_reference[coarse_index] =
              per_solid_angle_delta;
          coarse_radial_delta_has_reference[coarse_index] = true;
        }
        coarse_radial_delta_per_solid_angle_weighted_delta_sum[coarse_index] =
            Add(
                coarse_radial_delta_per_solid_angle_weighted_delta_sum[coarse_index],
                Scale(
                    Subtract(
                        per_solid_angle_delta,
                        coarse_radial_delta_per_solid_angle_reference[coarse_index]),
                    angular_factor));
        coarse_radial_delta_extensive[coarse_index] =
            Add(coarse_radial_delta_extensive[coarse_index], extensive_delta);
        coarse_radial_delta_solid_angle_weight[coarse_index] += angular_factor;
      }
    }

    if (collect_radial_flux_diagnostics) {
      append_primitive_spread_report(
          "a4_radial_ppm_input_cell_center",
          "macro_radial_line_cell_centers",
          radial_input_rho_min,
          radial_input_rho_max,
          radial_input_vn_min,
          radial_input_vn_max,
          radial_input_pressure_min,
          radial_input_pressure_max,
          radial_input_count);
      append_primitive_spread_report(
          "a4_radial_ppm_left_edge_before_trace",
          "macro_radial_interfaces_left_parabolic_edge",
          edge_left_rho_min,
          edge_left_rho_max,
          edge_left_vn_min,
          edge_left_vn_max,
          edge_left_pressure_min,
          edge_left_pressure_max,
          edge_interface_count);
      append_primitive_spread_report(
          "a4_radial_ppm_right_edge_before_trace",
          "macro_radial_interfaces_right_parabolic_edge",
          edge_right_rho_min,
          edge_right_rho_max,
          edge_right_vn_min,
          edge_right_vn_max,
          edge_right_pressure_min,
          edge_right_pressure_max,
          edge_interface_count);
      append_primitive_spread_report(
          "a4_radial_ppm_left_traced_interface",
          "macro_radial_interfaces_left_state",
          traced_left_rho_min,
          traced_left_rho_max,
          traced_left_vn_min,
          traced_left_vn_max,
          traced_left_pressure_min,
          traced_left_pressure_max,
          traced_interface_count);
      append_primitive_spread_report(
          "a4_radial_ppm_right_traced_interface",
          "macro_radial_interfaces_right_state",
          traced_right_rho_min,
          traced_right_rho_max,
          traced_right_vn_min,
          traced_right_vn_max,
          traced_right_pressure_min,
          traced_right_pressure_max,
          traced_interface_count);
      struct FluxSpreadSummary {
        double relative_spread{0.0};
        std::size_t radial_face{0u};
        double min_flux{0.0};
        double max_flux{0.0};
        double mean_flux{0.0};
        double absolute_spread{0.0};
        std::size_t absolute_radial_face{0u};
        double absolute_min_flux{0.0};
        double absolute_max_flux{0.0};
        double absolute_mean_flux{0.0};
      };
      const auto summarize_flux_spread =
          [&](const std::vector<double>& min_flux,
              const std::vector<double>& max_flux,
              const std::vector<double>& sum_flux,
              const std::vector<std::size_t>& count_flux) {
            FluxSpreadSummary result;
            for (std::size_t radial_face = 0; radial_face <= map.fine_radial_cells; ++radial_face) {
              if (count_flux[radial_face] == 0u) {
                continue;
              }
              const double mean_flux =
                  sum_flux[radial_face] /
                  static_cast<double>(count_flux[radial_face]);
              const double relative_spread =
                  std::abs(mean_flux) > 0.0
                      ? (max_flux[radial_face] - min_flux[radial_face]) /
                            std::abs(mean_flux)
                      : 0.0;
              const double absolute_spread =
                  max_flux[radial_face] - min_flux[radial_face];
              if (relative_spread > result.relative_spread) {
                result.relative_spread = relative_spread;
                result.radial_face = radial_face;
                result.min_flux = min_flux[radial_face];
                result.max_flux = max_flux[radial_face];
                result.mean_flux = mean_flux;
              }
              if (absolute_spread > result.absolute_spread) {
                result.absolute_spread = absolute_spread;
                result.absolute_radial_face = radial_face;
                result.absolute_min_flux = min_flux[radial_face];
                result.absolute_max_flux = max_flux[radial_face];
                result.absolute_mean_flux = mean_flux;
              }
            }
            return result;
          };
      const auto mass_flux_spread = summarize_flux_spread(
          radial_mass_flux_min,
          radial_mass_flux_max,
          radial_mass_flux_sum,
          radial_mass_flux_count);
      const auto mom_r_flux_spread = summarize_flux_spread(
          radial_mom_r_flux_min,
          radial_mom_r_flux_max,
          radial_mom_r_flux_sum,
          radial_mom_r_flux_count);
      const auto e_total_flux_spread = summarize_flux_spread(
          radial_e_total_flux_min,
          radial_e_total_flux_max,
          radial_e_total_flux_sum,
          radial_e_total_flux_count);

      std::ostringstream flux_report;
      flux_report << std::setprecision(17)
                  << "stage=macro_radial_ppm_face_flux"
                  << "; representation=macro_radial_interfaces"
                  << "; valid=true"
                  << "; max_mass_flux_relative_spread="
                  << mass_flux_spread.relative_spread
                  << "; max_mass_flux_spread_radial_face="
                  << mass_flux_spread.radial_face
                  << "; min_mass_flux=" << mass_flux_spread.min_flux
                  << "; max_mass_flux=" << mass_flux_spread.max_flux
                  << "; mean_mass_flux=" << mass_flux_spread.mean_flux
                  << "; max_mass_flux_absolute_spread="
                  << mass_flux_spread.absolute_spread
                  << "; max_mass_flux_absolute_spread_radial_face="
                  << mass_flux_spread.absolute_radial_face
                  << "; min_mass_flux_absolute_face="
                  << mass_flux_spread.absolute_min_flux
                  << "; max_mass_flux_absolute_face="
                  << mass_flux_spread.absolute_max_flux
                  << "; mean_mass_flux_absolute_face="
                  << mass_flux_spread.absolute_mean_flux
                  << "; max_mom_r_flux_relative_spread="
                  << mom_r_flux_spread.relative_spread
                  << "; max_mom_r_flux_spread_radial_face="
                  << mom_r_flux_spread.radial_face
                  << "; min_mom_r_flux=" << mom_r_flux_spread.min_flux
                  << "; max_mom_r_flux=" << mom_r_flux_spread.max_flux
                  << "; mean_mom_r_flux=" << mom_r_flux_spread.mean_flux
                  << "; max_mom_r_flux_absolute_spread="
                  << mom_r_flux_spread.absolute_spread
                  << "; max_mom_r_flux_absolute_spread_radial_face="
                  << mom_r_flux_spread.absolute_radial_face
                  << "; min_mom_r_flux_absolute_face="
                  << mom_r_flux_spread.absolute_min_flux
                  << "; max_mom_r_flux_absolute_face="
                  << mom_r_flux_spread.absolute_max_flux
                  << "; mean_mom_r_flux_absolute_face="
                  << mom_r_flux_spread.absolute_mean_flux
                  << "; max_e_total_flux_relative_spread="
                  << e_total_flux_spread.relative_spread
                  << "; max_e_total_flux_spread_radial_face="
                  << e_total_flux_spread.radial_face
                  << "; min_e_total_flux=" << e_total_flux_spread.min_flux
                  << "; max_e_total_flux=" << e_total_flux_spread.max_flux
                  << "; mean_e_total_flux=" << e_total_flux_spread.mean_flux
                  << "; max_e_total_flux_absolute_spread="
                  << e_total_flux_spread.absolute_spread
                  << "; max_e_total_flux_absolute_spread_radial_face="
                  << e_total_flux_spread.absolute_radial_face
                  << "; min_e_total_flux_absolute_face="
                  << e_total_flux_spread.absolute_min_flux
                  << "; max_e_total_flux_absolute_face="
                  << e_total_flux_spread.absolute_max_flux
                  << "; mean_e_total_flux_absolute_face="
                  << e_total_flux_spread.absolute_mean_flux;
      const auto append_detail =
          [&](std::string_view prefix, const MacroRadialFluxLineDetail& detail) {
            if (!detail.valid) {
              flux_report << "; " << prefix << "_detail_valid=false";
              return;
            }
            flux_report << "; " << prefix << "_detail_valid=true"
                        << "; " << prefix << "_theta=" << detail.theta
                        << "; " << prefix << "_phi=" << detail.phi
                        << "; " << prefix << "_mass_flux=" << detail.mass_flux
                        << "; " << prefix << "_mom_r_flux=" << detail.mom_r_flux
                        << "; " << prefix << "_e_total_flux=" << detail.e_total_flux
                        << "; " << prefix << "_face_speed=" << detail.face_speed
                        << "; " << prefix << "_s_left=" << detail.s_left
                        << "; " << prefix << "_s_star=" << detail.s_star
                        << "; " << prefix << "_s_right=" << detail.s_right
                        << "; " << prefix << "_moving_branch="
                        << HllcActiveRegionName(detail.moving_branch)
                        << "; " << prefix << "_downgraded_to_first_order="
                        << (detail.downgraded_to_first_order ? "true" : "false")
                        << "; " << prefix << "_shared_shell_fallback="
                        << (detail.shared_shell_fallback ? "true" : "false")
                        << "; " << prefix << "_left_profile_troubled="
                        << (detail.left_profile_troubled ? "true" : "false")
                        << "; " << prefix << "_right_profile_troubled="
                        << (detail.right_profile_troubled ? "true" : "false")
                        << "; " << prefix << "_left_trace_failed="
                        << (detail.left_trace_failed ? "true" : "false")
                        << "; " << prefix << "_right_trace_failed="
                        << (detail.right_trace_failed ? "true" : "false")
                        << "; " << prefix << "_left_troubled_component="
                        << detail.left_troubled_component
                        << "; " << prefix << "_right_troubled_component="
                        << detail.right_troubled_component
                        << "; " << prefix << "_left_trouble_q_im1="
                        << detail.left_trouble_q_im1
                        << "; " << prefix << "_left_trouble_q_i="
                        << detail.left_trouble_q_i
                        << "; " << prefix << "_left_trouble_q_ip1="
                        << detail.left_trouble_q_ip1
                        << "; " << prefix << "_left_trouble_total_variation="
                        << detail.left_trouble_total_variation
                        << "; " << prefix << "_left_trouble_end_to_end="
                        << detail.left_trouble_end_to_end
                        << "; " << prefix << "_right_trouble_q_im1="
                        << detail.right_trouble_q_im1
                        << "; " << prefix << "_right_trouble_q_i="
                        << detail.right_trouble_q_i
                        << "; " << prefix << "_right_trouble_q_ip1="
                        << detail.right_trouble_q_ip1
                        << "; " << prefix << "_right_trouble_total_variation="
                        << detail.right_trouble_total_variation
                        << "; " << prefix << "_right_trouble_end_to_end="
                        << detail.right_trouble_end_to_end
                        << "; " << prefix << "_right_component1_delta_left="
                        << detail.right_component1_smoothness.delta_left
                        << "; " << prefix << "_right_component1_delta_right="
                        << detail.right_component1_smoothness.delta_right
                        << "; " << prefix << "_right_component1_reversed_slope="
                        << detail.right_component1_smoothness.reversed_slope
                        << "; " << prefix << "_right_component1_tolerance="
                        << detail.right_component1_smoothness.tolerance
                        << "; " << prefix << "_left_width=" << detail.left_width
                        << "; " << prefix << "_right_width=" << detail.right_width
                        << "; " << prefix << "_left_coarse_index=" << detail.left_coarse_index
                        << "; " << prefix << "_right_coarse_index=" << detail.right_coarse_index;
            AppendPrimitiveReport(flux_report, std::string(prefix) + "_left_traced", detail.left_traced);
            AppendPrimitiveReport(flux_report, std::string(prefix) + "_right_traced", detail.right_traced);
            AppendPrimitiveReport(flux_report, std::string(prefix) + "_left_cell", detail.left_cell_center);
            AppendPrimitiveReport(flux_report, std::string(prefix) + "_right_cell", detail.right_cell_center);
            flux_report << "; " << prefix << "_stencil_logical_radial=";
            for (std::size_t index = 0; index < detail.stencil_logical_radial.size(); ++index) {
              if (index > 0u) {
                flux_report << ",";
              }
              flux_report << detail.stencil_logical_radial[index];
            }
            flux_report << "; " << prefix << "_stencil_rho=";
            for (std::size_t index = 0; index < detail.stencil_rho.size(); ++index) {
              if (index > 0u) {
                flux_report << ",";
              }
              flux_report << detail.stencil_rho[index];
            }
            flux_report << "; " << prefix << "_stencil_v_n=";
            for (std::size_t index = 0; index < detail.stencil_v_n.size(); ++index) {
              if (index > 0u) {
                flux_report << ",";
              }
              flux_report << detail.stencil_v_n[index];
            }
            flux_report << "; " << prefix << "_stencil_pressure=";
            for (std::size_t index = 0; index < detail.stencil_pressure.size(); ++index) {
              if (index > 0u) {
                flux_report << ",";
              }
              flux_report << detail.stencil_pressure[index];
            }
            if (detail.left_coarse_index < map.coarse_cells.size()) {
              AppendMacroCellReport(
                  flux_report,
                  std::string(prefix) + "_left_cell_meta",
                  map.coarse_cells[detail.left_coarse_index]);
            }
            if (detail.right_coarse_index < map.coarse_cells.size()) {
              AppendMacroCellReport(
                  flux_report,
                  std::string(prefix) + "_right_cell_meta",
                  map.coarse_cells[detail.right_coarse_index]);
            }
          };
      append_detail(
          "min_mass_flux_line",
          radial_mass_flux_min_detail[mass_flux_spread.radial_face]);
      append_detail(
          "max_mass_flux_line",
          radial_mass_flux_max_detail[mass_flux_spread.radial_face]);
      append_detail(
          "min_abs_mass_flux_line",
          radial_mass_flux_min_detail[mass_flux_spread.absolute_radial_face]);
      append_detail(
          "max_abs_mass_flux_line",
          radial_mass_flux_max_detail[mass_flux_spread.absolute_radial_face]);
      append_detail(
          "min_mom_r_flux_line",
          radial_mom_r_flux_min_detail[mom_r_flux_spread.radial_face]);
      append_detail(
          "max_mom_r_flux_line",
          radial_mom_r_flux_max_detail[mom_r_flux_spread.radial_face]);
      append_detail(
          "min_abs_mom_r_flux_line",
          radial_mom_r_flux_min_detail[mom_r_flux_spread.absolute_radial_face]);
      append_detail(
          "max_abs_mom_r_flux_line",
          radial_mom_r_flux_max_detail[mom_r_flux_spread.absolute_radial_face]);
      append_detail(
          "min_e_total_flux_line",
          radial_e_total_flux_min_detail[e_total_flux_spread.radial_face]);
      append_detail(
          "max_e_total_flux_line",
          radial_e_total_flux_max_detail[e_total_flux_spread.radial_face]);
      append_detail(
          "min_abs_e_total_flux_line",
          radial_e_total_flux_min_detail[e_total_flux_spread.absolute_radial_face]);
      append_detail(
          "max_abs_e_total_flux_line",
          radial_e_total_flux_max_detail[e_total_flux_spread.absolute_radial_face]);
      summary->angular_stage_report_lines.push_back(flux_report.str());

      for (const auto& detail : debug_face_smoothness_details) {
        std::ostringstream smoothness_report;
        smoothness_report << std::setprecision(17)
                          << "stage=macro_radial_ppm_component1_smoothness"
                          << "; representation=macro_radial_interfaces"
                          << "; valid=true"
                          << "; radial_face=" << options.debug_macro_radial_ppm_face
                          << "; theta=" << detail.theta
                          << "; phi=" << detail.phi
                          << "; mass_flux=" << detail.mass_flux
                          << "; downgraded_to_first_order="
                          << (detail.downgraded_to_first_order ? "true" : "false")
                          << "; shared_shell_fallback="
                          << (detail.shared_shell_fallback ? "true" : "false")
                          << "; right_profile_troubled="
                          << (detail.right_profile_troubled ? "true" : "false")
                          << "; right_trace_failed="
                          << (detail.right_trace_failed ? "true" : "false")
                          << "; q_im1=" << detail.right_component1_smoothness.q_im1
                          << "; q_i=" << detail.right_component1_smoothness.q_i
                          << "; q_ip1=" << detail.right_component1_smoothness.q_ip1
                          << "; delta_left=" << detail.right_component1_smoothness.delta_left
                          << "; delta_right=" << detail.right_component1_smoothness.delta_right
                          << "; sign_reversal="
                          << (detail.right_component1_smoothness.sign_reversal ? "true" : "false")
                          << "; primary_slope=" << detail.right_component1_smoothness.primary_slope
                          << "; reversed_slope=" << detail.right_component1_smoothness.reversed_slope
                          << "; tolerance=" << detail.right_component1_smoothness.tolerance;
        summary->angular_stage_report_lines.push_back(smoothness_report.str());
      }
    }
    if (summary != nullptr) {
      summary->radial_update_wall_s += ElapsedHydroSecondsSince(macro_radial_start);
    }
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "post_macro_radial_coarse_delta",
        "macro_coarse_state_plus_extensive_delta",
        BuildCoarseAngularStageMetrics(map, coarse_states, coarse_extensive_delta)));
    if (macro_ale_direct_hllc_mode) {
      summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
          "a4_after_radial_budget_vq_prediction",
          "macro_coarse_vq_prediction_new_geometry",
          BuildCoarseA4ExtensiveStageMetrics(
              map,
              coarse_states,
              coarse_extensive_delta,
              coarse_radial_delta_per_solid_angle_reference,
              coarse_radial_delta_per_solid_angle_weighted_delta_sum,
              coarse_radial_delta_extensive,
              coarse_radial_delta_solid_angle_weight,
              geometry,
              *ale_updated_geometry)));
    }
  }

  std::vector<std::vector<dec3d::mesh::MacroZoneThetaBand>> theta_bands_by_radial(
      map.fine_radial_cells);
  for (const auto& band : map.theta_bands) {
    if (band.radial_index < theta_bands_by_radial.size()) {
      theta_bands_by_radial[band.radial_index].push_back(band);
    }
  }

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    const auto& radial_bands = theta_bands_by_radial[radial];
    std::vector<std::vector<std::size_t>> radial_band_cells;
    radial_band_cells.reserve(radial_bands.size());
    for (const auto& band : radial_bands) {
      radial_band_cells.push_back(CoarseCellIndicesForBand(
          map,
          radial,
          band.theta_begin,
          band.theta_end));
    }

    if (options.apply_phi_sweep) {
      const auto macro_phi_start = std::chrono::steady_clock::now();
      std::vector<HydroConservativeState> phi_ghosted_line;
      std::vector<double> phi_effective_widths;
      std::vector<HydroConservativeState> phi_interface_fluxes;
      for (std::size_t band_index = 0; band_index < radial_bands.size(); ++band_index) {
        const auto& band = radial_bands[band_index];
        const auto& band_cells = radial_band_cells[band_index];
        if (band_cells.empty()) {
          failure_reason = "macro-zoning phi sweep could not find coarse band cells";
          return false;
        }

        const double face_area = CoarsePhiFaceArea(
            geometry,
            radial,
            band.theta_begin,
            band.theta_end);
        phi_interface_fluxes.resize(band_cells.size() + 1u);
        if (options.use_ppm_reconstruction &&
            options.macro_zoning_use_phi_ppm) {
          std::size_t downgraded_interface_count = 0u;
          if (!BuildMacroPhiPpmLine(
                  band_cells,
                  map,
                  coarse_states,
                  geometry,
                  ppm_ghost_layers,
                  phi_ghosted_line,
                  phi_effective_widths,
                  failure_reason)) {
            return false;
          }

          if (!BuildPpmInterfaceFluxes(
                  phi_ghosted_line,
                  SweepDirection::phi,
                  band_cells.size(),
                  ppm_ghost_layers,
                  dt_s,
                  phi_effective_widths,
                  phi_interface_fluxes,
                  downgraded_interface_count,
                  failure_reason)) {
            return false;
          }
          if (summary != nullptr) {
            summary->phi_ppm_summary.executed = true;
            summary->phi_ppm_summary.ghost_layers = ppm_ghost_layers;
            summary->phi_ppm_summary.line_count += 1u;
            summary->phi_ppm_summary.downgraded_interface_count += downgraded_interface_count;
          }
        } else {
          for (std::size_t local = 0; local < band_cells.size(); ++local) {
            const std::size_t left_index = band_cells[local];
            const std::size_t right_index = band_cells[(local + 1u) % band_cells.size()];
            phi_interface_fluxes[local + 1u] = ComputeDirectionalInterfaceFlux(
                coarse_states[left_index],
                coarse_states[right_index],
                SweepDirection::phi);
            if (!phi_interface_fluxes[local + 1u].is_finite()) {
              failure_reason = "macro-zoning phi coarse interface flux is non-finite";
              return false;
            }
          }
          phi_interface_fluxes[0] = phi_interface_fluxes.back();
        }

        for (std::size_t local = 0; local < band_cells.size(); ++local) {
          const std::size_t coarse_index = band_cells[local];
          const auto face_minus = Scale(phi_interface_fluxes[local], dt_s * face_area);
          const auto face_plus = Scale(phi_interface_fluxes[local + 1u], dt_s * face_area);
          coarse_extensive_delta[coarse_index] =
              Add(coarse_extensive_delta[coarse_index], Subtract(face_minus, face_plus));
        }
      }
      if (summary != nullptr) {
        summary->phi_update_wall_s += ElapsedHydroSecondsSince(macro_phi_start);
      }
    }

    if (summary != nullptr &&
        options.debug_angular_stage_diagnostics &&
        macro_ale_direct_hllc_mode &&
        options.apply_phi_sweep) {
      std::ostringstream stage;
      stage << "a4_after_phi_sweep_radial_" << radial
            << "_budget_vq_prediction";
      summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
          stage.str(),
          "macro_coarse_vq_prediction_new_geometry",
          BuildCoarseA4ExtensiveStageMetrics(
              map,
              coarse_states,
              coarse_extensive_delta,
              coarse_radial_delta_per_solid_angle_reference,
              coarse_radial_delta_per_solid_angle_weighted_delta_sum,
              coarse_radial_delta_extensive,
              coarse_radial_delta_solid_angle_weight,
              geometry,
              *ale_updated_geometry)));
    }

    if (options.apply_theta_sweep && radial_bands.size() > 1u) {
      const auto macro_theta_start = std::chrono::steady_clock::now();
      if (options.use_ppm_reconstruction &&
          options.macro_zoning_use_theta_ppm &&
          radial_bands.size() >= ppm_ghost_layers) {
        std::vector<std::size_t> theta_line_cell_indices;
        std::vector<HydroConservativeState> theta_ghosted_line;
        std::vector<double> theta_effective_widths;
        std::vector<HydroConservativeState> theta_interface_fluxes;
        for (std::size_t fine_phi = 0; fine_phi < map.fine_phi_cells; ++fine_phi) {
          theta_interface_fluxes.resize(radial_bands.size() + 1u);
          std::size_t downgraded_interface_count = 0u;
          if (!BuildMacroThetaPpmLine(
                  radial_bands,
                  radial_band_cells,
                  map,
                  coarse_states,
                  geometry,
                  fine_phi,
                  ppm_ghost_layers,
                  theta_line_cell_indices,
                  theta_ghosted_line,
                  theta_effective_widths,
                  failure_reason)) {
            return false;
          }

          if (!BuildPpmInterfaceFluxes(
                  theta_ghosted_line,
                  SweepDirection::theta,
                  radial_bands.size(),
                  ppm_ghost_layers,
                  dt_s,
                  theta_effective_widths,
                  theta_interface_fluxes,
                  downgraded_interface_count,
                  failure_reason)) {
            return false;
          }
          if (summary != nullptr) {
            summary->theta_ppm_summary.executed = true;
            summary->theta_ppm_summary.ghost_layers = ppm_ghost_layers;
            summary->theta_ppm_summary.line_count += 1u;
            summary->theta_ppm_summary.downgraded_interface_count += downgraded_interface_count;
          }

          for (std::size_t band_index = 0; band_index + 1u < radial_bands.size(); ++band_index) {
            const double face_area = CoarseThetaFaceArea(
                geometry,
                radial,
                radial_bands[band_index].theta_end,
                fine_phi,
                fine_phi + 1u);
            auto lower_flux = theta_interface_fluxes[band_index + 1u];
            auto upper_flux = theta_interface_fluxes[band_index + 1u];
            if (macro_ale_direct_hllc_mode && options.apply_geometric_source) {
              const auto lower_primitive =
                  RecoverPrimitiveState(coarse_states[theta_line_cell_indices[band_index]]);
              const auto upper_primitive =
                  RecoverPrimitiveState(coarse_states[theta_line_cell_indices[band_index + 1u]]);
              if (!lower_primitive.is_physical() || !upper_primitive.is_physical()) {
                failure_reason =
                    "direct moving-face HLLC pressure-balanced theta flux encountered non-physical coarse state";
                macro_ale_hllc_diagnostics->failure_class =
                    MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
                return false;
              }
              lower_flux.mom_theta -= lower_primitive.pressure;
              upper_flux.mom_theta -= upper_primitive.pressure;
            }
            const auto lower_extensive_flux = Scale(lower_flux, dt_s * face_area);
            const auto upper_extensive_flux = Scale(upper_flux, dt_s * face_area);
            coarse_extensive_delta[theta_line_cell_indices[band_index]] =
                Subtract(coarse_extensive_delta[theta_line_cell_indices[band_index]], lower_extensive_flux);
            coarse_extensive_delta[theta_line_cell_indices[band_index + 1u]] =
                Add(coarse_extensive_delta[theta_line_cell_indices[band_index + 1u]], upper_extensive_flux);
          }
        }
      } else {
        for (std::size_t band_index = 0; band_index + 1u < radial_bands.size(); ++band_index) {
          const auto& lower_band = radial_bands[band_index];
          const auto& lower_cells = radial_band_cells[band_index];
          const auto& upper_cells = radial_band_cells[band_index + 1u];
          if (lower_cells.empty() || upper_cells.empty()) {
            failure_reason = "macro-zoning theta sweep could not find adjacent coarse band cells";
            return false;
          }

          std::size_t lower_local = 0u;
          std::size_t upper_local = 0u;
          while (lower_local < lower_cells.size() && upper_local < upper_cells.size()) {
            const auto& lower_cell = map.coarse_cells[lower_cells[lower_local]];
            const auto& upper_cell = map.coarse_cells[upper_cells[upper_local]];
            const std::size_t overlap_begin = std::max(lower_cell.phi_begin, upper_cell.phi_begin);
            const std::size_t overlap_end = std::min(lower_cell.phi_end, upper_cell.phi_end);

            if (overlap_begin < overlap_end) {
              const auto interface_flux = ComputeDirectionalInterfaceFlux(
                  coarse_states[lower_cells[lower_local]],
                  coarse_states[upper_cells[upper_local]],
                  SweepDirection::theta);
              if (!interface_flux.is_finite()) {
                failure_reason = "macro-zoning theta coarse interface flux is non-finite";
                return false;
              }

              const double face_area = CoarseThetaFaceArea(
                  geometry,
                  radial,
                  lower_band.theta_end,
                  overlap_begin,
                  overlap_end);
              auto lower_flux = interface_flux;
              auto upper_flux = interface_flux;
              if (macro_ale_direct_hllc_mode && options.apply_geometric_source) {
                const auto lower_primitive =
                    RecoverPrimitiveState(coarse_states[lower_cells[lower_local]]);
                const auto upper_primitive =
                    RecoverPrimitiveState(coarse_states[upper_cells[upper_local]]);
                if (!lower_primitive.is_physical() || !upper_primitive.is_physical()) {
                  failure_reason =
                      "direct moving-face HLLC pressure-balanced theta flux encountered non-physical coarse state";
                  macro_ale_hllc_diagnostics->failure_class =
                      MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
                  return false;
                }
                lower_flux.mom_theta -= lower_primitive.pressure;
                upper_flux.mom_theta -= upper_primitive.pressure;
              }
              const auto lower_extensive_flux = Scale(lower_flux, dt_s * face_area);
              const auto upper_extensive_flux = Scale(upper_flux, dt_s * face_area);
              coarse_extensive_delta[lower_cells[lower_local]] =
                  Subtract(coarse_extensive_delta[lower_cells[lower_local]], lower_extensive_flux);
              coarse_extensive_delta[upper_cells[upper_local]] =
                  Add(coarse_extensive_delta[upper_cells[upper_local]], upper_extensive_flux);
            }

            if (lower_cell.phi_end <= upper_cell.phi_end) {
              ++lower_local;
            }
            if (upper_cell.phi_end <= lower_cell.phi_end) {
              ++upper_local;
            }
          }
        }
      }
      if (summary != nullptr) {
        summary->theta_update_wall_s += ElapsedHydroSecondsSince(macro_theta_start);
      }
    }

    if (summary != nullptr &&
        options.debug_angular_stage_diagnostics &&
        macro_ale_direct_hllc_mode &&
        options.apply_theta_sweep) {
      std::ostringstream stage;
      stage << "a4_after_theta_sweep_radial_" << radial
            << "_budget_vq_prediction";
      summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
          stage.str(),
          "macro_coarse_vq_prediction_new_geometry",
          BuildCoarseA4ExtensiveStageMetrics(
              map,
              coarse_states,
              coarse_extensive_delta,
              coarse_radial_delta_per_solid_angle_reference,
              coarse_radial_delta_per_solid_angle_weighted_delta_sum,
              coarse_radial_delta_extensive,
              coarse_radial_delta_solid_angle_weight,
              geometry,
              *ale_updated_geometry)));
    }
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "post_macro_angular_coarse_delta",
        "macro_coarse_state_plus_extensive_delta",
        BuildCoarseAngularStageMetrics(map, coarse_states, coarse_extensive_delta)));
    if (macro_ale_direct_hllc_mode) {
      summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
          "a4_after_angular_budget_vq_prediction",
          "macro_coarse_vq_prediction_new_geometry",
          BuildCoarseA4ExtensiveStageMetrics(
              map,
              coarse_states,
              coarse_extensive_delta,
              coarse_radial_delta_per_solid_angle_reference,
              coarse_radial_delta_per_solid_angle_weighted_delta_sum,
              coarse_radial_delta_extensive,
              coarse_radial_delta_solid_angle_weight,
              geometry,
              *ale_updated_geometry)));
    }
    auto coarse_flux_source_delta = coarse_extensive_delta;
    auto debug_source_extensive_delta = coarse_source_extensive_delta;
    if (!coarse_source_budget_integrated) {
      debug_source_extensive_delta =
          BuildCoarseGeometricSourceExtensiveDelta(
              map,
              coarse_states,
              geometry,
              dt_s,
              false);
      for (std::size_t coarse_index = 0; coarse_index < coarse_flux_source_delta.size(); ++coarse_index) {
        coarse_flux_source_delta[coarse_index] = Add(
            coarse_flux_source_delta[coarse_index],
            debug_source_extensive_delta[coarse_index]);
      }
    }
    auto debug_flux_only_extensive_delta = coarse_extensive_delta;
    if (coarse_source_budget_integrated) {
      for (std::size_t coarse_index = 0; coarse_index < debug_flux_only_extensive_delta.size(); ++coarse_index) {
        debug_flux_only_extensive_delta[coarse_index] = Subtract(
            debug_flux_only_extensive_delta[coarse_index],
            debug_source_extensive_delta[coarse_index]);
      }
    }
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "post_macro_geometric_source_coarse_prediction",
        "macro_coarse_state_plus_flux_and_source_extensive_delta",
        BuildCoarseAngularStageMetrics(map, coarse_states, coarse_flux_source_delta)));
    summary->angular_stage_report_lines.push_back(FormatCoarseGeometricSourceBalance(
        map,
        coarse_states,
        debug_flux_only_extensive_delta,
        debug_source_extensive_delta));
  }

  const auto macro_state_update_start = std::chrono::steady_clock::now();
  auto updated_states = coarse_states;
  for (std::size_t coarse_index = 0; coarse_index < updated_states.size(); ++coarse_index) {
    if (macro_ale_direct_hllc_mode) {
      updated_states[coarse_index] = BuildMacroRadialFactorizedCoarseState(
          coarse_states[coarse_index],
          coarse_extensive_delta[coarse_index],
          coarse_radial_delta_per_solid_angle_reference[coarse_index],
          coarse_radial_delta_per_solid_angle_weighted_delta_sum[coarse_index],
          coarse_radial_delta_extensive[coarse_index],
          coarse_radial_delta_solid_angle_weight[coarse_index],
          geometry,
          *ale_updated_geometry,
          map.coarse_cells[coarse_index]);
    } else {
      updated_states[coarse_index] = BuildMacroRadialFactorizedCoarseState(
          coarse_states[coarse_index],
          coarse_extensive_delta[coarse_index],
          coarse_radial_delta_per_solid_angle_reference[coarse_index],
          coarse_radial_delta_per_solid_angle_weighted_delta_sum[coarse_index],
          coarse_radial_delta_extensive[coarse_index],
          coarse_radial_delta_solid_angle_weight[coarse_index],
          geometry,
          geometry,
          map.coarse_cells[coarse_index]);
    }
    if (!IsPhysicalCellState(updated_states[coarse_index])) {
      failure_reason = "macro-zoning coarse update produced a non-physical state";
      return false;
    }
  }
  if (summary != nullptr) {
    summary->state_update_wall_s += ElapsedHydroSecondsSince(macro_state_update_start);
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        macro_ale_direct_hllc_mode
            ? "a4_after_coarse_commit_state"
            : "post_macro_coarse_commit_state",
        macro_ale_direct_hllc_mode
            ? "macro_coarse_state_new_geometry"
            : "macro_coarse_state_old_geometry",
        BuildCoarseStateAngularStageMetrics(
            map,
            updated_states,
            macro_ale_direct_hllc_mode ? *ale_updated_geometry : geometry)));
  }

  StoreMacroZoneStates(updated_states, updated_coarse);
  if (summary != nullptr) {
    summary->update_wall_s += ElapsedHydroSecondsSince(macro_update_start);
  }

  const auto prolong_start = std::chrono::steady_clock::now();
  const auto prolonged = dec3d::mesh::ProlongMacroZoneHydroPackage(updated_coarse);
  if (!prolonged.is_complete()) {
    failure_reason = prolonged.failure_reason.empty()
                         ? "macro-zoning prolong failed"
                         : prolonged.failure_reason;
    return false;
  }

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(snapshot.shape, radial, theta, phi);
        const HydroConservativeState prolonged_state{
            prolonged.rho(radial, theta, phi),
            prolonged.mom_r(radial, theta, phi),
              prolonged.mom_theta(radial, theta, phi),
              prolonged.mom_phi(radial, theta, phi),
              prolonged.e_fluid_total(radial, theta, phi),
              prolonged.chi_e(radial, theta, phi),
              snapshot.cells[index].alpha_chi,
              snapshot.cells[index].radiation_chi};
        accumulated_delta[index] =
            Add(accumulated_delta[index], Subtract(prolonged_state, snapshot.cells[index]));
      }
    }
  }
  if (summary != nullptr) {
    summary->prolong_wall_s += ElapsedHydroSecondsSince(prolong_start);
  }

  if (summary != nullptr && options.debug_angular_stage_diagnostics) {
    summary->angular_stage_report_lines.push_back(FormatAngularStageMetrics(
        "post_macro_prolong_fine_delta",
        "fine_snapshot_plus_accumulated_delta",
        BuildFineAngularStageMetrics(snapshot, accumulated_delta, geometry)));
  }

  if (summary != nullptr) {
    summary->executed = true;
    summary->radial_executed =
        options.apply_radial_sweep &&
        options.use_ppm_reconstruction &&
        options.macro_zoning_use_radial_ppm;
    summary->theta_executed = options.apply_theta_sweep;
    summary->phi_executed = options.apply_phi_sweep;
    summary->coarse_cells = map.coarse_cells.size();
    summary->theta_bands = map.theta_bands.size();
    summary->detect_report_line = map.report_line;
    summary->restrict_report_line = coarse_package.report_line;

    std::ostringstream coarse_update;
    coarse_update << "coarse_cells=" << summary->coarse_cells
                  << "; theta_bands=" << summary->theta_bands
                  << "; radial_executed=" << (summary->radial_executed ? "true" : "false")
                  << "; theta_executed=" << (summary->theta_executed ? "true" : "false")
                  << "; phi_executed=" << (summary->phi_executed ? "true" : "false")
                  << "; radial_ppm=" << (options.macro_zoning_use_radial_ppm ? "true" : "false")
                  << "; theta_ppm=" << (options.macro_zoning_use_theta_ppm ? "true" : "false")
                  << "; phi_ppm=" << (options.macro_zoning_use_phi_ppm ? "true" : "false")
                  << "; coarse_factor=" << options.macro_zoning_coarse_factor;
    summary->coarse_update_report_line = coarse_update.str();

    std::ostringstream prolong_report;
    prolong_report << "coarse_cells=" << summary->coarse_cells
                   << "; prolongation=zeroth_order_constant_fill";
    summary->prolong_report_line = prolong_report.str();

    if (summary->radial_ppm_summary.executed) {
      std::ostringstream radial_report;
      radial_report << "direction=radial"
                    << "; ghost_layers=" << summary->radial_ppm_summary.ghost_layers
                    << "; line_count=" << summary->radial_ppm_summary.line_count
                    << "; downgraded_interface_count=" << summary->radial_ppm_summary.downgraded_interface_count
                    << "; macro_zoning=true";
      summary->radial_ppm_summary.report_line = radial_report.str();
    }
    if (summary->theta_ppm_summary.executed) {
      std::ostringstream theta_report;
      theta_report << "direction=theta"
                   << "; ghost_layers=" << summary->theta_ppm_summary.ghost_layers
                   << "; line_count=" << summary->theta_ppm_summary.line_count
                   << "; downgraded_interface_count=" << summary->theta_ppm_summary.downgraded_interface_count
                   << "; macro_zoning=true";
      summary->theta_ppm_summary.report_line = theta_report.str();
    }
    if (summary->phi_ppm_summary.executed) {
      std::ostringstream phi_report;
      phi_report << "direction=phi"
                 << "; ghost_layers=" << summary->phi_ppm_summary.ghost_layers
                 << "; line_count=" << summary->phi_ppm_summary.line_count
                 << "; downgraded_interface_count=" << summary->phi_ppm_summary.downgraded_interface_count
                 << "; macro_zoning=true";
      summary->phi_ppm_summary.report_line = phi_report.str();
    }
  }

  return true;
}

[[nodiscard]] dec3d::core::Array3D<HydroConservativeState> BuildSnapshotArray(
    const HydroStateSnapshot& snapshot) {
  dec3d::core::Array3D<HydroConservativeState> state_array(
      snapshot.shape.radial_cells,
      snapshot.shape.theta_cells,
      snapshot.shape.phi_cells);

  for (std::size_t radial = 0; radial < snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.shape.phi_cells; ++phi) {
        state_array(radial, theta, phi) =
            snapshot.cells[LinearIndex(snapshot.shape, radial, theta, phi)];
      }
    }
  }

  return state_array;
}

[[nodiscard]] bool AccumulateDirectionalSweepDelta(
    SweepDirection direction,
    const HydroStateSnapshot& snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const StaticGridHydroOptions& options,
    bool use_extensive_update,
    std::vector<HydroConservativeState>& accumulated_delta,
    std::string& failure_reason,
    const RadialGhostOverride* radial_ghost_override = nullptr,
    std::string* ghost_report_line = nullptr,
    DirectionalPpmSummary* ppm_summary = nullptr,
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal = nullptr,
    RadialAleFluxSummary* radial_ale_summary = nullptr,
    std::vector<HydroConservativeState>* radial_delta_per_solid_angle = nullptr,
    std::vector<HydroConservativeState>* radial_delta_extensive = nullptr) {
  const GridShape& shape = snapshot.shape;
  if (!snapshot.is_complete() ||
      accumulated_delta.size() != snapshot.cells.size()) {
    failure_reason = "hydrodynamic snapshot is incomplete";
    return false;
  }
  if ((radial_delta_per_solid_angle != nullptr &&
       radial_delta_per_solid_angle->size() != snapshot.cells.size()) ||
      (radial_delta_extensive != nullptr &&
       radial_delta_extensive->size() != snapshot.cells.size())) {
    failure_reason = "radial ALE factorized delta buffers must match snapshot size";
    return false;
  }

  dec3d::mesh::RadialAleLocalProposalWindow radial_ale_window;
  if (direction == SweepDirection::radial && options.apply_radial_ale_flux_correction) {
    if (radial_ale_proposal == nullptr) {
      failure_reason = "radial ALE flux correction requires a mesh update proposal";
      return false;
    }
    if (!radial_ale_proposal->radial_ale ||
        std::abs(radial_ale_proposal->dt_s - dt_s) > 1.0e-12) {
      failure_reason = "radial ALE flux correction requires a current-step proposal";
      return false;
    }
    radial_ale_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
        *radial_ale_proposal,
        options.radial_ale_global_face_begin_index,
        shape.radial_cells + 1u);
    if (!radial_ale_window.is_complete()) {
      failure_reason = radial_ale_window.failure_reason.empty()
                           ? "radial ALE flux correction requires a complete local proposal window"
                           : radial_ale_window.failure_reason;
      return false;
    }
  }

  const std::size_t ghost_layers =
      options.use_ppm_reconstruction ? options.reconstruction_ghost_layers : 1u;
  const auto snapshot_array = BuildSnapshotArray(snapshot);
  const auto prepared_ghosts = [&]() {
    if (direction == SweepDirection::radial) {
      const auto lower_boundary_mode = InferRadialLowerBoundaryMode(geometry);
      if (radial_ghost_override != nullptr) {
        return FillHydroRadialGhosts(
            snapshot_array,
            ghost_layers,
            *radial_ghost_override,
            lower_boundary_mode);
      }
      return FillHydroRadialGhosts(
          snapshot_array,
          ghost_layers,
          lower_boundary_mode);
    }

    return PrepareDirectionalGhosts(
        snapshot_array,
        ToHydroDirection(direction),
        ghost_layers);
  }();
  if (!prepared_ghosts.success) {
    failure_reason = prepared_ghosts.failure_reason.empty()
                         ? "directional ghost preparation failed"
                         : prepared_ghosts.failure_reason;
    return false;
  }
  if (ghost_report_line != nullptr) {
    *ghost_report_line = prepared_ghosts.report_line;
  }

  if (direction == SweepDirection::radial) {
    std::vector<bool> zero_area_radial_faces(shape.radial_cells + 1u, false);
    for (std::size_t radial_face = 0; radial_face <= shape.radial_cells; ++radial_face) {
      zero_area_radial_faces[radial_face] =
          RadialFaceAreaPerSolidAngle(geometry, radial_face) <= 0.0;
    }
    std::vector<RadialAleFaceFluxBucket> radial_ale_face_flux_debug;
    const bool collect_radial_ale_face_flux_debug =
        options.debug_angular_stage_diagnostics &&
        options.use_ppm_reconstruction &&
        options.apply_radial_ale_flux_correction &&
        radial_ale_proposal != nullptr;
    if (collect_radial_ale_face_flux_debug) {
      radial_ale_face_flux_debug.assign(
          shape.radial_cells + 1u,
          RadialAleFaceFluxBucket{});
    }
    std::size_t downgraded_interfaces = 0u;

    const auto accumulate_radial_delta =
        [&](std::size_t theta,
            std::size_t phi,
            const std::vector<HydroConservativeState>& interface_fluxes) {
          for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
            HydroConservativeState delta;
            if (use_extensive_update) {
              const auto face_plus = Scale(
                  interface_fluxes[radial + 1u],
                  RadialFaceArea(geometry, radial + 1u, theta, phi));
              const auto face_minus = Scale(
                  interface_fluxes[radial],
                  RadialFaceArea(geometry, radial, theta, phi));
              delta = Scale(Subtract(face_minus, face_plus), dt_s);
              if (radial_delta_per_solid_angle != nullptr &&
                  radial_delta_extensive != nullptr) {
                const auto face_plus_per_solid_angle = Scale(
                    interface_fluxes[radial + 1u],
                    RadialFaceAreaPerSolidAngle(geometry, radial + 1u));
                const auto face_minus_per_solid_angle = Scale(
                    interface_fluxes[radial],
                    RadialFaceAreaPerSolidAngle(geometry, radial));
                const auto per_solid_angle_delta = Scale(
                    Subtract(face_minus_per_solid_angle, face_plus_per_solid_angle),
                    dt_s);
                const auto index = LinearIndex(shape, radial, theta, phi);
                (*radial_delta_per_solid_angle)[index] = Add(
                    (*radial_delta_per_solid_angle)[index],
                    per_solid_angle_delta);
                (*radial_delta_extensive)[index] = Add(
                    (*radial_delta_extensive)[index],
                    delta);
              }
            } else {
              const double radial_volume_factor = RadialShellVolumeFactor(geometry, radial);
              if (!(radial_volume_factor > 0.0) || !std::isfinite(radial_volume_factor)) {
                failure_reason = "radial shell volume factor must remain positive";
                return false;
              }
              const auto face_plus = Scale(
                  interface_fluxes[radial + 1u],
                  RadialFaceAreaPerSolidAngle(geometry, radial + 1u));
              const auto face_minus = Scale(
                  interface_fluxes[radial],
                  RadialFaceAreaPerSolidAngle(geometry, radial));
              delta = Scale(Subtract(face_minus, face_plus), dt_s / radial_volume_factor);
            }
            accumulated_delta[LinearIndex(shape, radial, theta, phi)] =
                Add(accumulated_delta[LinearIndex(shape, radial, theta, phi)], delta);
          }
          return true;
        };

    if (options.use_ppm_reconstruction) {
      struct FineRadialPpmLineWork {
        std::size_t theta{0u};
        std::size_t phi{0u};
        std::vector<HydroConservativeState> ghosted_line;
        std::vector<double> effective_widths;
        PpmLineReconstructionResult line_reconstruction;
        std::vector<HydroConservativeState> interface_fluxes;
      };

      std::vector<FineRadialPpmLineWork> radial_ppm_lines;
      radial_ppm_lines.reserve(shape.theta_cells * shape.phi_cells);
      std::vector<bool> shell_fallback_mask(shape.radial_cells + 1u, false);

      for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
          FineRadialPpmLineWork line_work;
          line_work.theta = theta;
          line_work.phi = phi;
          line_work.ghosted_line =
              ExtractGhostedLine(prepared_ghosts.ghosted_states, direction, theta, phi);
          line_work.effective_widths =
              ExtractLineEffectiveWidths(geometry, shape, direction, theta, phi, ghost_layers);
          line_work.line_reconstruction = ReconstructPpmLine(
              line_work.ghosted_line,
              ToHydroDirection(direction),
              shape.radial_cells,
              ghost_layers,
              dt_s,
              line_work.effective_widths);
          if (!line_work.line_reconstruction.success) {
            failure_reason = line_work.line_reconstruction.failure_reason.empty()
                                 ? "radial PPM line reconstruction failed"
                                 : line_work.line_reconstruction.failure_reason;
            return false;
          }
          AccumulateShellWidePpmFallbackMask(
              line_work.line_reconstruction,
              shell_fallback_mask);
          radial_ppm_lines.push_back(std::move(line_work));
        }
      }

      for (auto& line_work : radial_ppm_lines) {
        line_work.interface_fluxes.assign(
            shape.radial_cells + 1u,
            HydroConservativeState{});
        if (!BuildPpmInterfaceFluxesFromReconstruction(
                line_work.line_reconstruction,
                line_work.ghosted_line,
                direction,
                ghost_layers,
                shell_fallback_mask,
                line_work.interface_fluxes,
                downgraded_interfaces,
                failure_reason,
                options.apply_radial_ale_flux_correction
                    ? &radial_ale_window.radial_face_velocities
                    : nullptr,
                &zero_area_radial_faces)) {
          return false;
        }

        if (collect_radial_ale_face_flux_debug) {
          for (std::size_t radial_face = 0;
               radial_face < line_work.line_reconstruction.interfaces.size() &&
               radial_face < line_work.interface_fluxes.size();
               ++radial_face) {
            if (radial_face < zero_area_radial_faces.size() &&
                zero_area_radial_faces[radial_face]) {
              continue;
            }
            const auto& interface_state =
                line_work.line_reconstruction.interfaces[radial_face];
            AccumulateRadialAleFaceFluxDebug(
                radial_ale_face_flux_debug,
                radial_face,
                line_work.theta,
                line_work.phi,
                radial_ale_window.radial_face_velocities[radial_face],
                interface_state,
                interface_state.downgraded_to_first_order ||
                    shell_fallback_mask[radial_face],
                line_work.interface_fluxes[radial_face]);
          }
        }

        if (!accumulate_radial_delta(
                line_work.theta,
                line_work.phi,
                line_work.interface_fluxes)) {
          return false;
        }
      }
    } else {
      for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
        for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
          std::vector<HydroConservativeState> interface_fluxes(shape.radial_cells + 1u);
          for (std::size_t radial_face = 0; radial_face <= shape.radial_cells; ++radial_face) {
            if (zero_area_radial_faces[radial_face]) {
              interface_fluxes[radial_face] = HydroConservativeState{};
              continue;
            }

            const auto left_state =
                ToDirectionalState(prepared_ghosts.ghosted_states(radial_face, theta, phi), direction);
            const auto right_state =
                ToDirectionalState(prepared_ghosts.ghosted_states(radial_face + 1u, theta, phi), direction);
            const auto hllc = SolveHllcRiemann(left_state, right_state);
            if (!hllc.is_complete()) {
              failure_reason = hllc.failure_reason.empty()
                                   ? "radial HLLC solve failed"
                                   : hllc.failure_reason;
              return false;
            }
            auto interface_flux = hllc.interface_flux;
            if (options.apply_radial_ale_flux_correction) {
              interface_flux = ComputeAleCorrectedFlux(
                  hllc,
                  left_state,
                  right_state,
                  radial_ale_window.radial_face_velocities[radial_face]);
            }
            interface_fluxes[radial_face] = FromDirectionalFlux(interface_flux, direction);
            if (!interface_fluxes[radial_face].is_finite()) {
              failure_reason = "radial interface flux is non-finite";
              return false;
            }
          }

          if (!accumulate_radial_delta(theta, phi, interface_fluxes)) {
            return false;
          }
        }
      }
    }

    if (ppm_summary != nullptr && options.use_ppm_reconstruction) {
      ppm_summary->executed = true;
      ppm_summary->ghost_layers = ghost_layers;
      ppm_summary->line_count = shape.theta_cells * shape.phi_cells;
      ppm_summary->downgraded_interface_count = downgraded_interfaces;
      std::ostringstream report;
      report << "direction=radial"
             << "; ghost_layers=" << ghost_layers
             << "; line_count=" << ppm_summary->line_count
             << "; downgraded_interface_count=" << downgraded_interfaces;
      ppm_summary->report_line = report.str();
    }
    if (radial_ale_summary != nullptr && options.apply_radial_ale_flux_correction) {
      radial_ale_summary->executed = true;
      double max_face_speed = 0.0;
      for (const double face_speed : radial_ale_window.radial_face_velocities) {
        max_face_speed = std::max(max_face_speed, std::abs(face_speed));
      }
      radial_ale_summary->max_face_speed = max_face_speed;
      std::ostringstream report;
      report << "implementation_id=" << radial_ale_proposal->implementation_id
             << "; radial_face_indexing=global"
             << "; global_face_begin=" << options.radial_ale_global_face_begin_index
             << "; local_radial_faces=" << radial_ale_window.radial_face_velocities.size()
             << "; max_face_speed=" << max_face_speed;
      radial_ale_summary->report_line = report.str();
      if (collect_radial_ale_face_flux_debug) {
        radial_ale_summary->debug_report_lines.push_back(
            BuildRadialAleFaceFluxDebugReport(
                radial_ale_face_flux_debug,
                radial_ale_window.radial_face_velocities));
      }
    }
    return true;
  }

  if (direction == SweepDirection::theta) {
    std::vector<HydroConservativeState> interface_fluxes(shape.theta_cells + 1u);
    std::size_t downgraded_interfaces = 0u;
    for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
      for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
        if (options.use_ppm_reconstruction) {
          const auto ghosted_line =
              ExtractGhostedLine(prepared_ghosts.ghosted_states, direction, radial, phi);
          const auto effective_widths =
              ExtractLineEffectiveWidths(geometry, shape, direction, radial, phi, ghost_layers);
          if (!BuildPpmInterfaceFluxes(
                  ghosted_line,
                  direction,
                  shape.theta_cells,
                  ghost_layers,
                  dt_s,
                  effective_widths,
                  interface_fluxes,
                  downgraded_interfaces,
                  failure_reason)) {
            return false;
          }
        } else {
          for (std::size_t theta_face = 0; theta_face <= shape.theta_cells; ++theta_face) {
            interface_fluxes[theta_face] = ComputeDirectionalInterfaceFlux(
                prepared_ghosts.ghosted_states(radial, theta_face, phi),
                prepared_ghosts.ghosted_states(radial, theta_face + 1u, phi),
                direction);
            if (!interface_fluxes[theta_face].is_finite()) {
              failure_reason = "theta interface flux is non-finite";
              return false;
            }
          }
        }

        for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
          const auto face_plus = Scale(interface_fluxes[theta + 1u], ThetaFaceArea(geometry, radial, theta + 1u, phi));
          const auto face_minus = Scale(interface_fluxes[theta], ThetaFaceArea(geometry, radial, theta, phi));
          const double update_scale =
              use_extensive_update
                  ? dt_s
                  : dt_s / CellVolume(geometry, shape, radial, theta, phi);
          const auto delta = Scale(Subtract(face_minus, face_plus), update_scale);
          accumulated_delta[LinearIndex(shape, radial, theta, phi)] =
              Add(accumulated_delta[LinearIndex(shape, radial, theta, phi)], delta);
        }
      }
    }

    if (ppm_summary != nullptr && options.use_ppm_reconstruction) {
      ppm_summary->executed = true;
      ppm_summary->ghost_layers = ghost_layers;
      ppm_summary->line_count = shape.radial_cells * shape.phi_cells;
      ppm_summary->downgraded_interface_count = downgraded_interfaces;
      std::ostringstream report;
      report << "direction=theta"
             << "; ghost_layers=" << ghost_layers
             << "; line_count=" << ppm_summary->line_count
             << "; downgraded_interface_count=" << downgraded_interfaces;
      ppm_summary->report_line = report.str();
    }
    return true;
  }

  std::size_t downgraded_interfaces = 0u;
  for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
      std::vector<HydroConservativeState> interface_fluxes(shape.phi_cells + 1u);
      if (options.use_ppm_reconstruction) {
        const auto ghosted_line =
            ExtractGhostedLine(prepared_ghosts.ghosted_states, direction, radial, theta);
        const auto effective_widths =
            ExtractLineEffectiveWidths(geometry, shape, direction, radial, theta, ghost_layers);
        if (!BuildPpmInterfaceFluxes(
                ghosted_line,
                direction,
                shape.phi_cells,
                ghost_layers,
                dt_s,
                effective_widths,
                interface_fluxes,
                downgraded_interfaces,
                failure_reason)) {
          return false;
        }
      } else {
        for (std::size_t phi_face = 0; phi_face <= shape.phi_cells; ++phi_face) {
          interface_fluxes[phi_face] = ComputeDirectionalInterfaceFlux(
              prepared_ghosts.ghosted_states(radial, theta, phi_face),
              prepared_ghosts.ghosted_states(radial, theta, phi_face + 1u),
              direction);
          if (!interface_fluxes[phi_face].is_finite()) {
            failure_reason = "phi interface flux is non-finite";
            return false;
          }
        }
      }

      for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
        const auto face_plus = Scale(interface_fluxes[phi + 1u], PhiFaceArea(geometry, radial, theta));
        const auto face_minus = Scale(interface_fluxes[phi], PhiFaceArea(geometry, radial, theta));
        const double update_scale =
            use_extensive_update
                ? dt_s
                : dt_s / CellVolume(geometry, shape, radial, theta, phi);
        const auto delta = Scale(Subtract(face_minus, face_plus), update_scale);
        accumulated_delta[LinearIndex(shape, radial, theta, phi)] =
            Add(accumulated_delta[LinearIndex(shape, radial, theta, phi)], delta);
      }
    }
  }

  if (ppm_summary != nullptr && options.use_ppm_reconstruction) {
    ppm_summary->executed = true;
    ppm_summary->ghost_layers = ghost_layers;
    ppm_summary->line_count = shape.radial_cells * shape.theta_cells;
    ppm_summary->downgraded_interface_count = downgraded_interfaces;
    std::ostringstream report;
    report << "direction=phi"
           << "; ghost_layers=" << ghost_layers
           << "; line_count=" << ppm_summary->line_count
           << "; downgraded_interface_count=" << downgraded_interfaces;
    ppm_summary->report_line = report.str();
  }
  return true;
}

[[nodiscard]] HydroBudgetResidualSummary BuildBudgetResidualSummary(
    const HydroStateSnapshot& old_snapshot,
    const std::vector<HydroConservativeState>& accumulated_delta,
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  HydroBudgetResidualSummary summary;
  if (!old_snapshot.is_complete() ||
      accumulated_delta.size() != old_snapshot.cells.size() ||
      !hydro_view.is_complete() ||
      !geometry.is_valid()) {
    return summary;
  }

  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double volume = CellVolume(geometry, old_snapshot.shape, radial, theta, phi);
        const auto& old_state = old_snapshot.cells[index];
        const auto& flux_delta = accumulated_delta[index];
        const auto flux_state = Add(old_state, flux_delta);
        const auto new_state = LoadCellState(hydro_view, radial, theta, phi);

        summary.old_mass += CellMass(old_state, volume);
        summary.flux_mass_delta += CellMass(flux_delta, volume);
        summary.source_mass_delta += CellMass(Subtract(new_state, flux_state), volume);
        summary.new_mass += CellMass(new_state, volume);

        summary.old_mom_r += CellRadialMomentum(old_state, volume);
        summary.flux_mom_r_delta += CellRadialMomentum(flux_delta, volume);
        summary.source_mom_r_delta +=
            CellRadialMomentum(Subtract(new_state, flux_state), volume);
        summary.new_mom_r += CellRadialMomentum(new_state, volume);

        summary.old_e_fluid_total += CellFluidEnergy(old_state, volume);
        summary.flux_e_fluid_total_delta += CellFluidEnergy(flux_delta, volume);
        summary.source_e_fluid_total_delta +=
            CellFluidEnergy(Subtract(new_state, flux_state), volume);
        summary.new_e_fluid_total += CellFluidEnergy(new_state, volume);
      }
    }
  }

  summary.mass_residual =
      summary.old_mass + summary.flux_mass_delta + summary.source_mass_delta - summary.new_mass;
  summary.mom_r_residual =
      summary.old_mom_r +
      summary.flux_mom_r_delta +
      summary.source_mom_r_delta - summary.new_mom_r;
  summary.e_fluid_total_residual =
      summary.old_e_fluid_total +
      summary.flux_e_fluid_total_delta +
      summary.source_e_fluid_total_delta -
      summary.new_e_fluid_total;

  std::ostringstream report;
  report << std::setprecision(17);
  report << "budget_mass_old=" << summary.old_mass
         << "; budget_mass_flux_delta=" << summary.flux_mass_delta
         << "; budget_mass_source_delta=" << summary.source_mass_delta
         << "; budget_mass_new=" << summary.new_mass
         << "; budget_mass_residual=" << summary.mass_residual
         << "; budget_mom_r_old=" << summary.old_mom_r
         << "; budget_mom_r_flux_delta=" << summary.flux_mom_r_delta
         << "; budget_mom_r_source_delta=" << summary.source_mom_r_delta
         << "; budget_mom_r_new=" << summary.new_mom_r
         << "; budget_mom_r_residual=" << summary.mom_r_residual
         << "; budget_e_old=" << summary.old_e_fluid_total
         << "; budget_e_flux_delta=" << summary.flux_e_fluid_total_delta
         << "; budget_e_source_delta=" << summary.source_e_fluid_total_delta
         << "; budget_e_new=" << summary.new_e_fluid_total
         << "; budget_e_residual=" << summary.e_fluid_total_residual;
  summary.report_line = report.str();

  return summary;
}

[[nodiscard]] HydroBudgetResidualSummary BuildAleBudgetResidualSummary(
    const HydroStateSnapshot& old_snapshot,
    const std::vector<HydroConservativeState>& flux_delta_extensive,
    const std::vector<HydroConservativeState>& source_delta_extensive,
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) noexcept {
  HydroBudgetResidualSummary summary;
  if (!old_snapshot.is_complete() ||
      flux_delta_extensive.size() != old_snapshot.cells.size() ||
      source_delta_extensive.size() != old_snapshot.cells.size() ||
      !hydro_view.is_complete() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return summary;
  }

  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double old_volume = CellVolume(old_geometry, old_snapshot.shape, radial, theta, phi);
        const double new_volume = CellVolume(new_geometry, old_snapshot.shape, radial, theta, phi);
        const auto old_extensive = ToExtensiveState(old_snapshot.cells[index], old_volume);
        const auto new_state = LoadCellState(hydro_view, radial, theta, phi);
        const auto new_extensive = ToExtensiveState(new_state, new_volume);

        summary.old_mass += ExtensiveMass(old_extensive);
        summary.flux_mass_delta += ExtensiveMass(flux_delta_extensive[index]);
        summary.source_mass_delta += ExtensiveMass(source_delta_extensive[index]);
        summary.new_mass += ExtensiveMass(new_extensive);

        summary.old_mom_r += ExtensiveRadialMomentum(old_extensive);
        summary.flux_mom_r_delta += ExtensiveRadialMomentum(flux_delta_extensive[index]);
        summary.source_mom_r_delta += ExtensiveRadialMomentum(source_delta_extensive[index]);
        summary.new_mom_r += ExtensiveRadialMomentum(new_extensive);

        summary.old_e_fluid_total += ExtensiveFluidEnergy(old_extensive);
        summary.flux_e_fluid_total_delta += ExtensiveFluidEnergy(flux_delta_extensive[index]);
        summary.source_e_fluid_total_delta += ExtensiveFluidEnergy(source_delta_extensive[index]);
        summary.new_e_fluid_total += ExtensiveFluidEnergy(new_extensive);
      }
    }
  }

  summary.mass_residual =
      summary.old_mass + summary.flux_mass_delta + summary.source_mass_delta - summary.new_mass;
  summary.mom_r_residual =
      summary.old_mom_r +
      summary.flux_mom_r_delta +
      summary.source_mom_r_delta - summary.new_mom_r;
  summary.e_fluid_total_residual =
      summary.old_e_fluid_total +
      summary.flux_e_fluid_total_delta +
      summary.source_e_fluid_total_delta -
      summary.new_e_fluid_total;

  std::ostringstream report;
  report << std::setprecision(17);
  report << "budget_mass_old=" << summary.old_mass
         << "; budget_mass_flux_delta=" << summary.flux_mass_delta
         << "; budget_mass_source_delta=" << summary.source_mass_delta
         << "; budget_mass_new=" << summary.new_mass
         << "; budget_mass_residual=" << summary.mass_residual
         << "; budget_mom_r_old=" << summary.old_mom_r
         << "; budget_mom_r_flux_delta=" << summary.flux_mom_r_delta
         << "; budget_mom_r_source_delta=" << summary.source_mom_r_delta
         << "; budget_mom_r_new=" << summary.new_mom_r
         << "; budget_mom_r_residual=" << summary.mom_r_residual
         << "; budget_e_old=" << summary.old_e_fluid_total
         << "; budget_e_flux_delta=" << summary.flux_e_fluid_total_delta
         << "; budget_e_source_delta=" << summary.source_e_fluid_total_delta
         << "; budget_e_new=" << summary.new_e_fluid_total
         << "; budget_e_residual=" << summary.e_fluid_total_residual;
  summary.report_line = report.str();

  return summary;
}

[[nodiscard]] std::vector<HydroConservativeState> BuildOldMeshExtensiveDelta(
    const HydroStateSnapshot& old_snapshot,
    const dec3d::core::Array3D<HydroConservativeState>& staged_cells,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  std::vector<HydroConservativeState> delta;
  if (!old_snapshot.is_complete() ||
      staged_cells.extent_r() != old_snapshot.shape.radial_cells ||
      staged_cells.extent_theta() != old_snapshot.shape.theta_cells ||
      staged_cells.extent_phi() != old_snapshot.shape.phi_cells ||
      !geometry.is_valid()) {
    return delta;
  }

  delta.assign(old_snapshot.cells.size(), HydroConservativeState{});
  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double volume = CellVolume(geometry, old_snapshot.shape, radial, theta, phi);
        const auto old_extensive = ToExtensiveState(old_snapshot.cells[index], volume);
        const auto staged_extensive =
            ToExtensiveState(staged_cells(radial, theta, phi), volume);
        delta[index] = Subtract(staged_extensive, old_extensive);
      }
    }
  }
  return delta;
}

[[nodiscard]] std::vector<HydroConservativeState> BuildCommittedAleExtensiveDelta(
    const HydroStateSnapshot& old_snapshot,
    const dec3d::core::Array3D<HydroConservativeState>& staged_cells,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) noexcept {
  std::vector<HydroConservativeState> delta;
  if (!old_snapshot.is_complete() ||
      staged_cells.extent_r() != old_snapshot.shape.radial_cells ||
      staged_cells.extent_theta() != old_snapshot.shape.theta_cells ||
      staged_cells.extent_phi() != old_snapshot.shape.phi_cells ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return delta;
  }

  delta.assign(old_snapshot.cells.size(), HydroConservativeState{});
  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double old_volume =
            CellVolume(old_geometry, old_snapshot.shape, radial, theta, phi);
        const double new_volume =
            CellVolume(new_geometry, old_snapshot.shape, radial, theta, phi);
        const auto old_extensive =
            ToExtensiveState(old_snapshot.cells[index], old_volume);
        const auto staged_extensive =
            ToExtensiveState(staged_cells(radial, theta, phi), new_volume);
        delta[index] = Subtract(staged_extensive, old_extensive);
      }
    }
  }
  return delta;
}

[[nodiscard]] std::vector<HydroConservativeState> BuildAleExtensiveDeltaFromAverageDelta(
    const HydroStateSnapshot& old_snapshot,
    const std::vector<HydroConservativeState>& average_delta,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) noexcept {
  std::vector<HydroConservativeState> delta;
  if (!old_snapshot.is_complete() ||
      average_delta.size() != old_snapshot.cells.size() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return delta;
  }

  delta.assign(old_snapshot.cells.size(), HydroConservativeState{});
  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double old_volume =
            CellVolume(old_geometry, old_snapshot.shape, radial, theta, phi);
        const double new_volume =
            CellVolume(new_geometry, old_snapshot.shape, radial, theta, phi);
        const auto old_extensive =
            ToExtensiveState(old_snapshot.cells[index], old_volume);
        const auto projected_new_state =
            Add(old_snapshot.cells[index], average_delta[index]);
        const auto projected_new_extensive =
            ToExtensiveState(projected_new_state, new_volume);
        delta[index] = Subtract(projected_new_extensive, old_extensive);
      }
    }
  }
  return delta;
}

[[nodiscard]] HydroConservativeState BuildAleFullTransportBudgetResidual(
    const HydroStateSnapshot& old_snapshot,
    const std::vector<HydroConservativeState>& extensive_delta,
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry) noexcept {
  HydroConservativeState old_total;
  HydroConservativeState delta_total;
  HydroConservativeState new_total;
  if (!old_snapshot.is_complete() ||
      extensive_delta.size() != old_snapshot.cells.size() ||
      !hydro_view.is_complete() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }

  for (std::size_t radial = 0; radial < old_snapshot.shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < old_snapshot.shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < old_snapshot.shape.phi_cells; ++phi) {
        const auto index = LinearIndex(old_snapshot.shape, radial, theta, phi);
        const double old_volume =
            CellVolume(old_geometry, old_snapshot.shape, radial, theta, phi);
        const double new_volume =
            CellVolume(new_geometry, old_snapshot.shape, radial, theta, phi);
        old_total = Add(
            old_total,
            ToExtensiveState(old_snapshot.cells[index], old_volume));
        delta_total = Add(delta_total, extensive_delta[index]);
        new_total = Add(
            new_total,
            ToExtensiveState(LoadCellState(hydro_view, radial, theta, phi), new_volume));
      }
    }
  }

  return Subtract(Add(old_total, delta_total), new_total);
}

[[nodiscard]] bool AccumulateAleSourceExtensiveDelta(
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    std::vector<HydroConservativeState>& source_delta_extensive,
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string& failure_reason) {
  const auto source_terms = ComputeGeometricSourceTerms(source_snapshot, geometry);
  if (!source_terms.success) {
    failure_reason = source_terms.failure_reason.empty()
                         ? "geometric source term evaluation failed"
                         : source_terms.failure_reason;
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.source.geometric_step.failed",
        failure_reason);
    return false;
  }
  if (!source_terms.is_complete(source_snapshot.states.size()) ||
      source_delta_extensive.size() != source_snapshot.states.size()) {
    failure_reason = "geometric source terms are incomplete for ALE extensive update";
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.source.geometric_step.failed",
        failure_reason);
    return false;
  }

  const GridShape shape{
      source_snapshot.radial_cells,
      source_snapshot.theta_cells,
      source_snapshot.phi_cells};
  for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
        const auto index = LinearIndex(shape, radial, theta, phi);
        const double volume = CellVolume(geometry, shape, radial, theta, phi);
        source_delta_extensive[index] = Add(
            source_delta_extensive[index],
            Scale(source_terms.source_terms[index], dt_s * volume));
      }
    }
  }

  AppendDiagnostic(
      diagnostics,
      "p1.hydro.source.geometric_step",
      "hydro geometric source step updated momentum through extensive old-time source accumulation without modifying rho, E_fluid_total, or chi_e");
  return true;
}

}  // namespace

bool HydroBudgetResidualSummary::is_complete() const noexcept {
  return !report_line.empty();
}

bool ParseHydroBudgetResidualSummary(
    std::string_view report_line,
    HydroBudgetResidualSummary* summary) noexcept {
  if (summary == nullptr || report_line.empty()) {
    return false;
  }

  HydroBudgetResidualSummary parsed{};
  if (!ParseDoubleField(report_line, "budget_mass_old", &parsed.old_mass) ||
      !ParseDoubleField(report_line, "budget_mass_flux_delta", &parsed.flux_mass_delta) ||
      !ParseDoubleField(report_line, "budget_mass_source_delta", &parsed.source_mass_delta) ||
      !ParseDoubleField(report_line, "budget_mass_new", &parsed.new_mass) ||
      !ParseDoubleField(report_line, "budget_mass_residual", &parsed.mass_residual) ||
      !ParseDoubleField(report_line, "budget_mom_r_old", &parsed.old_mom_r) ||
      !ParseDoubleField(report_line, "budget_mom_r_flux_delta", &parsed.flux_mom_r_delta) ||
      !ParseDoubleField(report_line, "budget_mom_r_source_delta", &parsed.source_mom_r_delta) ||
      !ParseDoubleField(report_line, "budget_mom_r_new", &parsed.new_mom_r) ||
      !ParseDoubleField(report_line, "budget_mom_r_residual", &parsed.mom_r_residual) ||
      !ParseDoubleField(report_line, "budget_e_old", &parsed.old_e_fluid_total) ||
      !ParseDoubleField(report_line, "budget_e_flux_delta", &parsed.flux_e_fluid_total_delta) ||
      !ParseDoubleField(report_line, "budget_e_source_delta", &parsed.source_e_fluid_total_delta) ||
      !ParseDoubleField(report_line, "budget_e_new", &parsed.new_e_fluid_total) ||
      !ParseDoubleField(report_line, "budget_e_residual", &parsed.e_fluid_total_residual)) {
    return false;
  }

  parsed.report_line = std::string(report_line);
  *summary = std::move(parsed);
  return summary->is_complete();
}

bool TryExtractHydroBudgetResidualSummary(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    HydroBudgetResidualSummary* summary) noexcept {
  if (summary == nullptr) {
    return false;
  }

  for (const auto& entry : diagnostics.entries) {
    if ((entry.code == "p1.hydro.budget.mass" ||
         entry.code == "p1.hydro.budget.mom_r" ||
         entry.code == "p1.hydro.budget.e_fluid_total") &&
        ParseHydroBudgetResidualSummary(entry.message, summary)) {
      return true;
    }
  }

  return false;
}

bool StaticGridHydroResult::is_complete() const noexcept {
  return success && diagnostics.has_entries() && !report_line.empty() && budget.is_complete();
}

StaticGridHydroResult AdvanceStaticGridHydro(
    dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const StaticGridHydroOptions& options,
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal) noexcept {
  const RadialGhostOverride empty_override{};
  return AdvanceStaticGridHydro(
      hydro_view,
      geometry,
      dt_s,
      empty_override,
      options,
      radial_ale_proposal);
}

StaticGridHydroResult AdvanceStaticGridHydro(
    dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s,
    const RadialGhostOverride& radial_ghost_override,
    const StaticGridHydroOptions& options,
    const dec3d::core::MeshUpdateProposal* radial_ale_proposal) noexcept {
  StaticGridHydroResult result;

  if (!hydro_view.is_complete()) {
    result.failure_reason = "hydro view is incomplete for static-grid hydro";
  } else if (!geometry.is_valid()) {
    result.failure_reason = "spherical geometry metadata is invalid for static-grid hydro";
  } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
    result.failure_reason = "time step must be finite and positive for static-grid hydro";
  } else if (!(options.apply_radial_sweep || options.apply_theta_sweep || options.apply_phi_sweep)) {
    result.failure_reason = "static-grid hydro requires at least one directional sweep";
  } else if (options.enable_radiation_hydro_terms && hydro_view.radiation_chi.empty()) {
    result.failure_reason =
        "radiation hydro terms require a radiation pressure scalar bundle in the hydro view";
  } else if (options.enable_alpha_hydro_terms &&
             (!hydro_view.operator_local_alpha_chi || hydro_view.alpha_chi.empty())) {
    result.failure_reason =
        "alpha hydro terms require alpha_chi in the hydro scalar bundle";
  } else {
    auto timing_start = std::chrono::steady_clock::now();
    const auto hydrodynamic_snapshot = CaptureSnapshot(hydro_view);
    result.timing_snapshot_wall_s += ElapsedHydroSecondsSince(timing_start);
    const bool macro_ale_requested =
        options.use_macro_zoning && options.apply_radial_ale_flux_correction;
    const bool macro_ale_direct_hllc_mode =
        macro_ale_requested && options.use_macro_ale_direct_moving_face_hllc;
    const bool macro_ale_compatibility_mode =
        macro_ale_requested && !macro_ale_direct_hllc_mode;
    const bool defer_macro_ale_writeback =
        macro_ale_compatibility_mode &&
        options.defer_macro_ale_compatibility_writeback;
    const bool use_ale_extensive_update =
        options.apply_radial_ale_flux_correction &&
        !macro_ale_compatibility_mode &&
        !macro_ale_direct_hllc_mode;
    const bool macro_handles_radial_sweep =
        options.use_macro_zoning &&
        options.use_ppm_reconstruction &&
        options.apply_radial_sweep &&
        options.macro_zoning_use_radial_ppm;
    bool macro_ale_whole_domain_valid = false;
    timing_start = std::chrono::steady_clock::now();
    std::vector<HydroConservativeState> accumulated_delta(
        hydrodynamic_snapshot.cells.size(),
        HydroConservativeState{});
    std::vector<HydroConservativeState> radial_delta_per_solid_angle(
        hydrodynamic_snapshot.cells.size(),
        HydroConservativeState{});
    std::vector<HydroConservativeState> radial_delta_extensive(
        hydrodynamic_snapshot.cells.size(),
        HydroConservativeState{});
    std::vector<HydroConservativeState> source_delta(
        hydrodynamic_snapshot.cells.size(),
        HydroConservativeState{});
    dec3d::mesh::SphericalGeometryMetadata ale_updated_geometry = geometry;
    std::vector<HydroConservativeState> macro_ale_extensive_delta(
        hydrodynamic_snapshot.cells.size(),
        HydroConservativeState{});
    dec3d::state::HydroStateView* update_view = &hydro_view;
    dec3d::state::HydroStateView staged_work_view;
    result.timing_scratch_wall_s += ElapsedHydroSecondsSince(timing_start);
    timing_start = std::chrono::steady_clock::now();
    const auto geometric_source_snapshot =
        options.apply_geometric_source ? CaptureGeometricSourceSnapshot(hydro_view)
                                       : GeometricSourceSnapshot{};
    result.timing_snapshot_wall_s += ElapsedHydroSecondsSince(timing_start);

    if (defer_macro_ale_writeback) {
      staged_work_view = BuildOwnedHydroWorkViewFromSnapshot(hydrodynamic_snapshot);
      if (!staged_work_view.is_complete()) {
        result.failure_reason = "macro ALE staged working view is incomplete";
      } else {
        update_view = &staged_work_view;
      }
    }

    if (result.failure_reason.empty() && options.apply_radial_ale_flux_correction) {
      if (radial_ale_proposal == nullptr) {
        if (macro_ale_direct_hllc_mode) {
          result.failure_reason =
              "direct moving-face HLLC macro ALE requires a mesh update proposal";
          result.macro_ale_hllc_diagnostics.failure_class =
              MacroAleHllcFailureClass::proposal_window_missing;
        } else {
          result.failure_reason = "ALE update requires a mesh update proposal";
        }
      } else if (macro_ale_direct_hllc_mode &&
                 !radial_ale_proposal->radial_face_indexing_is_global) {
        result.failure_reason =
            "direct moving-face HLLC macro ALE requires global radial face indexing";
        result.macro_ale_hllc_diagnostics.failure_class =
            MacroAleHllcFailureClass::proposal_not_global_indexed;
      } else if (macro_ale_direct_hllc_mode &&
                 !radial_ale_proposal->has_radial_face_window(
                     options.radial_ale_global_face_begin_index,
                     geometry.radial_faces.size())) {
        result.failure_reason =
            "direct moving-face HLLC macro ALE requires a complete local global-face proposal window";
        result.macro_ale_hllc_diagnostics.failure_class =
            MacroAleHllcFailureClass::proposal_window_missing;
      } else if (macro_ale_compatibility_mode &&
                 !radial_ale_proposal->radial_face_indexing_is_global) {
        result.failure_reason =
            "macro ALE compatibility mode requires global radial face indexing";
      } else if (macro_ale_compatibility_mode &&
                 !radial_ale_proposal->has_radial_face_window(
                     options.radial_ale_global_face_begin_index,
                     geometry.radial_faces.size())) {
        result.failure_reason =
            "macro ALE compatibility mode requires a complete local global-face proposal window";
      } else {
        const auto ale_commit_preview =
            dec3d::mesh::ApplyRadialAleMeshUpdateProposal(
                *radial_ale_proposal,
                ale_updated_geometry,
                options.radial_ale_global_face_begin_index);
        if (!ale_commit_preview.success) {
          result.failure_reason =
              ale_commit_preview.failure_reason.empty()
                  ? "ALE extensive update could not construct proposal geometry"
                  : ale_commit_preview.failure_reason;
          if (macro_ale_direct_hllc_mode) {
            result.macro_ale_hllc_diagnostics.failure_class =
                MacroAleHllcFailureClass::geometry_preview_failed;
          }
        } else if (macro_ale_compatibility_mode) {
          macro_ale_whole_domain_valid = true;
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.macro_ale.geometry_preview",
              ale_commit_preview.report_line +
                  "; radial_face_indexing=global; global_face_begin=" +
                  std::to_string(options.radial_ale_global_face_begin_index));
        }
      }
    }

    if (result.failure_reason.empty() && options.apply_radial_sweep && !macro_handles_radial_sweep) {
      std::string ghost_report_line;
      DirectionalPpmSummary ppm_summary;
      RadialAleFluxSummary radial_ale_summary;
      const auto sweep_start = std::chrono::steady_clock::now();
      const bool sweep_ok = AccumulateDirectionalSweepDelta(
              SweepDirection::radial,
              hydrodynamic_snapshot,
              geometry,
              dt_s,
              options,
              use_ale_extensive_update,
              accumulated_delta,
              result.failure_reason,
              radial_ghost_override.ghost_layers == 0u ? nullptr : &radial_ghost_override,
              &ghost_report_line,
              &ppm_summary,
              radial_ale_proposal,
              &radial_ale_summary,
              use_ale_extensive_update ? &radial_delta_per_solid_angle : nullptr,
              use_ale_extensive_update ? &radial_delta_extensive : nullptr);
      result.timing_radial_sweep_wall_s += ElapsedHydroSecondsSince(sweep_start);
      if (!sweep_ok) {
        if (result.failure_reason.empty()) {
          result.failure_reason = "radial static-grid hydro sweep failed";
        }
      } else {
        result.radial_executed = true;
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.direction.radial",
            "executed radial directional hydro update");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.radial_boundary_contract.executed",
            "radial hydro path executed ghost-state boundary contract instead of boundary-cell physical-flux substitution");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.ghost.direction.radial",
            ghost_report_line.empty() ? "prepared radial directional hydro ghosts"
                                      : ghost_report_line);
        if (radial_ale_summary.executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.ale.radial_flux_correction.executed",
              radial_ale_summary.report_line.empty()
                  ? "executed thesis-style radial ALE flux correction"
                  : radial_ale_summary.report_line);
          if (options.debug_angular_stage_diagnostics) {
            for (const auto& debug_report : radial_ale_summary.debug_report_lines) {
              AppendDiagnostic(
                  result.diagnostics,
                  "p1.hydro.debug.radial_ale.face_flux",
                  debug_report);
            }
          }
        }
        if (ppm_summary.executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.executed",
              "executed direction-local characteristic traced-interface PPM reconstruction before radial HLLC flux assembly");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.ng3",
              "PPM reconstruction consumed ng=3 directional ghosts");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.direction.radial",
              ppm_summary.report_line);
        }
        if (options.debug_angular_stage_diagnostics) {
          AppendFineAngularStageDiagnostic(
              result.diagnostics,
              "post_fine_radial_delta",
              hydrodynamic_snapshot,
              accumulated_delta,
              geometry);
        }
      }
    }

    if (result.failure_reason.empty() &&
        options.use_macro_zoning &&
        (macro_handles_radial_sweep || options.apply_theta_sweep || options.apply_phi_sweep)) {
      MacroZonedHydroSummary macro_summary;
      if (!AccumulateMacroZonedHydroDelta(
              hydrodynamic_snapshot,
              *update_view,
              geometry,
              dt_s,
              options,
              radial_ghost_override.ghost_layers == 0u ? nullptr : &radial_ghost_override,
              accumulated_delta,
              result.failure_reason,
              macro_ale_whole_domain_valid,
              &macro_summary,
              macro_ale_direct_hllc_mode ? radial_ale_proposal : nullptr,
              macro_ale_direct_hllc_mode ? &ale_updated_geometry : nullptr,
              macro_ale_direct_hllc_mode ? &result.macro_ale_hllc_diagnostics : nullptr,
              macro_ale_direct_hllc_mode ? &result.macro_ale_hllc_local_face_window : nullptr)) {
        if (result.failure_reason.empty()) {
          result.failure_reason = "macro-zoned coarse hydro update failed";
        }
      } else {
        result.timing_macro_detect_wall_s += macro_summary.detect_wall_s;
        result.timing_macro_restrict_wall_s += macro_summary.restrict_wall_s;
        result.timing_macro_update_wall_s += macro_summary.update_wall_s;
        result.timing_macro_radial_update_wall_s += macro_summary.radial_update_wall_s;
        result.timing_macro_theta_update_wall_s += macro_summary.theta_update_wall_s;
        result.timing_macro_phi_update_wall_s += macro_summary.phi_update_wall_s;
        result.timing_macro_state_update_wall_s += macro_summary.state_update_wall_s;
        result.timing_macro_prolong_wall_s += macro_summary.prolong_wall_s;
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.macro_zoning.detected",
            macro_summary.detect_report_line.empty()
                ? "detected thesis-style macro-zoning coarse map"
                : macro_summary.detect_report_line);
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.macro_zoning.restrict.executed",
            macro_summary.restrict_report_line.empty()
                ? "restricted fine hydro package onto macro-zoned coarse work view"
                : macro_summary.restrict_report_line);
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.macro_zoning.coarse_update.executed",
            macro_summary.coarse_update_report_line.empty()
                ? "executed single-rank macro-zoned coarse angular hydro update"
                : macro_summary.coarse_update_report_line);
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.macro_zoning.prolong.executed",
            macro_summary.prolong_report_line.empty()
                ? "prolonged macro-zoned coarse work view back onto fine hydro state"
                : macro_summary.prolong_report_line);
        if (options.debug_angular_stage_diagnostics) {
          for (const auto& angular_report : macro_summary.angular_stage_report_lines) {
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.debug.angular.stage",
                angular_report);
          }
        }
        if (macro_summary.coarse_ghost_bootstrap_executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.macro_zoning.coarse_ghost_bootstrap.executed",
              macro_summary.coarse_ghost_bootstrap_report_line.empty()
                  ? "restricted fine radial halo override onto macro-zoned coarse ghost work view"
                  : macro_summary.coarse_ghost_bootstrap_report_line);
        }
        if (macro_summary.coarse_source_budget_executed) {
          result.geometric_source_executed = true;
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.macro_ale_hllc.coarse_source_budget",
              macro_summary.coarse_source_budget_report_line.empty()
                  ? "direct moving-face HLLC macro ALE integrated geometric source in the coarse extensive budget"
                  : macro_summary.coarse_source_budget_report_line);
        }
        if (macro_summary.radial_ppm_summary.executed ||
            macro_summary.theta_ppm_summary.executed ||
            macro_summary.phi_ppm_summary.executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.executed",
              "executed direction-local characteristic traced-interface PPM reconstruction on macro-zoned coarse work lines before HLLC flux assembly");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.ng3",
              "PPM reconstruction consumed ng=3 directional ghosts");
          if (macro_summary.radial_ppm_summary.executed) {
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.reconstruction.ppm.direction.radial",
                macro_summary.radial_ppm_summary.report_line);
          }
          if (macro_summary.theta_ppm_summary.executed) {
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.reconstruction.ppm.direction.theta",
                macro_summary.theta_ppm_summary.report_line);
          }
          if (macro_summary.phi_ppm_summary.executed) {
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.reconstruction.ppm.direction.phi",
                macro_summary.phi_ppm_summary.report_line);
          }
        }

        if (macro_summary.radial_executed) {
          result.radial_executed = true;
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.direction.radial",
              "executed radial directional hydro update through macro-zoned coarse/fine coupling");
        }
        if (options.apply_theta_sweep) {
          result.theta_executed = true;
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.direction.theta",
              "executed theta directional hydro update through macro-zoned coarse/fine coupling");
        }
        if (options.apply_phi_sweep) {
          result.phi_executed = true;
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.direction.phi",
              "executed phi directional hydro update through macro-zoned coarse/fine coupling");
        }
      }
    }

    if (result.failure_reason.empty() && !options.use_macro_zoning && options.apply_theta_sweep) {
      std::string ghost_report_line;
      DirectionalPpmSummary ppm_summary;
      const auto sweep_start = std::chrono::steady_clock::now();
      const bool sweep_ok = AccumulateDirectionalSweepDelta(
              SweepDirection::theta,
              hydrodynamic_snapshot,
              geometry,
              dt_s,
              options,
              use_ale_extensive_update,
              accumulated_delta,
              result.failure_reason,
              nullptr,
              &ghost_report_line,
              &ppm_summary);
      result.timing_theta_sweep_wall_s += ElapsedHydroSecondsSince(sweep_start);
      if (!sweep_ok) {
        if (result.failure_reason.empty()) {
          result.failure_reason = "theta static-grid hydro sweep failed";
        }
      } else {
        result.theta_executed = true;
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.direction.theta",
            "executed theta directional hydro update");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.theta_pole_contract.executed",
            "static-grid hydro executed theta pole remap instead of boundary-cell physical-flux substitution");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.ghost.direction.theta",
            ghost_report_line.empty() ? "prepared theta directional hydro ghosts"
                                      : ghost_report_line);
        if (ppm_summary.executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.executed",
              "executed direction-local characteristic traced-interface PPM reconstruction before theta HLLC flux assembly");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.ng3",
              "PPM reconstruction consumed ng=3 directional ghosts");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.direction.theta",
              ppm_summary.report_line);
        }
        if (options.debug_angular_stage_diagnostics) {
          AppendFineAngularStageDiagnostic(
              result.diagnostics,
              "post_fine_theta_delta",
              hydrodynamic_snapshot,
              accumulated_delta,
              geometry);
        }
      }
    }

    if (result.failure_reason.empty() && !options.use_macro_zoning && options.apply_phi_sweep) {
      std::string ghost_report_line;
      DirectionalPpmSummary ppm_summary;
      const auto sweep_start = std::chrono::steady_clock::now();
      const bool sweep_ok = AccumulateDirectionalSweepDelta(
              SweepDirection::phi,
              hydrodynamic_snapshot,
              geometry,
              dt_s,
              options,
              use_ale_extensive_update,
              accumulated_delta,
              result.failure_reason,
              nullptr,
              &ghost_report_line,
              &ppm_summary);
      result.timing_phi_sweep_wall_s += ElapsedHydroSecondsSince(sweep_start);
      if (!sweep_ok) {
        if (result.failure_reason.empty()) {
          result.failure_reason = "phi static-grid hydro sweep failed";
        }
      } else {
        result.phi_executed = true;
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.direction.phi",
            "executed phi directional hydro update");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.phi_periodic_contract.executed",
            "phi hydro path executed periodic boundary contract through a unified directional boundary entry");
        AppendDiagnostic(
            result.diagnostics,
            "p1.hydro.ghost.direction.phi",
            ghost_report_line.empty() ? "prepared phi directional hydro ghosts"
                                      : ghost_report_line);
        if (ppm_summary.executed) {
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.executed",
              "executed direction-local characteristic traced-interface PPM reconstruction before phi HLLC flux assembly");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.ng3",
              "PPM reconstruction consumed ng=3 directional ghosts");
          AppendDiagnostic(
              result.diagnostics,
              "p1.hydro.reconstruction.ppm.direction.phi",
              ppm_summary.report_line);
        }
        if (options.debug_angular_stage_diagnostics) {
          AppendFineAngularStageDiagnostic(
              result.diagnostics,
              "post_fine_phi_delta",
              hydrodynamic_snapshot,
              accumulated_delta,
              geometry);
        }
      }
    }

    if (result.failure_reason.empty()) {
      if (use_ale_extensive_update) {
        if (options.apply_geometric_source) {
          const auto source_start = std::chrono::steady_clock::now();
          const bool source_ok = AccumulateAleSourceExtensiveDelta(
              geometric_source_snapshot,
              geometry,
              dt_s,
              source_delta,
              result.diagnostics,
              result.failure_reason);
          result.timing_source_wall_s += ElapsedHydroSecondsSince(source_start);
          if (!source_ok) {
            if (result.failure_reason.empty()) {
              result.failure_reason = "ALE geometric source accumulation failed";
            }
          } else {
            result.geometric_source_executed = true;
          }
        }

        if (result.failure_reason.empty()) {
          const auto commit_start = std::chrono::steady_clock::now();
          for (std::size_t radial = 0; radial < hydrodynamic_snapshot.shape.radial_cells; ++radial) {
            for (std::size_t theta = 0; theta < hydrodynamic_snapshot.shape.theta_cells; ++theta) {
              for (std::size_t phi = 0; phi < hydrodynamic_snapshot.shape.phi_cells; ++phi) {
                const auto index = LinearIndex(hydrodynamic_snapshot.shape, radial, theta, phi);
                const auto updated_state = BuildAleFactorizedUpdatedState(
                    hydrodynamic_snapshot.cells[index],
                    accumulated_delta[index],
                    source_delta[index],
                    radial_delta_per_solid_angle[index],
                    radial_delta_extensive[index],
                    geometry,
                    ale_updated_geometry,
                    radial,
                    theta,
                    phi);
                if (!IsPhysicalCellState(updated_state)) {
                  result.failure_reason =
                      "ALE extensive hydro update produced a non-physical state";
                  break;
                }
                StoreCellState(updated_state, *update_view, radial, theta, phi);
              }
              if (!result.failure_reason.empty()) {
                break;
              }
            }
            if (!result.failure_reason.empty()) {
              break;
            }
          }
          result.timing_commit_wall_s += ElapsedHydroSecondsSince(commit_start);
        }
        if (result.failure_reason.empty() && options.debug_angular_stage_diagnostics) {
          AppendAleCommitDensityBalanceDiagnostic(
              result.diagnostics,
              hydrodynamic_snapshot,
              accumulated_delta,
              source_delta,
              radial_delta_per_solid_angle,
              radial_delta_extensive,
              geometry,
              ale_updated_geometry);
          AppendFineAngularStageDiagnosticFromView(
              result.diagnostics,
              "post_ale_extensive_commit",
              *update_view,
              hydrodynamic_snapshot.shape,
              ale_updated_geometry);
        }
      } else {
        const auto commit_start = std::chrono::steady_clock::now();
        for (std::size_t radial = 0; radial < hydrodynamic_snapshot.shape.radial_cells; ++radial) {
          for (std::size_t theta = 0; theta < hydrodynamic_snapshot.shape.theta_cells; ++theta) {
            for (std::size_t phi = 0; phi < hydrodynamic_snapshot.shape.phi_cells; ++phi) {
              const auto index = LinearIndex(hydrodynamic_snapshot.shape, radial, theta, phi);
              const auto updated_state =
                  Add(hydrodynamic_snapshot.cells[index], accumulated_delta[index]);
              if (!IsPhysicalCellState(updated_state)) {
                result.failure_reason = "combined static-grid hydro update produced a non-physical state";
                break;
              }
              StoreCellState(updated_state, *update_view, radial, theta, phi);
            }
            if (!result.failure_reason.empty()) {
              break;
            }
          }
          if (!result.failure_reason.empty()) {
            break;
          }
        }
        result.timing_commit_wall_s += ElapsedHydroSecondsSince(commit_start);
        if (result.failure_reason.empty() && options.debug_angular_stage_diagnostics) {
          AppendFineAngularStageDiagnosticFromView(
              result.diagnostics,
              "post_flux_commit",
              *update_view,
              hydrodynamic_snapshot.shape,
              geometry);
        }

        if (result.failure_reason.empty() && options.apply_geometric_source &&
            !macro_ale_direct_hllc_mode) {
          const bool use_macro_zoned_geometric_source =
              options.use_macro_zoning &&
              (options.apply_theta_sweep || options.apply_phi_sweep);
          const auto source_start = std::chrono::steady_clock::now();
          const auto source_step =
              use_macro_zoned_geometric_source
                  ? ApplyMacroZonedGeometricSourceStep(
                        *update_view,
                        geometric_source_snapshot,
                        geometry,
                        options.macro_zoning_coarse_factor,
                        dt_s)
                  : ApplyGeometricSourceStep(
                        *update_view,
                        geometric_source_snapshot,
                        geometry,
                        dt_s);
          result.timing_source_wall_s += ElapsedHydroSecondsSince(source_start);
          result.diagnostics.entries.insert(
              result.diagnostics.entries.end(),
              source_step.diagnostics.entries.begin(),
              source_step.diagnostics.entries.end());
          if (!source_step.success) {
            result.failure_reason = source_step.failure_reason;
          } else {
            result.geometric_source_executed = true;
            if (options.debug_angular_stage_diagnostics) {
              AppendFineAngularStageDiagnosticFromView(
                  result.diagnostics,
                  "post_geometric_source",
                  *update_view,
                  hydrodynamic_snapshot.shape,
                  geometry);
            }
          }
        }
        if (result.failure_reason.empty() &&
            options.apply_geometric_source &&
            macro_ale_direct_hllc_mode &&
            !result.geometric_source_executed) {
          result.failure_reason =
              "direct moving-face HLLC macro ALE did not integrate geometric source into the coarse extensive budget";
          result.macro_ale_hllc_diagnostics.failure_class =
              MacroAleHllcFailureClass::required_diagnostic_missing;
        }
      }
    }

    if (result.failure_reason.empty() && macro_ale_compatibility_mode) {
      const auto staged_old_mesh_cells = CaptureCellArray(*update_view);
      macro_ale_extensive_delta =
          BuildOldMeshExtensiveDelta(hydrodynamic_snapshot, staged_old_mesh_cells, geometry);
      if (macro_ale_extensive_delta.size() != hydrodynamic_snapshot.cells.size()) {
        result.failure_reason =
            "macro ALE compatibility mode could not stage old-mesh extensive delta";
      } else {
        const auto local_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
            *radial_ale_proposal,
            options.radial_ale_global_face_begin_index,
            geometry.radial_faces.size());
        if (!local_window.is_complete()) {
          result.failure_reason =
              local_window.failure_reason.empty()
                  ? "macro ALE compatibility mode could not build local proposal window"
                  : local_window.failure_reason;
        } else {
          const auto remap = RemapHydroStateRadiallyConservative(
              staged_old_mesh_cells,
              geometry,
              ale_updated_geometry,
              local_window.proposed_radial_faces);
          if (!remap.success || !remap.is_complete()) {
            result.failure_reason =
                remap.failure_reason.empty()
                    ? "macro ALE compatibility radial remap failed"
                    : remap.failure_reason;
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.macro_ale.remap.failed",
                remap.diagnostics.report_line.empty()
                    ? result.failure_reason
                    : remap.diagnostics.report_line);
          } else {
            result.macro_ale_staged_cells = remap.remapped_cells;
            result.macro_ale_staged_geometry = ale_updated_geometry;
            result.macro_ale_remap_diagnostics = remap.diagnostics;
            result.macro_ale_staged_writeback_available = true;
            result.macro_ale_staged_geometry_available = true;
            StoreCellArray(remap.remapped_cells, *update_view);
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.macro_ale.compatibility.executed",
                "mode=post_hydro_radial_remap; radial_face_indexing=global; global_face_begin=" +
                    std::to_string(options.radial_ale_global_face_begin_index) +
                    "; remap_order=first_order_proposal_mapped_overlap");
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.macro_ale.remap.executed",
                remap.diagnostics.report_line);
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.macro_ale.remap.conservation",
                remap.diagnostics.report_line);
            if (defer_macro_ale_writeback) {
              AppendDiagnostic(
                  result.diagnostics,
                  "p1.hydro.macro_ale.transaction.staged",
                  "macro ALE compatibility staged hydro and previewed geometry without canonical writeback; published_from_staged_geometry=false");
            } else {
              AppendDiagnostic(
                  result.diagnostics,
                  "p1.hydro.macro_ale.transaction.clean",
                  "macro ALE compatibility remap published only after proposal preview, macro hydro, geometric source, and remap validation succeeded");
            }
            if (options.debug_angular_stage_diagnostics) {
              AppendFineAngularStageDiagnosticFromView(
                  result.diagnostics,
                  "post_macro_ale_remap",
                  *update_view,
                  hydrodynamic_snapshot.shape,
                  ale_updated_geometry);
            }
          }
        }
      }
    }

    if (result.failure_reason.empty()) {
      if (macro_ale_direct_hllc_mode) {
        macro_ale_extensive_delta = BuildAleExtensiveDeltaFromAverageDelta(
            hydrodynamic_snapshot,
            accumulated_delta,
            geometry,
            ale_updated_geometry);
        if (macro_ale_extensive_delta.size() != hydrodynamic_snapshot.cells.size()) {
          result.failure_reason =
              "direct moving-face HLLC macro ALE could not build projected extensive delta from coarse budget";
          result.macro_ale_hllc_diagnostics.failure_class =
              MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
        }
      }
    }

    if (result.failure_reason.empty()) {
      const auto budget_start = std::chrono::steady_clock::now();
      if (use_ale_extensive_update) {
        result.budget = BuildAleBudgetResidualSummary(
            hydrodynamic_snapshot,
            accumulated_delta,
            source_delta,
            *update_view,
            geometry,
            ale_updated_geometry);
      } else if (macro_ale_direct_hllc_mode) {
        result.budget = BuildAleBudgetResidualSummary(
            hydrodynamic_snapshot,
            macro_ale_extensive_delta,
            source_delta,
            *update_view,
            geometry,
            ale_updated_geometry);
      } else if (macro_ale_compatibility_mode) {
        result.budget = BuildAleBudgetResidualSummary(
            hydrodynamic_snapshot,
            macro_ale_extensive_delta,
            source_delta,
            *update_view,
            geometry,
            ale_updated_geometry);
      } else {
        result.budget = BuildBudgetResidualSummary(
            hydrodynamic_snapshot,
            accumulated_delta,
            *update_view,
            geometry);
      }
      if (!result.budget.is_complete()) {
        result.failure_reason = "static-grid hydro budget residual summary is incomplete";
      }
      if (result.failure_reason.empty() && macro_ale_direct_hllc_mode) {
        const auto a4_full_residual = BuildAleFullTransportBudgetResidual(
            hydrodynamic_snapshot,
            macro_ale_extensive_delta,
            *update_view,
            geometry,
            ale_updated_geometry);
        if (!a4_full_residual.is_finite()) {
          result.failure_reason =
              "direct moving-face HLLC macro ALE full transport budget residual is incomplete";
          result.macro_ale_hllc_diagnostics.failure_class =
              MacroAleHllcFailureClass::extensive_budget_residual_exceeded;
        } else {
          result.macro_ale_hllc_staged_cells = CaptureCellArray(*update_view);
          result.macro_ale_hllc_staged_geometry = ale_updated_geometry;
          result.macro_ale_hllc_staged_writeback_available = true;
          result.macro_ale_hllc_staged_geometry_available = true;
          result.macro_ale_hllc_diagnostics.executed = true;
          result.macro_ale_hllc_diagnostics.global_face_begin =
              options.radial_ale_global_face_begin_index;
          result.macro_ale_hllc_diagnostics.local_face_count =
              geometry.radial_faces.size();
          result.macro_ale_hllc_diagnostics.radiation_group_count =
              update_view->radiation_chi.size();
          result.macro_ale_hllc_diagnostics.global_mass_residual =
              a4_full_residual.rho;
          result.macro_ale_hllc_diagnostics.global_mom_r_residual =
              a4_full_residual.mom_r;
          result.macro_ale_hllc_diagnostics.global_mom_theta_residual =
              a4_full_residual.mom_theta;
          result.macro_ale_hllc_diagnostics.global_mom_phi_residual =
              a4_full_residual.mom_phi;
          result.macro_ale_hllc_diagnostics.global_e_fluid_total_residual =
              a4_full_residual.e_fluid_total;
          result.macro_ale_hllc_diagnostics.global_chi_e_residual =
              a4_full_residual.chi_e;
          result.macro_ale_hllc_diagnostics.global_alpha_chi_residual =
              a4_full_residual.alpha_chi;
          double radiation_chi_residual = 0.0;
          for (const double residual : a4_full_residual.radiation_chi) {
            radiation_chi_residual =
                std::max(radiation_chi_residual, std::abs(residual));
          }
          result.macro_ale_hllc_diagnostics.global_radiation_chi_residual =
              radiation_chi_residual;
          result.macro_ale_hllc_diagnostics.local_max_mass_residual =
              std::abs(a4_full_residual.rho);
          result.macro_ale_hllc_diagnostics.local_max_mom_r_residual =
              std::abs(a4_full_residual.mom_r);
          result.macro_ale_hllc_diagnostics.local_max_mom_theta_residual =
              std::abs(a4_full_residual.mom_theta);
          result.macro_ale_hllc_diagnostics.local_max_mom_phi_residual =
              std::abs(a4_full_residual.mom_phi);
          result.macro_ale_hllc_diagnostics.local_max_e_fluid_total_residual =
              std::abs(a4_full_residual.e_fluid_total);
          result.macro_ale_hllc_diagnostics.local_max_chi_e_residual =
              std::abs(a4_full_residual.chi_e);
          result.macro_ale_hllc_diagnostics.local_max_alpha_chi_residual =
              std::abs(a4_full_residual.alpha_chi);
          result.macro_ale_hllc_diagnostics.local_max_radiation_chi_residual =
              radiation_chi_residual;
          result.macro_ale_hllc_diagnostics.global_stage_ok = true;
          result.macro_ale_hllc_diagnostics.global_publish_ok = true;
          result.macro_ale_hllc_diagnostics.no_partial_canonical_writeback = true;
          result.macro_ale_hllc_diagnostics.report_line =
              BuildMacroAleHllcReportLine(result.macro_ale_hllc_diagnostics);
          AppendMacroAleHllcDiagnostics(
              result.diagnostics,
              result.macro_ale_hllc_diagnostics);
          result.macro_ale_hllc_executed = true;
        }
      }
      result.timing_budget_wall_s += ElapsedHydroSecondsSince(budget_start);
    }
  }

  const auto diagnostics_start = std::chrono::steady_clock::now();
  result.success = result.failure_reason.empty();
  if (result.success) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.static_grid.executed",
        "static-grid hydro helper executed the requested directional and source path");
    AppendDiagnostic(
        result.diagnostics,
        "p3.radiation.hydro_terms.stage_H",
        std::string("stage_id=H; radiation_hydro_terms_enabled=") +
            (options.enable_radiation_hydro_terms ? "true" : "false") +
            "; radiation_advection=" +
            (options.enable_radiation_hydro_terms ? "enabled" : "disabled") +
            "; radiation_pressure_work=" +
            (options.enable_radiation_hydro_terms ? "enabled" : "disabled") +
            "; advected_radiation_scalar=P_g_power_3_over_4; passive_Ug_advection=false; radiation_group_count=" +
            std::to_string(hydro_view.radiation_chi.size()));
    if (options.enable_alpha_hydro_terms) {
      AppendDiagnostic(
          result.diagnostics,
          "p4.alpha.hydro_terms.stage_H",
          "stage_id=H; alpha_hydro_terms_enabled=true"
          "; advected_alpha_scalar=P_alpha_power_3_over_5"
          "; passive_epsilon_alpha_advection=false"
          "; alpha_hydro_terms_share_hydro_scalar_bundle=true");
    }
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.budget.mass",
        result.budget.report_line);
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.budget.mom_r",
        result.budget.report_line);
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.budget.e_fluid_total",
        result.budget.report_line);
  } else if (!result.diagnostics.has_entries()) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.static_grid.failed",
        result.failure_reason.empty() ? "static-grid hydro failed" : result.failure_reason);
  }
  result.timing_diagnostics_wall_s += ElapsedHydroSecondsSince(diagnostics_start);

  std::ostringstream report;
  report << std::setprecision(17)
         << "static_grid_hydro_success=" << (result.success ? "true" : "false")
         << "; source=" << (result.geometric_source_executed ? "true" : "false")
         << "; radial=" << (result.radial_executed ? "true" : "false")
         << "; theta=" << (result.theta_executed ? "true" : "false")
         << "; phi=" << (result.phi_executed ? "true" : "false")
         << "; radiation_hydro_terms_enabled="
         << (options.enable_radiation_hydro_terms ? "true" : "false")
         << "; advected_radiation_scalar=P_g_power_3_over_4"
         << "; passive_Ug_advection=false"
         << "; alpha_hydro_terms_enabled="
         << (options.enable_alpha_hydro_terms ? "true" : "false")
         << "; h_hydro_snapshot_wall_s=" << result.timing_snapshot_wall_s
         << "; h_hydro_scratch_wall_s=" << result.timing_scratch_wall_s
         << "; h_hydro_radial_sweep_wall_s=" << result.timing_radial_sweep_wall_s
         << "; h_hydro_macro_detect_wall_s=" << result.timing_macro_detect_wall_s
         << "; h_hydro_macro_restrict_wall_s=" << result.timing_macro_restrict_wall_s
         << "; h_hydro_macro_update_wall_s=" << result.timing_macro_update_wall_s
         << "; h_hydro_macro_radial_update_wall_s="
         << result.timing_macro_radial_update_wall_s
         << "; h_hydro_macro_theta_update_wall_s="
         << result.timing_macro_theta_update_wall_s
         << "; h_hydro_macro_phi_update_wall_s="
         << result.timing_macro_phi_update_wall_s
         << "; h_hydro_macro_state_update_wall_s="
         << result.timing_macro_state_update_wall_s
         << "; h_hydro_macro_prolong_wall_s=" << result.timing_macro_prolong_wall_s
         << "; h_hydro_theta_sweep_wall_s=" << result.timing_theta_sweep_wall_s
         << "; h_hydro_phi_sweep_wall_s=" << result.timing_phi_sweep_wall_s
         << "; h_hydro_commit_wall_s=" << result.timing_commit_wall_s
         << "; h_hydro_source_wall_s=" << result.timing_source_wall_s
         << "; h_hydro_budget_wall_s=" << result.timing_budget_wall_s
         << "; h_hydro_diagnostics_wall_s=" << result.timing_diagnostics_wall_s;
  if (options.enable_alpha_hydro_terms) {
    report << "; advected_alpha_scalar=P_alpha_power_3_over_5"
           << "; passive_epsilon_alpha_advection=false"
           << "; alpha_hydro_terms_share_hydro_scalar_bundle=true";
  }
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();
  AppendDiagnostic(
      result.diagnostics,
      "p1.hydro.static_grid.timing",
      result.report_line);

  return result;
}

bool PublishMacroAleStagedHydroState(
    dec3d::state::HydroStateView& hydro_view,
    const StaticGridHydroResult& staged_result,
    std::string* failure_reason) noexcept {
  if (!hydro_view.is_complete()) {
    if (failure_reason != nullptr) {
      *failure_reason =
          "cannot publish staged macro ALE state into incomplete hydro view";
    }
    return false;
  }
  const dec3d::core::Array3D<HydroConservativeState>* staged_cells = nullptr;
  if (staged_result.macro_ale_staged_writeback_available &&
      !staged_result.macro_ale_staged_cells.empty()) {
    staged_cells = &staged_result.macro_ale_staged_cells;
  } else if (staged_result.macro_ale_hllc_staged_writeback_available &&
             !staged_result.macro_ale_hllc_staged_cells.empty()) {
    staged_cells = &staged_result.macro_ale_hllc_staged_cells;
  }
  if (!staged_result.success || staged_cells == nullptr) {
    if (failure_reason != nullptr) {
      *failure_reason = "staged macro ALE result is not publishable";
    }
    return false;
  }
  if (staged_cells->extent_r() != hydro_view.rho->extent_r() ||
      staged_cells->extent_theta() != hydro_view.rho->extent_theta() ||
      staged_cells->extent_phi() != hydro_view.rho->extent_phi()) {
    if (failure_reason != nullptr) {
      *failure_reason = "staged macro ALE result shape does not match hydro view";
    }
    return false;
  }

  StoreCellArray(*staged_cells, hydro_view);
  if (failure_reason != nullptr) {
    failure_reason->clear();
  }
  return true;
}

}  // namespace dec3d::hydro
