#include "hydro/reconstruction/ppm_reconstruction.hpp"

#include "state/hydro_state/hydro_view.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

namespace dec3d::hydro {

namespace {

constexpr std::size_t kComponentCount = 7u;
constexpr double kSmall = 1.0e-12;
constexpr double kSigmaTolerance = 1.0e-8;

using PrimitiveVector = std::array<double, kComponentCount>;
using CharacteristicMatrix = std::array<PrimitiveVector, kComponentCount>;

struct CharacteristicBasis {
  CharacteristicMatrix R{};
  PrimitiveVector lambda{};
  double sound_speed_squared{0.0};
  double density_over_sound_speed{0.0};
  bool success{false};
};

struct PpmCellProfile {
  CharacteristicBasis basis;
  PrimitiveVector c_center{};
  PrimitiveVector c_left{};
  PrimitiveVector c_right{};
  PrimitiveVector c6{};
  DirectionalPrimitiveState w_center{};
  DirectionalPrimitiveState w_left_edge{};
  DirectionalPrimitiveState w_right_edge{};
  bool troubled{false};
  int troubled_component{-1};
  double trouble_q_im1{0.0};
  double trouble_q_i{0.0};
  double trouble_q_ip1{0.0};
  double trouble_total_variation{0.0};
  double trouble_end_to_end{0.0};
  PpmSmoothnessProbe component1_smoothness;
};

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] bool IsFinite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] bool RadiationBundleFinite(const std::vector<double>& values) noexcept {
  for (const double value : values) {
    if (!IsFinite(value) || value < 0.0) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] CharacteristicBasis BuildCharacteristicBasis(
    const DirectionalPrimitiveState& primitive) noexcept {
  CharacteristicBasis basis;
  if (!primitive.is_physical()) {
    return basis;
  }

  const double gamma = dec3d::state::HydroIdealGasGamma();
  const double sound_speed =
      std::sqrt(gamma * primitive.pressure / primitive.rho);
  if (!(sound_speed > 0.0) || !IsFinite(sound_speed)) {
    return basis;
  }

  const double c2 = sound_speed * sound_speed;
  const double c_over_rho = sound_speed / primitive.rho;
  basis.sound_speed_squared = c2;
  basis.density_over_sound_speed = primitive.rho / sound_speed;
  basis.lambda = {
      primitive.v_n - sound_speed,
      primitive.v_n,
      primitive.v_n,
      primitive.v_n,
      primitive.v_n,
      primitive.v_n,
      primitive.v_n + sound_speed};
  basis.R = {{
      {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0},
      {-c_over_rho, 0.0, 0.0, 0.0, 0.0, 0.0, c_over_rho},
      {0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
      {c2, 0.0, 0.0, 0.0, 0.0, 0.0, c2},
      {0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0}}};

  basis.success = true;
  return basis;
}

[[nodiscard]] PrimitiveVector PrimitiveToVector(
    const DirectionalPrimitiveState& primitive) noexcept {
  return {
      primitive.rho,
      primitive.v_n,
      primitive.v_t1,
      primitive.v_t2,
      primitive.pressure,
      primitive.chi_e,
      primitive.alpha_chi};
}

[[nodiscard]] DirectionalPrimitiveState VectorToPrimitive(
    const PrimitiveVector& vector,
    const std::vector<double>& radiation_chi) noexcept {
  return {
      vector[0],
      vector[1],
      vector[2],
      vector[3],
      vector[4],
      vector[5],
      vector[6],
      radiation_chi};
}

[[nodiscard]] PrimitiveVector Multiply(
    const CharacteristicMatrix& matrix,
    const PrimitiveVector& vector) noexcept {
  PrimitiveVector result{};
  for (std::size_t row = 0; row < kComponentCount; ++row) {
    for (std::size_t col = 0; col < kComponentCount; ++col) {
      result[row] += matrix[row][col] * vector[col];
    }
  }
  return result;
}

[[nodiscard]] PrimitiveVector Add(
    const PrimitiveVector& lhs,
    const PrimitiveVector& rhs) noexcept {
  PrimitiveVector result{};
  for (std::size_t component = 0; component < kComponentCount; ++component) {
    result[component] = lhs[component] + rhs[component];
  }
  return result;
}

[[nodiscard]] PrimitiveVector ProjectPrimitiveDelta(
    const DirectionalPrimitiveState& center,
    const DirectionalPrimitiveState& sample,
    const CharacteristicBasis& basis) noexcept {
  PrimitiveVector projected{};
  if (!basis.success) {
    return projected;
  }
  const double delta_rho = sample.rho - center.rho;
  const double delta_vn = sample.v_n - center.v_n;
  const double delta_pressure = sample.pressure - center.pressure;
  const double pressure_characteristic = delta_pressure / basis.sound_speed_squared;
  const double velocity_characteristic = basis.density_over_sound_speed * delta_vn;

  projected[0] = 0.5 * (pressure_characteristic - velocity_characteristic);
  projected[1] = delta_rho - pressure_characteristic;
  projected[2] = sample.v_t1 - center.v_t1;
  projected[3] = sample.v_t2 - center.v_t2;
  projected[4] = sample.chi_e - center.chi_e;
  projected[5] = sample.alpha_chi - center.alpha_chi;
  projected[6] = 0.5 * (pressure_characteristic + velocity_characteristic);
  return projected;
}

[[nodiscard]] double TvdSlope(
    double q_im1,
    double q_i,
    double q_ip1) noexcept {
  const double delta_left = q_i - q_im1;
  const double delta_right = q_ip1 - q_i;
  if (delta_left * delta_right <= 0.0) {
    return 0.0;
  }

  const double centered = 0.5 * (q_ip1 - q_im1);
  const double magnitude = std::min({
      std::abs(centered),
      2.0 * std::abs(delta_left),
      2.0 * std::abs(delta_right)});
  return std::copysign(magnitude, centered);
}

[[nodiscard]] PpmSmoothnessProbe BuildSmoothnessProbe(
    int component,
    double q_im1,
    double q_i,
    double q_ip1) noexcept {
  PpmSmoothnessProbe probe;
  probe.component = component;
  probe.q_im1 = q_im1;
  probe.q_i = q_i;
  probe.q_ip1 = q_ip1;

  const double delta_left = q_i - q_im1;
  const double delta_right = q_ip1 - q_i;
  probe.delta_left = delta_left;
  probe.delta_right = delta_right;
  if (delta_left * delta_right >= 0.0) {
    return probe;
  }

  probe.sign_reversal = true;
  const double primary_slope = std::max(std::abs(delta_left), std::abs(delta_right));
  const double reversed_slope = std::min(std::abs(delta_left), std::abs(delta_right));
  const double tolerance = std::max(1.0e-10, 1.0e-3 * primary_slope);
  probe.primary_slope = primary_slope;
  probe.reversed_slope = reversed_slope;
  probe.tolerance = tolerance;
  return probe;
}

[[nodiscard]] bool HasNonSmoothLocalVariation(
    const PpmSmoothnessProbe& probe) noexcept {
  return probe.sign_reversal && probe.reversed_slope > probe.tolerance;
}

void ApplyMonotoneReset(
    double center,
    double& left_edge,
    double& right_edge) noexcept {
  if ((right_edge - center) * (center - left_edge) <= 0.0) {
    left_edge = center;
    right_edge = center;
    return;
  }

  const double left_offset = center - left_edge;
  const double right_offset = right_edge - center;
  if (std::abs(right_offset) >= 2.0 * std::abs(left_offset)) {
    right_edge = center + 2.0 * left_offset;
  }
  if (std::abs(left_offset) >= 2.0 * std::abs(right_offset)) {
    left_edge = center - 2.0 * right_offset;
  }

  if ((right_edge - center) * (center - left_edge) <= 0.0) {
    left_edge = center;
    right_edge = center;
  }
}

[[nodiscard]] PpmCellProfile BuildCellProfile(
    const std::vector<DirectionalPrimitiveState>& primitives,
    std::size_t cell_index) noexcept {
  PpmCellProfile profile;
  profile.w_center = primitives[cell_index];
  profile.basis = BuildCharacteristicBasis(profile.w_center);
  if (!profile.basis.success) {
    profile.troubled = true;
    profile.troubled_component = -2;
    profile.w_left_edge = primitives[cell_index];
    profile.w_right_edge = primitives[cell_index];
    return profile;
  }

  const auto project =
      [&](std::size_t index) noexcept {
        return ProjectPrimitiveDelta(profile.w_center, primitives[index], profile.basis);
      };

  const auto c_im2 = project(cell_index - 2u);
  const auto c_im1 = project(cell_index - 1u);
  const auto c_i = project(cell_index);
  const auto c_ip1 = project(cell_index + 1u);
  const auto c_ip2 = project(cell_index + 2u);

  profile.c_center = c_i;
  const auto w_center_vector = PrimitiveToVector(profile.w_center);
  for (std::size_t component = 0; component < kComponentCount; ++component) {
    const auto smoothness = BuildSmoothnessProbe(
        static_cast<int>(component),
        c_im1[component],
        c_i[component],
        c_ip1[component]);
    if (component == 1u) {
      profile.component1_smoothness = smoothness;
    }
    if (HasNonSmoothLocalVariation(smoothness)) {
      profile.troubled = true;
      if (profile.troubled_component < 0) {
        profile.troubled_component = static_cast<int>(component);
        profile.trouble_q_im1 = smoothness.q_im1;
        profile.trouble_q_i = smoothness.q_i;
        profile.trouble_q_ip1 = smoothness.q_ip1;
        profile.trouble_total_variation =
            std::abs(smoothness.delta_right) +
            std::abs(smoothness.delta_left);
        profile.trouble_end_to_end =
            std::abs(smoothness.q_ip1 - smoothness.q_im1);
      }
    }

    const double delta_im1 = TvdSlope(c_im2[component], c_im1[component], c_i[component]);
    const double delta_i = TvdSlope(c_im1[component], c_i[component], c_ip1[component]);
    const double delta_ip1 = TvdSlope(c_i[component], c_ip1[component], c_ip2[component]);

    double left_edge =
        0.5 * (c_im1[component] + c_i[component]) -
        (delta_i - delta_im1) / 6.0;
    double right_edge =
        0.5 * (c_i[component] + c_ip1[component]) -
        (delta_ip1 - delta_i) / 6.0;
    ApplyMonotoneReset(c_i[component], left_edge, right_edge);

    profile.c_left[component] = left_edge;
    profile.c_right[component] = right_edge;
    profile.c6[component] = 6.0 * c_i[component] - 3.0 * (left_edge + right_edge);
  }

  profile.w_left_edge = VectorToPrimitive(
      Add(w_center_vector, Multiply(profile.basis.R, profile.c_left)),
      profile.w_center.radiation_chi);
  profile.w_right_edge = VectorToPrimitive(
      Add(w_center_vector, Multiply(profile.basis.R, profile.c_right)),
      profile.w_center.radiation_chi);
  if (!profile.w_left_edge.is_physical() || !profile.w_right_edge.is_physical()) {
    profile.troubled = true;
    profile.troubled_component = -3;
    profile.w_left_edge = primitives[cell_index];
    profile.w_right_edge = primitives[cell_index];
  }

  return profile;
}

[[nodiscard]] bool TraceFromCell(
    const PpmCellProfile& profile,
    double dt_s,
    double dx_eff,
    double face_speed,
    bool right_face,
    DirectionalPrimitiveState* traced_primitive) noexcept {
  if (traced_primitive == nullptr || profile.troubled) {
    return false;
  }
  if (!(dt_s > 0.0) ||
      !(dx_eff > 0.0) ||
      !IsFinite(dt_s) ||
      !IsFinite(dx_eff) ||
      !IsFinite(face_speed)) {
    return false;
  }

  PrimitiveVector traced_characteristic{};
  for (std::size_t component = 0; component < kComponentCount; ++component) {
    const double relative_characteristic_speed =
        profile.basis.lambda[component] - face_speed;
    double sigma = std::abs(relative_characteristic_speed) * dt_s / dx_eff;
    if (!IsFinite(sigma) || sigma > 1.0 + kSigmaTolerance) {
      return false;
    }
    sigma = std::clamp(sigma, 0.0, 1.0);

    const double delta = profile.c_right[component] - profile.c_left[component];
    const double correction =
        1.0 - ((2.0 * sigma) / 3.0);
    if (right_face) {
      traced_characteristic[component] =
          relative_characteristic_speed >= 0.0
              ? profile.c_right[component] -
                    0.5 * sigma * (delta - correction * profile.c6[component])
              : profile.c_right[component];
    } else {
      traced_characteristic[component] =
          relative_characteristic_speed <= 0.0
              ? profile.c_left[component] +
                    0.5 * sigma * (delta + correction * profile.c6[component])
              : profile.c_left[component];
    }
  }

  *traced_primitive = VectorToPrimitive(
      Add(
          PrimitiveToVector(profile.w_center),
          Multiply(profile.basis.R, traced_characteristic)),
      profile.w_center.radiation_chi);
  return traced_primitive->is_physical();
}

[[nodiscard]] PpmInterfaceState BuildFirstOrderInterfaceState(
    const DirectionalPrimitiveState& left,
    const DirectionalPrimitiveState& right) noexcept {
  return {left, right, true};
}

}  // namespace

bool DirectionalPrimitiveState::is_physical() const noexcept {
  return IsFinite(rho) &&
         IsFinite(v_n) &&
         IsFinite(v_t1) &&
         IsFinite(v_t2) &&
         IsFinite(pressure) &&
         IsFinite(chi_e) &&
         IsFinite(alpha_chi) &&
         RadiationBundleFinite(radiation_chi) &&
         rho > 0.0 &&
         pressure > 0.0 &&
         chi_e >= 0.0 &&
         alpha_chi >= 0.0;
}

bool PpmInterfaceState::is_complete() const noexcept {
  return left.is_physical() && right.is_physical();
}

bool PpmLineReconstructionResult::is_complete(std::size_t expected_interface_count) const noexcept {
  return success &&
         ghost_layers_consumed >= 3u &&
         interface_count == expected_interface_count &&
         interfaces.size() == expected_interface_count &&
         diagnostics.has_entries() &&
         !report_line.empty();
}

DirectionalPrimitiveState ToDirectionalPrimitive(
    const HydroConservativeState& canonical_state,
    HydroDirection direction) noexcept {
  const auto canonical_primitive = RecoverPrimitiveState(canonical_state);
  switch (direction) {
    case HydroDirection::radial:
      return {
          canonical_primitive.rho,
          canonical_primitive.v_r,
          canonical_primitive.v_theta,
          canonical_primitive.v_phi,
          canonical_primitive.pressure,
          canonical_primitive.chi_e,
          canonical_primitive.alpha_chi,
          canonical_primitive.radiation_chi};
    case HydroDirection::theta:
      return {
          canonical_primitive.rho,
          canonical_primitive.v_theta,
          canonical_primitive.v_r,
          canonical_primitive.v_phi,
          canonical_primitive.pressure,
          canonical_primitive.chi_e,
          canonical_primitive.alpha_chi,
          canonical_primitive.radiation_chi};
    case HydroDirection::phi:
      return {
          canonical_primitive.rho,
          canonical_primitive.v_phi,
          canonical_primitive.v_r,
          canonical_primitive.v_theta,
          canonical_primitive.pressure,
          canonical_primitive.chi_e,
          canonical_primitive.alpha_chi,
          canonical_primitive.radiation_chi};
  }

  return {};
}

HydroConservativeState MakeDirectionalConservativeState(
    const DirectionalPrimitiveState& directional_primitive) noexcept {
  return MakeConservativeState({
      directional_primitive.rho,
      directional_primitive.v_n,
      directional_primitive.v_t1,
      directional_primitive.v_t2,
      directional_primitive.pressure,
      directional_primitive.chi_e,
      directional_primitive.alpha_chi,
      directional_primitive.radiation_chi});
}

PpmLineReconstructionResult ReconstructPpmLine(
    const std::vector<HydroConservativeState>& ghosted_line,
    HydroDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers,
    double dt_s,
    const std::vector<double>& effective_cell_widths,
    const std::vector<double>* moving_face_speeds) noexcept {
  PpmLineReconstructionResult result;
  result.direction = direction;
  result.ghost_layers_consumed = ghost_layers;

  if (ghost_layers < 3u) {
    result.failure_reason = "PPM reconstruction requires at least three ghost layers";
  } else if (interior_cells == 0u) {
    result.failure_reason = "PPM reconstruction requires at least one interior cell";
  } else if (ghosted_line.size() != interior_cells + (2u * ghost_layers)) {
    result.failure_reason = "ghosted PPM line size does not match interior cells plus ghost layers";
  } else if (effective_cell_widths.size() != ghosted_line.size()) {
    result.failure_reason = "effective cell widths must match the ghosted PPM line size";
  } else if (moving_face_speeds != nullptr &&
             moving_face_speeds->size() != interior_cells + 1u) {
    result.failure_reason = "moving face speeds must match the PPM interface count";
  } else if (!(dt_s > 0.0) || !IsFinite(dt_s)) {
    result.failure_reason = "PPM reconstruction requires a finite positive dt";
  } else {
    std::vector<DirectionalPrimitiveState> primitives;
    primitives.reserve(ghosted_line.size());
    for (const auto& conservative_state : ghosted_line) {
      const auto primitive = ToDirectionalPrimitive(conservative_state, direction);
      if (!primitive.is_physical()) {
        result.failure_reason = "ghosted PPM line contains a non-physical primitive state";
        break;
      }
      primitives.push_back(primitive);
    }

    if (result.failure_reason.empty()) {
      std::vector<PpmCellProfile> profiles(primitives.size());
      std::size_t troubled_cell_count = 0u;
      for (std::size_t cell = ghost_layers - 1u;
           cell <= ghost_layers + interior_cells;
           ++cell) {
        profiles[cell] = BuildCellProfile(primitives, cell);
        if (profiles[cell].troubled) {
          ++troubled_cell_count;
        }
      }

      result.interfaces.reserve(interior_cells + 1u);
      for (std::size_t face = 0; face <= interior_cells; ++face) {
        const std::size_t left_cell = ghost_layers + face - 1u;
        const std::size_t right_cell = ghost_layers + face;
        const auto& left_profile = profiles[left_cell];
        const auto& right_profile = profiles[right_cell];

        DirectionalPrimitiveState left_traced{};
        DirectionalPrimitiveState right_traced{};
        const double face_speed =
            moving_face_speeds == nullptr ? 0.0 : (*moving_face_speeds)[face];
        const bool left_trace_failed =
            !TraceFromCell(
                left_profile,
                dt_s,
                effective_cell_widths[left_cell],
                face_speed,
                true,
                &left_traced);
        const bool right_trace_failed =
            !TraceFromCell(
                right_profile,
                dt_s,
                effective_cell_widths[right_cell],
                face_speed,
                false,
                &right_traced);
        const bool fallback =
            left_profile.troubled ||
            right_profile.troubled ||
            left_trace_failed ||
            right_trace_failed;

        if (fallback ||
            !left_traced.is_physical() ||
            !right_traced.is_physical()) {
          auto interface_state = BuildFirstOrderInterfaceState(
              primitives[left_cell],
              primitives[right_cell]);
          interface_state.left_edge_before_trace = left_profile.w_right_edge;
          interface_state.right_edge_before_trace = right_profile.w_left_edge;
          interface_state.edge_state_diagnostics_available = true;
          interface_state.left_profile_troubled = left_profile.troubled;
          interface_state.right_profile_troubled = right_profile.troubled;
          interface_state.left_trace_failed = left_trace_failed;
          interface_state.right_trace_failed = right_trace_failed;
          interface_state.left_troubled_component = left_profile.troubled_component;
          interface_state.right_troubled_component = right_profile.troubled_component;
          interface_state.left_trouble_q_im1 = left_profile.trouble_q_im1;
          interface_state.left_trouble_q_i = left_profile.trouble_q_i;
          interface_state.left_trouble_q_ip1 = left_profile.trouble_q_ip1;
          interface_state.left_trouble_total_variation = left_profile.trouble_total_variation;
          interface_state.left_trouble_end_to_end = left_profile.trouble_end_to_end;
          interface_state.right_trouble_q_im1 = right_profile.trouble_q_im1;
          interface_state.right_trouble_q_i = right_profile.trouble_q_i;
          interface_state.right_trouble_q_ip1 = right_profile.trouble_q_ip1;
          interface_state.right_trouble_total_variation = right_profile.trouble_total_variation;
          interface_state.right_trouble_end_to_end = right_profile.trouble_end_to_end;
          interface_state.left_component1_smoothness = left_profile.component1_smoothness;
          interface_state.right_component1_smoothness = right_profile.component1_smoothness;
          result.interfaces.push_back(interface_state);
          ++result.downgraded_interface_count;
        } else {
          PpmInterfaceState interface_state{left_traced, right_traced, false};
          interface_state.left_edge_before_trace = left_profile.w_right_edge;
          interface_state.right_edge_before_trace = right_profile.w_left_edge;
          interface_state.edge_state_diagnostics_available = true;
          interface_state.left_profile_troubled = left_profile.troubled;
          interface_state.right_profile_troubled = right_profile.troubled;
          interface_state.left_trace_failed = left_trace_failed;
          interface_state.right_trace_failed = right_trace_failed;
          interface_state.left_troubled_component = left_profile.troubled_component;
          interface_state.right_troubled_component = right_profile.troubled_component;
          interface_state.left_trouble_q_im1 = left_profile.trouble_q_im1;
          interface_state.left_trouble_q_i = left_profile.trouble_q_i;
          interface_state.left_trouble_q_ip1 = left_profile.trouble_q_ip1;
          interface_state.left_trouble_total_variation = left_profile.trouble_total_variation;
          interface_state.left_trouble_end_to_end = left_profile.trouble_end_to_end;
          interface_state.right_trouble_q_im1 = right_profile.trouble_q_im1;
          interface_state.right_trouble_q_i = right_profile.trouble_q_i;
          interface_state.right_trouble_q_ip1 = right_profile.trouble_q_ip1;
          interface_state.right_trouble_total_variation = right_profile.trouble_total_variation;
          interface_state.right_trouble_end_to_end = right_profile.trouble_end_to_end;
          interface_state.left_component1_smoothness = left_profile.component1_smoothness;
          interface_state.right_component1_smoothness = right_profile.component1_smoothness;
          result.interfaces.push_back(interface_state);
        }
      }

      result.success = true;
      result.interface_count = result.interfaces.size();
      std::ostringstream report;
      report << "ppm_reconstruction_success=true"
             << "; direction="
             << (direction == HydroDirection::radial ? "radial"
                 : direction == HydroDirection::theta ? "theta"
                 : "phi")
             << "; ghost_layers=" << ghost_layers
             << "; interface_count=" << result.interface_count
             << "; troubled_cell_count=" << troubled_cell_count
             << "; downgraded_interface_count=" << result.downgraded_interface_count;
      result.report_line = report.str();
      AppendDiagnostic(
          result.diagnostics,
          "p1.hydro.reconstruction.ppm.line",
          result.report_line);
    }
  }

  if (!result.success) {
    if (result.failure_reason.empty()) {
      result.failure_reason = "PPM reconstruction failed";
    }
    std::ostringstream report;
    report << "ppm_reconstruction_success=false"
           << "; ghost_layers=" << ghost_layers
           << "; failure_reason=" << result.failure_reason;
    result.report_line = report.str();
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.reconstruction.ppm.failed",
        result.report_line);
  }

  return result;
}

void AccumulateShellWidePpmFallbackMask(
    const PpmLineReconstructionResult& line_reconstruction,
    std::vector<bool>& shell_face_mask) noexcept {
  if (shell_face_mask.empty()) {
    return;
  }

  if (!line_reconstruction.success ||
      line_reconstruction.interfaces.size() < shell_face_mask.size()) {
    std::fill(shell_face_mask.begin(), shell_face_mask.end(), true);
    return;
  }

  for (std::size_t face = 0; face < shell_face_mask.size(); ++face) {
    const auto& interface_state = line_reconstruction.interfaces[face];
    if (interface_state.downgraded_to_first_order ||
        interface_state.left_profile_troubled ||
        interface_state.right_profile_troubled ||
        interface_state.left_trace_failed ||
        interface_state.right_trace_failed) {
      shell_face_mask[face] = true;
    }
  }
}

std::vector<bool> BuildShellWidePpmFallbackMask(
    const std::vector<PpmLineReconstructionResult>& angular_line_reconstructions,
    std::size_t interface_count) noexcept {
  std::vector<bool> shell_face_mask(interface_count, false);
  for (const auto& line_reconstruction : angular_line_reconstructions) {
    AccumulateShellWidePpmFallbackMask(line_reconstruction, shell_face_mask);
  }
  return shell_face_mask;
}

}  // namespace dec3d::hydro
