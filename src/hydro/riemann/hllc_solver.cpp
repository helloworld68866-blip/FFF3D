#include "hydro/riemann/hllc_solver.hpp"

#include "state/hydro_state/hydro_view.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace dec3d::hydro {

namespace {

[[nodiscard]] bool IsFinite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] bool RadiationBundleFinite(
    const std::vector<double>& values) noexcept {
  return std::all_of(values.begin(), values.end(), [](double value) {
    return std::isfinite(value) && value >= 0.0;
  });
}

[[nodiscard]] bool SameRadiationBundleSize(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  return lhs.radiation_chi.size() == rhs.radiation_chi.size();
}

[[nodiscard]] std::vector<double> ScaleBundle(
    const std::vector<double>& values,
    double factor) {
  std::vector<double> result(values.size(), 0.0);
  for (std::size_t index = 0; index < values.size(); ++index) {
    result[index] = values[index] * factor;
  }
  return result;
}

[[nodiscard]] std::vector<double> AddBundle(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) {
  std::vector<double> result(lhs.size(), 0.0);
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    result[index] = lhs[index] + rhs[index];
  }
  return result;
}

[[nodiscard]] std::vector<double> SubtractBundle(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) {
  std::vector<double> result(lhs.size(), 0.0);
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    result[index] = lhs[index] - rhs[index];
  }
  return result;
}

[[nodiscard]] double SoundSpeed(const HydroPrimitiveState& primitive) noexcept {
  return std::sqrt(dec3d::state::HydroIdealGasGamma() * primitive.pressure / primitive.rho);
}

[[nodiscard]] std::string FormatPrimitiveForDiagnostics(
    const char* label,
    const HydroPrimitiveState& primitive) {
  std::ostringstream output;
  output << std::setprecision(17)
         << label << "{rho=" << primitive.rho
         << ",v_r=" << primitive.v_r
         << ",v_theta=" << primitive.v_theta
         << ",v_phi=" << primitive.v_phi
         << ",pressure=" << primitive.pressure
         << ",chi_e=" << primitive.chi_e
         << ",alpha_chi=" << primitive.alpha_chi
         << '}';
  return output.str();
}

[[nodiscard]] HydroConservativeState Scale(
    const HydroConservativeState& state,
    double factor) noexcept {
  HydroConservativeState result{
      state.rho * factor,
      state.mom_r * factor,
      state.mom_theta * factor,
      state.mom_phi * factor,
      state.e_fluid_total * factor,
      state.chi_e * factor,
      state.alpha_chi * factor};
  result.radiation_chi = ScaleBundle(state.radiation_chi, factor);
  return result;
}

[[nodiscard]] HydroConservativeState Add(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  HydroConservativeState result{
      lhs.rho + rhs.rho,
      lhs.mom_r + rhs.mom_r,
      lhs.mom_theta + rhs.mom_theta,
      lhs.mom_phi + rhs.mom_phi,
      lhs.e_fluid_total + rhs.e_fluid_total,
      lhs.chi_e + rhs.chi_e,
      lhs.alpha_chi + rhs.alpha_chi};
  result.radiation_chi = AddBundle(lhs.radiation_chi, rhs.radiation_chi);
  return result;
}

[[nodiscard]] HydroConservativeState Subtract(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  HydroConservativeState result{
      lhs.rho - rhs.rho,
      lhs.mom_r - rhs.mom_r,
      lhs.mom_theta - rhs.mom_theta,
      lhs.mom_phi - rhs.mom_phi,
      lhs.e_fluid_total - rhs.e_fluid_total,
      lhs.chi_e - rhs.chi_e,
      lhs.alpha_chi - rhs.alpha_chi};
  result.radiation_chi = SubtractBundle(lhs.radiation_chi, rhs.radiation_chi);
  return result;
}

[[nodiscard]] HydroConservativeState BuildStarState(
    const HydroPrimitiveState& primitive,
    const HydroConservativeState& conservative,
    double s_k,
    double s_star,
    double p_star) noexcept {
  const double scale = (s_k - primitive.v_r) / (s_k - s_star);
  const double rho_star = conservative.rho * scale;

  HydroConservativeState star;
  star.rho = rho_star;
  star.mom_r = rho_star * s_star;
  star.mom_theta = rho_star * primitive.v_theta;
  star.mom_phi = rho_star * primitive.v_phi;
  star.e_fluid_total =
      ((s_k - primitive.v_r) * conservative.e_fluid_total -
       primitive.pressure * primitive.v_r +
       p_star * s_star) /
      (s_k - s_star);
  star.chi_e = conservative.chi_e * scale;
  star.alpha_chi = conservative.alpha_chi * scale;
  star.radiation_chi = ScaleBundle(conservative.radiation_chi, scale);
  return star;
}

}  // namespace

bool HydroPrimitiveState::is_physical() const noexcept {
  return IsFinite(rho) &&
         IsFinite(v_r) &&
         IsFinite(v_theta) &&
         IsFinite(v_phi) &&
         IsFinite(pressure) &&
         IsFinite(chi_e) &&
         IsFinite(alpha_chi) &&
         RadiationBundleFinite(radiation_chi) &&
         rho > 0.0 &&
         pressure > 0.0 &&
         chi_e >= 0.0 &&
         alpha_chi >= 0.0;
}

bool HydroConservativeState::is_finite() const noexcept {
  return IsFinite(rho) &&
         IsFinite(mom_r) &&
         IsFinite(mom_theta) &&
         IsFinite(mom_phi) &&
         IsFinite(e_fluid_total) &&
         IsFinite(chi_e) &&
         IsFinite(alpha_chi) &&
         RadiationBundleFinite(radiation_chi);
}

bool HllcWaveStructure::is_complete() const noexcept {
  return IsFinite(s_left) &&
         IsFinite(s_star) &&
         IsFinite(s_right) &&
         s_left <= s_star &&
         s_star <= s_right;
}

bool HllcResult::is_complete() const noexcept {
  return success &&
         waves.is_complete() &&
         left_star_state.is_finite() &&
         right_star_state.is_finite() &&
         interface_flux.is_finite() &&
         active_region != HllcActiveRegion::invalid;
}

bool MovingInterfaceHllcFluxResult::is_complete() const noexcept {
  return success &&
         waves.is_complete() &&
         active_region != HllcActiveRegion::invalid &&
         static_zero_region != HllcActiveRegion::invalid &&
         selected_state.is_finite() &&
         selected_flux.is_finite() &&
         flux.is_finite() &&
         IsFinite(face_speed) &&
         mode == "moving_interface_hllc" &&
         remap_order == "none" &&
         failure_reason.empty();
}

HydroConservativeState MakeConservativeState(
    const HydroPrimitiveState& primitive) noexcept {
  HydroConservativeState conservative;
  if (!primitive.is_physical()) {
    return conservative;
  }

  conservative.rho = primitive.rho;
  conservative.mom_r = primitive.rho * primitive.v_r;
  conservative.mom_theta = primitive.rho * primitive.v_theta;
  conservative.mom_phi = primitive.rho * primitive.v_phi;
  conservative.e_fluid_total =
      primitive.pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) +
      0.5 * primitive.rho *
          (primitive.v_r * primitive.v_r +
           primitive.v_theta * primitive.v_theta +
           primitive.v_phi * primitive.v_phi);
  conservative.chi_e = primitive.chi_e;
  conservative.alpha_chi = primitive.alpha_chi;
  conservative.radiation_chi = primitive.radiation_chi;
  return conservative;
}

HydroPrimitiveState RecoverPrimitiveState(
    const HydroConservativeState& conservative) noexcept {
  HydroPrimitiveState primitive;
  if (!conservative.is_finite() || conservative.rho <= 0.0) {
    return primitive;
  }

  primitive.rho = conservative.rho;
  primitive.v_r = conservative.mom_r / conservative.rho;
  primitive.v_theta = conservative.mom_theta / conservative.rho;
  primitive.v_phi = conservative.mom_phi / conservative.rho;
  const double kinetic =
      0.5 * conservative.rho *
      (primitive.v_r * primitive.v_r +
       primitive.v_theta * primitive.v_theta +
       primitive.v_phi * primitive.v_phi);
  primitive.pressure =
      (dec3d::state::HydroIdealGasGamma() - 1.0) *
      (conservative.e_fluid_total - kinetic);
  primitive.chi_e = conservative.chi_e;
  primitive.alpha_chi = conservative.alpha_chi;
  primitive.radiation_chi = conservative.radiation_chi;
  return primitive;
}

HydroConservativeState ComputePhysicalFlux(
    const HydroConservativeState& conservative) noexcept {
  const auto primitive = RecoverPrimitiveState(conservative);
  HydroConservativeState flux;
  if (!primitive.is_physical()) {
    return flux;
  }

  flux.rho = conservative.mom_r;
  flux.mom_r = conservative.mom_r * primitive.v_r + primitive.pressure;
  flux.mom_theta = conservative.mom_theta * primitive.v_r;
  flux.mom_phi = conservative.mom_phi * primitive.v_r;
  flux.e_fluid_total =
      primitive.v_r * (conservative.e_fluid_total + primitive.pressure);
  flux.chi_e = conservative.chi_e * primitive.v_r;
  flux.alpha_chi = conservative.alpha_chi * primitive.v_r;
  flux.radiation_chi = ScaleBundle(conservative.radiation_chi, primitive.v_r);
  return flux;
}

HllcResult SolveHllcRiemann(
    const HydroConservativeState& left,
    const HydroConservativeState& right) noexcept {
  HllcResult result;

  if (!SameRadiationBundleSize(left, right)) {
    result.failure_reason = "radiation scalar bundle size mismatch";
    return result;
  }

  const auto left_primitive = RecoverPrimitiveState(left);
  const auto right_primitive = RecoverPrimitiveState(right);
  if (!left_primitive.is_physical() || !right_primitive.is_physical()) {
    result.failure_reason = "left or right state is not physical";
    return result;
  }

  const double c_left = SoundSpeed(left_primitive);
  const double c_right = SoundSpeed(right_primitive);
  result.waves.s_left =
      std::min(left_primitive.v_r - c_left, right_primitive.v_r - c_right);
  result.waves.s_right =
      std::max(left_primitive.v_r + c_left, right_primitive.v_r + c_right);

  const double denominator =
      left_primitive.rho * (result.waves.s_left - left_primitive.v_r) -
      right_primitive.rho * (result.waves.s_right - right_primitive.v_r);
  if (!IsFinite(denominator) || std::abs(denominator) < 1.0e-12) {
    result.failure_reason = "HLLC contact-wave denominator is ill-conditioned";
    return result;
  }

  result.waves.s_star =
      (right_primitive.pressure - left_primitive.pressure +
       left_primitive.rho * left_primitive.v_r * (result.waves.s_left - left_primitive.v_r) -
       right_primitive.rho * right_primitive.v_r * (result.waves.s_right - right_primitive.v_r)) /
      denominator;

  if (!result.waves.is_complete()) {
    result.failure_reason = "HLLC wave speeds are incomplete";
    return result;
  }

  const double p_star_left =
      left_primitive.pressure +
      left_primitive.rho * (result.waves.s_left - left_primitive.v_r) *
          (result.waves.s_star - left_primitive.v_r);
  const double p_star_right =
      right_primitive.pressure +
      right_primitive.rho * (result.waves.s_right - right_primitive.v_r) *
          (result.waves.s_star - right_primitive.v_r);
  const double p_star = 0.5 * (p_star_left + p_star_right);
  if (!IsFinite(p_star) || p_star <= 0.0) {
    std::ostringstream reason;
    reason << std::setprecision(17)
           << "HLLC star pressure is not physical"
           << "; p_star=" << p_star
           << "; p_star_left=" << p_star_left
           << "; p_star_right=" << p_star_right
           << "; s_left=" << result.waves.s_left
           << "; s_star=" << result.waves.s_star
           << "; s_right=" << result.waves.s_right
           << "; " << FormatPrimitiveForDiagnostics("left", left_primitive)
           << "; " << FormatPrimitiveForDiagnostics("right", right_primitive);
    result.failure_reason = reason.str();
    return result;
  }

  result.left_star_state = BuildStarState(
      left_primitive,
      left,
      result.waves.s_left,
      result.waves.s_star,
      p_star);
  result.right_star_state = BuildStarState(
      right_primitive,
      right,
      result.waves.s_right,
      result.waves.s_star,
      p_star);

  const auto left_flux = ComputePhysicalFlux(left);
  const auto right_flux = ComputePhysicalFlux(right);

  if (0.0 <= result.waves.s_left) {
    result.interface_flux = left_flux;
    result.active_region = HllcActiveRegion::left_flux;
  } else if (0.0 <= result.waves.s_star) {
    result.interface_flux = Add(
        left_flux,
        Scale(Subtract(result.left_star_state, left), result.waves.s_left));
    result.active_region = HllcActiveRegion::left_star;
  } else if (0.0 <= result.waves.s_right) {
    result.interface_flux = Add(
        right_flux,
        Scale(Subtract(result.right_star_state, right), result.waves.s_right));
    result.active_region = HllcActiveRegion::right_star;
  } else {
    result.interface_flux = right_flux;
    result.active_region = HllcActiveRegion::right_flux;
  }

  result.success = true;
  return result;
}

bool MovingInterfaceBranch::is_complete() const noexcept {
  return active_region != HllcActiveRegion::invalid &&
         upwind_state.is_finite() &&
         upwind_flux.is_finite();
}

MovingInterfaceBranch SelectMovingInterfaceBranch(
    const HllcResult& result,
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept {
  MovingInterfaceBranch branch;
  if (!result.is_complete() || !IsFinite(face_speed)) {
    return branch;
  }

  const auto left_flux = ComputePhysicalFlux(left);
  const auto right_flux = ComputePhysicalFlux(right);
  const auto left_star_flux = Add(
      left_flux,
      Scale(Subtract(result.left_star_state, left), result.waves.s_left));
  const auto right_star_flux = Add(
      right_flux,
      Scale(Subtract(result.right_star_state, right), result.waves.s_right));

  if (face_speed <= result.waves.s_left) {
    branch.active_region = HllcActiveRegion::left_flux;
    branch.upwind_state = left;
    branch.upwind_flux = left_flux;
  } else if (face_speed <= result.waves.s_star) {
    branch.active_region = HllcActiveRegion::left_star;
    branch.upwind_state = result.left_star_state;
    branch.upwind_flux = left_star_flux;
  } else if (face_speed <= result.waves.s_right) {
    branch.active_region = HllcActiveRegion::right_star;
    branch.upwind_state = result.right_star_state;
    branch.upwind_flux = right_star_flux;
  } else {
    branch.active_region = HllcActiveRegion::right_flux;
    branch.upwind_state = right;
    branch.upwind_flux = right_flux;
  }

  return branch;
}

HydroConservativeState ComputeAleCorrectedFlux(
    const HllcResult& result,
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept {
  const auto branch =
      SelectMovingInterfaceBranch(result, left, right, face_speed);
  if (!branch.is_complete()) {
    return {};
  }

  return Subtract(
      branch.upwind_flux,
      Scale(branch.upwind_state, face_speed));
}

MovingInterfaceHllcFluxResult SolveMovingInterfaceHllcFlux(
    const HydroConservativeState& left,
    const HydroConservativeState& right,
    double face_speed) noexcept {
  MovingInterfaceHllcFluxResult moving;
  moving.face_speed = face_speed;
  moving.mode = "moving_interface_hllc";
  moving.remap_order = "none";

  const auto hllc = SolveHllcRiemann(left, right);
  if (!hllc.is_complete()) {
    moving.failure_reason =
        hllc.failure_reason.empty() ? "HLLC solve failed" : hllc.failure_reason;
    return moving;
  }
  moving.waves = hllc.waves;
  moving.static_zero_region = hllc.active_region;

  const auto branch = SelectMovingInterfaceBranch(
      hllc,
      left,
      right,
      face_speed);
  if (!branch.is_complete()) {
    moving.failure_reason =
        "moving-interface HLLC branch classification failed";
    return moving;
  }

  moving.active_region = branch.active_region;
  moving.selected_state = branch.upwind_state;
  moving.selected_flux = branch.upwind_flux;
  moving.flux = Subtract(
      branch.upwind_flux,
      Scale(branch.upwind_state, face_speed));
  if (!moving.flux.is_finite()) {
    moving.failure_reason = "moving-interface HLLC flux is not finite";
    return moving;
  }

  moving.success = true;
  return moving;
}

}  // namespace dec3d::hydro
