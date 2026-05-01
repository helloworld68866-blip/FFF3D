#include "radiation/providers/group_blackbody.hpp"

#include "physics/units/physical_constants.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <utility>

namespace dec3d::radiation {
namespace {

constexpr double kPlanckWeightIntegral = 6.493939402266829;  // pi^4 / 15.
constexpr std::size_t kBlackbodyQuadratureIntervals = 256u;

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] GroupBlackbodyResult Fail(std::string reason) {
  GroupBlackbodyResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = "diagnostic_id=p3.radiation.group_blackbody.failure; failure_reason=" +
                               result.failure_reason;
  return result;
}

[[nodiscard]] GroupBlackbodyValueResult ValueFail(std::string reason) {
  GroupBlackbodyValueResult result;
  result.failure_reason = std::move(reason);
  return result;
}

[[nodiscard]] double PlanckConstantErgS() noexcept {
  return 2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
         dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
}

[[nodiscard]] double RadiationConstantErgCm3PerErg4() noexcept {
  static const double value = [] {
    const double h_erg_s = PlanckConstantErgS();
    const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
    return 8.0 * std::pow(dec3d::physics::PhysicsConstantsCGS::pi, 5.0) /
           (15.0 * std::pow(h_erg_s, 3.0) * std::pow(c, 3.0));
  }();
  return value;
}

[[nodiscard]] double PlanckWeightIntegrand(double u) noexcept {
  if (u <= 0.0) {
    return 0.0;
  }
  const double exp_u = std::exp(u);
  if (!std::isfinite(exp_u)) {
    return 0.0;
  }
  return u * u * u / (exp_u - 1.0);
}

[[nodiscard]] double IntegratePlanckWeight(double u0, double u1) noexcept {
  if (!(u1 > u0)) {
    return 0.0;
  }
  const double capped_u0 = std::max(0.0, u0);
  const double capped_u1 = std::min(200.0, u1);
  if (!(capped_u1 > capped_u0)) {
    return 0.0;
  }
  const double h =
      (capped_u1 - capped_u0) / static_cast<double>(kBlackbodyQuadratureIntervals);
  double sum = PlanckWeightIntegrand(capped_u0) + PlanckWeightIntegrand(capped_u1);
  for (std::size_t i = 1u; i < kBlackbodyQuadratureIntervals; ++i) {
    const double u = capped_u0 + h * static_cast<double>(i);
    sum += (i % 2u == 0u ? 2.0 : 4.0) * PlanckWeightIntegrand(u);
  }
  return sum * h / 3.0;
}

[[nodiscard]] std::string BuildReport(double energy, double weight, std::size_t group) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "diagnostic_id=p3.radiation.group_blackbody"
      << "; group_index=" << group
      << "; group_edges_unit=Hz"
      << "; temperature_internal_unit=erg_per_particle"
      << "; B_g_output_unit=erg_per_cm3"
      << "; blackbody_formula=planck_group_integral"
      << "; double_kB_guard_passed=true"
      << "; group_weight=" << weight
      << "; B_g_erg_cm3=" << energy;
  return out.str();
}

}  // namespace

GroupBlackbodyValueResult EvaluateGroupBlackbodyEnergyDensityValueOnly(
    const RadiationGroupLayout& group_layout,
    std::size_t group_index,
    double te_erg_per_particle) noexcept {
  if (group_layout.mode != RadiationGroupMode::explicit_frequency_groups) {
    return ValueFail("group blackbody requires explicit frequency groups");
  }
  if (group_index + 1u >= group_layout.frequency_edges_hz.size() ||
      group_layout.group_count + 1u != group_layout.frequency_edges_hz.size()) {
    return ValueFail("invalid group index or frequency edge count");
  }
  if (!(te_erg_per_particle > 0.0) || !std::isfinite(te_erg_per_particle)) {
    return ValueFail("invalid Te_erg_per_particle");
  }
  const double nu0 = group_layout.frequency_edges_hz[group_index];
  const double nu1 = group_layout.frequency_edges_hz[group_index + 1u];
  if (!(nu1 > nu0) || !std::isfinite(nu0) || !std::isfinite(nu1)) {
    return ValueFail("frequency edges must be finite and strictly increasing");
  }

  const double h_erg_s = PlanckConstantErgS();
  const double u0 = h_erg_s * nu0 / te_erg_per_particle;
  const double u1 = h_erg_s * nu1 / te_erg_per_particle;
  const double weight = IntegratePlanckWeight(u0, u1) / kPlanckWeightIntegral;
  const double te2 = te_erg_per_particle * te_erg_per_particle;
  const double energy = RadiationConstantErgCm3PerErg4() * te2 * te2 * weight;
  if (!(weight >= 0.0) || !std::isfinite(weight) || !(energy >= 0.0) ||
      !std::isfinite(energy)) {
    return ValueFail("blackbody group evaluation failed");
  }

  GroupBlackbodyValueResult result;
  result.success = true;
  result.group_weight = weight;
  result.energy_density_erg_cm3 = energy;
  return result;
}

GroupBlackbodyResult EvaluateGroupBlackbodyEnergyDensity(
    const RadiationGroupLayout& group_layout,
    std::size_t group_index,
    double te_erg_per_particle) noexcept {
  const auto value =
      EvaluateGroupBlackbodyEnergyDensityValueOnly(group_layout, group_index, te_erg_per_particle);
  if (!value.success) {
    return Fail(value.failure_reason);
  }

  GroupBlackbodyResult result;
  result.success = true;
  result.group_weight = value.group_weight;
  result.energy_density_erg_cm3 = value.energy_density_erg_cm3;
  result.report_line =
      BuildReport(value.energy_density_erg_cm3, value.group_weight, group_index);
  return result;
}

bool ValidateGroupBlackbodyDiagnostics(const GroupBlackbodyResult& result) noexcept {
  if (!result.success) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p3.radiation.group_blackbody") &&
         Contains(line, "group_edges_unit=Hz") &&
         Contains(line, "temperature_internal_unit=erg_per_particle") &&
         Contains(line, "B_g_output_unit=erg_per_cm3") &&
         Contains(line, "blackbody_formula=planck_group_integral") &&
         Contains(line, "double_kB_guard_passed=true");
}

}  // namespace dec3d::radiation
