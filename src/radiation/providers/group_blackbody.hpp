#pragma once

#include "radiation/groups/radiation_group_layout.hpp"

#include <cstddef>
#include <string>

namespace dec3d::radiation {

struct GroupBlackbodyResult {
  bool success{false};
  double energy_density_erg_cm3{0.0};
  double group_weight{0.0};
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct GroupBlackbodyValueResult {
  bool success{false};
  double energy_density_erg_cm3{0.0};
  double group_weight{0.0};
  std::string failure_reason;
};

[[nodiscard]] GroupBlackbodyValueResult EvaluateGroupBlackbodyEnergyDensityValueOnly(
    const RadiationGroupLayout& group_layout,
    std::size_t group_index,
    double te_erg_per_particle) noexcept;

[[nodiscard]] GroupBlackbodyResult EvaluateGroupBlackbodyEnergyDensity(
    const RadiationGroupLayout& group_layout,
    std::size_t group_index,
    double te_erg_per_particle) noexcept;

[[nodiscard]] bool ValidateGroupBlackbodyDiagnostics(
    const GroupBlackbodyResult& result) noexcept;

}  // namespace dec3d::radiation
