#include "mesh/boundary/spherical_scalar_remap.hpp"

#include <sstream>

namespace dec3d::mesh {

namespace {

[[nodiscard]] ScalarRemapValidation ValidationFailure(const char* reason) {
  ScalarRemapValidation validation;
  validation.failure_reason = reason;
  std::ostringstream report;
  report << "scalar_remap_topology_valid=false"
         << "; phi_half_turn_available=false"
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false";
  validation.report_line = report.str();
  return validation;
}

[[nodiscard]] ScalarRemapResult RemapFailure(
    const char* kind,
    std::size_t ghost_layer,
    const char* reason) {
  ScalarRemapResult result;
  result.failure_reason = reason;
  std::ostringstream report;
  report << "scalar_remap_success=false"
         << "; remap_kind=" << kind
         << "; ghost_layer=" << ghost_layer
         << "; failure_reason=" << reason
         << "; canonical_state_mutated=false";
  result.report_line = report.str();
  return result;
}

[[nodiscard]] ScalarRemapResult RemapSuccess(
    const char* kind,
    std::size_t ghost_layer,
    ScalarRemapIndex mapped_index) {
  ScalarRemapResult result;
  result.success = true;
  result.mapped_index = mapped_index;
  std::ostringstream report;
  report << "scalar_remap_success=true"
         << "; remap_kind=" << kind
         << "; ghost_layer=" << ghost_layer
         << "; mapped_index=" << mapped_index.radial << "," << mapped_index.theta << ","
         << mapped_index.phi
         << "; value_transform=identity";
  result.report_line = report.str();
  return result;
}

}  // namespace

bool ScalarRemapLayout::is_valid() const noexcept {
  return radial_cells > 0 && theta_cells > 0 && phi_cells > 0;
}

ScalarRemapValidation ValidateScalarHalfTurnTopology(
    const ScalarRemapLayout& layout) noexcept {
  if (!layout.is_valid()) {
    return ValidationFailure("scalar remap topology requires positive cell counts");
  }
  if (layout.phi_cells != 1u && (layout.phi_cells % 2u) != 0u) {
    return ValidationFailure("scalar remap requires one or an even phi cell count");
  }

  ScalarRemapValidation validation;
  validation.success = true;
  validation.phi_half_turn_available = true;
  std::ostringstream report;
  report << "scalar_remap_topology_valid=true"
         << "; phi_half_turn_available=true"
         << "; phi_cells=" << layout.phi_cells
         << "; value_transform=identity";
  validation.report_line = report.str();
  return validation;
}

std::size_t MapPhiHalfTurn(std::size_t phi, std::size_t phi_cells) noexcept {
  if (phi_cells == 0u || phi >= phi_cells) {
    return phi_cells;
  }
  if (phi_cells == 1u) {
    return 0u;
  }
  if ((phi_cells % 2u) != 0u) {
    return phi_cells;
  }
  return (phi + (phi_cells / 2u)) % phi_cells;
}

ScalarRemapResult MapScalarOriginNeighbor(
    std::size_t ghost_layer,
    std::size_t theta,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept {
  if (ghost_layer == 0u) {
    return RemapFailure("origin", ghost_layer, "scalar origin remap requires a positive ghost layer");
  }
  const auto validation = ValidateScalarHalfTurnTopology(layout);
  if (!validation.success) {
    return RemapFailure("origin", ghost_layer, validation.failure_reason.c_str());
  }
  if (theta >= layout.theta_cells || phi >= layout.phi_cells) {
    return RemapFailure("origin", ghost_layer, "scalar origin remap input index is out of range");
  }
  if (ghost_layer > layout.radial_cells) {
    return RemapFailure("origin", ghost_layer, "scalar origin remap ghost layer exceeds radial cells");
  }

  return RemapSuccess(
      "origin",
      ghost_layer,
      ScalarRemapIndex{
          ghost_layer - 1u,
          layout.theta_cells - 1u - theta,
          MapPhiHalfTurn(phi, layout.phi_cells)});
}

ScalarRemapResult MapScalarLowerPoleNeighbor(
    std::size_t ghost_layer,
    std::size_t radial,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept {
  if (ghost_layer == 0u) {
    return RemapFailure("lower_pole", ghost_layer, "scalar pole remap requires a positive ghost layer");
  }
  const auto validation = ValidateScalarHalfTurnTopology(layout);
  if (!validation.success) {
    return RemapFailure("lower_pole", ghost_layer, validation.failure_reason.c_str());
  }
  if (radial >= layout.radial_cells || phi >= layout.phi_cells) {
    return RemapFailure("lower_pole", ghost_layer, "scalar lower pole remap input index is out of range");
  }
  if (ghost_layer > layout.theta_cells) {
    return RemapFailure("lower_pole", ghost_layer, "scalar lower pole remap ghost layer exceeds theta cells");
  }

  return RemapSuccess(
      "lower_pole",
      ghost_layer,
      ScalarRemapIndex{radial, ghost_layer - 1u, MapPhiHalfTurn(phi, layout.phi_cells)});
}

ScalarRemapResult MapScalarUpperPoleNeighbor(
    std::size_t ghost_layer,
    std::size_t radial,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept {
  if (ghost_layer == 0u) {
    return RemapFailure("upper_pole", ghost_layer, "scalar pole remap requires a positive ghost layer");
  }
  const auto validation = ValidateScalarHalfTurnTopology(layout);
  if (!validation.success) {
    return RemapFailure("upper_pole", ghost_layer, validation.failure_reason.c_str());
  }
  if (radial >= layout.radial_cells || phi >= layout.phi_cells) {
    return RemapFailure("upper_pole", ghost_layer, "scalar upper pole remap input index is out of range");
  }
  if (ghost_layer > layout.theta_cells) {
    return RemapFailure("upper_pole", ghost_layer, "scalar upper pole remap ghost layer exceeds theta cells");
  }

  return RemapSuccess(
      "upper_pole",
      ghost_layer,
      ScalarRemapIndex{
          radial,
          layout.theta_cells - ghost_layer,
          MapPhiHalfTurn(phi, layout.phi_cells)});
}

}  // namespace dec3d::mesh
