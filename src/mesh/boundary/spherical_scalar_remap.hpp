#pragma once

#include <cstddef>
#include <string>

namespace dec3d::mesh {

struct ScalarRemapLayout {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};

  [[nodiscard]] bool is_valid() const noexcept;
};

struct ScalarRemapIndex {
  std::size_t radial{0};
  std::size_t theta{0};
  std::size_t phi{0};
};

struct ScalarRemapValidation {
  bool success{false};
  bool phi_half_turn_available{false};
  std::string failure_reason;
  std::string report_line;
};

struct ScalarRemapResult {
  bool success{false};
  ScalarRemapIndex mapped_index;
  std::string failure_reason;
  std::string report_line;
};

[[nodiscard]] ScalarRemapValidation ValidateScalarHalfTurnTopology(
    const ScalarRemapLayout& layout) noexcept;

[[nodiscard]] std::size_t MapPhiHalfTurn(
    std::size_t phi,
    std::size_t phi_cells) noexcept;

[[nodiscard]] ScalarRemapResult MapScalarOriginNeighbor(
    std::size_t ghost_layer,
    std::size_t theta,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept;

[[nodiscard]] ScalarRemapResult MapScalarLowerPoleNeighbor(
    std::size_t ghost_layer,
    std::size_t radial,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept;

[[nodiscard]] ScalarRemapResult MapScalarUpperPoleNeighbor(
    std::size_t ghost_layer,
    std::size_t radial,
    std::size_t phi,
    const ScalarRemapLayout& layout) noexcept;

}  // namespace dec3d::mesh
