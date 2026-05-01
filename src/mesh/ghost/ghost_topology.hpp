#pragma once

#include "mesh/ownership/radial_ownership.hpp"

namespace dec3d::mesh {

struct GhostTopology {
  bool initialized{false};
  bool has_radial_inner_neighbor{false};
  bool has_radial_outer_neighbor{false};

  [[nodiscard]] bool is_initialized() const noexcept { return initialized; }
  [[nodiscard]] bool has_neighbors() const noexcept {
    return has_radial_inner_neighbor || has_radial_outer_neighbor;
  }
};

[[nodiscard]] GhostTopology BuildGhostTopology(
    const RadialOwnershipSlice& local_slice,
    std::size_t rank_index,
    std::size_t rank_count) noexcept;

}  // namespace dec3d::mesh
