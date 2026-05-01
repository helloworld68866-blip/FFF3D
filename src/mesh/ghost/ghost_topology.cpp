#include "mesh/ghost/ghost_topology.hpp"

namespace dec3d::mesh {

[[maybe_unused]] static constexpr const char* kGhostLayerRole = "generic.ghost.topology";

GhostTopology BuildGhostTopology(
    const RadialOwnershipSlice& local_slice,
    std::size_t rank_index,
    std::size_t rank_count) noexcept {
  GhostTopology topology;
  topology.initialized = local_slice.initialized && rank_count > 0 && rank_index < rank_count;

  if (!topology.initialized) {
    return topology;
  }

  topology.has_radial_inner_neighbor = rank_index > 0;
  topology.has_radial_outer_neighbor = rank_index + 1 < rank_count;
  return topology;
}

}  // namespace dec3d::mesh
