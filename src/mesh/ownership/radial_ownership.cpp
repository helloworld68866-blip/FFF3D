#include "mesh/ownership/radial_ownership.hpp"

namespace dec3d::mesh {

RadialOwnershipDecomposition BuildRadialOwnership(
    std::size_t global_radial_cells,
    std::size_t rank_count) noexcept {
  RadialOwnershipDecomposition decomposition;
  decomposition.global_radial_cells = global_radial_cells;

  if (rank_count == 0) {
    decomposition.failure_reason = "rank_count must be positive";
    return decomposition;
  }

  if (global_radial_cells == 0) {
    decomposition.failure_reason = "global_radial_cells must be positive";
    return decomposition;
  }

  if (rank_count > global_radial_cells) {
    decomposition.failure_reason = "rank_count may not exceed global_radial_cells in P0";
    return decomposition;
  }

  decomposition.valid = true;
  decomposition.slices.resize(rank_count);

  const std::size_t base_cells = global_radial_cells / rank_count;
  const std::size_t remainder = global_radial_cells % rank_count;

  std::size_t cursor = 0;
  for (std::size_t rank_index = 0; rank_index < rank_count; ++rank_index) {
    const std::size_t local_cells = base_cells + (rank_index < remainder ? 1 : 0);
    decomposition.slices[rank_index].begin_index = cursor;
    decomposition.slices[rank_index].end_index = cursor + local_cells;
    decomposition.slices[rank_index].initialized = true;
    cursor += local_cells;
  }

  return decomposition;
}

}  // namespace dec3d::mesh
