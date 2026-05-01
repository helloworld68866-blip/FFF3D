#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::mesh {

struct RadialOwnershipSlice {
  std::size_t begin_index{0};
  std::size_t end_index{0};
  bool initialized{false};

  [[nodiscard]] std::size_t local_cell_count() const noexcept {
    return end_index >= begin_index ? end_index - begin_index : 0;
  }
};

struct RadialOwnershipDecomposition {
  bool valid{false};
  std::size_t global_radial_cells{0};
  std::string failure_reason;
  std::vector<RadialOwnershipSlice> slices;

  [[nodiscard]] bool is_valid() const noexcept { return valid; }
};

[[nodiscard]] RadialOwnershipDecomposition BuildRadialOwnership(
    std::size_t global_radial_cells,
    std::size_t rank_count) noexcept;

}  // namespace dec3d::mesh
