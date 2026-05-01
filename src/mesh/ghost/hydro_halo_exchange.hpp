#pragma once

#include "core/array/array3d.hpp"

#include <cstddef>
#include <string>

namespace dec3d::mesh {

struct RadialHaloFieldExchange {
  bool success{false};
  bool mpi_initialized{false};
  int rank{-1};
  int rank_count{0};
  bool has_inner_neighbor{false};
  bool has_outer_neighbor{false};
  std::size_t ghost_layers{0};
  dec3d::core::Array3D<double> inner_ghost_values;
  dec3d::core::Array3D<double> outer_ghost_values;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete(
      std::size_t theta_cells,
      std::size_t phi_cells,
      std::size_t required_ghost_layers) const noexcept;
};

[[nodiscard]] RadialHaloFieldExchange ExchangeRadialHaloField(
    const dec3d::core::Array3D<double>& local_field,
    std::size_t ghost_layers) noexcept;

}  // namespace dec3d::mesh
