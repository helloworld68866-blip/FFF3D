#include "mesh/ghost/hydro_halo_exchange.hpp"

#include <mpi.h>

#include <sstream>
#include <vector>

namespace dec3d::mesh {

namespace {

[[nodiscard]] RadialHaloFieldExchange FailedExchange(
    std::size_t ghost_layers,
    std::string failure_reason) noexcept {
  RadialHaloFieldExchange exchange;
  exchange.ghost_layers = ghost_layers;
  exchange.failure_reason = std::move(failure_reason);

  std::ostringstream report;
  report << "radial_halo_exchange_success=false"
         << "; ghost_layers=" << ghost_layers;
  if (!exchange.failure_reason.empty()) {
    report << "; failure_reason=" << exchange.failure_reason;
  }
  exchange.report_line = report.str();
  return exchange;
}

void CopyFlatBufferToGhostArray(
    const std::vector<double>& buffer,
    dec3d::core::Array3D<double>& ghost_values) {
  std::size_t cursor = 0;
  for (std::size_t radial = 0; radial < ghost_values.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < ghost_values.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < ghost_values.extent_phi(); ++phi) {
        ghost_values(radial, theta, phi) = buffer[cursor++];
      }
    }
  }
}

}  // namespace

bool RadialHaloFieldExchange::is_complete(
    std::size_t theta_cells,
    std::size_t phi_cells,
    std::size_t required_ghost_layers) const noexcept {
  if (!success ||
      !mpi_initialized ||
      rank < 0 ||
      rank_count <= 0 ||
      ghost_layers != required_ghost_layers ||
      report_line.empty()) {
    return false;
  }

  const auto has_expected_shape =
      [theta_cells, phi_cells, required_ghost_layers](
          const dec3d::core::Array3D<double>& values) noexcept {
        return values.extent_r() == required_ghost_layers &&
               values.extent_theta() == theta_cells &&
               values.extent_phi() == phi_cells;
      };

  if (has_inner_neighbor && !has_expected_shape(inner_ghost_values)) {
    return false;
  }
  if (has_outer_neighbor && !has_expected_shape(outer_ghost_values)) {
    return false;
  }

  return true;
}

RadialHaloFieldExchange ExchangeRadialHaloField(
    const dec3d::core::Array3D<double>& local_field,
    std::size_t ghost_layers) noexcept {
  if (ghost_layers == 0u) {
    return FailedExchange(ghost_layers, "radial halo exchange requires at least one ghost layer");
  }
  if (local_field.extent_r() < ghost_layers ||
      local_field.extent_theta() == 0u ||
      local_field.extent_phi() == 0u) {
    return FailedExchange(
        ghost_layers,
        "radial halo exchange requires local field extents >= ghost layers and non-empty angular dimensions");
  }

  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  if (mpi_initialized == 0) {
    return FailedExchange(ghost_layers, "radial halo exchange requires MPI to be initialized");
  }

  RadialHaloFieldExchange exchange;
  exchange.mpi_initialized = true;
  exchange.ghost_layers = ghost_layers;
  MPI_Comm_rank(MPI_COMM_WORLD, &exchange.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &exchange.rank_count);

  const std::size_t plane_size = local_field.extent_theta() * local_field.extent_phi();
  const std::size_t halo_value_count = ghost_layers * plane_size;

  if (exchange.rank > 0) {
    exchange.has_inner_neighbor = true;
    exchange.inner_ghost_values = dec3d::core::Array3D<double>(
        ghost_layers,
        local_field.extent_theta(),
        local_field.extent_phi());
    std::vector<double> receive_buffer(halo_value_count, 0.0);
    const int result = MPI_Sendrecv(
        const_cast<double*>(local_field.storage().data()),
        static_cast<int>(halo_value_count),
        MPI_DOUBLE,
        exchange.rank - 1,
        101,
        receive_buffer.data(),
        static_cast<int>(halo_value_count),
        MPI_DOUBLE,
        exchange.rank - 1,
        202,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    if (result != MPI_SUCCESS) {
      return FailedExchange(ghost_layers, "inner radial halo exchange sendrecv failed");
    }
    CopyFlatBufferToGhostArray(receive_buffer, exchange.inner_ghost_values);
  }

  if (exchange.rank + 1 < exchange.rank_count) {
    exchange.has_outer_neighbor = true;
    exchange.outer_ghost_values = dec3d::core::Array3D<double>(
        ghost_layers,
        local_field.extent_theta(),
        local_field.extent_phi());
    std::vector<double> receive_buffer(halo_value_count, 0.0);
    const auto send_offset =
        (local_field.extent_r() - ghost_layers) * plane_size;
    const int result = MPI_Sendrecv(
        const_cast<double*>(local_field.storage().data() + static_cast<std::ptrdiff_t>(send_offset)),
        static_cast<int>(halo_value_count),
        MPI_DOUBLE,
        exchange.rank + 1,
        202,
        receive_buffer.data(),
        static_cast<int>(halo_value_count),
        MPI_DOUBLE,
        exchange.rank + 1,
        101,
        MPI_COMM_WORLD,
        MPI_STATUS_IGNORE);
    if (result != MPI_SUCCESS) {
      return FailedExchange(ghost_layers, "outer radial halo exchange sendrecv failed");
    }
    CopyFlatBufferToGhostArray(receive_buffer, exchange.outer_ghost_values);
  }

  exchange.success = true;
  std::ostringstream report;
  report << "radial_halo_exchange_success=true"
         << "; rank=" << exchange.rank
         << "; rank_count=" << exchange.rank_count
         << "; ghost_layers=" << exchange.ghost_layers
         << "; inner_neighbor=" << (exchange.has_inner_neighbor ? "true" : "false")
         << "; outer_neighbor=" << (exchange.has_outer_neighbor ? "true" : "false");
  exchange.report_line = report.str();
  return exchange;
}

}  // namespace dec3d::mesh
