#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
  int mpi_initialized = 0;

  try {
    MPI_Init(&argc, &argv);
    MPI_Initialized(&mpi_initialized);

    int rank = 0;
    int rank_count = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

    DEC3D_CHECK_EQ(rank_count, 3);

    dec3d::core::Array3D<double> local_field(2, 2, 2);
    for (std::size_t radial = 0; radial < local_field.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < local_field.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < local_field.extent_phi(); ++phi) {
          local_field(radial, theta, phi) =
              1000.0 * static_cast<double>(rank) +
              100.0 * static_cast<double>(radial) +
              10.0 * static_cast<double>(theta) +
              static_cast<double>(phi);
        }
      }
    }

    const auto exchange = dec3d::mesh::ExchangeRadialHaloField(local_field, 1u);
    DEC3D_CHECK(exchange.is_complete(2u, 2u, 1u));
    DEC3D_CHECK(exchange.report_line.find("rank=") != std::string::npos);

    if (rank > 0) {
      DEC3D_CHECK(exchange.has_inner_neighbor);
      for (std::size_t theta = 0; theta < 2u; ++theta) {
        for (std::size_t phi = 0; phi < 2u; ++phi) {
          const double expected =
              1000.0 * static_cast<double>(rank - 1) +
              100.0 * 1.0 +
              10.0 * static_cast<double>(theta) +
              static_cast<double>(phi);
          DEC3D_CHECK(std::abs(exchange.inner_ghost_values(0, theta, phi) - expected) < 1.0e-12);
        }
      }
    } else {
      DEC3D_CHECK(!exchange.has_inner_neighbor);
    }

    if (rank + 1 < rank_count) {
      DEC3D_CHECK(exchange.has_outer_neighbor);
      for (std::size_t theta = 0; theta < 2u; ++theta) {
        for (std::size_t phi = 0; phi < 2u; ++phi) {
          const double expected =
              1000.0 * static_cast<double>(rank + 1) +
              10.0 * static_cast<double>(theta) +
              static_cast<double>(phi);
          DEC3D_CHECK(std::abs(exchange.outer_ghost_values(0, theta, phi) - expected) < 1.0e-12);
        }
      }
    } else {
      DEC3D_CHECK(!exchange.has_outer_neighbor);
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return 1;
  }
}
