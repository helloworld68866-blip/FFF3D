#include "mesh/ghost/ghost_topology.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

void CheckCoverage(const dec3d::mesh::RadialOwnershipDecomposition& decomposition) {
  std::size_t expected_begin = 0;

  for (const auto& slice : decomposition.slices) {
    DEC3D_CHECK(slice.initialized);
    DEC3D_CHECK_EQ(slice.begin_index, expected_begin);
    DEC3D_CHECK(slice.end_index >= slice.begin_index);
    expected_begin = slice.end_index;
  }

  DEC3D_CHECK_EQ(expected_begin, decomposition.global_radial_cells);
}

}  // namespace

int main() {
  try {
    using namespace dec3d::mesh;

    const auto one_rank = BuildRadialOwnership(8, 1);
    DEC3D_CHECK(one_rank.is_valid());
    DEC3D_CHECK_EQ(one_rank.slices.size(), static_cast<std::size_t>(1));
    CheckCoverage(one_rank);
    DEC3D_CHECK_EQ(one_rank.slices[0].local_cell_count(), static_cast<std::size_t>(8));

    const auto two_rank = BuildRadialOwnership(8, 2);
    DEC3D_CHECK(two_rank.is_valid());
    DEC3D_CHECK_EQ(two_rank.slices.size(), static_cast<std::size_t>(2));
    CheckCoverage(two_rank);

    const auto three_rank = BuildRadialOwnership(10, 3);
    DEC3D_CHECK(three_rank.is_valid());
    DEC3D_CHECK_EQ(three_rank.slices.size(), static_cast<std::size_t>(3));
    CheckCoverage(three_rank);

    const auto oversubscribed = BuildRadialOwnership(2, 4);
    DEC3D_CHECK(!oversubscribed.is_valid());
    DEC3D_CHECK(oversubscribed.slices.empty());

    const auto zero_rank = BuildRadialOwnership(8, 0);
    DEC3D_CHECK(!zero_rank.is_valid());
    DEC3D_CHECK(zero_rank.slices.empty());

    const auto single_rank_topology = BuildGhostTopology(one_rank.slices[0], 0, one_rank.slices.size());
    DEC3D_CHECK(single_rank_topology.is_initialized());
    DEC3D_CHECK(!single_rank_topology.has_neighbors());

    const auto middle_topology = BuildGhostTopology(three_rank.slices[1], 1, three_rank.slices.size());
    DEC3D_CHECK(middle_topology.is_initialized());
    DEC3D_CHECK(middle_topology.has_radial_inner_neighbor);
    DEC3D_CHECK(middle_topology.has_radial_outer_neighbor);
    DEC3D_CHECK(middle_topology.has_neighbors());

    const auto outer_topology = BuildGhostTopology(three_rank.slices[2], 2, three_rank.slices.size());
    DEC3D_CHECK(outer_topology.is_initialized());
    DEC3D_CHECK(outer_topology.has_radial_inner_neighbor);
    DEC3D_CHECK(!outer_topology.has_radial_outer_neighbor);

    RadialOwnershipSlice invalid_slice;
    const auto invalid_topology = BuildGhostTopology(invalid_slice, 99, three_rank.slices.size());
    DEC3D_CHECK(!invalid_topology.is_initialized());
    DEC3D_CHECK(!invalid_topology.has_neighbors());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
