#include "mesh/boundary/spherical_scalar_remap.hpp"
#include "test_assert.hpp"

#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

namespace {

dec3d::mesh::ScalarRemapLayout Layout(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) {
  return dec3d::mesh::ScalarRemapLayout{radial, theta, phi};
}

void CheckIndex(
    const dec3d::mesh::ScalarRemapIndex& actual,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    const std::string& label) {
  DEC3D_CHECK_EQ(actual.radial, radial);
  DEC3D_CHECK_EQ(actual.theta, theta);
  DEC3D_CHECK_EQ(actual.phi, phi);
  if (actual.radial != radial || actual.theta != theta || actual.phi != phi) {
    dec3d::test::Fail("scalar remap index", __FILE__, __LINE__, label);
  }
}

}  // namespace

int main() {
  try {
    using dec3d::mesh::MapPhiHalfTurn;
    using dec3d::mesh::MapScalarLowerPoleNeighbor;
    using dec3d::mesh::MapScalarOriginNeighbor;
    using dec3d::mesh::MapScalarUpperPoleNeighbor;
    using dec3d::mesh::ValidateScalarHalfTurnTopology;

    {
      const auto validation = ValidateScalarHalfTurnTopology(Layout(4, 4, 8));
      DEC3D_CHECK(validation.success);
      DEC3D_CHECK(validation.phi_half_turn_available);
      DEC3D_CHECK(validation.report_line.find("scalar_remap_topology_valid=true") !=
                  std::string::npos);
      DEC3D_CHECK_EQ(MapPhiHalfTurn(1, 8), std::size_t{5});
      DEC3D_CHECK_EQ(MapPhiHalfTurn(7, 8), std::size_t{3});
    }

    {
      const auto validation = ValidateScalarHalfTurnTopology(Layout(4, 4, 7));
      DEC3D_CHECK(!validation.success);
      DEC3D_CHECK(!validation.phi_half_turn_available);
      DEC3D_CHECK(validation.report_line.find("scalar_remap_topology_valid=false") !=
                  std::string::npos);
      DEC3D_CHECK(validation.failure_reason.find("even phi") != std::string::npos);
    }

    {
      const auto mapped = MapScalarOriginNeighbor(1, 1, 2, Layout(4, 4, 8));
      DEC3D_CHECK(mapped.success);
      CheckIndex(mapped.mapped_index, 0, 2, 6, "origin maps to antipodal first shell");
      DEC3D_CHECK(mapped.report_line.find("remap_kind=origin") != std::string::npos);
      DEC3D_CHECK(mapped.report_line.find("value_transform=identity") != std::string::npos);
    }

    {
      const auto mapped = MapScalarLowerPoleNeighbor(1, 2, 3, Layout(4, 4, 8));
      DEC3D_CHECK(mapped.success);
      CheckIndex(mapped.mapped_index, 2, 0, 7, "lower pole maps same band half-turn phi");
      DEC3D_CHECK(mapped.report_line.find("remap_kind=lower_pole") != std::string::npos);
      DEC3D_CHECK(mapped.report_line.find("value_transform=identity") != std::string::npos);
    }

    {
      const auto mapped = MapScalarUpperPoleNeighbor(1, 2, 3, Layout(4, 4, 8));
      DEC3D_CHECK(mapped.success);
      CheckIndex(mapped.mapped_index, 2, 3, 7, "upper pole maps same band half-turn phi");
      DEC3D_CHECK(mapped.report_line.find("remap_kind=upper_pole") != std::string::npos);
      DEC3D_CHECK(mapped.report_line.find("value_transform=identity") != std::string::npos);
    }

    {
      const auto mapped = MapScalarOriginNeighbor(0, 1, 2, Layout(4, 4, 8));
      DEC3D_CHECK(!mapped.success);
      DEC3D_CHECK(mapped.failure_reason.find("positive ghost layer") != std::string::npos);
      DEC3D_CHECK(mapped.report_line.find("canonical_state_mutated=false") !=
                  std::string::npos);
    }

    std::cout << "spherical scalar remap contract passed\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
