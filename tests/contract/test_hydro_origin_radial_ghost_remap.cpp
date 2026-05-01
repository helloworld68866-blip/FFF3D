#include "core/array/array3d.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "test_assert.hpp"

#include <iostream>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState MakeState(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double seed =
      static_cast<double>((radial + 1u) * 100u + (theta + 1u) * 10u + phi);
  return dec3d::hydro::MakeConservativeState({
      1.0 + (1.0e-3 * seed),
      0.10 + (1.0e-4 * seed),
      -0.20 + (1.0e-4 * seed),
      0.30 + (1.0e-4 * seed),
      10.0 + (1.0e-3 * seed),
      0.40 + (1.0e-4 * seed)});
}

[[nodiscard]] dec3d::hydro::HydroConservativeState BuildExpectedOriginGhost(
    const dec3d::hydro::HydroConservativeState& mapped) noexcept {
  return {
      mapped.rho,
      -mapped.mom_r,
      mapped.mom_theta,
      -mapped.mom_phi,
      mapped.e_fluid_total,
      mapped.chi_e};
}

void CheckStateEquals(
    const dec3d::hydro::HydroConservativeState& actual,
    const dec3d::hydro::HydroConservativeState& expected) {
  DEC3D_CHECK_EQ(actual.rho, expected.rho);
  DEC3D_CHECK_EQ(actual.mom_r, expected.mom_r);
  DEC3D_CHECK_EQ(actual.mom_theta, expected.mom_theta);
  DEC3D_CHECK_EQ(actual.mom_phi, expected.mom_phi);
  DEC3D_CHECK_EQ(actual.e_fluid_total, expected.e_fluid_total);
  DEC3D_CHECK_EQ(actual.chi_e, expected.chi_e);
}

}  // namespace

int main() {
  try {
    using dec3d::core::Array3D;
    using dec3d::hydro::FillHydroRadialGhosts;
    using dec3d::hydro::HydroDirection;
    using dec3d::hydro::HydroConservativeState;

    constexpr std::size_t kGhostLayers = 3u;
    constexpr std::size_t kRadialCells = 4u;
    constexpr std::size_t kThetaCells = 4u;
    constexpr std::size_t kPhiCells = 6u;

    Array3D<HydroConservativeState> interior(kRadialCells, kThetaCells, kPhiCells);
    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          interior(radial, theta, phi) = MakeState(radial, theta, phi);
        }
      }
    }

    const auto prepared = FillHydroRadialGhosts(interior, kGhostLayers);
    DEC3D_CHECK(prepared.is_complete(kRadialCells, kThetaCells, kPhiCells, kGhostLayers));
    DEC3D_CHECK(prepared.direction == HydroDirection::radial);

    for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
      for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
        const std::size_t mapped_theta = kThetaCells - 1u - theta;
        const std::size_t mapped_phi = (phi + (kPhiCells / 2u)) % kPhiCells;
        for (std::size_t ghost = 0; ghost < kGhostLayers; ++ghost) {
          const std::size_t interior_radial = kGhostLayers - 1u - ghost;
          const auto expected =
              BuildExpectedOriginGhost(interior(interior_radial, mapped_theta, mapped_phi));
          CheckStateEquals(prepared.ghosted_states(ghost, theta, phi), expected);
        }
      }
    }

    Array3D<HydroConservativeState> odd_phi(kRadialCells, kThetaCells, 5u);
    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < odd_phi.extent_phi(); ++phi) {
          odd_phi(radial, theta, phi) = MakeState(radial, theta, phi);
        }
      }
    }
    const auto invalid_prepared = FillHydroRadialGhosts(odd_phi, kGhostLayers);
    DEC3D_CHECK(!invalid_prepared.success);
    DEC3D_CHECK(!invalid_prepared.failure_reason.empty());

    Array3D<HydroConservativeState> too_few_radial_cells(2u, kThetaCells, kPhiCells);
    for (std::size_t radial = 0; radial < too_few_radial_cells.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          too_few_radial_cells(radial, theta, phi) = MakeState(radial, theta, phi);
        }
      }
    }
    const auto underresolved_prepared =
        FillHydroRadialGhosts(too_few_radial_cells, kGhostLayers);
    DEC3D_CHECK(!underresolved_prepared.success);
    DEC3D_CHECK(!underresolved_prepared.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
