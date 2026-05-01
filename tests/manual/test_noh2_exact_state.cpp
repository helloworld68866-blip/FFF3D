#define main dec3d_noh_base_main_for_noh2_exact_state_test
#include "benchmark_noh_ppm_nr256_mpi24.cpp"
#undef main

#include "test_assert.hpp"

#include <cmath>

int main() {
  constexpr double kRho0 = 1.0;
  constexpr double kSpecificInternalEnergy0 = 1.0;
  constexpr double kElectronEnergyFraction = 0.5;
  constexpr double kTime = 0.2;
  constexpr double kRadius = 0.75;

  const auto state = Noh2ExactState(
      kRho0,
      kSpecificInternalEnergy0,
      kElectronEnergyFraction,
      kTime,
      kRadius);

  constexpr double kGamma = 5.0 / 3.0;
  const double compression = 1.0 / (1.0 - kTime);
  const double expected_rho = kRho0 * std::pow(compression, 3.0);
  const double expected_velocity = -kRadius * compression;
  const double expected_specific_internal_energy =
      kSpecificInternalEnergy0 * std::pow(compression, 2.0);
  const double expected_pressure =
      expected_rho * expected_specific_internal_energy * (kGamma - 1.0);
  const double expected_e_total =
      expected_rho * expected_specific_internal_energy +
      0.5 * expected_rho * expected_velocity * expected_velocity;

  const auto primitive = dec3d::hydro::RecoverPrimitiveState(state);
  DEC3D_CHECK(std::abs(state.rho - expected_rho) <= 1.0e-12);
  DEC3D_CHECK(std::abs(state.mom_r - expected_rho * expected_velocity) <= 1.0e-12);
  DEC3D_CHECK(std::abs(state.mom_theta) <= 1.0e-15);
  DEC3D_CHECK(std::abs(state.mom_phi) <= 1.0e-15);
  DEC3D_CHECK(std::abs(state.e_fluid_total - expected_e_total) <= 1.0e-12);
  DEC3D_CHECK(primitive.is_physical());
  DEC3D_CHECK(std::abs(primitive.pressure - expected_pressure) <= 1.0e-12);
  DEC3D_CHECK(std::abs(primitive.v_r - expected_velocity) <= 1.0e-12);
  return 0;
}
