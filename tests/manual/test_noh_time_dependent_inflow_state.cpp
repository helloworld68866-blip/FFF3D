#define main dec3d_noh_base_main_for_inflow_state_test
#include "benchmark_noh_ppm_nr256_mpi24.cpp"
#undef main

#include "test_assert.hpp"

#include <cmath>

int main() {
  constexpr double kRho0 = 1.0;
  constexpr double kPressure0 = 1.0e-6;
  constexpr double kRadialVelocity0 = -1.0;
  constexpr double kElectronEnergyFraction = 0.5;
  constexpr double kTime = 0.4;
  constexpr double kRadius = 1.0;

  const auto state = NohExactpackUpstreamInflowGhostState(
      kRho0,
      kPressure0,
      kRadialVelocity0,
      kElectronEnergyFraction,
      kTime,
      kRadius);

  const double expected_rho = kRho0 * std::pow(1.0 + kTime / kRadius, 2.0);
  const double expected_e_total =
      kPressure0 / (dec3d::state::HydroIdealGasGamma() - 1.0) +
      0.5 * expected_rho * kRadialVelocity0 * kRadialVelocity0;
  const auto primitive = dec3d::hydro::RecoverPrimitiveState(state);

  DEC3D_CHECK(std::abs(state.rho - expected_rho) <= 1.0e-12);
  DEC3D_CHECK(std::abs(state.mom_r - expected_rho * kRadialVelocity0) <= 1.0e-12);
  DEC3D_CHECK(std::abs(state.mom_theta) <= 1.0e-15);
  DEC3D_CHECK(std::abs(state.mom_phi) <= 1.0e-15);
  DEC3D_CHECK(std::abs(state.e_fluid_total - expected_e_total) <= 1.0e-12);
  DEC3D_CHECK(primitive.is_physical());
  DEC3D_CHECK(std::abs(primitive.pressure - kPressure0) <= 1.0e-12);
  DEC3D_CHECK(std::abs(primitive.v_r - kRadialVelocity0) <= 1.0e-12);
  return 0;
}
