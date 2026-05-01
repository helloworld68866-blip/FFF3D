#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ThermodynamicFeasibilityThresholds;
    using dec3d::state::ValidateThermodynamicFeasibility;

    auto state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 1});
    state.rho(0, 0, 0) = 1.0;
    state.mom_r(0, 0, 0) = 0.0;
    state.mom_theta(0, 0, 0) = 0.0;
    state.mom_phi(0, 0, 0) = 0.0;
    state.e_fluid_total(0, 0, 0) = 2.0;
    state.e_electron(0, 0, 0) = 3.0;

    const auto result = ValidateThermodynamicFeasibility(
        state,
        ThermodynamicFeasibilityThresholds{0.0, 0.0});

    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(!result.failure_reason.empty());
    DEC3D_CHECK(result.minimum_ion_energy < 0.0);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
