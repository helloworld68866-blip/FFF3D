#include "state/canonical_state/canonical_state.hpp"
#include "state/canonical_state/state_utilities.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ValidateCanonicalStateLayout;

    const CanonicalStateLayout no_radiation_layout{2, 3, 4, 0};
    const auto validation = ValidateCanonicalStateLayout(no_radiation_layout);
    DEC3D_CHECK(validation.success);

    const auto state = CanonicalState::Create(no_radiation_layout);
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::rho));
    DEC3D_CHECK(!state.HasAuthoritativeStorage(AuthoritativeField::radiation_groups));

    const auto invalid_validation = ValidateCanonicalStateLayout(CanonicalStateLayout{0, 1, 1, 0});
    DEC3D_CHECK(!invalid_validation.success);
    DEC3D_CHECK(!invalid_validation.failure_reason.empty());

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
