#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::CachedField;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    state.rho(0, 0, 0) = 5.0;
    state.MarkCachedFieldValid(CachedField::chi_e);

    DEC3D_CHECK(!state.TryCommitCachedField(CachedField::chi_e));
    DEC3D_CHECK_EQ(state.rho(0, 0, 0), 5.0);
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::chi_e));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
