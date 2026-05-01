#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    const auto state = CanonicalState::Create(CanonicalStateLayout{2, 3, 4, 2});

    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::rho));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::mom_r));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::mom_theta));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::mom_phi));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::e_fluid_total));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::e_electron));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::radiation_groups));
    DEC3D_CHECK(state.HasAuthoritativeStorage(AuthoritativeField::alpha_state));

    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::electron_temperature));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::ion_temperature));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::electron_pressure));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::ion_pressure));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::chi_e));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::conductivity));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::opacity));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::sound_speed));
    DEC3D_CHECK(state.HasCachedFieldMetadata(CachedField::velocity));

    DEC3D_CHECK(!state.HasCachedFieldStorage(CachedField::electron_temperature));
    DEC3D_CHECK(!state.HasCachedFieldStorage(CachedField::chi_e));
    DEC3D_CHECK(!state.HasCachedFieldStorage(CachedField::opacity));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
