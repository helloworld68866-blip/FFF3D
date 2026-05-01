#include "core/diagnostics/stage_contracts.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::core::Combine;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    state.MarkCachedFieldValid(CachedField::electron_temperature);
    state.MarkCachedFieldValid(CachedField::ion_temperature);
    state.MarkCachedFieldValid(CachedField::electron_pressure);
    state.MarkCachedFieldValid(CachedField::ion_pressure);
    state.MarkCachedFieldValid(CachedField::chi_e);
    state.MarkCachedFieldValid(CachedField::conductivity);
    state.MarkCachedFieldValid(CachedField::opacity);
    state.MarkCachedFieldValid(CachedField::sound_speed);
    state.MarkCachedFieldValid(CachedField::velocity);

    state.ApplyAuthoritativeWrite(Combine({AuthoritativeField::e_electron}));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::electron_temperature));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::ion_temperature));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::electron_pressure));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::ion_pressure));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::chi_e));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::conductivity));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::opacity));

    state.MarkCachedFieldValid(CachedField::velocity);
    state.MarkCachedFieldValid(CachedField::sound_speed);
    state.ApplyAuthoritativeWrite(Combine({AuthoritativeField::mom_r}));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::velocity));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::sound_speed));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
