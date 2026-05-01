#include "core/diagnostics/stage_contracts.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    using namespace dec3d::core;

    const CachedFieldMask electron_invalidated =
        InvalidatedCachesForWriteMask(Combine({AuthoritativeField::e_electron}));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::electron_temperature));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::ion_temperature));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::electron_pressure));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::ion_pressure));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::chi_e));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::conductivity));
    DEC3D_CHECK(MaskContains(electron_invalidated, CachedField::opacity));

    const CachedFieldMask momentum_invalidated =
        InvalidatedCachesForWriteMask(Combine({AuthoritativeField::mom_r}));
    DEC3D_CHECK(MaskContains(momentum_invalidated, CachedField::velocity));
    DEC3D_CHECK(MaskContains(momentum_invalidated, CachedField::sound_speed));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
