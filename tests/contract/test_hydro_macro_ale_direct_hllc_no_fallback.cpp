#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <iostream>
#include <string>

namespace {

[[nodiscard]] bool HasDiagnostic(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& code) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }
  return false;
}

}  // namespace

int main() {
  try {
    dec3d::hydro::StaticGridHydroOptions options{};
    options.use_macro_zoning = true;
    options.use_ppm_reconstruction = true;
    options.apply_radial_ale_flux_correction = true;
    options.use_macro_ale_direct_moving_face_hllc = true;

    DEC3D_CHECK(options.use_macro_ale_direct_moving_face_hllc);
    DEC3D_CHECK_EQ(
        dec3d::hydro::MacroAleHllcModeName(
            dec3d::hydro::MacroAleHllcMode::direct_moving_face_hllc),
        std::string("direct_moving_face_hllc"));
    DEC3D_CHECK_EQ(
        dec3d::hydro::MacroAleHllcRemapOrder(),
        std::string("none"));

    auto state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{4u, 4u, 4u, 0u});
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          state.rho(radial, theta, phi) = 1.0;
          state.mom_r(radial, theta, phi) = 0.01;
          state.mom_theta(radial, theta, phi) = 0.0;
          state.mom_phi(radial, theta, phi) = 0.0;
          state.e_fluid_total(radial, theta, phi) = 2.5;
          state.e_electron(radial, theta, phi) = 0.5;
        }
      }
    }

    auto view = dec3d::state::BuildHydroWorkView(state);
    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{4u, 4u, 4u, 0.1, 0.5});

    auto result = dec3d::hydro::AdvanceStaticGridHydro(
        view,
        geometry,
        1.0e-5,
        options,
        nullptr);

    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(
        result.failure_reason.find("direct moving-face HLLC") !=
        std::string::npos);
    DEC3D_CHECK(!HasDiagnostic(
        result.diagnostics,
        "p1.hydro.macro_ale.compatibility.executed"));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
