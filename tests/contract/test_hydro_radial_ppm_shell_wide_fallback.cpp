#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <charconv>
#include <cmath>
#include <iostream>
#include <string>
#include <string_view>

namespace {

void SeedAngularlyMixedRadialProfile(dec3d::state::CanonicalState& state) {
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const double radius_fraction =
        static_cast<double>(radial) /
        static_cast<double>(state.layout.radial_cells - 1u);
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        double rho = 1.0 + (0.2 * radius_fraction);
        if (theta == 0u && radial == 5u) {
          rho = 1.8;
        }
        if (theta == 0u && radial == 6u) {
          rho = 0.95;
        }

        constexpr double kRadialVelocity = -0.2;
        constexpr double kPressure = 0.1;
        const auto conservative = dec3d::hydro::MakeConservativeState(
            dec3d::hydro::HydroPrimitiveState{
                rho,
                kRadialVelocity,
                0.0,
                0.0,
                kPressure,
                std::pow(0.04, 3.0 / 5.0)});
        state.rho(radial, theta, phi) = conservative.rho;
        state.mom_r(radial, theta, phi) = conservative.mom_r;
        state.mom_theta(radial, theta, phi) = conservative.mom_theta;
        state.mom_phi(radial, theta, phi) = conservative.mom_phi;
        state.e_fluid_total(radial, theta, phi) =
            conservative.e_fluid_total;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(0.04);
      }
    }
  }
}

[[nodiscard]] bool TryParseSizeField(
    std::string_view message,
    std::string_view key,
    std::size_t* value) noexcept {
  const std::string token = std::string(key) + "=";
  const std::size_t begin = message.find(token);
  if (begin == std::string_view::npos) {
    return false;
  }

  const std::size_t value_begin = begin + token.size();
  std::size_t value_end = message.find(';', value_begin);
  if (value_end == std::string_view::npos) {
    value_end = message.size();
  }
  while (value_end > value_begin && message[value_end - 1u] == ' ') {
    --value_end;
  }

  const auto field = message.substr(value_begin, value_end - value_begin);
  const auto* first = field.data();
  const auto* last = first + field.size();
  const auto parsed = std::from_chars(first, last, *value);
  return parsed.ec == std::errc{} && parsed.ptr == last;
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::RadialGhostOverride;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    constexpr std::size_t kRadialCells = 12u;
    constexpr std::size_t kThetaCells = 2u;
    constexpr std::size_t kPhiCells = 2u;
    constexpr double kDt = 1.0e-4;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto state = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    SeedAngularlyMixedRadialProfile(state);

    const auto proposal = BuildRadialAleMeshUpdateProposal(
        state.rho,
        state.mom_r,
        geometry,
        kDt);
    DEC3D_CHECK(proposal.is_complete(geometry.radial_faces.size()));

    StaticGridHydroOptions options;
    options.apply_geometric_source = false;
    options.apply_theta_sweep = false;
    options.apply_phi_sweep = false;
    options.use_ppm_reconstruction = true;
    options.reconstruction_ghost_layers = 3u;
    options.apply_radial_ale_flux_correction = true;
    options.debug_angular_stage_diagnostics = true;

    auto hydro_view = BuildHydroWorkView(state);
    DEC3D_CHECK(hydro_view.is_complete());
    const auto result = AdvanceStaticGridHydro(
        hydro_view,
        geometry,
        kDt,
        RadialGhostOverride{},
        options,
        &proposal);
    if (!result.is_complete()) {
      std::cerr << result.failure_reason << '\n';
    }
    DEC3D_CHECK(result.is_complete());

    bool saw_face_debug = false;
    for (const auto& entry : result.diagnostics.entries) {
      if (entry.code != "p1.hydro.debug.radial_ale.face_flux") {
        continue;
      }

      saw_face_debug = true;
      std::size_t line_count = 0u;
      std::size_t downgraded_count = 0u;
      DEC3D_CHECK(TryParseSizeField(entry.message, "count", &line_count));
      DEC3D_CHECK(TryParseSizeField(
          entry.message,
          "downgraded_count",
          &downgraded_count));
      DEC3D_CHECK_EQ(line_count, kThetaCells * kPhiCells);
      DEC3D_CHECK_EQ(downgraded_count, line_count);
    }
    DEC3D_CHECK(saw_face_debug);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
