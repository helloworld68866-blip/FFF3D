#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>

namespace {

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
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
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    constexpr std::size_t kRadialCells = 8u;
    constexpr std::size_t kThetaCells = 4u;
    constexpr std::size_t kPhiCells = 4u;

    auto state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      const bool left = radial < kRadialCells / 2u;
      const double rho = left ? 1.0 : 0.125;
      const double pressure = left ? 1.0 : 0.1;
      const double electron_pressure = left ? 0.4 : 0.04;
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          state.rho(radial, theta, phi) = rho;
          state.mom_r(radial, theta, phi) = 0.0;
          state.mom_theta(radial, theta, phi) = 0.0;
          state.mom_phi(radial, theta, phi) = 0.0;
          state.e_fluid_total(radial, theta, phi) =
              pressure / (dec3d::state::HydroIdealGasGamma() - 1.0);
          state.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(electron_pressure);
        }
      }
    }

    HydroOperator hydro;
    hydro.SetStaticGridOptions(StaticGridHydroOptions{
        true,
        true,
        true,
        true,
        true,
        3u});

    const StageContext context{
        0.0,
        5.0e-4,
        1,
        PhaseId::p1,
        "p1-v0.1-ppm",
        "mesh:p1.ppm",
        "ownership:rank0",
        "diagnostics:p1.hydro.ppm"};
    DEC3D_CHECK(hydro.bind(context, geometry, state));

    const auto result = hydro.advance();
    DEC3D_CHECK(result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.reconstruction.ppm.direction.radial"));
    DEC3D_CHECK(result.execution_evidence.has_value());
    DEC3D_CHECK(result.execution_evidence->implementation_id == "p1.hydro.operator.ppm_static_grid");

    bool saw_full_ppm_message = false;
    bool saw_narrow_ppm_message = false;
    for (const auto& entry : result.diagnostics.entries) {
      if (entry.code == "p1.hydro.stage.ppm_reconstruction") {
        if (entry.message.find("characteristic") != std::string::npos &&
            entry.message.find("traced-interface") != std::string::npos) {
          saw_full_ppm_message = true;
        }
        if (entry.message.find("narrow") != std::string::npos) {
          saw_narrow_ppm_message = true;
        }
      }
    }
    DEC3D_CHECK(saw_full_ppm_message);
    DEC3D_CHECK(!saw_narrow_ppm_message);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
