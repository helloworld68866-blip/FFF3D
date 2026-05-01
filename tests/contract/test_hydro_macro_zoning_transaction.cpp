#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace {

[[nodiscard]] bool SameStorage(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) noexcept {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    if (lhs[index] != rhs[index]) {
      return false;
    }
  }
  return true;
}

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

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{2, 6, 8, 0.05, 0.25});
    DEC3D_CHECK(geometry.is_valid());

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 6, 8, 0});
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          const double density = 1.0 + 0.05 * static_cast<double>(radial);
          state.rho(radial, theta, phi) = density;
          state.mom_r(radial, theta, phi) = 0.01 * density;
          state.mom_theta(radial, theta, phi) = 0.02 * std::sin(static_cast<double>(theta + 1u));
          state.mom_phi(radial, theta, phi) = -0.01 * std::cos(static_cast<double>(phi + 1u));
          state.e_fluid_total(radial, theta, phi) = 2.5 + 0.1 * density;
          state.e_electron(radial, theta, phi) = 0.4;
        }
      }
    }

    const auto snapshot_rho = state.rho.storage();
    const auto snapshot_mom_r = state.mom_r.storage();
    const auto snapshot_mom_theta = state.mom_theta.storage();
    const auto snapshot_mom_phi = state.mom_phi.storage();
    const auto snapshot_e_fluid_total = state.e_fluid_total.storage();
    const auto snapshot_e_electron = state.e_electron.storage();

    const StageContext context{
        0.0,
        1.0e-3,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.macro-zoning",
        "ownership:rank0",
        "diagnostics:p1.hydro"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, state));

    StaticGridHydroOptions options{};
    options.apply_geometric_source = false;
    options.apply_radial_sweep = false;
    options.apply_theta_sweep = true;
    options.apply_phi_sweep = true;
    options.use_macro_zoning = true;
    hydro.SetStaticGridOptions(options);

    const auto result = hydro.advance();
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(!result.failure_reason.empty());
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.stage.static_grid_failed"));
    DEC3D_CHECK(SameStorage(state.rho.storage(), snapshot_rho));
    DEC3D_CHECK(SameStorage(state.mom_r.storage(), snapshot_mom_r));
    DEC3D_CHECK(SameStorage(state.mom_theta.storage(), snapshot_mom_theta));
    DEC3D_CHECK(SameStorage(state.mom_phi.storage(), snapshot_mom_phi));
    DEC3D_CHECK(SameStorage(state.e_fluid_total.storage(), snapshot_e_fluid_total));
    DEC3D_CHECK(SameStorage(state.e_electron.storage(), snapshot_e_electron));

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
