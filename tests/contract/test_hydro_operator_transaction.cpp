#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
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
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    const auto geometry = BuildSphericalGeometry(SphericalMeshDescriptor{2, 4, 4, 1.0, 2.0});
    DEC3D_CHECK(geometry.is_valid());

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 4, 4, 0});
    const dec3d::hydro::HydroPrimitiveState left{
        1.0,
        0.0,
        0.30,
        -0.20,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const dec3d::hydro::HydroPrimitiveState right{
        0.125,
        0.0,
        -0.10,
        0.25,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};
    const auto left_state = dec3d::hydro::MakeConservativeState(left);
    const auto right_state = dec3d::hydro::MakeConservativeState(right);

    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        state.rho(0, theta, phi) = left_state.rho;
        state.mom_r(0, theta, phi) = left_state.mom_r;
        state.mom_theta(0, theta, phi) = left_state.mom_theta;
        state.mom_phi(0, theta, phi) = left_state.mom_phi;
        state.e_fluid_total(0, theta, phi) = left_state.e_fluid_total;
        state.e_electron(0, theta, phi) = ElectronEnergyDensityFromPressure(0.4);

        state.rho(1, theta, phi) = right_state.rho;
        state.mom_r(1, theta, phi) = right_state.mom_r;
        state.mom_theta(1, theta, phi) = right_state.mom_theta;
        state.mom_phi(1, theta, phi) = right_state.mom_phi;
        state.e_fluid_total(1, theta, phi) = right_state.e_fluid_total;
        state.e_electron(1, theta, phi) = ElectronEnergyDensityFromPressure(0.04);
      }
    }

    const dec3d::hydro::HydroPrimitiveState fragile{
        1.0,
        0.0,
        0.7,
        0.7,
        1.0e-6,
        std::pow(1.0e-6, 3.0 / 5.0)};
    const auto fragile_state = dec3d::hydro::MakeConservativeState(fragile);
    const auto last_theta = state.rho.extent_theta() - 1u;
    const auto last_phi = state.rho.extent_phi() - 1u;
    state.rho(1, last_theta, last_phi) = fragile_state.rho;
    state.mom_r(1, last_theta, last_phi) = fragile_state.mom_r;
    state.mom_theta(1, last_theta, last_phi) = fragile_state.mom_theta;
    state.mom_phi(1, last_theta, last_phi) = fragile_state.mom_phi;
    state.e_fluid_total(1, last_theta, last_phi) = fragile_state.e_fluid_total;
    state.e_electron(1, last_theta, last_phi) = ElectronEnergyDensityFromPressure(1.0e-6);

    const auto snapshot_rho = state.rho.storage();
    const auto snapshot_mom_r = state.mom_r.storage();
    const auto snapshot_mom_theta = state.mom_theta.storage();
    const auto snapshot_mom_phi = state.mom_phi.storage();
    const auto snapshot_e_fluid_total = state.e_fluid_total.storage();
    const auto snapshot_e_electron = state.e_electron.storage();

    const StageContext context{
        0.0,
        1.0,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1",
        "ownership:rank0",
        "diagnostics:p1.hydro"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, state));

    const auto result = hydro.advance();
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(!result.is_semantically_complete());
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(!result.failure_reason.empty());
    DEC3D_CHECK(result.diagnostics.has_entries());
    DEC3D_CHECK(result.execution_evidence.has_value());
    DEC3D_CHECK(result.execution_evidence->entered_stage);
    DEC3D_CHECK(HasDiagnosticCode(result.diagnostics, "p1.hydro.stage.failed"));
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
