#include "core/diagnostics/stage_contracts.hpp"
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
    using dec3d::core::AuthoritativeField;
    using dec3d::core::CachedField;
    using dec3d::state::BuildHydroStateView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::CommitHydroWriteback;
    using dec3d::state::BuildHydroElectronWritebackMask;
    using dec3d::state::ElectronEnergyDensityFromPressure;
    using dec3d::state::ElectronPressureFromChiE;
    using dec3d::state::HydroAlphaViewOptions;
    using dec3d::state::HydroIdealGasGamma;

    auto state = CanonicalState::Create(CanonicalStateLayout{2, 2, 2, 1});
    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          state.rho(radial, theta, phi) = 1.0;
          state.mom_r(radial, theta, phi) = 0.0;
          state.mom_theta(radial, theta, phi) = 0.0;
          state.mom_phi(radial, theta, phi) = 0.0;
          state.e_fluid_total(radial, theta, phi) = 20.0;
        }
      }
    }
    state.e_electron(0, 0, 0) = 6.0;
    state.e_electron(0, 0, 1) = 9.0;
    state.radiation_groups[0](0, 0, 0) = 77.0;
    state.alpha_state.storage(0, 0, 0) = 88.0;

    const auto initial_radiation = state.radiation_groups[0](0, 0, 0);
    const auto initial_alpha = state.alpha_state.storage(0, 0, 0);

    auto hydro_view = BuildHydroStateView(state);
    DEC3D_CHECK(hydro_view.is_complete());
    DEC3D_CHECK(hydro_view.operator_local_chi_e);
    DEC3D_CHECK(!state.HasCachedFieldStorage(CachedField::chi_e));
    DEC3D_CHECK(!state.IsCachedFieldValid(CachedField::chi_e));
    DEC3D_CHECK_EQ(hydro_view.chi_e.size(), state.rho.size());
    DEC3D_CHECK_EQ(hydro_view.radiation_chi.size(), 1u);

    const double expected_pressure = (HydroIdealGasGamma() - 1.0) * state.e_electron(0, 0, 0);
    const double expected_chi_e = std::pow(expected_pressure, 3.0 / 5.0);
    DEC3D_CHECK(std::abs(hydro_view.chi_e(0, 0, 0) - expected_chi_e) < 1.0e-12);
    const double expected_chi_rad = std::pow(state.radiation_groups[0](0, 0, 0) / 3.0, 3.0 / 4.0);
    DEC3D_CHECK(std::abs(hydro_view.radiation_chi[0](0, 0, 0) - expected_chi_rad) < 1.0e-12);

    const double updated_chi_e = 3.5;
    hydro_view.chi_e(0, 0, 0) = updated_chi_e;
    const auto writeback = CommitHydroWriteback(state, hydro_view);
    DEC3D_CHECK(writeback.success);
    DEC3D_CHECK(writeback.electron_channel_updated);
    DEC3D_CHECK_EQ(writeback.updated_fields, BuildHydroElectronWritebackMask());
    DEC3D_CHECK(HasDiagnosticCode(writeback.diagnostics, "p1.hydro.writeback.electron_channel"));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask, BuildHydroElectronWritebackMask());
    DEC3D_CHECK(!dec3d::core::MaskContains(writeback.updated_fields, AuthoritativeField::radiation_groups));
    DEC3D_CHECK(!dec3d::core::MaskContains(writeback.updated_fields, AuthoritativeField::alpha_state));
    DEC3D_CHECK(std::abs(
        state.e_electron(0, 0, 0) -
        ElectronEnergyDensityFromPressure(ElectronPressureFromChiE(updated_chi_e))) < 1.0e-12);
    DEC3D_CHECK_EQ(state.radiation_groups[0](0, 0, 0), initial_radiation);
    DEC3D_CHECK_EQ(state.alpha_state.storage(0, 0, 0), initial_alpha);
    DEC3D_CHECK(dec3d::core::MaskContains(state.last_invalidated_cached_mask, CachedField::electron_temperature));
    DEC3D_CHECK(dec3d::core::MaskContains(state.last_invalidated_cached_mask, CachedField::chi_e));

    auto radiation_state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 1});
    radiation_state.rho(0, 0, 0) = 1.0;
    radiation_state.mom_r(0, 0, 0) = 0.0;
    radiation_state.mom_theta(0, 0, 0) = 0.0;
    radiation_state.mom_phi(0, 0, 0) = 0.0;
    radiation_state.e_fluid_total(0, 0, 0) = 20.0;
    radiation_state.e_electron(0, 0, 0) = 2.0;
    radiation_state.radiation_groups[0](0, 0, 0) = 12.0;
    auto radiation_view = BuildHydroStateView(radiation_state);
    DEC3D_CHECK(radiation_view.is_complete());
    DEC3D_CHECK_EQ(radiation_view.radiation_chi.size(), 1u);
    radiation_view.radiation_chi[0](0, 0, 0) = std::pow(81.0 / 3.0, 3.0 / 4.0);
    const auto radiation_writeback = CommitHydroWriteback(
        radiation_state,
        radiation_view,
        dec3d::core::ToMask(AuthoritativeField::radiation_groups));
    DEC3D_CHECK(radiation_writeback.success);
    DEC3D_CHECK(HasDiagnosticCode(
        radiation_writeback.diagnostics,
        "p3.radiation.hydro_terms.writeback"));
    DEC3D_CHECK(dec3d::core::MaskContains(
        radiation_writeback.updated_fields,
        AuthoritativeField::radiation_groups));
    DEC3D_CHECK_EQ(
        radiation_state.last_authoritative_write_mask,
        dec3d::core::ToMask(AuthoritativeField::radiation_groups));
    DEC3D_CHECK(std::abs(radiation_state.radiation_groups[0](0, 0, 0) - 81.0) < 1.0e-11);
    DEC3D_CHECK_EQ(radiation_state.e_electron(0, 0, 0), 2.0);

    auto alpha_state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 0});
    alpha_state.rho(0, 0, 0) = 1.0;
    alpha_state.mom_r(0, 0, 0) = 0.0;
    alpha_state.mom_theta(0, 0, 0) = 0.0;
    alpha_state.mom_phi(0, 0, 0) = 0.0;
    alpha_state.e_fluid_total(0, 0, 0) = 20.0;
    alpha_state.e_electron(0, 0, 0) = 2.0;
    alpha_state.alpha_state.storage(0, 0, 0) = 150.0;

    HydroAlphaViewOptions alpha_options;
    alpha_options.enabled = true;
    auto alpha_view = BuildHydroStateView(alpha_state, alpha_options);
    DEC3D_CHECK(alpha_view.is_complete());
    DEC3D_CHECK(alpha_view.operator_local_alpha_chi);
    const double expected_alpha_chi = std::pow((2.0 / 3.0) * 150.0, 3.0 / 5.0);
    DEC3D_CHECK(std::abs(alpha_view.alpha_chi(0, 0, 0) - expected_alpha_chi) < 1.0e-12);
    alpha_view.alpha_chi(0, 0, 0) = std::pow((2.0 / 3.0) * 300.0, 3.0 / 5.0);
    const auto alpha_writeback = CommitHydroWriteback(
        alpha_state,
        alpha_view,
        dec3d::core::ToMask(AuthoritativeField::alpha_state));
    DEC3D_CHECK(alpha_writeback.success);
    DEC3D_CHECK(HasDiagnosticCode(alpha_writeback.diagnostics, "p4.alpha.hydro_terms.writeback"));
    DEC3D_CHECK(dec3d::core::MaskContains(alpha_writeback.updated_fields, AuthoritativeField::alpha_state));
    DEC3D_CHECK_EQ(alpha_state.last_authoritative_write_mask,
                   dec3d::core::ToMask(AuthoritativeField::alpha_state));
    DEC3D_CHECK(std::abs(alpha_state.alpha_state.storage(0, 0, 0) - 300.0) < 1.0e-10);
    DEC3D_CHECK_EQ(alpha_state.e_electron(0, 0, 0), 2.0);

    auto disabled_alpha_state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 0});
    disabled_alpha_state.rho(0, 0, 0) = 1.0;
    disabled_alpha_state.mom_r(0, 0, 0) = 0.0;
    disabled_alpha_state.mom_theta(0, 0, 0) = 0.0;
    disabled_alpha_state.mom_phi(0, 0, 0) = 0.0;
    disabled_alpha_state.e_fluid_total(0, 0, 0) = 20.0;
    disabled_alpha_state.e_electron(0, 0, 0) = 2.0;
    disabled_alpha_state.alpha_state.storage(0, 0, 0) = -1.0;
    auto disabled_view = BuildHydroStateView(disabled_alpha_state);
    DEC3D_CHECK(disabled_view.is_complete());
    DEC3D_CHECK(disabled_view.alpha_chi.empty());
    DEC3D_CHECK(!disabled_view.operator_local_alpha_chi);

    auto infeasible_state = CanonicalState::Create(CanonicalStateLayout{1, 1, 1, 0});
    infeasible_state.rho(0, 0, 0) = 1.0;
    infeasible_state.mom_r(0, 0, 0) = 0.0;
    infeasible_state.mom_theta(0, 0, 0) = 0.0;
    infeasible_state.mom_phi(0, 0, 0) = 0.0;
    infeasible_state.e_fluid_total(0, 0, 0) = 1.0;
    infeasible_state.e_electron(0, 0, 0) = 0.1;

    auto infeasible_view = BuildHydroStateView(infeasible_state);
    DEC3D_CHECK(infeasible_view.is_complete());
    infeasible_view.chi_e(0, 0, 0) = std::pow(1.0, 3.0 / 5.0);

    const auto original_e_electron = infeasible_state.e_electron(0, 0, 0);
    const auto original_write_mask = infeasible_state.last_authoritative_write_mask;
    const auto infeasible_writeback = CommitHydroWriteback(infeasible_state, infeasible_view);
    DEC3D_CHECK(!infeasible_writeback.success);
    DEC3D_CHECK(HasDiagnosticCode(infeasible_writeback.diagnostics, "p1.hydro.writeback.failed"));
    DEC3D_CHECK_EQ(infeasible_state.e_electron(0, 0, 0), original_e_electron);
    DEC3D_CHECK_EQ(infeasible_state.last_authoritative_write_mask, original_write_mask);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
