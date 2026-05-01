#include "alpha/alpha_operator.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>

namespace {

dec3d::mesh::SphericalGeometryMetadata MakePatchGeometry(std::size_t radial_cells) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  for (std::size_t r = 0u; r <= radial_cells; ++r) {
    geometry.radial_faces[r] = 0.1 + 0.1 * static_cast<double>(r);
  }
  geometry.theta_faces = {0.5, 1.2};
  geometry.phi_faces = {0.0, 2.0 * dec3d::physics::PhysicsConstantsCGS::pi};
  geometry.cell_volumes.assign(radial_cells, 0.0);
  geometry.global_volume = 0.0;
  const double polar_factor =
      std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1]);
  const double azimuthal_factor = geometry.phi_faces[1] - geometry.phi_faces[0];
  for (std::size_t r = 0u; r < radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) /
        3.0;
    geometry.cell_volumes[r] = radial_factor * polar_factor * azimuthal_factor;
    geometry.global_volume += geometry.cell_volumes[r];
  }
  return geometry;
}

dec3d::transport::GenericDiffusionBoundaryPolicy MakePatchBoundary() {
  dec3d::transport::GenericDiffusionBoundaryPolicy policy;
  policy.inner_radial = dec3d::transport::DiffusionBoundaryKind::interior_patch_no_origin;
  policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::neumann_zero_flux;
  policy.theta_lower = dec3d::transport::DiffusionBoundaryKind::interior_patch_no_pole;
  policy.theta_upper = dec3d::transport::DiffusionBoundaryKind::interior_patch_no_pole;
  policy.phi = dec3d::transport::DiffusionBoundaryKind::periodic;
  return policy;
}

dec3d::state::CanonicalState BuildAlphaOperatorState(
    std::size_t radial_cells,
    double alpha0,
    double alpha1 = 0.0,
    double rho = 10.0) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{radial_cells, 1u, 1u, 0u});
  const double Te = dec3d::physics::ErgFromKeV(4.0);
  const double Ti = dec3d::physics::ErgFromKeV(3.0);
  const double ni = rho / dec3d::physics::DefaultMeanDTIonMassG();
  const double ne = ni;
  constexpr double gamma_minus_one = 2.0 / 3.0;
  for (std::size_t r = 0u; r < radial_cells; ++r) {
    state.rho(r, 0u, 0u) = rho;
    state.mom_r(r, 0u, 0u) = 0.0;
    state.mom_theta(r, 0u, 0u) = 0.0;
    state.mom_phi(r, 0u, 0u) = 0.0;
    state.e_electron(r, 0u, 0u) = ne * Te / gamma_minus_one;
    state.e_fluid_total(r, 0u, 0u) = (ne * Te + ni * Ti) / gamma_minus_one;
    state.alpha_state.storage(r, 0u, 0u) =
        (r == 0u || radial_cells == 1u) ? alpha0 : alpha1;
  }
  return state;
}

dec3d::alpha::OneGroupAlphaTransportOptions MakeOptions(
    double dt_s,
    dec3d::alpha::AlphaReactivityModel model =
        dec3d::alpha::AlphaReactivityModel::bosch_hale_dt,
    double constant_reactivity = 1.0e-40) {
  dec3d::alpha::OneGroupAlphaTransportOptions options;
  options.dt_s = dt_s;
  options.boundary_policy = MakePatchBoundary();
  options.provider_options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.provider_options.tau_model =
      dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.provider_options.reactivity_model = model;
  options.provider_options.constant_dt_reactivity_cm3_s = constant_reactivity;
  options.provider_options.recovery_options.mean_ion_mass_g =
      dec3d::physics::DefaultMeanDTIonMassG();
  options.provider_options.recovery_options.zbar = dec3d::physics::DefaultZbar();
  options.provider_options.alpha_energy_floor_erg_cm3 = 0.0;
  options.electron_energy_floor_erg_per_cm3 = 0.0;
  options.ion_energy_floor_erg_per_cm3 = 0.0;
  options.alpha_energy_floor_erg_per_cm3 = 0.0;
  return options;
}

void CheckNear(double actual, double expected, double tolerance, const std::string& label) {
  if (std::abs(actual - expected) > tolerance) {
    std::ostringstream detail;
    detail << label << ": expected " << expected << ", got " << actual
           << ", tolerance " << tolerance;
    dec3d::test::Fail("near", __FILE__, __LINE__, detail.str());
  }
}

void CheckSuccess(const dec3d::alpha::OneGroupAlphaTransportResult& result) {
  if (!result.success) {
    dec3d::test::Fail("alpha transport success", __FILE__, __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

void CheckStateUnchanged(
    const dec3d::state::CanonicalState& state,
    const dec3d::core::Array3D<double>& alpha_before,
    const dec3d::core::Array3D<double>& electron_before,
    const dec3d::core::Array3D<double>& total_before) {
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    CheckNear(state.alpha_state.storage(r, 0u, 0u), alpha_before(r, 0u, 0u), 0.0,
              "alpha unchanged");
    CheckNear(state.e_electron(r, 0u, 0u), electron_before(r, 0u, 0u), 0.0,
              "electron unchanged");
    CheckNear(state.e_fluid_total(r, 0u, 0u), total_before(r, 0u, 0u), 0.0,
              "fluid total unchanged");
  }
}

}  // namespace

int RunTests() {
  {
    auto state = BuildAlphaOperatorState(1u, 3.0e9);
    const auto before_alpha = state.alpha_state.storage;
    const auto before_electron = state.e_electron;
    const auto before_total = state.e_fluid_total;
    const auto result = dec3d::alpha::ApplyOneGroupAlphaTransport(
        state, MakePatchGeometry(1u), MakeOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(dec3d::alpha::ValidateAlphaOneGroupOperatorDiagnostics(result));
    DEC3D_CHECK(result.report_line.find("diagnostic_id=p4.alpha.one_group_operator") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("assembly_report_present=true") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("min_epsilon_alpha_before=") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("max_epsilon_alpha_before=") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("min_epsilon_alpha_after=") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("max_epsilon_alpha_after=") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("alpha_equation_budget_residual=") !=
                std::string::npos);
    DEC3D_CHECK(result.alpha_equation_budget_residual <=
                1.0e-10 * std::max(1.0, std::abs(result.birth_source_total)));
    DEC3D_CHECK_EQ(result.updated_fields,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("global_budget_semantics=single_rank_whole_domain_volume_integral") !=
                std::string::npos);
    CheckStateUnchanged(state, before_alpha, before_electron, before_total);
  }

  {
    auto state = BuildAlphaOperatorState(1u, 2.0e9);
    const auto provider =
        dec3d::alpha::BuildAlphaCoefficientArrays(state, MakeOptions(1.0e-13).provider_options);
    DEC3D_CHECK(provider.success);
    const double tau = provider.coefficients.tau_alphae_s(0u, 0u, 0u);
    const double birth = provider.coefficients.birth_source_erg_cm3_s(0u, 0u, 0u);
    const double alpha_old = state.alpha_state.storage(0u, 0u, 0u);
    const double electron_old = state.e_electron(0u, 0u, 0u);
    const double total_old = state.e_fluid_total(0u, 0u, 0u);
    const double dt = 1.0e-13;
    const double alpha_expected = (alpha_old + dt * birth) / (1.0 + dt / tau);
    const double electron_expected = electron_old + dt * alpha_expected / tau;

    const auto result = dec3d::alpha::ApplyOneGroupAlphaTransport(
        state, MakePatchGeometry(1u), MakeOptions(dt));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::alpha::ValidateAlphaOneGroupOperatorDiagnostics(result));
    CheckNear(state.alpha_state.storage(0u, 0u, 0u),
              alpha_expected,
              1.0e-11 * std::max(1.0, std::abs(alpha_expected)),
              "one-cell alpha implicit source oracle");
    CheckNear(state.e_electron(0u, 0u, 0u),
              electron_expected,
              1.0e-11 * std::max(1.0, std::abs(electron_expected)),
              "one-cell electron deposition oracle");
    CheckNear(state.e_fluid_total(0u, 0u, 0u),
              total_old + (electron_expected - electron_old),
              1.0e-11 * std::max(1.0, std::abs(total_old)),
              "one-cell total fluid oracle");
    DEC3D_CHECK_EQ(result.updated_fields, dec3d::alpha::AlphaOperatorWriteMask());
    DEC3D_CHECK(result.report_line.find("reactivity_model_executed=bosch_hale_dt") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("thesis_reactivity_claim_allowed=true") !=
                std::string::npos);
  }

  {
    auto state = BuildAlphaOperatorState(2u, 0.0, 1.0, 1.0e-6);
    const double low_before = state.alpha_state.storage(0u, 0u, 0u);
    const double high_before = state.alpha_state.storage(1u, 0u, 0u);
    auto options = MakeOptions(
        1.0e-11,
        dec3d::alpha::AlphaReactivityModel::constant_user_supplied,
        1.0e-80);
    // This stress case intentionally uses low-density alpha coefficients to make
    // two-cell coupling visible; keep the tolerance local to the oracle.
    options.serial_reference_options.residual_tolerance = 1.0e-2;
    const auto result =
        dec3d::alpha::ApplyOneGroupAlphaTransport(state, MakePatchGeometry(2u), options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::alpha::ValidateAlphaOneGroupOperatorDiagnostics(result));
    DEC3D_CHECK(state.alpha_state.storage(0u, 0u, 0u) > low_before);
    DEC3D_CHECK(state.alpha_state.storage(1u, 0u, 0u) < high_before);
    DEC3D_CHECK(result.assembly_report.find("nonzero_count=") != std::string::npos);
  }

  {
    auto state = BuildAlphaOperatorState(1u, 2.0e9);
    const auto before_alpha = state.alpha_state.storage;
    const auto before_electron = state.e_electron;
    const auto before_total = state.e_fluid_total;
    auto options = MakeOptions(1.0e-13);
    options.provider_options.reactivity_model = dec3d::alpha::AlphaReactivityModel::missing;
    const auto result =
        dec3d::alpha::ApplyOneGroupAlphaTransport(state, MakePatchGeometry(1u), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.updated_fields,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);
    CheckStateUnchanged(state, before_alpha, before_electron, before_total);
  }

  {
    auto state = BuildAlphaOperatorState(1u, 2.0e9);
    const auto before_alpha = state.alpha_state.storage;
    const auto before_electron = state.e_electron;
    const auto before_total = state.e_fluid_total;
    auto options = MakeOptions(1.0e-13);
    options.electron_energy_floor_erg_per_cm3 = 1.0e99;
    const auto result =
        dec3d::alpha::ApplyOneGroupAlphaTransport(state, MakePatchGeometry(1u), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.updated_fields,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);
    CheckStateUnchanged(state, before_alpha, before_electron, before_total);
  }

  {
    auto state = BuildAlphaOperatorState(1u, 2.0e9);
    const auto before_alpha = state.alpha_state.storage;
    const auto before_electron = state.e_electron;
    const auto before_total = state.e_fluid_total;
    auto options = MakeOptions(1.0e-13);
    options.ion_energy_floor_erg_per_cm3 = 1.0e99;
    const auto result =
        dec3d::alpha::ApplyOneGroupAlphaTransport(state, MakePatchGeometry(1u), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.updated_fields,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);
    CheckStateUnchanged(state, before_alpha, before_electron, before_total);
  }

  {
    auto state = BuildAlphaOperatorState(1u, 2.0e9);
    const double rho_before = state.rho(0u, 0u, 0u);
    const double mom_r_before = state.mom_r(0u, 0u, 0u);
    const double mom_theta_before = state.mom_theta(0u, 0u, 0u);
    const double mom_phi_before = state.mom_phi(0u, 0u, 0u);
    const auto result = dec3d::alpha::ApplyOneGroupAlphaTransport(
        state, MakePatchGeometry(1u), MakeOptions(1.0e-13));
    CheckSuccess(result);
    DEC3D_CHECK_EQ(result.updated_fields, dec3d::alpha::AlphaOperatorWriteMask());
    CheckNear(state.rho(0u, 0u, 0u), rho_before, 0.0, "rho unchanged");
    CheckNear(state.mom_r(0u, 0u, 0u), mom_r_before, 0.0, "mom_r unchanged");
    CheckNear(state.mom_theta(0u, 0u, 0u), mom_theta_before, 0.0, "mom_theta unchanged");
    CheckNear(state.mom_phi(0u, 0u, 0u), mom_phi_before, 0.0, "mom_phi unchanged");
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask, dec3d::alpha::AlphaOperatorWriteMask());
  }

  {
    dec3d::alpha::OneGroupAlphaTransportResult missing;
    missing.success = true;
    DEC3D_CHECK(!dec3d::alpha::ValidateAlphaOneGroupOperatorDiagnostics(missing));
  }

  return 0;
}

int main() {
  try {
    return RunTests();
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
