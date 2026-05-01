#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "radiation/transport/one_group_gray_diffusion.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

namespace {

dec3d::mesh::SphericalGeometryMetadata MakeOneCellPatchGeometry() {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces = {0.1, 0.3};
  geometry.theta_faces = {0.5, 1.2};
  geometry.phi_faces = {0.0, 2.0 * std::acos(-1.0)};
  const double radial_factor =
      (std::pow(geometry.radial_faces[1], 3) -
       std::pow(geometry.radial_faces[0], 3)) /
      3.0;
  const double polar_factor =
      std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1]);
  const double azimuthal_factor = geometry.phi_faces[1] - geometry.phi_faces[0];
  geometry.cell_volumes = {radial_factor * polar_factor * azimuthal_factor};
  geometry.global_volume = geometry.cell_volumes[0];
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

dec3d::state::CanonicalState MakeTwoGroupOneCellState(
    double u0,
    double u1,
    double e_electron,
    double e_ion) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 2u});
  state.rho(0u, 0u, 0u) = 1.0;
  state.mom_r(0u, 0u, 0u) = 0.0;
  state.mom_theta(0u, 0u, 0u) = 0.0;
  state.mom_phi(0u, 0u, 0u) = 0.0;
  state.radiation_groups[0](0u, 0u, 0u) = u0;
  state.radiation_groups[1](0u, 0u, 0u) = u1;
  state.e_electron(0u, 0u, 0u) = e_electron;
  state.e_fluid_total(0u, 0u, 0u) = e_electron + e_ion;
  return state;
}

dec3d::core::Array3D<double> OneCellArray(double value) {
  return dec3d::core::Array3D<double>(1u, 1u, 1u, value);
}

dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions MakeTwoGroupOptions(
    double dt_s) {
  dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions options;
  options.dt_s = dt_s;
  options.group_layout =
      dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
          {1.0e15, 2.0e15, 4.0e15});
  options.boundary_policy = MakePatchBoundary();
  options.radiation_boundary_model =
      dec3d::radiation::OneGroupGrayRadiationBoundaryModel::
          contract_zero_flux_or_scalar_remap;
  options.radiation_energy_floor_erg_per_cm3 = 0.0;
  options.electron_energy_floor_erg_per_cm3 = 0.0;
  options.ion_energy_floor_erg_per_cm3 = 0.0;
  options.coefficients.Dbar_cm2_per_s = {OneCellArray(0.0), OneCellArray(0.0)};
  options.coefficients.kappaP_cm_inv = {OneCellArray(1.0e-4), OneCellArray(2.0e-4)};
  options.coefficients.B_erg_per_cm3 = {OneCellArray(5.0), OneCellArray(7.0)};
  return options;
}

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    dec3d::test::Fail(label, __FILE__, __LINE__);
  }
}

double ExactImplicitSourceUpdate(double u_old, double dt_s, double kappa, double b) {
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double rate = dt_s * c * kappa;
  return (u_old + rate * b) / (1.0 + rate);
}

void CheckSuccess(
    const dec3d::radiation::MultigroupGrayRadiationMatterCouplingResult& result) {
  if (!result.success) {
    dec3d::test::Fail("multigroup radiation matter coupling success",
                      __FILE__,
                      __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

}  // namespace

int RunTests() {
  {
    const double dt_s = 1.0e-6;
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    auto options = MakeTwoGroupOptions(dt_s);

    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(result));

    const double u0 = ExactImplicitSourceUpdate(2.0, dt_s, 1.0e-4, 5.0);
    const double u1 = ExactImplicitSourceUpdate(3.0, dt_s, 2.0e-4, 7.0);
    const double expected_delta_e = -((u0 - 2.0) + (u1 - 3.0));
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), u0, 1.0e-12,
              "group 0 exact source oracle");
    CheckNear(state.radiation_groups[1](0u, 0u, 0u), u1, 1.0e-12,
              "group 1 exact source oracle");
    CheckNear(state.e_electron(0u, 0u, 0u), 40.0 + expected_delta_e, 1.0e-12,
              "two-group accumulated electron update");
    CheckNear(state.e_fluid_total(0u, 0u, 0u), 90.0 + expected_delta_e, 1.0e-12,
              "two-group accumulated total-fluid update");
    DEC3D_CHECK_EQ(result.per_group_solve_count, 2u);
    DEC3D_CHECK(result.all_groups_updated);
  }

  {
    const double dt_s = 1.0e-6;
    auto multigroup_state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 1u});
    multigroup_state.rho(0u, 0u, 0u) = 1.0;
    multigroup_state.radiation_groups[0](0u, 0u, 0u) = 2.0;
    multigroup_state.e_electron(0u, 0u, 0u) = 20.0;
    multigroup_state.e_fluid_total(0u, 0u, 0u) = 50.0;

    auto one_group_state = multigroup_state;

    dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions multi_options;
    multi_options.dt_s = dt_s;
    multi_options.group_layout =
        dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout({1.0e15, 2.0e15});
    multi_options.boundary_policy = MakePatchBoundary();
    multi_options.coefficients.Dbar_cm2_per_s = {OneCellArray(0.0)};
    multi_options.coefficients.kappaP_cm_inv = {OneCellArray(1.0e-4)};
    multi_options.coefficients.B_erg_per_cm3 = {OneCellArray(5.0)};

    dec3d::radiation::OneGroupGrayRadiationMatterCouplingOptions one_options;
    one_options.radiation_options.dt_s = dt_s;
    one_options.radiation_options.group_layout =
        dec3d::radiation::MakeGrayFullSpectrumRadiationGroupLayout();
    one_options.radiation_options.boundary_policy = MakePatchBoundary();
    one_options.radiation_options.coefficients.Dbar_cm2_per_s = 0.0;
    one_options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    one_options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;

    const auto multi_result =
        dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
            multigroup_state, MakeOneCellPatchGeometry(), multi_options);
    const auto one_result =
        dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
            one_group_state, MakeOneCellPatchGeometry(), one_options);
    CheckSuccess(multi_result);
    if (!one_result.success) {
      dec3d::test::Fail("P3-2 one-group baseline success", __FILE__, __LINE__);
    }
    CheckNear(multigroup_state.radiation_groups[0](0u, 0u, 0u),
              one_group_state.radiation_groups[0](0u, 0u, 0u),
              1.0e-12,
              "group_count one matches P3-2 U");
    CheckNear(multigroup_state.e_electron(0u, 0u, 0u),
              one_group_state.e_electron(0u, 0u, 0u),
              1.0e-12,
              "group_count one matches P3-2 electron");
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    auto options = MakeTwoGroupOptions(1.0e-6);
    const double before_group0 = state.radiation_groups[0](0u, 0u, 0u);
    const double before_group1 = state.radiation_groups[1](0u, 0u, 0u);
    options.coefficients.Dbar_cm2_per_s[1](0u, 0u, 0u) = -1.0;

    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.first_failing_group_index, 1u);
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), before_group0, 0.0,
              "group 0 not published after group 1 failure");
    CheckNear(state.radiation_groups[1](0u, 0u, 0u), before_group1, 0.0,
              "group 1 not published after group 1 failure");
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeTwoGroupOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(result.per_group_solve_count, 0u);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    auto options = MakeTwoGroupOptions(1.0e-6);
    options.radiation_boundary_model =
        dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;
    options.coefficients.Dbar_cm2_per_s = {OneCellArray(4.0), OneCellArray(3.0)};
    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), options);
    CheckSuccess(result);
    DEC3D_CHECK(result.boundary_leak_total_all_groups >= -1.0e-18);
    CheckNear(result.delta_radiation_total_all_groups + result.delta_electron_total +
                  result.boundary_leak_total_all_groups,
              0.0,
              1.0e-10,
              "all-group Marshak budget closes");
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeTwoGroupOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(result));
    const std::string token = "per_group_solve_count=";
    result.report_line.erase(result.report_line.find(token), token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(result));
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    auto options = MakeTwoGroupOptions(1.0e-6);
    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), options);
    CheckSuccess(result);
    DEC3D_CHECK(state.radiation_groups[0](0u, 0u, 0u) != 2.0);
    DEC3D_CHECK(state.radiation_groups[1](0u, 0u, 0u) != 3.0);
    DEC3D_CHECK(result.report_line.find("all_groups_updated=true") != std::string::npos);
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeTwoGroupOptions(1.0e-6));
    CheckSuccess(result);
    CheckNear(result.boundary_leak_total_all_groups, 0.0, 1.0e-12,
              "zero-flux all-group boundary leak");
    CheckNear(result.global_radiation_plus_electron_residual, 0.0, 1.0e-12,
              "zero-flux all-group radiation plus electron residual");
  }

  {
    auto state = MakeTwoGroupOneCellState(2.0, 3.0, 40.0, 50.0);
    const auto result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeTwoGroupOptions(1.0e-6));
    CheckSuccess(result);
    DEC3D_CHECK(result.report_line.find("frequency_edges_unit=Hz") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("group_layout_report_present=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("per_group_matrix_report_present=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("per_group_solve_report_present=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("first_failing_group_index=none") != std::string::npos);
    DEC3D_CHECK_EQ(result.per_group_radiation_reports.size(), 2u);
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
