#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/transport/one_group_gray_diffusion.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

namespace {

dec3d::mesh::SphericalGeometryMetadata MakePatchGeometry() {
  return dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{2u, 2u, 2u, 0.1, 0.3});
}

dec3d::mesh::SphericalGeometryMetadata MakeOneCellPatchGeometry() {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces = {0.1, 0.3};
  geometry.theta_faces = {0.5, 1.2};
  geometry.phi_faces = {0.0, 2.0 * std::acos(-1.0)};
  const double radial_factor =
      (std::pow(geometry.radial_faces[1], 3) - std::pow(geometry.radial_faces[0], 3)) / 3.0;
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

dec3d::state::CanonicalState MakeRadiationState() {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{2u, 2u, 2u, 1u});
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        state.radiation_groups[0](r, theta, phi) = 4.0 + static_cast<double>(r + theta + phi);
        state.e_electron(r, theta, phi) = 11.0;
        state.e_fluid_total(r, theta, phi) = 17.0;
      }
    }
  }
  return state;
}

dec3d::state::CanonicalState MakeOneCellRadiationState(double ug) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 1u});
  state.radiation_groups[0](0u, 0u, 0u) = ug;
  state.e_electron(0u, 0u, 0u) = 11.0;
  state.e_fluid_total(0u, 0u, 0u) = 17.0;
  return state;
}

dec3d::state::CanonicalState MakeOneCellRadiationMatterState(
    double ug,
    double e_electron,
    double e_ion) {
  auto state = MakeOneCellRadiationState(ug);
  state.rho(0u, 0u, 0u) = 1.0;
  state.mom_r(0u, 0u, 0u) = 0.0;
  state.mom_theta(0u, 0u, 0u) = 0.0;
  state.mom_phi(0u, 0u, 0u) = 0.0;
  state.e_electron(0u, 0u, 0u) = e_electron;
  state.e_fluid_total(0u, 0u, 0u) = e_electron + e_ion;
  return state;
}

dec3d::radiation::OneGroupGrayRadiationOptions MakeOptions(double dt_s) {
  dec3d::radiation::OneGroupGrayRadiationOptions options;
  options.dt_s = dt_s;
  options.group_layout = dec3d::radiation::MakeGrayFullSpectrumRadiationGroupLayout();
  options.boundary_policy = MakePatchBoundary();
  options.coefficients.Dbar_cm2_per_s = 1.0;
  options.coefficients.kappaP_cm_inv = 1.0e-4;
  options.coefficients.Bgray_erg_per_cm3 = 5.0;
  return options;
}

dec3d::radiation::OneGroupGrayRadiationOptions MakeMarshakOptions(double dt_s) {
  auto options = MakeOptions(dt_s);
  options.radiation_boundary_model =
      dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;
  return options;
}

dec3d::radiation::OneGroupGrayRadiationMatterCouplingOptions MakeCouplingOptions(
    double dt_s) {
  dec3d::radiation::OneGroupGrayRadiationMatterCouplingOptions options;
  options.radiation_options = MakeOptions(dt_s);
  options.electron_energy_floor_erg_per_cm3 = 0.0;
  options.ion_energy_floor_erg_per_cm3 = 0.0;
  return options;
}

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    dec3d::test::Fail(label, __FILE__, __LINE__);
  }
}

double OneCellOuterArea(const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const double radius = geometry.radial_faces.back();
  const double polar_factor = std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1]);
  const double dphi = geometry.phi_faces[1] - geometry.phi_faces[0];
  return radius * radius * polar_factor * dphi;
}

double OneCellMarshakGFace(double dbar, const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  const double center = 0.5 * (geometry.radial_faces[0] + geometry.radial_faces[1]);
  const double distance = geometry.radial_faces[1] - center;
  const double h = 0.5 * dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double d_over_distance = dbar / distance;
  return (h * d_over_distance) / (h + d_over_distance);
}

double VolumeWeightedRadiationTotal(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double total = 0.0;
  std::size_t linear = 0u;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        total += state.radiation_groups[0](r, theta, phi) * geometry.cell_volumes[linear++];
      }
    }
  }
  return total;
}

double VolumeWeightedElectronTotal(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  double total = 0.0;
  std::size_t linear = 0u;
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        total += state.e_electron(r, theta, phi) * geometry.cell_volumes[linear++];
      }
    }
  }
  return total;
}

void CheckSuccess(const dec3d::radiation::OneGroupGrayRadiationResult& result) {
  if (!result.success) {
    dec3d::test::Fail("radiation result success", __FILE__, __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

void CheckSuccess(
    const dec3d::radiation::OneGroupGrayRadiationMatterCouplingResult& result) {
  if (!result.success) {
    dec3d::test::Fail("radiation matter coupling result success", __FILE__, __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

void CheckRadiationUnchanged(
    const dec3d::state::CanonicalState& state,
    const dec3d::core::Array3D<double>& before) {
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        CheckNear(state.radiation_groups[0](r, theta, phi), before(r, theta, phi), 0.0,
                  "radiation group unchanged");
        CheckNear(state.e_electron(r, theta, phi), 11.0, 0.0, "electron energy unchanged");
        CheckNear(state.e_fluid_total(r, theta, phi), 17.0, 0.0, "fluid energy unchanged");
      }
    }
  }
}

}  // namespace

int RunTests() {
  {
    auto state = MakeRadiationState();
    const auto before = state.radiation_groups[0];
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), MakeOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("metadata_written=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("marshak_enabled=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("parity_claim_allowed=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("matrix_report_present=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("solve_report_present=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("min_Dbar_cm2_per_s=") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("max_Dbar_cm2_per_s=") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("min_kappaP_cm_inv=") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("max_kappaP_cm_inv=") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("min_Bgray_erg_per_cm3=") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("max_Bgray_erg_per_cm3=") != std::string::npos);
    CheckRadiationUnchanged(state, before);
  }

  {
    auto state = MakeOneCellRadiationState(3.0);
    const auto before = state.radiation_groups[0];
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakeOneCellPatchGeometry(), MakeMarshakOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.marshak_enabled);
    DEC3D_CHECK(result.report_line.find("radiation_boundary_model=thesis_marshak_vacuum") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("marshak_enabled=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("metadata_written=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), before(0u, 0u, 0u), 0.0,
              "dt zero Marshak radiation unchanged");
  }

  {
    auto state = MakeRadiationState();
    const auto before = state.radiation_groups[0];
    auto options = MakeOptions(1.0e-6);
    options.coefficients.Dbar_cm2_per_s = 0.0;
    options.coefficients.kappaP_cm_inv = 0.0;
    options.coefficients.Bgray_erg_per_cm3 = 0.0;

    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(state, MakePatchGeometry(), options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("metadata_written=false") != std::string::npos);
    CheckRadiationUnchanged(state, before);
  }

  {
    const auto geometry = MakeOneCellPatchGeometry();
    auto zero_flux_state = MakeOneCellRadiationState(8.0);
    auto marshak_state = MakeOneCellRadiationState(8.0);
    auto zero_flux_options = MakeOptions(1.0e-6);
    zero_flux_options.coefficients.Dbar_cm2_per_s = 4.0;
    zero_flux_options.coefficients.kappaP_cm_inv = 0.0;
    zero_flux_options.coefficients.Bgray_erg_per_cm3 = 0.0;
    auto marshak_options = zero_flux_options;
    marshak_options.radiation_boundary_model =
        dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;

    const auto zero_flux_result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        zero_flux_state, geometry, zero_flux_options);
    const auto marshak_result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        marshak_state, geometry, marshak_options);
    CheckSuccess(zero_flux_result);
    CheckSuccess(marshak_result);
    DEC3D_CHECK(marshak_result.marshak_enabled);
    DEC3D_CHECK(marshak_result.marshak_nonzero_diagonal_loss_count > 0u);
    DEC3D_CHECK(zero_flux_state.radiation_groups[0](0u, 0u, 0u) >
                marshak_state.radiation_groups[0](0u, 0u, 0u));
    CheckNear(marshak_state.e_electron(0u, 0u, 0u), 11.0, 0.0,
              "Marshak leaves electron energy unchanged");
    CheckNear(marshak_state.e_fluid_total(0u, 0u, 0u), 17.0, 0.0,
              "Marshak leaves fluid energy unchanged");
  }

  {
    const double dt_s = 1.0e-6;
    const double u_old = 2.0;
    const auto geometry = MakeOneCellPatchGeometry();
    auto state = MakeOneCellRadiationState(u_old);
    auto options = MakeMarshakOptions(dt_s);
    options.coefficients.Dbar_cm2_per_s = 3.0;
    options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.coefficients.Bgray_erg_per_cm3 = 5.0;

    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(state, geometry, options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
    DEC3D_CHECK(result.marshak_enabled);

    const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
    const double outer_area = OneCellOuterArea(geometry);
    const double volume = geometry.cell_volumes[0];
    const double g_face = OneCellMarshakGFace(options.coefficients.Dbar_cm2_per_s, geometry);
    const double expected =
        (u_old + dt_s * c * options.coefficients.kappaP_cm_inv *
                     options.coefficients.Bgray_erg_per_cm3) /
        (1.0 + dt_s * (c * options.coefficients.kappaP_cm_inv +
                       (outer_area / volume) * g_face));
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), expected, 1.0e-12,
              "one-cell Marshak plus absorption/emission oracle");
  }

  {
    auto state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{1u, 2u, 2u, 1u});
    for (std::size_t theta = 0u; theta < 2u; ++theta) {
      for (std::size_t phi = 0u; phi < 2u; ++phi) {
        state.radiation_groups[0](0u, theta, phi) = 2.0;
        state.e_electron(0u, theta, phi) = 11.0;
        state.e_fluid_total(0u, theta, phi) = 17.0;
      }
    }

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{1u, 2u, 2u, 0.1, 0.2});

    auto options = MakeOptions(1.0e-6);
    options.coefficients.Dbar_cm2_per_s = 0.0;
    options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.coefficients.Bgray_erg_per_cm3 = 5.0;
    options.boundary_policy.inner_radial =
        dec3d::transport::DiffusionBoundaryKind::interior_patch_no_origin;
    options.boundary_policy.theta_lower =
        dec3d::transport::DiffusionBoundaryKind::scalar_pole_remap_required;
    options.boundary_policy.theta_upper =
        dec3d::transport::DiffusionBoundaryKind::scalar_pole_remap_required;

    const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
    const double rate = options.dt_s * c * options.coefficients.kappaP_cm_inv;
    const double expected = (2.0 + rate * options.coefficients.Bgray_erg_per_cm3) / (1.0 + rate);

    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, geometry, options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
    for (std::size_t theta = 0u; theta < 2u; ++theta) {
      for (std::size_t phi = 0u; phi < 2u; ++phi) {
        CheckNear(state.radiation_groups[0](0u, theta, phi), expected, 1.0e-12,
                  "gray radiation source oracle");
      }
    }
    DEC3D_CHECK_EQ(result.updated_fields,
                  dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask, result.updated_fields);
    DEC3D_CHECK(result.report_line.find("metadata_written=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=true") != std::string::npos);
    for (std::size_t theta = 0u; theta < 2u; ++theta) {
      for (std::size_t phi = 0u; phi < 2u; ++phi) {
        CheckNear(state.e_electron(0u, theta, phi), 11.0, 0.0, "electron energy unchanged");
        CheckNear(state.e_fluid_total(0u, theta, phi), 17.0, 0.0, "fluid energy unchanged");
      }
    }
  }

  {
    auto state = MakeOneCellRadiationState(2.0);
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakeOneCellPatchGeometry(), MakeMarshakOptions(1.0e-12));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));

    const auto token = std::string{"marshak_nonzero_diagonal_loss_count="};
    const auto pos = result.report_line.find(token);
    DEC3D_CHECK(pos != std::string::npos);
    result.report_line.erase(pos, token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
  }

  {
    auto state = MakeRadiationState();
    state.radiation_groups[0](0u, 0u, 0u) = -1.0e-30;
    const auto before = state.radiation_groups[0];
    auto options = MakeOptions(1.0e-12);
    options.radiation_energy_floor_erg_per_cm3 = 0.0;
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.failure_diagnostics.find(
                    "diagnostic_id=p3.radiation.one_group_gray_diffusion.failure") !=
                std::string::npos);
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), before(0u, 0u, 0u), 0.0,
              "floor failure leaves radiation state unchanged");
  }

  {
    auto state = MakeRadiationState();
    auto options = MakeOptions(1.0e-12);
    options.coefficients.Dbar_cm2_per_s = -1.0;
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("gray radiation coefficients must be finite and nonnegative") !=
                std::string::npos);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
  }

  {
    auto state = MakeRadiationState();
    auto options = MakeOptions(1.0e-12);
    options.group_layout = dec3d::radiation::RadiationGroupLayout{};
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("radiation group mode is missing") != std::string::npos);
  }

  {
    auto state = MakeRadiationState();
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), MakeOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));

    result.report_line.erase(
        result.report_line.find("marshak_enabled=false"),
        std::string{"marshak_enabled=false"}.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
  }

  {
    auto state = MakeRadiationState();
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), MakeOptions(0.0));
    CheckSuccess(result);
    const auto token = std::string{"radiation_energy_floor_erg_per_cm3="};
    result.report_line.erase(result.report_line.find(token), token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
  }

  {
    auto state = MakeRadiationState();
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), MakeOptions(0.0));
    CheckSuccess(result);
    const auto token = std::string{"matrix_report_present=true"};
    result.report_line.erase(result.report_line.find(token), token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
  }

  {
    auto state = MakeRadiationState();
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationDiffusion(
        state, MakePatchGeometry(), MakeOptions(0.0));
    CheckSuccess(result);
    const auto token = std::string{"min_Dbar_cm2_per_s="};
    result.report_line.erase(result.report_line.find(token), token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationDiagnostics(result));
  }

  {
    const double dt_s = 1.0e-6;
    const double u_old = 2.0;
    const double e_e_old = 20.0;
    const double e_i_old = 30.0;
    auto state = MakeOneCellRadiationMatterState(u_old, e_e_old, e_i_old);
    auto options = MakeCouplingOptions(dt_s);
    options.radiation_options.coefficients.Dbar_cm2_per_s = 0.0;
    options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;

    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), options);
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationMatterCouplingDiagnostics(result));

    const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
    const double rate = dt_s * c * options.radiation_options.coefficients.kappaP_cm_inv;
    const double expected_u =
        (u_old + rate * options.radiation_options.coefficients.Bgray_erg_per_cm3) /
        (1.0 + rate);
    const double expected_exchange = -(expected_u - u_old);
    CheckNear(state.radiation_groups[0](0u, 0u, 0u), expected_u, 1.0e-12,
              "P3-2 one-cell U oracle");
    CheckNear(state.e_electron(0u, 0u, 0u), e_e_old + expected_exchange, 1.0e-12,
              "P3-2 one-cell electron exchange");
    CheckNear(state.e_fluid_total(0u, 0u, 0u), e_e_old + e_i_old + expected_exchange, 1.0e-12,
              "P3-2 one-cell total fluid exchange");
    CheckNear(
        state.e_fluid_total(0u, 0u, 0u) - state.e_electron(0u, 0u, 0u),
        e_i_old,
        1.0e-12,
        "P3-2 ion reservoir remains derived and unchanged");
  }

  {
    const auto geometry = MakeOneCellPatchGeometry();
    auto state = MakeOneCellRadiationMatterState(2.0, 20.0, 30.0);
    const double before =
        VolumeWeightedRadiationTotal(state, geometry) +
        VolumeWeightedElectronTotal(state, geometry);
    auto options = MakeCouplingOptions(1.0e-6);
    options.radiation_options.coefficients.Dbar_cm2_per_s = 0.0;
    options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;

    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(state, geometry, options);
    CheckSuccess(result);
    const double after =
        VolumeWeightedRadiationTotal(state, geometry) +
        VolumeWeightedElectronTotal(state, geometry);
    CheckNear(after, before, 1.0e-12, "P3-2 zero-flux radiation plus electron conservation");
    CheckNear(result.boundary_leak_total, 0.0, 1.0e-12, "P3-2 zero-flux boundary leak");
    CheckNear(result.global_radiation_plus_electron_residual, 0.0, 1.0e-12,
              "P3-2 zero-flux residual");
  }

  {
    const auto geometry = MakeOneCellPatchGeometry();
    auto state = MakeOneCellRadiationMatterState(8.0, 20.0, 30.0);
    auto options = MakeCouplingOptions(1.0e-6);
    options.radiation_options = MakeMarshakOptions(1.0e-6);
    options.radiation_options.coefficients.Dbar_cm2_per_s = 4.0;
    options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;

    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(state, geometry, options);
    CheckSuccess(result);
    DEC3D_CHECK(result.marshak_enabled);
    DEC3D_CHECK(result.boundary_leak_total >= -1.0e-18);
    CheckNear(
        result.delta_radiation_total + result.delta_electron_total + result.boundary_leak_total,
        0.0,
        1.0e-10,
        "P3-2 Marshak leakage closes radiation plus electron budget");
  }

  {
    auto state = MakeOneCellRadiationMatterState(2.0, 20.0, 30.0);
    const auto result = dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeCouplingOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
  }

  {
    auto state = MakeOneCellRadiationMatterState(2.0, 1.0e-12, 30.0);
    auto options = MakeCouplingOptions(1.0e-6);
    options.radiation_options.coefficients.Dbar_cm2_per_s = 0.0;
    options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;
    options.electron_energy_floor_erg_per_cm3 = 1.0e-9;
    const auto before_e = state.e_electron(0u, 0u, 0u);
    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
            state, MakeOneCellPatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    CheckNear(state.e_electron(0u, 0u, 0u), before_e, 0.0,
              "P3-2 electron floor failure no publish");
    DEC3D_CHECK(result.failure_diagnostics.find(
                    "diagnostic_id=p3.radiation.one_group_gray_matter_coupling.failure") !=
                std::string::npos);
  }

  {
    auto state = MakeOneCellRadiationMatterState(2.0, 20.0, 1.0e-12);
    auto options = MakeCouplingOptions(1.0e-6);
    options.radiation_options.coefficients.Dbar_cm2_per_s = 0.0;
    options.radiation_options.coefficients.kappaP_cm_inv = 1.0e-4;
    options.radiation_options.coefficients.Bgray_erg_per_cm3 = 5.0;
    options.ion_energy_floor_erg_per_cm3 = 1.0e-9;
    const auto before_total = state.e_fluid_total(0u, 0u, 0u);
    const auto result =
        dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
            state, MakeOneCellPatchGeometry(), options);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    CheckNear(state.e_fluid_total(0u, 0u, 0u), before_total, 0.0,
              "P3-2 ion floor failure no publish");
  }

  {
    auto state = MakeOneCellRadiationMatterState(2.0, 20.0, 30.0);
    auto result = dec3d::radiation::ApplyOneGroupGrayRadiationMatterCoupling(
        state, MakeOneCellPatchGeometry(), MakeCouplingOptions(0.0));
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateOneGroupGrayRadiationMatterCouplingDiagnostics(result));
    const std::string token = "boundary_leak_total=";
    result.report_line.erase(result.report_line.find(token), token.size());
    DEC3D_CHECK(!dec3d::radiation::ValidateOneGroupGrayRadiationMatterCouplingDiagnostics(result));
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
