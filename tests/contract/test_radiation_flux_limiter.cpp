#include "radiation/transport/radiation_flux_limiter.hpp"

#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace {

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    dec3d::test::Fail(
        label,
        __FILE__,
        __LINE__,
        std::to_string(actual) + " vs " + std::to_string(expected));
  }
}

void CheckComplete(const dec3d::radiation::RadiationFluxLimiterFaceResult& result) {
  DEC3D_CHECK(result.success);
  DEC3D_CHECK(dec3d::radiation::ValidateRadiationFluxLimiterDiagnostics(result));
  DEC3D_CHECK(result.report_line.find("face_limiter_time_level=old_time_lagged") !=
              std::string::npos);
  DEC3D_CHECK(result.report_line.find("fallback_used=false") != std::string::npos);
}

void CheckValueParity(
    const dec3d::radiation::RadiationFluxLimiterFaceInput& input,
    const char* label) {
  const auto detailed = dec3d::radiation::ApplyRadiationFluxLimiter(input);
  const auto values = dec3d::radiation::ApplyRadiationFluxLimiterValues(input);
  CheckComplete(detailed);
  DEC3D_CHECK(values.success);
  DEC3D_CHECK_EQ(values.limited, detailed.limited);
  CheckNear(
      values.effective_D_cm2_per_s,
      detailed.effective_D_cm2_per_s,
      std::max(1.0, std::abs(detailed.effective_D_cm2_per_s)) * 1.0e-15,
      label);
  CheckNear(
      values.limiter_scale,
      detailed.limiter_scale,
      1.0e-15,
      label);
  CheckNear(
      values.flux_ratio_before_limit,
      detailed.flux_ratio_before_limit,
      std::max(1.0, std::abs(detailed.flux_ratio_before_limit)) * 1.0e-15,
      label);
  DEC3D_CHECK(values.report_line.empty());
}

}  // namespace

int main() {
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;

  {
    dec3d::radiation::RadiationFluxLimiterFaceInput input;
    input.model = dec3d::radiation::RadiationFluxLimiterModel::disabled;
    input.unlimited_D_cm2_per_s = 4.0;
    input.U_face_erg_per_cm3 = 2.0;
    input.abs_grad_U_erg_per_cm4 = 7.0;
    input.U_floor_erg_per_cm3 = 1.0e-40;
    const auto result = dec3d::radiation::ApplyRadiationFluxLimiter(input);
    CheckValueParity(input, "disabled value-only parity");
    CheckComplete(result);
    DEC3D_CHECK(!result.limited);
    CheckNear(result.effective_D_cm2_per_s, 4.0, 0.0, "disabled D");
    CheckNear(result.limiter_scale, 1.0, 0.0, "disabled scale");
  }

  {
    dec3d::radiation::RadiationFluxLimiterFaceInput input;
    input.model = dec3d::radiation::RadiationFluxLimiterModel::harmonic_eq_5_209;
    input.unlimited_D_cm2_per_s = 4.0;
    input.U_face_erg_per_cm3 = 2.0;
    input.abs_grad_U_erg_per_cm4 = 8.0;
    input.U_floor_erg_per_cm3 = 1.0e-40;
    const auto result = dec3d::radiation::ApplyRadiationFluxLimiter(input);
    CheckValueParity(input, "harmonic value-only parity");
    CheckComplete(result);
    const double sigma = c / input.unlimited_D_cm2_per_s;
    const double r_face = input.abs_grad_U_erg_per_cm4 / input.U_face_erg_per_cm3;
    CheckNear(result.effective_D_cm2_per_s,
              c / (sigma + r_face),
              1.0e-12,
              "harmonic oracle");
    DEC3D_CHECK(result.effective_D_cm2_per_s <= input.unlimited_D_cm2_per_s);
  }

  {
    dec3d::radiation::RadiationFluxLimiterFaceInput input;
    input.model = dec3d::radiation::RadiationFluxLimiterModel::larsen_eq_5_210_n2;
    input.unlimited_D_cm2_per_s = 3.0;
    input.U_face_erg_per_cm3 = 4.0;
    input.abs_grad_U_erg_per_cm4 = 10.0;
    input.U_floor_erg_per_cm3 = 1.0e-40;
    const auto result = dec3d::radiation::ApplyRadiationFluxLimiter(input);
    CheckValueParity(input, "Larsen value-only parity");
    CheckComplete(result);
    const double sigma = c / input.unlimited_D_cm2_per_s;
    const double r_face = input.abs_grad_U_erg_per_cm4 / input.U_face_erg_per_cm3;
    CheckNear(result.effective_D_cm2_per_s,
              c / std::sqrt(sigma * sigma + r_face * r_face),
              1.0e-12,
              "Larsen n=2 oracle");
  }

  {
    dec3d::radiation::RadiationFluxLimiterFaceInput input;
    input.model = dec3d::radiation::RadiationFluxLimiterModel::minmax_eq_5_211;
    input.unlimited_D_cm2_per_s = 5.0;
    input.U_face_erg_per_cm3 = 2.0;
    input.abs_grad_U_erg_per_cm4 = 4.0 * c;
    input.U_floor_erg_per_cm3 = 1.0e-40;
    const auto result = dec3d::radiation::ApplyRadiationFluxLimiter(input);
    CheckValueParity(input, "minmax value-only parity");
    CheckComplete(result);
    CheckNear(result.effective_D_cm2_per_s, 0.5, 1.0e-12, "minmax oracle");
    DEC3D_CHECK(result.limited);
    DEC3D_CHECK(result.flux_ratio_before_limit > 1.0);
  }

  {
    dec3d::radiation::RadiationFluxLimiterFaceInput input;
    input.model = dec3d::radiation::RadiationFluxLimiterModel::harmonic_eq_5_209;
    input.unlimited_D_cm2_per_s = -1.0;
    input.U_face_erg_per_cm3 = 1.0;
    input.abs_grad_U_erg_per_cm4 = 1.0;
    input.U_floor_erg_per_cm3 = 1.0e-40;
    const auto result = dec3d::radiation::ApplyRadiationFluxLimiter(input);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_diagnostics.find(
                    "diagnostic_id=p3.radiation.flux_limiter.failure") !=
                std::string::npos);
  }

  {
    dec3d::transport::DistributedDiffusionRowOwnership ownership;
    ownership.success = true;
    ownership.rank = 0;
    ownership.rank_count = 1;
    ownership.global_radial_cells = 2;
    ownership.global_theta_cells = 2;
    ownership.global_phi_cells = 2;
    ownership.global_radial_begin = 0;
    ownership.global_radial_end = 2;
    ownership.local_row_count = 8;

    dec3d::state::CanonicalState state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{2, 2, 2, 1});
    dec3d::core::Array3D<double> unlimited_D(2, 2, 2, 1.0e12);
    for (std::size_t r = 0; r < 2; ++r) {
      for (std::size_t t = 0; t < 2; ++t) {
        for (std::size_t p = 0; p < 2; ++p) {
          state.radiation_groups[0](r, t, p) = r == 0 ? 100.0 : 1.0;
        }
      }
    }

    dec3d::mesh::SphericalGeometryMetadata geometry;
    geometry.valid = true;
    geometry.radial_faces = {1.0, 2.0, 3.0};
    geometry.theta_faces = {0.5, 1.0, 1.5};
    geometry.phi_faces = {0.0, 3.14159265358979323846, 6.28318530717958647692};
    geometry.cell_volumes.assign(8, 1.0);

    dec3d::radiation::RadiationFluxLimiterOptions options;
    options.model = dec3d::radiation::RadiationFluxLimiterModel::minmax_eq_5_211;
    options.U_floor_erg_per_cm3 = 1.0e-40;

    const auto arrays = dec3d::radiation::BuildRadiationFluxLimitedFaceCoefficients(
        ownership, geometry, state.radiation_groups[0], unlimited_D, options);
    DEC3D_CHECK(arrays.success);
    DEC3D_CHECK(arrays.face_coefficients.enabled);
    DEC3D_CHECK(arrays.limited_face_count_local > 0u);
    DEC3D_CHECK(arrays.min_limiter_scale < 1.0);
    DEC3D_CHECK(arrays.max_flux_ratio_before_limit > 1.0);
    DEC3D_CHECK(arrays.report_line.find(
                    "radiation_flux_limiter_model=minmax_eq_5_211") !=
                std::string::npos);
    DEC3D_CHECK(dec3d::radiation::ValidateRadiationFluxLimitedFaceDiagnostics(arrays));
  }

  {
    dec3d::transport::DistributedDiffusionRowOwnership ownership;
    ownership.success = true;
    ownership.rank = 0;
    ownership.rank_count = 1;
    ownership.global_radial_cells = 1;
    ownership.global_theta_cells = 1;
    ownership.global_phi_cells = 1;
    ownership.global_radial_begin = 0;
    ownership.global_radial_end = 1;
    ownership.local_row_count = 1;
    dec3d::state::CanonicalState state = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{1, 1, 1, 1});
    state.radiation_groups[0](0, 0, 0) = 1.0;
    dec3d::core::Array3D<double> unlimited_D(1, 1, 1, 7.0);
    dec3d::mesh::SphericalGeometryMetadata geometry;
    geometry.valid = true;
    geometry.radial_faces = {1.0, 2.0};
    geometry.theta_faces = {0.5, 1.5};
    geometry.phi_faces = {0.0, 6.28318530717958647692};
    geometry.cell_volumes.assign(1, 1.0);

    dec3d::radiation::RadiationFluxLimiterOptions options;
    options.model = dec3d::radiation::RadiationFluxLimiterModel::disabled;
    options.U_floor_erg_per_cm3 = 0.0;
    const auto arrays = dec3d::radiation::BuildRadiationFluxLimitedFaceCoefficients(
        ownership, geometry, state.radiation_groups[0], unlimited_D, options);
    DEC3D_CHECK(arrays.success);
    DEC3D_CHECK(!arrays.face_coefficients.enabled);
    DEC3D_CHECK_EQ(arrays.limited_face_count_local, std::size_t{0});
    DEC3D_CHECK(arrays.report_line.find("radiation_flux_limiter_enabled=false") !=
                std::string::npos);
  }

  return 0;
}
