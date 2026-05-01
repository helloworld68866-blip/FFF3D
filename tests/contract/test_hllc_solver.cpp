#include "hydro/riemann/hllc_solver.hpp"
#include "radiation/hydro_terms/radiation_hydro_terms.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <iostream>
#include <string>

namespace {

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << label << ": actual=" << actual << " expected=" << expected << '\n';
    DEC3D_CHECK(false);
  }
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::ComputePhysicalFlux;
    using dec3d::hydro::HllcActiveRegion;
    using dec3d::hydro::HydroConservativeState;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::hydro::SolveHllcRiemann;

    const HydroPrimitiveState left{
        1.0,
        0.0,
        0.35,
        -0.20,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const HydroPrimitiveState right{
        0.125,
        0.0,
        -0.15,
        0.25,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};

    const auto left_state = MakeConservativeState(left);
    const auto right_state = MakeConservativeState(right);
    const auto left_flux = ComputePhysicalFlux(left_state);
    const auto right_flux = ComputePhysicalFlux(right_state);
    const auto result = SolveHllcRiemann(left_state, right_state);

    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.waves.is_complete());
    DEC3D_CHECK(result.waves.s_left < result.waves.s_star);
    DEC3D_CHECK(result.waves.s_star < result.waves.s_right);
    DEC3D_CHECK(
        result.active_region == HllcActiveRegion::left_star ||
        result.active_region == HllcActiveRegion::right_star);
    DEC3D_CHECK(std::isfinite(result.left_star_state.rho));
    DEC3D_CHECK(std::isfinite(result.right_star_state.rho));
    DEC3D_CHECK(result.left_star_state.rho > 0.0);
    DEC3D_CHECK(result.right_star_state.rho > 0.0);
    DEC3D_CHECK(std::abs(result.interface_flux.rho - 0.5 * (left_flux.rho + right_flux.rho)) > 1.0e-8);
    DEC3D_CHECK(std::isfinite(result.interface_flux.chi_e));

    {
      HydroPrimitiveState primitive;
      primitive.rho = 2.0;
      primitive.v_r = 1.5;
      primitive.v_theta = 0.0;
      primitive.v_phi = 0.0;
      primitive.pressure = 10.0;
      primitive.chi_e = 1.2;
      primitive.radiation_chi = {2.0, 4.0};

      const auto conservative = MakeConservativeState(primitive);
      DEC3D_CHECK_EQ(conservative.radiation_chi.size(), std::size_t{2});
      DEC3D_CHECK(std::abs(conservative.radiation_chi[0] - 2.0) < 1.0e-14);
      DEC3D_CHECK(std::abs(conservative.radiation_chi[1] - 4.0) < 1.0e-14);

      const auto flux = ComputePhysicalFlux(conservative);
      DEC3D_CHECK_EQ(flux.radiation_chi.size(), std::size_t{2});
      DEC3D_CHECK(std::abs(flux.radiation_chi[0] - 3.0) < 1.0e-14);
      DEC3D_CHECK(std::abs(flux.radiation_chi[1] - 6.0) < 1.0e-14);
    }

    {
      HydroPrimitiveState left_rad;
      left_rad.rho = 1.0;
      left_rad.v_r = 0.4;
      left_rad.pressure = 4.0;
      left_rad.chi_e = 1.0;
      left_rad.radiation_chi = {
          dec3d::radiation::RadiationPressureScalarFromUg(9.0)};

      HydroPrimitiveState right_rad;
      right_rad.rho = 0.5;
      right_rad.v_r = -0.3;
      right_rad.pressure = 3.5;
      right_rad.chi_e = 1.0;
      right_rad.radiation_chi = {
          dec3d::radiation::RadiationPressureScalarFromUg(6.0)};

      const auto hllc = SolveHllcRiemann(
          MakeConservativeState(left_rad),
          MakeConservativeState(right_rad));
      DEC3D_CHECK(hllc.success);
      DEC3D_CHECK_EQ(hllc.interface_flux.radiation_chi.size(), std::size_t{1});
      DEC3D_CHECK(std::isfinite(hllc.interface_flux.radiation_chi[0]));
      DEC3D_CHECK_EQ(hllc.left_star_state.radiation_chi.size(), std::size_t{1});

      const double rho_ratio = hllc.left_star_state.rho / left_rad.rho;
      const double chi_ratio =
          hllc.left_star_state.radiation_chi[0] / left_rad.radiation_chi[0];
      CheckNear(chi_ratio, rho_ratio, 1.0e-12, "HLLC chi_rad scales with density");

      const double ug_left = 9.0;
      const double ug_star =
          dec3d::radiation::UgFromRadiationPressureScalar(
              hllc.left_star_state.radiation_chi[0]);
      CheckNear(
          ug_star / ug_left,
          std::pow(rho_ratio, 4.0 / 3.0),
          1.0e-12,
          "HLLC radiation Ug compression");

      const double passive_ug_ratio = rho_ratio;
      DEC3D_CHECK(std::abs((ug_star / ug_left) - passive_ug_ratio) > 1.0e-8);
    }

    {
      HydroPrimitiveState primitive;
      primitive.rho = 1.0;
      primitive.v_r = 0.0;
      primitive.pressure = 5.0;
      primitive.chi_e = 1.0;
      primitive.radiation_chi = {1.0};
      auto left_mismatch = MakeConservativeState(primitive);
      primitive.radiation_chi = {1.0, 2.0};
      auto right_mismatch = MakeConservativeState(primitive);
      const auto mismatch = SolveHllcRiemann(left_mismatch, right_mismatch);
      DEC3D_CHECK(!mismatch.success);
      DEC3D_CHECK(mismatch.failure_reason.find("radiation scalar bundle size mismatch") !=
                  std::string::npos);
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
