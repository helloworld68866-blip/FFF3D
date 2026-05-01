#include "physics/units/physical_constants.hpp"
#include "transport/thermal/electron_flux_limiter.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <string>

int main() {
  using dec3d::transport::ApplyElectronFluxLimiter;
  using dec3d::transport::ApplyElectronFluxLimiterValues;
  using dec3d::transport::ElectronFluxLimiterFaceInput;
  using dec3d::transport::ElectronFluxLimiterModel;
  using dec3d::transport::ValidateElectronFluxLimiterDiagnostics;

  {
    const ElectronFluxLimiterFaceInput input{
        ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa,
        0.65,
        dec3d::physics::ErgFromKeV(5.0),
        1.0e25,
        1.0e12,
        1.0e8};
    const auto result = ApplyElectronFluxLimiter(input);
    const auto values = ApplyElectronFluxLimiterValues(input);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(values.success);
    DEC3D_CHECK_EQ(values.limited, result.limited);
    DEC3D_CHECK(std::abs(values.scale - result.scale) < 1.0e-15);
    DEC3D_CHECK(std::abs(values.flux_ratio_before_limit - result.flux_ratio_before_limit) <
                std::max(1.0, std::abs(result.flux_ratio_before_limit)) * 1.0e-15);
    DEC3D_CHECK(std::abs(values.effective_kappa_cm_inv_s - result.effective_kappa_cm_inv_s) <
                std::max(1.0, std::abs(result.effective_kappa_cm_inv_s)) * 1.0e-15);
    DEC3D_CHECK(result.scale >= 0.0 && result.scale <= 1.0);
    DEC3D_CHECK(result.effective_kappa_cm_inv_s <= 1.0e12);
    DEC3D_CHECK(result.report_line.find("alpha_e_source=user_supplied") != std::string::npos);
    DEC3D_CHECK(ValidateElectronFluxLimiterDiagnostics(result));
  }

  {
    const auto result = ApplyElectronFluxLimiter(ElectronFluxLimiterFaceInput{
        ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa,
        0.65,
        dec3d::physics::ErgFromKeV(5.0),
        1.0e25,
        1.0,
        1.0e-20});
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(std::abs(result.scale - 1.0) < 1.0e-15);
    DEC3D_CHECK(!result.limited);
  }

  {
    const auto result = ApplyElectronFluxLimiter(ElectronFluxLimiterFaceInput{
        ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa,
        0.0,
        dec3d::physics::ErgFromKeV(1.0),
        1.0e25,
        1.0e12,
        1.0e8});
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_diagnostics.find("alpha_e must be user-supplied and positive") !=
                std::string::npos);
  }

  return 0;
}
