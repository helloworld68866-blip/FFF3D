#pragma once

#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <cstddef>
#include <string>

namespace dec3d::transport {
struct DistributedDiffusionRowOwnership;
}  // namespace dec3d::transport

namespace dec3d::radiation {

enum class RadiationFluxLimiterModel {
  disabled,
  harmonic_eq_5_209,
  larsen_eq_5_210_n2,
  minmax_eq_5_211,
  unsupported
};

struct RadiationFluxLimiterFaceInput {
  RadiationFluxLimiterModel model{RadiationFluxLimiterModel::disabled};
  double unlimited_D_cm2_per_s{0.0};
  double U_face_erg_per_cm3{0.0};
  double abs_grad_U_erg_per_cm4{0.0};
  double U_floor_erg_per_cm3{0.0};
  double grad_floor_erg_per_cm4{0.0};
};

struct RadiationFluxLimiterFaceResult {
  bool success{false};
  bool limited{false};
  double effective_D_cm2_per_s{0.0};
  double limiter_scale{1.0};
  double flux_ratio_before_limit{0.0};
  double R_face_cm_inv{0.0};
  double U_regularized_erg_per_cm3{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
};

struct RadiationFluxLimiterOptions {
  RadiationFluxLimiterModel model{RadiationFluxLimiterModel::disabled};
  double U_floor_erg_per_cm3{0.0};
  double grad_floor_erg_per_cm4{0.0};
};

struct RadiationFluxLimitedFaceCoefficientsResult {
  bool success{false};
  dec3d::transport::GenericDiffusionFaceEffectiveCoefficients face_coefficients;
  std::size_t limited_face_count_local{0u};
  double min_limiter_scale{1.0};
  double max_flux_ratio_before_limit{0.0};
  double min_effective_D_cm2_per_s{0.0};
  double max_effective_D_cm2_per_s{0.0};
  std::string report_line;
  std::string failure_diagnostics;
  std::string failure_reason;
};

[[nodiscard]] const char* RadiationFluxLimiterModelName(
    RadiationFluxLimiterModel model) noexcept;

[[nodiscard]] RadiationFluxLimiterFaceResult ApplyRadiationFluxLimiter(
    const RadiationFluxLimiterFaceInput& input) noexcept;

[[nodiscard]] RadiationFluxLimiterFaceResult ApplyRadiationFluxLimiterValues(
    const RadiationFluxLimiterFaceInput& input) noexcept;

[[nodiscard]] bool ValidateRadiationFluxLimiterDiagnostics(
    const RadiationFluxLimiterFaceResult& result) noexcept;

[[nodiscard]] RadiationFluxLimitedFaceCoefficientsResult
BuildRadiationFluxLimitedFaceCoefficients(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::Array3D<double>& U_old,
    const dec3d::core::Array3D<double>& unlimited_D,
    const RadiationFluxLimiterOptions& options) noexcept;

[[nodiscard]] bool ValidateRadiationFluxLimitedFaceDiagnostics(
    const RadiationFluxLimitedFaceCoefficientsResult& result) noexcept;

}  // namespace dec3d::radiation
