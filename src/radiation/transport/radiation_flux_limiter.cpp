#include "radiation/transport/radiation_flux_limiter.hpp"

#include "physics/units/physical_constants.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace dec3d::radiation {
namespace {

[[nodiscard]] const char* BoolToken(bool value) noexcept {
  return value ? "true" : "false";
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] double ArithmeticMean(double lhs, double rhs) noexcept {
  return 0.5 * (lhs + rhs);
}

[[nodiscard]] double RadialCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

[[nodiscard]] double ThetaCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] double PhiCenter(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t phi) noexcept {
  return 0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
}

[[nodiscard]] double CartesianDistance(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_a,
    std::size_t theta_a,
    std::size_t phi_a,
    std::size_t radial_b,
    std::size_t theta_b,
    std::size_t phi_b) noexcept {
  const double ra = RadialCenter(geometry, radial_a);
  const double rb = RadialCenter(geometry, radial_b);
  const double ta = ThetaCenter(geometry, theta_a);
  const double tb = ThetaCenter(geometry, theta_b);
  const double pa = PhiCenter(geometry, phi_a);
  const double pb = PhiCenter(geometry, phi_b);
  const double xa = ra * std::sin(ta) * std::cos(pa);
  const double ya = ra * std::sin(ta) * std::sin(pa);
  const double za = ra * std::cos(ta);
  const double xb = rb * std::sin(tb) * std::cos(pb);
  const double yb = rb * std::sin(tb) * std::sin(pb);
  const double zb = rb * std::cos(tb);
  const double dx = xa - xb;
  const double dy = ya - yb;
  const double dz = za - zb;
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

void BuildReport(
    RadiationFluxLimiterFaceResult& result,
    RadiationFluxLimiterModel model,
    const char* diagnostic_id) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=" << diagnostic_id
      << "; stage_id=R"
      << "; radiation_flux_limiter_model=" << RadiationFluxLimiterModelName(model)
      << "; radiation_flux_limiter_enabled="
      << BoolToken(model != RadiationFluxLimiterModel::disabled)
      << "; face_limiter_time_level=old_time_lagged"
      << "; U_floor_is_numerical_regularization=true"
      << "; effective_D_cm2_per_s=" << result.effective_D_cm2_per_s
      << "; limiter_scale=" << result.limiter_scale
      << "; limited=" << BoolToken(result.limited)
      << "; flux_ratio_before_limit=" << result.flux_ratio_before_limit
      << "; R_face_cm_inv=" << result.R_face_cm_inv
      << "; U_regularized_erg_per_cm3=" << result.U_regularized_erg_per_cm3
      << "; parity_claim_allowed=false"
      << "; fallback_used=false";
  result.report_line = out.str();
}

[[nodiscard]] RadiationFluxLimiterFaceResult Fail(
    RadiationFluxLimiterModel model,
    const char* reason) {
  RadiationFluxLimiterFaceResult result;
  result.success = false;
  result.failure_reason = reason;
  BuildReport(result, model, "p3.radiation.flux_limiter.failure");
  result.failure_diagnostics = result.report_line + "; failure_reason=" + reason;
  return result;
}

void AccumulateFaceStats(
    RadiationFluxLimitedFaceCoefficientsResult& result,
    const RadiationFluxLimiterFaceResult& face) {
  if (face.limited) {
    ++result.limited_face_count_local;
  }
  result.min_limiter_scale = std::min(result.min_limiter_scale, face.limiter_scale);
  result.max_flux_ratio_before_limit =
      std::max(result.max_flux_ratio_before_limit, face.flux_ratio_before_limit);
  if (result.min_effective_D_cm2_per_s == 0.0 &&
      result.max_effective_D_cm2_per_s == 0.0) {
    result.min_effective_D_cm2_per_s = face.effective_D_cm2_per_s;
    result.max_effective_D_cm2_per_s = face.effective_D_cm2_per_s;
  } else {
    result.min_effective_D_cm2_per_s =
        std::min(result.min_effective_D_cm2_per_s, face.effective_D_cm2_per_s);
    result.max_effective_D_cm2_per_s =
        std::max(result.max_effective_D_cm2_per_s, face.effective_D_cm2_per_s);
  }
}

[[nodiscard]] RadiationFluxLimitedFaceCoefficientsResult FaceArrayFailure(
    const char* reason,
    const RadiationFluxLimiterOptions& options) {
  RadiationFluxLimitedFaceCoefficientsResult result;
  result.success = false;
  result.failure_reason = reason;
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation.flux_limiter.faces.failure"
      << "; stage_id=R"
      << "; radiation_flux_limiter_model="
      << RadiationFluxLimiterModelName(options.model)
      << "; failure_reason=" << reason;
  result.failure_diagnostics = out.str();
  return result;
}

void BuildFaceArrayReport(
    RadiationFluxLimitedFaceCoefficientsResult& result,
    const RadiationFluxLimiterOptions& options) {
  std::ostringstream out;
  out << std::setprecision(17)
      << "diagnostic_id=p3.radiation.flux_limiter.faces"
      << "; stage_id=R"
      << "; radiation_flux_limiter_model="
      << RadiationFluxLimiterModelName(options.model)
      << "; radiation_flux_limiter_enabled="
      << BoolToken(options.model != RadiationFluxLimiterModel::disabled)
      << "; face_limiter_time_level=old_time_lagged"
      << "; U_floor_is_numerical_regularization=true"
      << "; U_floor_erg_per_cm3=" << options.U_floor_erg_per_cm3
      << "; limited_face_count_local=" << result.limited_face_count_local
      << "; min_limiter_scale=" << result.min_limiter_scale
      << "; max_flux_ratio_before_limit=" << result.max_flux_ratio_before_limit
      << "; min_effective_D_cm2_per_s=" << result.min_effective_D_cm2_per_s
      << "; max_effective_D_cm2_per_s=" << result.max_effective_D_cm2_per_s
      << "; face_effective_coefficients_enabled="
      << BoolToken(result.face_coefficients.enabled)
      << "; parity_claim_allowed=false"
      << "; fallback_used=false";
  result.report_line = out.str();
}

struct RadialLimiterHalo {
  bool lower_received{false};
  bool upper_received{false};
  std::vector<double> lower_U;
  std::vector<double> upper_U;
  std::vector<double> lower_D;
  std::vector<double> upper_D;
};

[[nodiscard]] RadialLimiterHalo ExchangeRadialLimiterHalo(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::core::Array3D<double>& U_old,
    const dec3d::core::Array3D<double>& unlimited_D,
    std::size_t local_radial_count) {
  RadialLimiterHalo halo;
  const std::size_t plane_size =
      ownership.global_theta_cells * ownership.global_phi_cells;
  halo.lower_U.assign(plane_size, 0.0);
  halo.upper_U.assign(plane_size, 0.0);
  halo.lower_D.assign(plane_size, 0.0);
  halo.upper_D.assign(plane_size, 0.0);
  if (ownership.rank_count <= 1) {
    return halo;
  }

  const int lower_rank = ownership.rank > 0 ? ownership.rank - 1 : MPI_PROC_NULL;
  const int upper_rank =
      ownership.rank + 1 < ownership.rank_count ? ownership.rank + 1 : MPI_PROC_NULL;
  std::vector<double> send_lower(2u * plane_size, 0.0);
  std::vector<double> send_upper(2u * plane_size, 0.0);
  std::vector<double> recv_lower(2u * plane_size, 0.0);
  std::vector<double> recv_upper(2u * plane_size, 0.0);

  for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
    for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
      const std::size_t slot = t * ownership.global_phi_cells + p;
      send_lower[slot] = U_old(0, t, p);
      send_lower[plane_size + slot] = unlimited_D(0, t, p);
      send_upper[slot] = U_old(local_radial_count - 1u, t, p);
      send_upper[plane_size + slot] = unlimited_D(local_radial_count - 1u, t, p);
    }
  }

  MPI_Status status{};
  MPI_Sendrecv(send_lower.data(),
               static_cast<int>(send_lower.size()),
               MPI_DOUBLE,
               lower_rank,
               2101,
               recv_upper.data(),
               static_cast<int>(recv_upper.size()),
               MPI_DOUBLE,
               upper_rank,
               2101,
               ownership.communicator,
               &status);
  MPI_Sendrecv(send_upper.data(),
               static_cast<int>(send_upper.size()),
               MPI_DOUBLE,
               upper_rank,
               2102,
               recv_lower.data(),
               static_cast<int>(recv_lower.size()),
               MPI_DOUBLE,
               lower_rank,
               2102,
               ownership.communicator,
               &status);

  halo.lower_received = lower_rank != MPI_PROC_NULL;
  halo.upper_received = upper_rank != MPI_PROC_NULL;
  for (std::size_t slot = 0; slot < plane_size; ++slot) {
    halo.lower_U[slot] = recv_lower[slot];
    halo.lower_D[slot] = recv_lower[plane_size + slot];
    halo.upper_U[slot] = recv_upper[slot];
    halo.upper_D[slot] = recv_upper[plane_size + slot];
  }
  return halo;
}

}  // namespace

const char* RadiationFluxLimiterModelName(RadiationFluxLimiterModel model) noexcept {
  switch (model) {
    case RadiationFluxLimiterModel::disabled:
      return "disabled";
    case RadiationFluxLimiterModel::harmonic_eq_5_209:
      return "harmonic_eq_5_209";
    case RadiationFluxLimiterModel::larsen_eq_5_210_n2:
      return "larsen_eq_5_210_n2";
    case RadiationFluxLimiterModel::minmax_eq_5_211:
      return "minmax_eq_5_211";
    case RadiationFluxLimiterModel::unsupported:
      return "unsupported";
  }
  return "unknown";
}

namespace {

RadiationFluxLimiterFaceResult ApplyRadiationFluxLimiterInternal(
    const RadiationFluxLimiterFaceInput& input,
    bool build_success_report) noexcept {
  if (input.model == RadiationFluxLimiterModel::unsupported) {
    return Fail(input.model, "unsupported radiation flux limiter model");
  }
  if (!std::isfinite(input.unlimited_D_cm2_per_s) ||
      input.unlimited_D_cm2_per_s < 0.0) {
    return Fail(input.model,
                "unlimited radiation diffusion coefficient must be non-negative");
  }
  if (!std::isfinite(input.U_face_erg_per_cm3) ||
      !std::isfinite(input.abs_grad_U_erg_per_cm4) ||
      !std::isfinite(input.U_floor_erg_per_cm3) ||
      !std::isfinite(input.grad_floor_erg_per_cm4) ||
      input.abs_grad_U_erg_per_cm4 < 0.0 ||
      input.U_floor_erg_per_cm3 < 0.0 ||
      input.grad_floor_erg_per_cm4 < 0.0) {
    return Fail(input.model,
                "radiation flux limiter input is not finite or non-negative");
  }

  RadiationFluxLimiterFaceResult result;
  result.success = true;
  const double c = dec3d::physics::PhysicsConstantsCGS::speed_of_light_cm_per_s;
  const double u_regularized =
      std::max(input.U_face_erg_per_cm3, input.U_floor_erg_per_cm3);
  result.U_regularized_erg_per_cm3 = u_regularized;
  result.effective_D_cm2_per_s = input.unlimited_D_cm2_per_s;

  if (input.model != RadiationFluxLimiterModel::disabled &&
      input.unlimited_D_cm2_per_s == 0.0) {
    if (build_success_report) {
      BuildReport(result, input.model, "p3.radiation.flux_limiter");
    }
    return result;
  }
  if (input.model != RadiationFluxLimiterModel::disabled && u_regularized <= 0.0) {
    return Fail(input.model,
                "radiation flux limiter U floor regularization is zero");
  }

  const double grad =
      std::max(input.abs_grad_U_erg_per_cm4, input.grad_floor_erg_per_cm4);
  result.R_face_cm_inv =
      input.model == RadiationFluxLimiterModel::disabled ? 0.0 : grad / u_regularized;
  const double free_stream_D =
      grad > 0.0 ? c * u_regularized / grad
                 : std::numeric_limits<double>::infinity();
  result.flux_ratio_before_limit =
      std::isfinite(free_stream_D) && free_stream_D > 0.0
          ? input.unlimited_D_cm2_per_s / free_stream_D
          : 0.0;

  switch (input.model) {
    case RadiationFluxLimiterModel::disabled:
      result.effective_D_cm2_per_s = input.unlimited_D_cm2_per_s;
      break;
    case RadiationFluxLimiterModel::harmonic_eq_5_209: {
      const double sigma = c / input.unlimited_D_cm2_per_s;
      result.effective_D_cm2_per_s = c / (sigma + result.R_face_cm_inv);
      break;
    }
    case RadiationFluxLimiterModel::larsen_eq_5_210_n2: {
      const double sigma = c / input.unlimited_D_cm2_per_s;
      result.effective_D_cm2_per_s =
          c / std::sqrt(sigma * sigma +
                        result.R_face_cm_inv * result.R_face_cm_inv);
      break;
    }
    case RadiationFluxLimiterModel::minmax_eq_5_211:
      result.effective_D_cm2_per_s =
          std::min(input.unlimited_D_cm2_per_s, free_stream_D);
      break;
    case RadiationFluxLimiterModel::unsupported:
      return Fail(input.model, "unsupported radiation flux limiter model");
  }

  if (!std::isfinite(result.effective_D_cm2_per_s) ||
      result.effective_D_cm2_per_s < 0.0) {
    return Fail(input.model,
                "radiation flux limiter produced invalid effective D");
  }
  // The analytic limiter formulas never increase D; enforce that invariant
  // after floating-point evaluation so roundoff cannot turn an unlimited face
  // into a spurious failure.
  result.effective_D_cm2_per_s =
      std::min(result.effective_D_cm2_per_s, input.unlimited_D_cm2_per_s);
  result.limiter_scale =
      input.unlimited_D_cm2_per_s > 0.0
          ? result.effective_D_cm2_per_s / input.unlimited_D_cm2_per_s
          : 1.0;
  result.limited = result.limiter_scale < 1.0 - 1.0e-14;
  if (build_success_report) {
    BuildReport(result, input.model, "p3.radiation.flux_limiter");
  }
  return result;
}

}  // namespace

RadiationFluxLimiterFaceResult ApplyRadiationFluxLimiter(
    const RadiationFluxLimiterFaceInput& input) noexcept {
  return ApplyRadiationFluxLimiterInternal(input, true);
}

RadiationFluxLimiterFaceResult ApplyRadiationFluxLimiterValues(
    const RadiationFluxLimiterFaceInput& input) noexcept {
  return ApplyRadiationFluxLimiterInternal(input, false);
}

bool ValidateRadiationFluxLimiterDiagnostics(
    const RadiationFluxLimiterFaceResult& result) noexcept {
  if (!result.success || result.report_line.empty()) {
    return false;
  }
  return Contains(result.report_line,
                  "diagnostic_id=p3.radiation.flux_limiter") &&
         Contains(result.report_line, "stage_id=R") &&
         Contains(result.report_line, "radiation_flux_limiter_model=") &&
         Contains(result.report_line, "radiation_flux_limiter_enabled=") &&
         Contains(result.report_line,
                  "face_limiter_time_level=old_time_lagged") &&
         Contains(result.report_line,
                  "U_floor_is_numerical_regularization=true") &&
         Contains(result.report_line, "effective_D_cm2_per_s=") &&
         Contains(result.report_line, "limiter_scale=") &&
         Contains(result.report_line, "limited=") &&
         Contains(result.report_line, "flux_ratio_before_limit=") &&
         Contains(result.report_line, "R_face_cm_inv=") &&
         Contains(result.report_line, "parity_claim_allowed=false") &&
         Contains(result.report_line, "fallback_used=false");
}

RadiationFluxLimitedFaceCoefficientsResult BuildRadiationFluxLimitedFaceCoefficients(
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::core::Array3D<double>& U_old,
    const dec3d::core::Array3D<double>& unlimited_D,
    const RadiationFluxLimiterOptions& options) noexcept {
  RadiationFluxLimitedFaceCoefficientsResult result;
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  if (!ownership.success || !geometry.is_valid() || local_radial_count == 0u) {
    return FaceArrayFailure(
        "radiation flux limiter ownership or geometry is invalid", options);
  }
  if (U_old.extent_r() != local_radial_count ||
      U_old.extent_theta() != ownership.global_theta_cells ||
      U_old.extent_phi() != ownership.global_phi_cells ||
      unlimited_D.extent_r() != local_radial_count ||
      unlimited_D.extent_theta() != ownership.global_theta_cells ||
      unlimited_D.extent_phi() != ownership.global_phi_cells) {
    return FaceArrayFailure("radiation flux limiter input shape mismatch", options);
  }
  if (options.model == RadiationFluxLimiterModel::disabled) {
    result.success = true;
    result.face_coefficients.enabled = false;
    BuildFaceArrayReport(result, options);
    return result;
  }

  result.face_coefficients.enabled = true;
  result.face_coefficients.radial_face_D = dec3d::core::Array3D<double>(
      local_radial_count + 1u,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      0.0);
  result.face_coefficients.theta_face_D = dec3d::core::Array3D<double>(
      local_radial_count,
      ownership.global_theta_cells + 1u,
      ownership.global_phi_cells,
      0.0);
  result.face_coefficients.phi_face_D = dec3d::core::Array3D<double>(
      local_radial_count,
      ownership.global_theta_cells,
      ownership.global_phi_cells,
      0.0);

  auto apply = [&](double D, double U_face, double grad) {
    RadiationFluxLimiterFaceInput input;
    input.model = options.model;
    input.unlimited_D_cm2_per_s = D;
    input.U_face_erg_per_cm3 = U_face;
    input.abs_grad_U_erg_per_cm4 = grad;
    input.U_floor_erg_per_cm3 = options.U_floor_erg_per_cm3;
    input.grad_floor_erg_per_cm4 = options.grad_floor_erg_per_cm4;
    return ApplyRadiationFluxLimiterValues(input);
  };

  const auto halo =
      ExchangeRadialLimiterHalo(ownership, U_old, unlimited_D, local_radial_count);

  for (std::size_t face = 0; face <= local_radial_count; ++face) {
    const std::size_t global_face = ownership.global_radial_begin + face;
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        double D = 0.0;
        double U_face = 0.0;
        double grad = 0.0;
        if (global_face == 0u) {
          D = unlimited_D(0, t, p);
          U_face = U_old(0, t, p);
        } else if (global_face == ownership.global_radial_cells) {
          D = unlimited_D(local_radial_count - 1u, t, p);
          U_face = U_old(local_radial_count - 1u, t, p);
        } else {
          const std::size_t slot = t * ownership.global_phi_cells + p;
          double lhs_U = 0.0;
          double rhs_U = 0.0;
          double lhs_D = 0.0;
          double rhs_D = 0.0;
          if (face == 0u) {
            if (!halo.lower_received) {
              return FaceArrayFailure(
                  "radiation flux limiter lower radial halo is missing", options);
            }
            lhs_U = halo.lower_U[slot];
            lhs_D = halo.lower_D[slot];
            rhs_U = U_old(0, t, p);
            rhs_D = unlimited_D(0, t, p);
          } else if (face == local_radial_count) {
            if (!halo.upper_received) {
              return FaceArrayFailure(
                  "radiation flux limiter upper radial halo is missing", options);
            }
            lhs_U = U_old(local_radial_count - 1u, t, p);
            lhs_D = unlimited_D(local_radial_count - 1u, t, p);
            rhs_U = halo.upper_U[slot];
            rhs_D = halo.upper_D[slot];
          } else {
            lhs_U = U_old(face - 1u, t, p);
            lhs_D = unlimited_D(face - 1u, t, p);
            rhs_U = U_old(face, t, p);
            rhs_D = unlimited_D(face, t, p);
          }
          D = ArithmeticMean(lhs_D, rhs_D);
          U_face = ArithmeticMean(lhs_U, rhs_U);
          const double distance =
              RadialCenter(geometry, global_face) -
              RadialCenter(geometry, global_face - 1u);
          grad = distance > 0.0 ? std::abs(rhs_U - lhs_U) / distance : 0.0;
        }
        const auto limited = apply(D, U_face, grad);
        if (!limited.success) {
          return FaceArrayFailure(limited.failure_reason.c_str(), options);
        }
        result.face_coefficients.radial_face_D(face, t, p) =
            limited.effective_D_cm2_per_s;
        AccumulateFaceStats(result, limited);
      }
    }
  }

  for (std::size_t r = 0; r < local_radial_count; ++r) {
    const std::size_t gr = ownership.global_radial_begin + r;
    for (std::size_t theta_face = 0; theta_face <= ownership.global_theta_cells;
         ++theta_face) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        double D = 0.0;
        double U_face = 0.0;
        double grad = 0.0;
        if (theta_face == 0u || theta_face == ownership.global_theta_cells) {
          const std::size_t theta_cell =
              theta_face == 0u ? 0u : ownership.global_theta_cells - 1u;
          D = unlimited_D(r, theta_cell, p);
          U_face = U_old(r, theta_cell, p);
        } else {
          D = ArithmeticMean(
              unlimited_D(r, theta_face - 1u, p), unlimited_D(r, theta_face, p));
          U_face = ArithmeticMean(
              U_old(r, theta_face - 1u, p), U_old(r, theta_face, p));
          const double distance =
              RadialCenter(geometry, gr) *
              (ThetaCenter(geometry, theta_face) -
               ThetaCenter(geometry, theta_face - 1u));
          grad = distance > 0.0
                     ? std::abs(U_old(r, theta_face, p) -
                                U_old(r, theta_face - 1u, p)) /
                           distance
                     : 0.0;
        }
        const auto limited = apply(D, U_face, grad);
        if (!limited.success) {
          return FaceArrayFailure(limited.failure_reason.c_str(), options);
        }
        result.face_coefficients.theta_face_D(r, theta_face, p) =
            limited.effective_D_cm2_per_s;
        AccumulateFaceStats(result, limited);
      }
    }
  }

  for (std::size_t r = 0; r < local_radial_count; ++r) {
    const std::size_t gr = ownership.global_radial_begin + r;
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const std::size_t next = (p + 1u) % ownership.global_phi_cells;
        const double D = ArithmeticMean(unlimited_D(r, t, p), unlimited_D(r, t, next));
        const double U_face = ArithmeticMean(U_old(r, t, p), U_old(r, t, next));
        const double distance = CartesianDistance(geometry, gr, t, p, gr, t, next);
        const double grad =
            distance > 0.0 ? std::abs(U_old(r, t, next) - U_old(r, t, p)) / distance
                           : 0.0;
        const auto limited = apply(D, U_face, grad);
        if (!limited.success) {
          return FaceArrayFailure(limited.failure_reason.c_str(), options);
        }
        result.face_coefficients.phi_face_D(r, t, p) =
            limited.effective_D_cm2_per_s;
        AccumulateFaceStats(result, limited);
      }
    }
  }

  result.success = true;
  BuildFaceArrayReport(result, options);
  return result;
}

bool ValidateRadiationFluxLimitedFaceDiagnostics(
    const RadiationFluxLimitedFaceCoefficientsResult& result) noexcept {
  if (!result.success || result.report_line.empty()) {
    return false;
  }
  return Contains(result.report_line,
                  "diagnostic_id=p3.radiation.flux_limiter.faces") &&
         Contains(result.report_line, "stage_id=R") &&
         Contains(result.report_line, "radiation_flux_limiter_model=") &&
         Contains(result.report_line, "radiation_flux_limiter_enabled=") &&
         Contains(result.report_line,
                  "face_limiter_time_level=old_time_lagged") &&
         Contains(result.report_line,
                  "U_floor_is_numerical_regularization=true") &&
         Contains(result.report_line, "limited_face_count_local=") &&
         Contains(result.report_line, "min_limiter_scale=") &&
         Contains(result.report_line, "max_flux_ratio_before_limit=") &&
         Contains(result.report_line, "face_effective_coefficients_enabled=") &&
         Contains(result.report_line, "parity_claim_allowed=false") &&
         Contains(result.report_line, "fallback_used=false");
}

}  // namespace dec3d::radiation
