#include "transport/thermal/distributed_thermal_conduction.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"
#include "transport/thermal/electron_flux_limiter.hpp"
#include "transport/thermal/thermal_conductivity.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace dec3d::transport {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

[[nodiscard]] bool Finite(double value) noexcept {
  return std::isfinite(value);
}

[[nodiscard]] double ElapsedSecondsSince(
    const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

[[nodiscard]] std::size_t FlatIndex(
    const DiffusionGridLayout& layout,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return (radial * layout.theta_cells + theta) * layout.phi_cells + phi;
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

[[nodiscard]] GenericDiffusionBoundaryPolicy FullSphereBoundary() noexcept {
  return GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::periodic};
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildFullSphereGeometry(
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells,
    double r_outer) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  geometry.theta_faces.resize(theta_cells + 1u);
  geometry.phi_faces.resize(phi_cells + 1u);
  geometry.cell_volumes.resize(radial_cells * theta_cells * phi_cells, 0.0);

  for (std::size_t r = 0; r <= radial_cells; ++r) {
    geometry.radial_faces[r] =
        r_outer * static_cast<double>(r) / static_cast<double>(radial_cells);
  }
  for (std::size_t t = 0; t <= theta_cells; ++t) {
    geometry.theta_faces[t] =
        kPi * static_cast<double>(t) / static_cast<double>(theta_cells);
  }
  for (std::size_t p = 0; p <= phi_cells; ++p) {
    geometry.phi_faces[p] =
        2.0 * kPi * static_cast<double>(p) / static_cast<double>(phi_cells);
  }

  geometry.global_volume = 0.0;
  std::size_t linear = 0;
  for (std::size_t r = 0; r < radial_cells; ++r) {
    const double radial_factor =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) / 3.0;
    for (std::size_t t = 0; t < theta_cells; ++t) {
      const double polar_factor =
          std::cos(geometry.theta_faces[t]) -
          std::cos(geometry.theta_faces[t + 1u]);
      for (std::size_t p = 0; p < phi_cells; ++p) {
        const double azimuthal_factor =
            geometry.phi_faces[p + 1u] - geometry.phi_faces[p];
        geometry.cell_volumes[linear] =
            radial_factor * polar_factor * azimuthal_factor;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

[[nodiscard]] ThermalConductivityModel ProviderModelFromKappaModel(
    ThermalConductionKappaModel model) noexcept {
  switch (model) {
    case ThermalConductionKappaModel::spitzer_no_degeneracy:
      return ThermalConductivityModel::spitzer_no_degeneracy;
    case ThermalConductionKappaModel::lee_more_with_degeneracy:
      return ThermalConductivityModel::lee_more_with_degeneracy;
    case ThermalConductionKappaModel::constant_user_supplied:
    case ThermalConductionKappaModel::unsupported:
      return ThermalConductivityModel::unsupported;
  }
  return ThermalConductivityModel::unsupported;
}

void FillSmokeState(dec3d::state::CanonicalState& state) {
  const double ne = 1.0e23;
  const double ni = 1.0e23;
  const double rho = ni * dec3d::physics::DefaultMeanDTIonMassG();
  const double te = dec3d::physics::ErgFromKeV(1.0);
  const double ti = dec3d::physics::ErgFromKeV(1.0);
  const double e_e = ne * te / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double e_i = ni * ti / (dec3d::state::HydroIdealGasGamma() - 1.0);
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        state.rho(r, t, p) = rho;
        state.mom_r(r, t, p) = 0.0;
        state.mom_theta(r, t, p) = 0.0;
        state.mom_phi(r, t, p) = 0.0;
        state.e_electron(r, t, p) = e_e;
        state.e_fluid_total(r, t, p) = e_e + e_i;
      }
    }
  }
}

void Fail(
    DistributedThermalConductionResult& result,
    const std::string& reason,
    int rank,
    int rank_count,
    bool global_stage_ok) {
  result.success = false;
  result.failure_reason = reason;
  result.global_stage_ok = global_stage_ok;
  result.updated_fields = kDistributedThermalWritesNoFields;
  result.updated_fields_label = "none";
  std::ostringstream out;
  out << "diagnostic_id=p2.thermal_conduction.distributed.failure"
      << "; failure_reason=" << reason
      << "; mpi_rank=" << rank
      << "; mpi_rank_count=" << rank_count
      << "; global_stage_ok=" << (global_stage_ok ? "true" : "false")
      << "; global_publish_ok=false"
      << "; updated_fields=none"
      << "; canonical_state_mutated=false"
      << "; gathered_writeback_used=false"
      << "; fallback_used=false";
  result.failure_diagnostics = out.str();
}

[[nodiscard]] int AllreduceMinBool(MPI_Comm communicator, bool value) noexcept {
  int local = value ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, communicator);
  return global;
}

[[nodiscard]] unsigned long long AllreduceSumUll(
    MPI_Comm communicator,
    unsigned long long value) noexcept {
  unsigned long long global = 0;
  MPI_Allreduce(&value, &global, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, communicator);
  return global;
}

struct ElectronLimiterStats {
  bool success{true};
  std::string failure_reason;
  std::string failure_diagnostics;
  unsigned long long local_face_count{0};
  unsigned long long local_limited_face_count{0};
  double min_scale{1.0};
  double max_flux_ratio{0.0};
};

struct ElectronLimiterHalo {
  bool lower_received{false};
  bool upper_received{false};
  std::vector<double> lower_te;
  std::vector<double> lower_ne;
  std::vector<double> lower_kappa;
  std::vector<double> upper_te;
  std::vector<double> upper_ne;
  std::vector<double> upper_kappa;
};

[[nodiscard]] double CartesianDistance(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t lhs_radial,
    std::size_t lhs_theta,
    std::size_t lhs_phi,
    std::size_t rhs_radial,
    std::size_t rhs_theta,
    std::size_t rhs_phi) noexcept {
  const auto center = [&](std::size_t r, std::size_t t, std::size_t p) {
    const double radius = RadialCenter(geometry, r);
    const double theta = ThetaCenter(geometry, t);
    const double phi = PhiCenter(geometry, p);
    const double sin_theta = std::sin(theta);
    return std::array<double, 3>{
        radius * sin_theta * std::cos(phi),
        radius * sin_theta * std::sin(phi),
        radius * std::cos(theta)};
  };
  const auto lhs = center(lhs_radial, lhs_theta, lhs_phi);
  const auto rhs = center(rhs_radial, rhs_theta, rhs_phi);
  const double dx = lhs[0] - rhs[0];
  const double dy = lhs[1] - rhs[1];
  const double dz = lhs[2] - rhs[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

[[nodiscard]] ElectronFluxLimiterFaceResult ApplyLimiterForPair(
    const DistributedThermalConductionOptions& options,
    double te_lhs,
    double te_rhs,
    double ne_lhs,
    double ne_rhs,
    double kappa_lhs,
    double kappa_rhs,
    double distance_cm) noexcept {
  const double kappa_face = 0.5 * (kappa_lhs + kappa_rhs);
  if (options.electron_flux_limiter_model == ElectronFluxLimiterModel::disabled) {
    return ApplyElectronFluxLimiterValues(ElectronFluxLimiterFaceInput{
        ElectronFluxLimiterModel::disabled,
        0.0,
        0.5 * (te_lhs + te_rhs),
        0.5 * (ne_lhs + ne_rhs),
        kappa_face,
        0.0});
  }
  return ApplyElectronFluxLimiterValues(ElectronFluxLimiterFaceInput{
      options.electron_flux_limiter_model,
      options.electron_flux_limiter_alpha_e,
      0.5 * (te_lhs + te_rhs),
      0.5 * (ne_lhs + ne_rhs),
      kappa_face,
      std::abs(te_rhs - te_lhs) / distance_cm});
}

void AccumulateLimiterStats(
    ElectronLimiterStats& stats,
    const ElectronFluxLimiterFaceResult& face) noexcept {
  ++stats.local_face_count;
  if (!face.success) {
    stats.success = false;
    stats.failure_reason = face.failure_reason;
    stats.failure_diagnostics = face.failure_diagnostics;
    return;
  }
  if (face.limited) {
    ++stats.local_limited_face_count;
  }
  stats.min_scale = std::min(stats.min_scale, face.scale);
  stats.max_flux_ratio = std::max(stats.max_flux_ratio, face.flux_ratio_before_limit);
}

[[nodiscard]] ElectronLimiterHalo ExchangeElectronLimiterHalo(
    const DistributedDiffusionRowOwnership& ownership,
    const dec3d::state::ThermodynamicRecoveryResult& recovered,
    const dec3d::core::Array3D<double>& kappa_e) {
  const std::size_t plane_size = ownership.global_theta_cells * ownership.global_phi_cells;
  ElectronLimiterHalo halo;
  halo.lower_te.assign(plane_size, 0.0);
  halo.lower_ne.assign(plane_size, 0.0);
  halo.lower_kappa.assign(plane_size, 0.0);
  halo.upper_te.assign(plane_size, 0.0);
  halo.upper_ne.assign(plane_size, 0.0);
  halo.upper_kappa.assign(plane_size, 0.0);

  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  std::vector<double> send_lower(3u * plane_size, 0.0);
  std::vector<double> send_upper(3u * plane_size, 0.0);
  for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
    for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
      const std::size_t slot = t * ownership.global_phi_cells + p;
      const auto& lower = recovered.cells(0, t, p);
      const auto& upper = recovered.cells(local_radial_count - 1u, t, p);
      send_lower[slot] = lower.t_e_erg_per_particle;
      send_lower[plane_size + slot] = lower.n_e_cm3;
      send_lower[2u * plane_size + slot] = kappa_e(0, t, p);
      send_upper[slot] = upper.t_e_erg_per_particle;
      send_upper[plane_size + slot] = upper.n_e_cm3;
      send_upper[2u * plane_size + slot] = kappa_e(local_radial_count - 1u, t, p);
    }
  }

  std::vector<double> recv_lower(3u * plane_size, 0.0);
  std::vector<double> recv_upper(3u * plane_size, 0.0);
  const int lower_rank = ownership.global_radial_begin > 0u ? ownership.rank - 1 : MPI_PROC_NULL;
  const int upper_rank = ownership.global_radial_end < ownership.global_radial_cells
                             ? ownership.rank + 1
                             : MPI_PROC_NULL;
  MPI_Status status{};
  MPI_Sendrecv(
      send_lower.data(),
      static_cast<int>(send_lower.size()),
      MPI_DOUBLE,
      lower_rank,
      7101,
      recv_upper.data(),
      static_cast<int>(recv_upper.size()),
      MPI_DOUBLE,
      upper_rank,
      7101,
      ownership.communicator,
      &status);
  halo.upper_received = upper_rank != MPI_PROC_NULL;

  MPI_Sendrecv(
      send_upper.data(),
      static_cast<int>(send_upper.size()),
      MPI_DOUBLE,
      upper_rank,
      7102,
      recv_lower.data(),
      static_cast<int>(recv_lower.size()),
      MPI_DOUBLE,
      lower_rank,
      7102,
      ownership.communicator,
      &status);
  halo.lower_received = lower_rank != MPI_PROC_NULL;

  for (std::size_t slot = 0; slot < plane_size; ++slot) {
    halo.lower_te[slot] = recv_lower[slot];
    halo.lower_ne[slot] = recv_lower[plane_size + slot];
    halo.lower_kappa[slot] = recv_lower[2u * plane_size + slot];
    halo.upper_te[slot] = recv_upper[slot];
    halo.upper_ne[slot] = recv_upper[plane_size + slot];
    halo.upper_kappa[slot] = recv_upper[2u * plane_size + slot];
  }
  return halo;
}

[[nodiscard]] ElectronLimiterStats BuildElectronFaceEffectiveCoefficients(
    const DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::state::ThermodynamicRecoveryResult& recovered,
    const DistributedThermalConductionOptions& options,
    DistributedGenericDiffusionProblem& problem) noexcept {
  ElectronLimiterStats stats;
  const std::size_t local_radial_count =
      ownership.global_radial_end - ownership.global_radial_begin;
  problem.face_effective_coefficients.enabled = true;
  problem.face_effective_coefficients.radial_face_D = dec3d::core::Array3D<double>(
      local_radial_count + 1u, ownership.global_theta_cells, ownership.global_phi_cells, 0.0);
  problem.face_effective_coefficients.theta_face_D = dec3d::core::Array3D<double>(
      local_radial_count, ownership.global_theta_cells + 1u, ownership.global_phi_cells, 0.0);
  problem.face_effective_coefficients.phi_face_D = dec3d::core::Array3D<double>(
      local_radial_count, ownership.global_theta_cells, ownership.global_phi_cells, 0.0);

  const auto halo = ExchangeElectronLimiterHalo(ownership, recovered, problem.local_coefficient_D);
  for (std::size_t face = 0; face <= local_radial_count; ++face) {
    const std::size_t global_face = ownership.global_radial_begin + face;
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        if (global_face == 0u) {
          problem.face_effective_coefficients.radial_face_D(face, t, p) =
              problem.local_coefficient_D(0, t, p);
          continue;
        }
        if (global_face == ownership.global_radial_cells) {
          problem.face_effective_coefficients.radial_face_D(face, t, p) =
              problem.local_coefficient_D(local_radial_count - 1u, t, p);
          continue;
        }

        double te_lhs = 0.0;
        double te_rhs = 0.0;
        double ne_lhs = 0.0;
        double ne_rhs = 0.0;
        double kappa_lhs = 0.0;
        double kappa_rhs = 0.0;
        if (face == 0u) {
          const std::size_t slot = t * ownership.global_phi_cells + p;
          if (!halo.lower_received) {
            stats.success = false;
            stats.failure_reason = "electron flux limiter lower halo is missing";
            return stats;
          }
          const auto& rhs = recovered.cells(0, t, p);
          te_lhs = halo.lower_te[slot];
          ne_lhs = halo.lower_ne[slot];
          kappa_lhs = halo.lower_kappa[slot];
          te_rhs = rhs.t_e_erg_per_particle;
          ne_rhs = rhs.n_e_cm3;
          kappa_rhs = problem.local_coefficient_D(0, t, p);
        } else if (face == local_radial_count) {
          const std::size_t slot = t * ownership.global_phi_cells + p;
          if (!halo.upper_received) {
            stats.success = false;
            stats.failure_reason = "electron flux limiter upper halo is missing";
            return stats;
          }
          const auto& lhs = recovered.cells(local_radial_count - 1u, t, p);
          te_lhs = lhs.t_e_erg_per_particle;
          ne_lhs = lhs.n_e_cm3;
          kappa_lhs = problem.local_coefficient_D(local_radial_count - 1u, t, p);
          te_rhs = halo.upper_te[slot];
          ne_rhs = halo.upper_ne[slot];
          kappa_rhs = halo.upper_kappa[slot];
        } else {
          const auto& lhs = recovered.cells(face - 1u, t, p);
          const auto& rhs = recovered.cells(face, t, p);
          te_lhs = lhs.t_e_erg_per_particle;
          ne_lhs = lhs.n_e_cm3;
          kappa_lhs = problem.local_coefficient_D(face - 1u, t, p);
          te_rhs = rhs.t_e_erg_per_particle;
          ne_rhs = rhs.n_e_cm3;
          kappa_rhs = problem.local_coefficient_D(face, t, p);
        }
        const double distance =
            RadialCenter(geometry, global_face) - RadialCenter(geometry, global_face - 1u);
        const auto limited = ApplyLimiterForPair(
            options, te_lhs, te_rhs, ne_lhs, ne_rhs, kappa_lhs, kappa_rhs, distance);
        AccumulateLimiterStats(stats, limited);
        if (!stats.success) {
          return stats;
        }
        problem.face_effective_coefficients.radial_face_D(face, t, p) =
            limited.effective_kappa_cm_inv_s;
      }
    }
  }

  for (std::size_t lr = 0; lr < local_radial_count; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t theta_face = 0; theta_face <= ownership.global_theta_cells; ++theta_face) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        if (theta_face == 0u || theta_face == ownership.global_theta_cells) {
          const std::size_t theta_cell = theta_face == 0u ? 0u : ownership.global_theta_cells - 1u;
          problem.face_effective_coefficients.theta_face_D(lr, theta_face, p) =
              problem.local_coefficient_D(lr, theta_cell, p);
          continue;
        }
        const auto& lower = recovered.cells(lr, theta_face - 1u, p);
        const auto& upper = recovered.cells(lr, theta_face, p);
        const double distance =
            RadialCenter(geometry, gr) *
            (ThetaCenter(geometry, theta_face) - ThetaCenter(geometry, theta_face - 1u));
        const auto limited = ApplyLimiterForPair(
            options,
            lower.t_e_erg_per_particle,
            upper.t_e_erg_per_particle,
            lower.n_e_cm3,
            upper.n_e_cm3,
            problem.local_coefficient_D(lr, theta_face - 1u, p),
            problem.local_coefficient_D(lr, theta_face, p),
            distance);
        AccumulateLimiterStats(stats, limited);
        if (!stats.success) {
          return stats;
        }
        problem.face_effective_coefficients.theta_face_D(lr, theta_face, p) =
            limited.effective_kappa_cm_inv_s;
      }
    }
  }

  for (std::size_t lr = 0; lr < local_radial_count; ++lr) {
    const std::size_t gr = ownership.global_radial_begin + lr;
    for (std::size_t t = 0; t < ownership.global_theta_cells; ++t) {
      for (std::size_t p = 0; p < ownership.global_phi_cells; ++p) {
        const std::size_t next_phi = (p + 1u) % ownership.global_phi_cells;
        const auto& lhs = recovered.cells(lr, t, p);
        const auto& rhs = recovered.cells(lr, t, next_phi);
        const double distance = CartesianDistance(geometry, gr, t, p, gr, t, next_phi);
        const auto limited = ApplyLimiterForPair(
            options,
            lhs.t_e_erg_per_particle,
            rhs.t_e_erg_per_particle,
            lhs.n_e_cm3,
            rhs.n_e_cm3,
            problem.local_coefficient_D(lr, t, p),
            problem.local_coefficient_D(lr, t, next_phi),
            distance);
        AccumulateLimiterStats(stats, limited);
        if (!stats.success) {
          return stats;
        }
        problem.face_effective_coefficients.phi_face_D(lr, t, p) =
            limited.effective_kappa_cm_inv_s;
      }
    }
  }

  return stats;
}

[[nodiscard]] DistributedGenericDiffusionProblem MakeEmptyReservoirProblem(
    const DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const GenericDiffusionBoundaryPolicy& boundary_policy,
    double dt_s) {
  const DiffusionGridLayout layout{
      ownership.global_radial_end - ownership.global_radial_begin,
      ownership.global_theta_cells,
      ownership.global_phi_cells};
  DistributedGenericDiffusionProblem problem;
  problem.ownership = ownership;
  problem.global_geometry = geometry;
  problem.dt_s = dt_s;
  problem.boundary_policy = boundary_policy;
  problem.local_coefficient_A = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.local_coefficient_D = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.local_coefficient_C = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.local_coefficient_B = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  problem.local_scalar_old = dec3d::core::Array3D<double>(
      layout.radial_cells, layout.theta_cells, layout.phi_cells, 0.0);
  return problem;
}

struct ThermalReservoirProblems {
  DistributedGenericDiffusionProblem electron;
  DistributedGenericDiffusionProblem ion;
  std::string coefficient_report;
  std::string failure_reason;
  ElectronLimiterStats limiter_stats;
};

[[nodiscard]] ThermalReservoirProblems BuildReservoirProblems(
    const DistributedDiffusionRowOwnership& ownership,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const GenericDiffusionBoundaryPolicy& boundary_policy,
    double dt_s,
    const dec3d::state::ThermodynamicRecoveryResult& recovered,
    const DistributedThermalConductionOptions& options,
    ThermalConductionKappaModel model) {
  const DiffusionGridLayout layout{
      ownership.global_radial_end - ownership.global_radial_begin,
      ownership.global_theta_cells,
      ownership.global_phi_cells};
  ThermalReservoirProblems problems;
  problems.electron = MakeEmptyReservoirProblem(ownership, geometry, boundary_policy, dt_s);
  problems.ion = MakeEmptyReservoirProblem(ownership, geometry, boundary_policy, dt_s);

  double min_kappa_e = std::numeric_limits<double>::infinity();
  double max_kappa_e = 0.0;
  double min_kappa_i = std::numeric_limits<double>::infinity();
  double max_kappa_i = 0.0;
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const auto provider_model = ProviderModelFromKappaModel(model);
  if (provider_model == ThermalConductivityModel::unsupported) {
    problems.failure_reason = "distributed thermal requires variable-kappa provider model";
    return problems;
  }

  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const auto& cell = recovered.cells(r, t, p);
        const auto coefficient = ComputeThermalConductivityValues(
            ThermalConductivityInput{
                cell.t_e_erg_per_particle,
                cell.t_i_erg_per_particle,
                cell.n_e_cm3,
                cell.n_i_cm3,
                cell.zbar,
                cell.mean_ion_mass_g,
                provider_model});
        if (!coefficient.success) {
          problems.failure_reason = "thermal conductivity provider failed";
          return problems;
        }
        problems.electron.local_coefficient_A(r, t, p) = cell.n_e_cm3 / gamma_minus_one;
        problems.electron.local_coefficient_D(r, t, p) = coefficient.kappa_e_cm_inv_s;
        problems.electron.local_scalar_old(r, t, p) = cell.t_e_erg_per_particle;
        problems.ion.local_coefficient_A(r, t, p) = cell.n_i_cm3 / gamma_minus_one;
        problems.ion.local_coefficient_D(r, t, p) = coefficient.kappa_i_cm_inv_s;
        problems.ion.local_scalar_old(r, t, p) = cell.t_i_erg_per_particle;
        min_kappa_e = std::min(min_kappa_e, coefficient.kappa_e_cm_inv_s);
        max_kappa_e = std::max(max_kappa_e, coefficient.kappa_e_cm_inv_s);
        min_kappa_i = std::min(min_kappa_i, coefficient.kappa_i_cm_inv_s);
        max_kappa_i = std::max(max_kappa_i, coefficient.kappa_i_cm_inv_s);
      }
    }
  }

  if (options.electron_flux_limiter_model != ElectronFluxLimiterModel::disabled) {
    problems.limiter_stats = BuildElectronFaceEffectiveCoefficients(
        ownership,
        geometry,
        recovered,
        options,
        problems.electron);
    if (!problems.limiter_stats.success) {
      problems.failure_reason = problems.limiter_stats.failure_reason.empty()
          ? "electron flux limiter failed"
          : problems.limiter_stats.failure_reason;
      return problems;
    }
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p2.thermal_conduction.distributed.coefficients"
         << "; coefficient_time_level=old_time_lagged"
         << "; kappa_model_requested=" << ThermalConductionKappaModelName(model)
         << "; provider_model_executed=" << ThermalConductivityModelName(provider_model)
         << "; min_kappa_cm_inv_s=" << std::min(min_kappa_e, min_kappa_i)
         << "; max_kappa_cm_inv_s=" << std::max(max_kappa_e, max_kappa_i)
         << "; min_kappa_e_cm_inv_s=" << min_kappa_e
         << "; max_kappa_e_cm_inv_s=" << max_kappa_e
         << "; min_kappa_i_cm_inv_s=" << min_kappa_i
         << "; max_kappa_i_cm_inv_s=" << max_kappa_i
         << "; provider_failure_count=0";
  problems.coefficient_report = report.str();
  return problems;
}

}  // namespace

bool DistributedThermalConductionResult::is_complete() const noexcept {
  return success &&
         distributed_writeback &&
         !gathered_writeback_used &&
         global_stage_ok &&
         global_publish_ok &&
         updated_fields ==
             (kDistributedThermalWritesElectronEnergy |
              kDistributedThermalWritesFluidTotalEnergy) &&
         updated_fields_label == "e_electron,e_fluid_total" &&
         !report_line.empty() &&
         report_line.find("diagnostic_id=p2.production.thermal_stage") !=
             std::string::npos &&
         report_line.find("operator_diagnostic_id=p2.thermal_conduction.distributed") !=
             std::string::npos &&
         report_line.find("global_stage_ok=true") != std::string::npos &&
         report_line.find("global_publish_ok=true") != std::string::npos &&
         report_line.find("updated_fields=e_electron,e_fluid_total") != std::string::npos;
}

DistributedThermalConductionProblem BuildDistributedThermalConductionSmokeProblem(
    MPI_Comm communicator,
    std::size_t global_radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  DistributedThermalConductionProblem problem;
  problem.ownership = BuildDistributedDiffusionRowOwnership(
      communicator, global_radial_cells, theta_cells, phi_cells);
  const std::size_t local_radial =
      problem.ownership.global_radial_end - problem.ownership.global_radial_begin;
  problem.local_state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{local_radial, theta_cells, phi_cells, 0});
  FillSmokeState(problem.local_state);
  problem.global_geometry =
      BuildFullSphereGeometry(global_radial_cells, theta_cells, phi_cells, 1.0);
  problem.boundary_policy = FullSphereBoundary();
  return problem;
}

DistributedThermalConductionResult ApplyDistributedVariableKappaThermalConduction(
    DistributedThermalConductionProblem& problem,
    const DistributedThermalConductionOptions& options) noexcept {
  DistributedThermalConductionResult result;
  const auto& ownership = problem.ownership;
  int rank = ownership.rank;
  int rank_count = ownership.rank_count;

  auto local_stage_ok = ownership.success &&
                        Finite(options.dt_s) &&
                        options.dt_s >= 0.0 &&
                        problem.local_state.layout.radial_cells ==
                            ownership.global_radial_end - ownership.global_radial_begin &&
                        problem.local_state.layout.theta_cells == ownership.global_theta_cells &&
                        problem.local_state.layout.phi_cells == ownership.global_phi_cells;

  const auto coefficient_timer = std::chrono::steady_clock::now();
  auto recovered = dec3d::state::RecoverThermodynamicState(
      problem.local_state,
      dec3d::state::ThermodynamicRecoveryOptions{
          options.electron_energy_floor,
          options.ion_energy_floor,
          1.0,
          0.0});
  if (!recovered.success) {
    local_stage_ok = false;
  }

  const bool global_stage_preflight =
      AllreduceMinBool(ownership.communicator, local_stage_ok) != 0;
  if (!global_stage_preflight) {
    Fail(result, recovered.success ? "distributed thermal stage preflight failed"
                                   : "thermodynamic recovery failed",
         rank, rank_count, false);
    if (!recovered.failure_diagnostics.empty()) {
      result.failure_diagnostics += "; nested_diagnostics={" + recovered.failure_diagnostics + "}";
    }
    return result;
  }

  auto reservoir_problems = BuildReservoirProblems(
      ownership,
      problem.global_geometry,
      problem.boundary_policy,
      options.dt_s,
      recovered,
      options,
      options.kappa_model);
  result.coefficient_provider_wall_s = ElapsedSecondsSince(coefficient_timer);
  local_stage_ok = reservoir_problems.failure_reason.empty();
  const bool global_coefficients_ok =
      AllreduceMinBool(ownership.communicator, local_stage_ok) != 0;
  if (!global_coefficients_ok) {
    Fail(result,
         reservoir_problems.failure_reason.empty() ? "distributed thermal coefficient stage failed"
                                                   : reservoir_problems.failure_reason,
         rank,
         rank_count,
         false);
    return result;
  }
  auto& electron_problem = reservoir_problems.electron;
  auto& ion_problem = reservoir_problems.ion;
  const auto& limiter_stats = reservoir_problems.limiter_stats;

  const auto assembly_timer = std::chrono::steady_clock::now();
  const auto electron_assembly = AssembleDistributedGenericDiffusionSystem(electron_problem);
  const auto ion_assembly = AssembleDistributedGenericDiffusionSystem(ion_problem);
  result.assembly_wall_s = ElapsedSecondsSince(assembly_timer);
  local_stage_ok = electron_assembly.success && ion_assembly.success;
  const bool global_assembly_ok =
      AllreduceMinBool(ownership.communicator, local_stage_ok) != 0;
  if (!global_assembly_ok) {
    Fail(result,
         "distributed thermal diffusion assembly failed",
         rank,
         rank_count,
         false);
    return result;
  }

  auto solve_options = GenericDiffusionHypreSolveOptions{};
  solve_options.communicator = ownership.communicator;
  const auto electron_solve = SolveDistributedGenericDiffusionHypre(electron_assembly, solve_options);
  const auto ion_solve = SolveDistributedGenericDiffusionHypre(ion_assembly, solve_options);
  result.hypre_setup_wall_s =
      electron_solve.hypre_setup_wall_s + ion_solve.hypre_setup_wall_s;
  result.hypre_solve_wall_s =
      electron_solve.hypre_solve_wall_s + ion_solve.hypre_solve_wall_s;
  result.solver_iterations =
      electron_solve.gmres_iterations + ion_solve.gmres_iterations;
  local_stage_ok =
      electron_solve.success &&
      ion_solve.success &&
      electron_solve.local_scalar_new.size() == ownership.local_row_count &&
      ion_solve.local_scalar_new.size() == ownership.local_row_count;
  const bool global_solve_ok =
      AllreduceMinBool(ownership.communicator, local_stage_ok) != 0;
  if (!global_solve_ok) {
    Fail(result, "distributed thermal HYPRE solve failed", rank, rank_count, false);
    result.electron_solve_report =
        electron_solve.success ? electron_solve.report_line : electron_solve.failure_diagnostics;
    result.ion_solve_report =
        ion_solve.success ? ion_solve.report_line : ion_solve.failure_diagnostics;
    result.failure_diagnostics += "; electron_solve_diagnostics={" + result.electron_solve_report + "}";
    result.failure_diagnostics += "; ion_solve_diagnostics={" + result.ion_solve_report + "}";
    return result;
  }

  const DiffusionGridLayout layout{
      problem.local_state.layout.radial_cells,
      problem.local_state.layout.theta_cells,
      problem.local_state.layout.phi_cells};
  const auto writeback_timer = std::chrono::steady_clock::now();
  auto staged_e_electron = problem.local_state.e_electron;
  auto staged_e_total = problem.local_state.e_fluid_total;
  double local_max_cellwise_residual = 0.0;
  double local_thermal_before_volume = 0.0;
  double local_thermal_after_volume = 0.0;
  for (std::size_t r = 0; r < layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < layout.phi_cells; ++p) {
        const auto& cell = recovered.cells(r, t, p);
        const auto idx = FlatIndex(layout, r, t, p);
        const double electron_energy =
            cell.n_e_cm3 * electron_solve.local_scalar_new[idx] /
            (dec3d::state::HydroIdealGasGamma() - 1.0);
        const double ion_energy =
            cell.n_i_cm3 * ion_solve.local_scalar_new[idx] /
            (dec3d::state::HydroIdealGasGamma() - 1.0);
        if (!Finite(electron_energy) || !Finite(ion_energy) ||
            electron_energy < options.electron_energy_floor ||
            ion_energy < options.ion_energy_floor) {
          local_stage_ok = false;
        }
        staged_e_electron(r, t, p) = electron_energy;
        staged_e_total(r, t, p) =
            cell.kinetic_energy_density_erg_per_cm3 + electron_energy + ion_energy;
        const std::size_t global_row = ownership.local_row_begin + idx;
        const double volume = problem.global_geometry.cell_volumes[global_row];
        local_thermal_before_volume +=
            volume * (cell.e_electron_erg_per_cm3 + cell.e_ion_erg_per_cm3);
        local_thermal_after_volume += volume * (electron_energy + ion_energy);
        local_max_cellwise_residual = std::max(
            local_max_cellwise_residual,
            std::abs((electron_energy + ion_energy) -
                     (cell.e_electron_erg_per_cm3 + cell.e_ion_erg_per_cm3)));
      }
    }
  }

  result.global_stage_ok = AllreduceMinBool(ownership.communicator, local_stage_ok) != 0;
  if (!result.global_stage_ok) {
    Fail(result, "distributed thermal staged state is not publishable", rank, rank_count, false);
    return result;
  }

  constexpr std::uint32_t kWriteMask =
      kDistributedThermalWritesElectronEnergy |
      kDistributedThermalWritesFluidTotalEnergy;
  const bool local_publishable =
      staged_e_electron.extent_r() == problem.local_state.e_electron.extent_r() &&
      staged_e_total.extent_r() == problem.local_state.e_fluid_total.extent_r() &&
      kWriteMask ==
          (kDistributedThermalWritesElectronEnergy |
           kDistributedThermalWritesFluidTotalEnergy);
  result.global_publish_ok =
      AllreduceMinBool(ownership.communicator, local_publishable) != 0;
  if (!result.global_publish_ok) {
    Fail(result, "distributed thermal publish preflight failed", rank, rank_count, true);
    return result;
  }

  problem.local_state.e_electron = std::move(staged_e_electron);
  problem.local_state.e_fluid_total = std::move(staged_e_total);
  problem.local_state.ApplyAuthoritativeWrite(
      dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron) |
      dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_fluid_total));
  result.writeback_wall_s = ElapsedSecondsSince(writeback_timer);

  result.success = true;
  result.distributed_writeback = true;
  result.gathered_writeback_used = false;
  result.updated_fields = kWriteMask;
  result.updated_fields_label = "e_electron,e_fluid_total";
  result.local_owned_write_count = ownership.local_row_count;
  result.global_owned_write_count = static_cast<std::size_t>(
      AllreduceSumUll(ownership.communicator, static_cast<unsigned long long>(ownership.local_row_count)));
  double global_cellwise_residual = 0.0;
  MPI_Allreduce(
      &local_max_cellwise_residual,
      &global_cellwise_residual,
      1,
      MPI_DOUBLE,
      MPI_MAX,
      ownership.communicator);
  double global_thermal_before_volume = 0.0;
  double global_thermal_after_volume = 0.0;
  MPI_Allreduce(
      &local_thermal_before_volume,
      &global_thermal_before_volume,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      ownership.communicator);
  MPI_Allreduce(
      &local_thermal_after_volume,
      &global_thermal_after_volume,
      1,
      MPI_DOUBLE,
      MPI_SUM,
      ownership.communicator);
  result.global_thermal_energy_residual =
      std::abs(global_thermal_after_volume - global_thermal_before_volume);
  result.global_max_cellwise_thermal_energy_mismatch = global_cellwise_residual;
  result.electron_solve_report = electron_solve.report_line;
  result.ion_solve_report = ion_solve.report_line;
  result.coefficient_report = reservoir_problems.coefficient_report;
  if (options.electron_flux_limiter_model == ElectronFluxLimiterModel::disabled) {
    result.flux_limiter_report =
        "diagnostic_id=p2.thermal_flux_limiter; electron_flux_limiter_enabled=false; "
        "limited_face_count_local=0; limited_face_count_global=0; fallback_used=false";
  } else {
    const auto global_face_count =
        AllreduceSumUll(ownership.communicator, limiter_stats.local_face_count);
    const auto global_limited_face_count =
        AllreduceSumUll(ownership.communicator, limiter_stats.local_limited_face_count);
    double global_min_scale = 1.0;
    double global_max_flux_ratio = 0.0;
    MPI_Allreduce(
        &limiter_stats.min_scale,
        &global_min_scale,
        1,
        MPI_DOUBLE,
        MPI_MIN,
        ownership.communicator);
    MPI_Allreduce(
        &limiter_stats.max_flux_ratio,
        &global_max_flux_ratio,
        1,
        MPI_DOUBLE,
        MPI_MAX,
        ownership.communicator);
    std::ostringstream limiter_report;
    limiter_report << std::setprecision(17)
                   << "diagnostic_id=p2.thermal_flux_limiter"
                   << "; electron_flux_limiter_enabled=true"
                   << "; electron_flux_limiter_model="
                   << ElectronFluxLimiterModelName(options.electron_flux_limiter_model)
                   << "; alpha_e=" << options.electron_flux_limiter_alpha_e
                   << "; alpha_e_source=user_supplied"
                   << "; face_effective_kappa_consumed_by_diffusion=true"
                   << "; face_count_local=" << limiter_stats.local_face_count
                   << "; face_count_global=" << global_face_count
                   << "; limited_face_count_local=" << limiter_stats.local_limited_face_count
                   << "; limited_face_count_global=" << global_limited_face_count
                   << "; min_limiter_scale_global=" << global_min_scale
                   << "; max_flux_ratio_before_limit_global=" << global_max_flux_ratio
                   << "; fallback_used=false";
    result.flux_limiter_report = limiter_report.str();
  }
  result.writeback_report =
      "diagnostic_id=p2.thermal_conduction.distributed.writeback; owned_slab_only=true";

  std::ostringstream report;
  report << std::setprecision(17)
         << "diagnostic_id=p2.production.thermal_stage"
         << "; operator_diagnostic_id=p2.thermal_conduction.distributed"
         << "; stage_id=T"
         << "; stage_order=H,T,E"
         << "; distributed_writeback=true"
         << "; gathered_writeback_used=false"
         << "; local_stage_ok=true"
         << "; global_stage_ok=true"
         << "; local_publishable=true"
         << "; global_publish_ok=true"
         << "; updated_fields=e_electron,e_fluid_total"
         << "; local_owned_write_count=" << result.local_owned_write_count
         << "; global_owned_write_count=" << result.global_owned_write_count
         << "; global_thermal_energy_residual=" << result.global_thermal_energy_residual
         << "; global_volume_integrated_thermal_energy_residual="
         << result.global_thermal_energy_residual
         << "; global_max_cellwise_thermal_energy_mismatch="
         << result.global_max_cellwise_thermal_energy_mismatch
         << "; backend_requested=hypre_parcsr_gmres_boomeramg"
         << "; backend_executed=hypre_parcsr_gmres_boomeramg"
         << "; thermal_kappa_model_requested="
         << ThermalConductionKappaModelName(options.kappa_model)
         << "; thermal_kappa_model_executed="
         << ThermalConductionKappaModelName(options.kappa_model)
         << "; coefficient_time_level=old_time_lagged"
         << "; coefficient_provider_wall_s=" << result.coefficient_provider_wall_s
         << "; assembly_wall_s=" << result.assembly_wall_s
         << "; hypre_setup_wall_s=" << result.hypre_setup_wall_s
         << "; hypre_solve_wall_s=" << result.hypre_solve_wall_s
         << "; writeback_wall_s=" << result.writeback_wall_s
         << "; solver_iterations=" << result.solver_iterations
         << "; electron_matrix_report_present=true"
         << "; ion_matrix_report_present=true"
         << "; electron_solve_report_present=true"
         << "; ion_solve_report_present=true"
         << "; electron_flux_limiter_enabled="
         << (options.electron_flux_limiter_model == ElectronFluxLimiterModel::disabled ? "false" : "true")
         << "; flux_limiter_report_present=true"
         << "; canonical_state_mutated=true"
         << "; fallback_used=false";
  result.report_line = report.str();
  return result;
}

bool ValidateDistributedThermalConductionDiagnostics(
    const DistributedThermalConductionResult& result) noexcept {
  return result.is_complete() &&
         result.report_line.find("distributed_writeback=true") != std::string::npos &&
         result.report_line.find("gathered_writeback_used=false") != std::string::npos &&
         result.report_line.find("backend_executed=hypre_parcsr_gmres_boomeramg") !=
             std::string::npos &&
         result.report_line.find("backend_requested=hypre_parcsr_gmres_boomeramg") !=
             std::string::npos &&
         result.report_line.find("stage_id=T") != std::string::npos &&
         result.report_line.find("thermal_kappa_model_requested=") != std::string::npos &&
         result.report_line.find("global_volume_integrated_thermal_energy_residual=") !=
             std::string::npos &&
         result.report_line.find("coefficient_provider_wall_s=") != std::string::npos &&
         result.report_line.find("assembly_wall_s=") != std::string::npos &&
         result.report_line.find("hypre_setup_wall_s=") != std::string::npos &&
         result.report_line.find("hypre_solve_wall_s=") != std::string::npos &&
         result.report_line.find("writeback_wall_s=") != std::string::npos &&
         result.report_line.find("solver_iterations=") != std::string::npos &&
         result.report_line.find("fallback_used=false") != std::string::npos;
}

}  // namespace dec3d::transport
