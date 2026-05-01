#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/macro_ale_hllc_mpi.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ale/radial_ale.hpp"
#include "mesh/ale/radial_ale_mpi.hpp"
#include "mesh/ghost/hydro_halo_exchange.hpp"
#include "mesh/ownership/radial_ownership.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "hydro_mpi_real_case_support.hpp"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#ifndef DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT
#define DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT 24
#endif

#ifndef DEC3D_NOH_MPI24_RADIAL_CELLS
#define DEC3D_NOH_MPI24_RADIAL_CELLS 256u
#endif

#ifndef DEC3D_NOH_MPI24_USE_MACRO_ZONING
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 0
#endif

#ifndef DEC3D_NOH_MPI24_USE_PPM
#define DEC3D_NOH_MPI24_USE_PPM 1
#endif

#ifndef DEC3D_NOH_MPI24_USE_ALE
#define DEC3D_NOH_MPI24_USE_ALE 0
#endif

#ifndef DEC3D_NOH_MPI24_USE_NOH2
#define DEC3D_NOH_MPI24_USE_NOH2 0
#endif

#ifndef DEC3D_NOH_MPI24_NOH2_SPECIFIC_INTERNAL_ENERGY0
#define DEC3D_NOH_MPI24_NOH2_SPECIFIC_INTERNAL_ENERGY0 1.0
#endif

#ifndef DEC3D_NOH_MPI24_APPLY_GEOMETRIC_SOURCE
#define DEC3D_NOH_MPI24_APPLY_GEOMETRIC_SOURCE 1
#endif

#ifndef DEC3D_NOH_MPI24_APPLY_THETA_SWEEP
#define DEC3D_NOH_MPI24_APPLY_THETA_SWEEP 1
#endif

#ifndef DEC3D_NOH_MPI24_APPLY_PHI_SWEEP
#define DEC3D_NOH_MPI24_APPLY_PHI_SWEEP 1
#endif

#ifndef DEC3D_NOH_MPI24_MACRO_RADIAL_PPM
#define DEC3D_NOH_MPI24_MACRO_RADIAL_PPM 1
#endif

#ifndef DEC3D_NOH_MPI24_MACRO_THETA_PPM
#define DEC3D_NOH_MPI24_MACRO_THETA_PPM 1
#endif

#ifndef DEC3D_NOH_MPI24_MACRO_PHI_PPM
#define DEC3D_NOH_MPI24_MACRO_PHI_PPM 1
#endif

#ifndef DEC3D_NOH_MPI24_FINAL_TIME
#define DEC3D_NOH_MPI24_FINAL_TIME 0.6
#endif

#ifndef DEC3D_NOH_MPI24_MAX_STEPS
#define DEC3D_NOH_MPI24_MAX_STEPS 1000000u
#endif

#ifndef DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG
#define DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG 0
#endif

#ifndef DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS
#define DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS 0u
#endif

#ifndef DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS 0
#endif

#ifndef DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE 0u
#endif

#ifndef DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_nr256_p1e6_mpi24"
#endif

#ifndef DEC3D_NOH_MPI24_CASE_NAME_LITERAL
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL "case_noh_spherical_ppm_nr256_p1e6_mpi24"
#endif

#ifndef DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL "mesh:p1.noh.nr256.mpi24"
#endif

#ifndef DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.nr256.mpi24"
#endif

namespace {

struct NohRadialProfileSample {
  std::size_t step_index{0};
  double time_s{0.0};
  double exact_shock_radius{0.0};
  double numerical_shock_radius{0.0};
  std::vector<double> radial_centers;
  std::vector<double> rho_avg;
  std::vector<double> radial_velocity_avg;
  std::vector<double> velocity_magnitude_avg;
  std::vector<double> mom_r_avg;
  std::vector<double> e_fluid_total_avg;
};

struct ShellSymmetryMetrics {
  double max_rho_angular_relative_spread{0.0};
  double tangential_to_radial_momentum_ratio{0.0};
  std::size_t max_spread_radial_index{0u};
  std::size_t min_rho_theta_index{0u};
  std::size_t min_rho_phi_index{0u};
  std::size_t max_rho_theta_index{0u};
  std::size_t max_rho_phi_index{0u};
};

struct Noh2ProfileErrorMetrics {
  double rho_relative_l1{0.0};
  double radial_velocity_relative_l1{0.0};
  double e_fluid_total_relative_l1{0.0};
  double rho_relative_linf{0.0};
  double radial_velocity_relative_linf{0.0};
  double e_fluid_total_relative_linf{0.0};
};

struct LocalShellSymmetryProfile {
  std::vector<double> spread;
  std::vector<int> min_theta;
  std::vector<int> min_phi;
  std::vector<int> max_theta;
  std::vector<int> max_phi;
  double tangential_sum{0.0};
  double radial_sum{0.0};
};

[[nodiscard]] double Noh2ExactDensity(double rho0, double time_s) noexcept;
[[nodiscard]] double Noh2ExactSpecificInternalEnergy(
    double specific_internal_energy0,
    double time_s) noexcept;
[[nodiscard]] double Noh2ExactRadialVelocity(double time_s, double radius) noexcept;
[[nodiscard]] dec3d::hydro::HydroConservativeState Noh2ExactState(
    double rho0,
    double specific_internal_energy0,
    double electron_energy_fraction,
    double time_s,
    double radius);

[[nodiscard]] std::vector<double> BuildLocalShellAverage(
    const dec3d::core::Array3D<double>& local_field) {
  std::vector<double> values(local_field.extent_r(), 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_field.extent_theta() * local_field.extent_phi());
  for (std::size_t radial = 0; radial < local_field.extent_r(); ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_field.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < local_field.extent_phi(); ++phi) {
        sum += local_field(radial, theta, phi);
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

[[nodiscard]] std::vector<double> BuildLocalRadialVelocityAverage(
    const dec3d::state::CanonicalState& local_state) {
  std::vector<double> values(local_state.layout.radial_cells, 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_state.layout.theta_cells * local_state.layout.phi_cells);
  for (std::size_t radial = 0; radial < local_state.layout.radial_cells; ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < local_state.layout.phi_cells; ++phi) {
        const double rho = local_state.rho(radial, theta, phi);
        sum += rho > 0.0 ? local_state.mom_r(radial, theta, phi) / rho : 0.0;
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

[[nodiscard]] std::vector<double> BuildLocalVelocityMagnitudeAverage(
    const dec3d::state::CanonicalState& local_state) {
  std::vector<double> values(local_state.layout.radial_cells, 0.0);
  const double inverse_angular =
      1.0 / static_cast<double>(local_state.layout.theta_cells * local_state.layout.phi_cells);
  for (std::size_t radial = 0; radial < local_state.layout.radial_cells; ++radial) {
    double sum = 0.0;
    for (std::size_t theta = 0; theta < local_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < local_state.layout.phi_cells; ++phi) {
        const double rho = local_state.rho(radial, theta, phi);
        if (rho > 0.0) {
          const double vr = local_state.mom_r(radial, theta, phi) / rho;
          const double vt = local_state.mom_theta(radial, theta, phi) / rho;
          const double vp = local_state.mom_phi(radial, theta, phi) / rho;
          sum += std::sqrt((vr * vr) + (vt * vt) + (vp * vp));
        }
      }
    }
    values[radial] = sum * inverse_angular;
  }
  return values;
}

[[nodiscard]] LocalShellSymmetryProfile BuildLocalShellSymmetryProfile(
    const dec3d::state::CanonicalState& local_state) {
  LocalShellSymmetryProfile profile;
  profile.spread.resize(local_state.layout.radial_cells, 0.0);
  profile.min_theta.resize(local_state.layout.radial_cells, 0);
  profile.min_phi.resize(local_state.layout.radial_cells, 0);
  profile.max_theta.resize(local_state.layout.radial_cells, 0);
  profile.max_phi.resize(local_state.layout.radial_cells, 0);

  for (std::size_t radial = 0; radial < local_state.layout.radial_cells; ++radial) {
    double min_rho = std::numeric_limits<double>::infinity();
    double max_rho = -std::numeric_limits<double>::infinity();
    double rho_sum = 0.0;
    std::size_t min_theta = 0u;
    std::size_t min_phi = 0u;
    std::size_t max_theta = 0u;
    std::size_t max_phi = 0u;
    for (std::size_t theta = 0; theta < local_state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < local_state.layout.phi_cells; ++phi) {
        const double rho = local_state.rho(radial, theta, phi);
        if (rho < min_rho) {
          min_rho = rho;
          min_theta = theta;
          min_phi = phi;
        }
        if (rho > max_rho) {
          max_rho = rho;
          max_theta = theta;
          max_phi = phi;
        }
        rho_sum += rho;
        profile.tangential_sum +=
            std::abs(local_state.mom_theta(radial, theta, phi)) +
            std::abs(local_state.mom_phi(radial, theta, phi));
        profile.radial_sum += std::abs(local_state.mom_r(radial, theta, phi));
      }
    }

    const double angular_count =
        static_cast<double>(local_state.layout.theta_cells * local_state.layout.phi_cells);
    const double mean_rho = rho_sum / angular_count;
    profile.spread[radial] =
        std::abs(mean_rho) > 0.0 ? (max_rho - min_rho) / std::abs(mean_rho) : 0.0;
    profile.min_theta[radial] = static_cast<int>(min_theta);
    profile.min_phi[radial] = static_cast<int>(min_phi);
    profile.max_theta[radial] = static_cast<int>(max_theta);
    profile.max_phi[radial] = static_cast<int>(max_phi);
  }

  return profile;
}

void GatherRadialAverageToRoot(
    const std::vector<double>& local_values,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    int rank,
    std::vector<double>* global_values) {
  std::vector<int> recvcounts;
  std::vector<int> displacements;
  if (rank == 0) {
    recvcounts.resize(decomposition.slices.size(), 0);
    displacements.resize(decomposition.slices.size(), 0);
    for (std::size_t index = 0; index < decomposition.slices.size(); ++index) {
      const auto& slice = decomposition.slices[index];
      recvcounts[index] = static_cast<int>(slice.local_cell_count());
      displacements[index] = static_cast<int>(slice.begin_index);
    }
  }

  MPI_Gatherv(
      const_cast<double*>(local_values.data()),
      static_cast<int>(local_values.size()),
      MPI_DOUBLE,
      rank == 0 ? global_values->data() : nullptr,
      rank == 0 ? recvcounts.data() : nullptr,
      rank == 0 ? displacements.data() : nullptr,
      MPI_DOUBLE,
      0,
      MPI_COMM_WORLD);
}

void GatherRadialIntToRoot(
    const std::vector<int>& local_values,
    const dec3d::mesh::RadialOwnershipDecomposition& decomposition,
    int rank,
    std::vector<int>* global_values) {
  std::vector<int> recvcounts;
  std::vector<int> displacements;
  if (rank == 0) {
    recvcounts.resize(decomposition.slices.size(), 0);
    displacements.resize(decomposition.slices.size(), 0);
    for (std::size_t index = 0; index < decomposition.slices.size(); ++index) {
      const auto& slice = decomposition.slices[index];
      recvcounts[index] = static_cast<int>(slice.local_cell_count());
      displacements[index] = static_cast<int>(slice.begin_index);
    }
  }

  MPI_Gatherv(
      const_cast<int*>(local_values.data()),
      static_cast<int>(local_values.size()),
      MPI_INT,
      rank == 0 ? global_values->data() : nullptr,
      rank == 0 ? recvcounts.data() : nullptr,
      rank == 0 ? displacements.data() : nullptr,
      MPI_INT,
      0,
      MPI_COMM_WORLD);
}

[[nodiscard]] double EstimateNohShockRadiusFromShellAverage(
    const std::vector<double>& rho_profile,
    const std::vector<double>& radial_faces,
    double exact_shock_radius) {
  if (rho_profile.size() < 2u || radial_faces.size() < rho_profile.size() + 1u) {
    return 0.0;
  }

  std::size_t shock_index = 0u;
  double strongest_fall = std::numeric_limits<double>::infinity();
  for (std::size_t radial = 0; radial + 1u < rho_profile.size(); ++radial) {
    const double left_center = 0.5 * (radial_faces[radial] + radial_faces[radial + 1u]);
    const double right_center = 0.5 * (radial_faces[radial + 1u] + radial_faces[radial + 2u]);
    const double interface_radius = 0.5 * (left_center + right_center);
    if (exact_shock_radius > 0.0 &&
        (interface_radius < 0.4 * exact_shock_radius ||
         interface_radius > 1.6 * exact_shock_radius)) {
      continue;
    }
    const double gradient =
        (rho_profile[radial + 1u] - rho_profile[radial]) /
        (right_center - left_center);
    if (gradient < strongest_fall) {
      strongest_fall = gradient;
      shock_index = radial;
    }
  }

  if (!std::isfinite(strongest_fall)) {
    return 0.0;
  }
  const double left_center = 0.5 * (radial_faces[shock_index] + radial_faces[shock_index + 1u]);
  const double right_center = 0.5 * (radial_faces[shock_index + 1u] + radial_faces[shock_index + 2u]);
  return 0.5 * (left_center + right_center);
}

[[nodiscard]] bool WriteNohRadialProfileSampleFile(
    const std::filesystem::path& path,
    const char* label,
    const NohRadialProfileSample& sample) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=noh_radial_profile\n";
  output << "# label=" << label << '\n';
  output << "# step_index=" << sample.step_index << '\n';
  output << std::setprecision(17);
  output << "# time_s=" << sample.time_s << '\n';
  output << "# exact_shock_radius=" << sample.exact_shock_radius << '\n';
  output << "# numerical_shock_radius=" << sample.numerical_shock_radius << '\n';
  output << "# radial_cells=" << sample.radial_centers.size() << '\n';
  output << "# columns=radial_index radial_center rho_avg radial_velocity_avg "
            "velocity_magnitude_avg mom_r_avg e_fluid_total_avg\n";
  for (std::size_t radial = 0; radial < sample.radial_centers.size(); ++radial) {
    output << radial << ' '
           << sample.radial_centers[radial] << ' '
           << sample.rho_avg[radial] << ' '
           << sample.radial_velocity_avg[radial] << ' '
           << sample.velocity_magnitude_avg[radial] << ' '
           << sample.mom_r_avg[radial] << ' '
           << sample.e_fluid_total_avg[radial] << '\n';
  }
  return true;
}

[[nodiscard]] bool WriteNohShockHistory(
    const std::filesystem::path& path,
    const std::vector<NohRadialProfileSample>& profile_history) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=noh_shock_radius_history\n";
  output << "# definition=numerical radius is strongest negative density gradient / outer shock front\n";
  output << "# columns=time_s numerical_radius exact_radius relative_error\n";
  output << std::setprecision(17);
  for (const auto& sample : profile_history) {
    const double relative_error =
        sample.exact_shock_radius > 0.0
            ? std::abs(sample.numerical_shock_radius - sample.exact_shock_radius) /
                  sample.exact_shock_radius
            : 0.0;
    output << sample.time_s << ' '
           << sample.numerical_shock_radius << ' '
           << sample.exact_shock_radius << ' '
           << relative_error << '\n';
  }
  return true;
}

[[nodiscard]] Noh2ProfileErrorMetrics MeasureNoh2ProfileError(
    const NohRadialProfileSample& sample,
    double rho0,
    double specific_internal_energy0) noexcept {
  Noh2ProfileErrorMetrics metrics;
  double rho_error_l1 = 0.0;
  double rho_exact_l1 = 0.0;
  double velocity_error_l1 = 0.0;
  double velocity_exact_l1 = 0.0;
  double energy_error_l1 = 0.0;
  double energy_exact_l1 = 0.0;

  for (std::size_t radial = 0; radial < sample.radial_centers.size(); ++radial) {
    const double radius = sample.radial_centers[radial];
    const double exact_rho = Noh2ExactDensity(rho0, sample.time_s);
    const double exact_velocity = Noh2ExactRadialVelocity(sample.time_s, radius);
    const double exact_sie =
        Noh2ExactSpecificInternalEnergy(specific_internal_energy0, sample.time_s);
    const double exact_e_total =
        (exact_rho * exact_sie) + 0.5 * exact_rho * exact_velocity * exact_velocity;

    const double rho_error = std::abs(sample.rho_avg[radial] - exact_rho);
    const double velocity_error =
        std::abs(sample.radial_velocity_avg[radial] - exact_velocity);
    const double energy_error =
        std::abs(sample.e_fluid_total_avg[radial] - exact_e_total);

    rho_error_l1 += rho_error;
    rho_exact_l1 += std::abs(exact_rho);
    velocity_error_l1 += velocity_error;
    velocity_exact_l1 += std::abs(exact_velocity);
    energy_error_l1 += energy_error;
    energy_exact_l1 += std::abs(exact_e_total);

    metrics.rho_relative_linf = std::max(
        metrics.rho_relative_linf,
        std::abs(exact_rho) > 0.0 ? rho_error / std::abs(exact_rho) : rho_error);
    metrics.radial_velocity_relative_linf = std::max(
        metrics.radial_velocity_relative_linf,
        std::abs(exact_velocity) > 0.0
            ? velocity_error / std::abs(exact_velocity)
            : velocity_error);
    metrics.e_fluid_total_relative_linf = std::max(
        metrics.e_fluid_total_relative_linf,
        std::abs(exact_e_total) > 0.0
            ? energy_error / std::abs(exact_e_total)
            : energy_error);
  }

  metrics.rho_relative_l1 =
      rho_exact_l1 > 0.0 ? rho_error_l1 / rho_exact_l1 : rho_error_l1;
  metrics.radial_velocity_relative_l1 =
      velocity_exact_l1 > 0.0 ? velocity_error_l1 / velocity_exact_l1 : velocity_error_l1;
  metrics.e_fluid_total_relative_l1 =
      energy_exact_l1 > 0.0 ? energy_error_l1 / energy_exact_l1 : energy_error_l1;
  return metrics;
}

[[nodiscard]] bool WriteNoh2ProfileErrorHistory(
    const std::filesystem::path& path,
    const std::vector<NohRadialProfileSample>& profile_history,
    double rho0,
    double specific_internal_energy0) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }

  output << "# kind=noh2_profile_error_history\n";
  output << "# columns=time_s rho_relative_l1 radial_velocity_relative_l1 "
            "e_fluid_total_relative_l1 rho_relative_linf "
            "radial_velocity_relative_linf e_fluid_total_relative_linf\n";
  output << std::setprecision(17);
  for (const auto& sample : profile_history) {
    const auto metrics = MeasureNoh2ProfileError(sample, rho0, specific_internal_energy0);
    output << sample.time_s << ' '
           << metrics.rho_relative_l1 << ' '
           << metrics.radial_velocity_relative_l1 << ' '
           << metrics.e_fluid_total_relative_l1 << ' '
           << metrics.rho_relative_linf << ' '
           << metrics.radial_velocity_relative_linf << ' '
           << metrics.e_fluid_total_relative_linf << '\n';
  }
  return true;
}

[[nodiscard]] ShellSymmetryMetrics MeasureShellSymmetry(
    const dec3d::state::CanonicalState& state) {
  ShellSymmetryMetrics metrics;
  double tangential_sum = 0.0;
  double radial_sum = 0.0;

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    double min_rho = std::numeric_limits<double>::infinity();
    double max_rho = -std::numeric_limits<double>::infinity();
    double rho_sum = 0.0;
    std::size_t min_theta = 0u;
    std::size_t min_phi = 0u;
    std::size_t max_theta = 0u;
    std::size_t max_phi = 0u;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        const double rho = state.rho(radial, theta, phi);
        if (rho < min_rho) {
          min_rho = rho;
          min_theta = theta;
          min_phi = phi;
        }
        if (rho > max_rho) {
          max_rho = rho;
          max_theta = theta;
          max_phi = phi;
        }
        rho_sum += rho;
        tangential_sum +=
            std::abs(state.mom_theta(radial, theta, phi)) +
            std::abs(state.mom_phi(radial, theta, phi));
        radial_sum += std::abs(state.mom_r(radial, theta, phi));
      }
    }

    const double angular_count =
        static_cast<double>(state.layout.theta_cells * state.layout.phi_cells);
    const double mean_rho = rho_sum / angular_count;
    if (std::abs(mean_rho) > 0.0) {
      const double spread = (max_rho - min_rho) / std::abs(mean_rho);
      if (spread > metrics.max_rho_angular_relative_spread) {
        metrics.max_rho_angular_relative_spread = spread;
        metrics.max_spread_radial_index = radial;
        metrics.min_rho_theta_index = min_theta;
        metrics.min_rho_phi_index = min_phi;
        metrics.max_rho_theta_index = max_theta;
        metrics.max_rho_phi_index = max_phi;
      }
    }
  }

  metrics.tangential_to_radial_momentum_ratio =
      radial_sum > 0.0 ? tangential_sum / radial_sum : 0.0;
  return metrics;
}

void SeedNohSphericalInflow(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before,
    double rho0,
    double pressure0,
    double radial_velocity0,
    double electron_energy_fraction) {
  const double e_internal = pressure0 / (dec3d::state::HydroIdealGasGamma() - 1.0);
  const double e_total =
      e_internal + 0.5 * rho0 * radial_velocity0 * radial_velocity0;
  const double e_electron = dec3d::state::ElectronEnergyDensityFromPressure(
      electron_energy_fraction * pressure0);

  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = rho0;
        state.mom_r(radial, theta, phi) = rho0 * radial_velocity0;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = e_total;
        state.e_electron(radial, theta, phi) = e_electron;

        before.rho(radial, theta, phi) = state.rho(radial, theta, phi);
        before.mom_r(radial, theta, phi) = state.mom_r(radial, theta, phi);
        before.mom_theta(radial, theta, phi) = state.mom_theta(radial, theta, phi);
        before.mom_phi(radial, theta, phi) = state.mom_phi(radial, theta, phi);
        before.e_fluid_total(radial, theta, phi) = state.e_fluid_total(radial, theta, phi);
        before.e_electron(radial, theta, phi) = state.e_electron(radial, theta, phi);
      }
    }
  }
}

[[nodiscard]] double Noh2Compression(double time_s) noexcept {
  return (time_s < 1.0) ? 1.0 / (1.0 - time_s) : std::numeric_limits<double>::infinity();
}

[[nodiscard]] double Noh2ExactDensity(
    double rho0,
    double time_s) noexcept {
  const double compression = Noh2Compression(time_s);
  return rho0 * compression * compression * compression;
}

[[nodiscard]] double Noh2ExactSpecificInternalEnergy(
    double specific_internal_energy0,
    double time_s) noexcept {
  const double compression = Noh2Compression(time_s);
  return specific_internal_energy0 * compression * compression;
}

[[nodiscard]] double Noh2ExactRadialVelocity(double time_s, double radius) noexcept {
  return -radius * Noh2Compression(time_s);
}

[[nodiscard]] dec3d::hydro::HydroConservativeState Noh2ExactState(
    double rho0,
    double specific_internal_energy0,
    double electron_energy_fraction,
    double time_s,
    double radius) {
  const double rho = Noh2ExactDensity(rho0, time_s);
  const double radial_velocity = Noh2ExactRadialVelocity(time_s, radius);
  const double specific_internal_energy =
      Noh2ExactSpecificInternalEnergy(specific_internal_energy0, time_s);
  const double pressure =
      rho * specific_internal_energy * (dec3d::state::HydroIdealGasGamma() - 1.0);
  return dec3d::hydro::HydroConservativeState{
      rho,
      rho * radial_velocity,
      0.0,
      0.0,
      (rho * specific_internal_energy) +
          0.5 * rho * radial_velocity * radial_velocity,
      dec3d::state::ChiEFromElectronPressure(electron_energy_fraction * pressure)};
}

void SeedNoh2SphericalFlow(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double rho0,
    double specific_internal_energy0,
    double electron_energy_fraction) {
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const double radial_center =
        0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
    const auto state_at_cell =
        Noh2ExactState(rho0, specific_internal_energy0, electron_energy_fraction, 0.0, radial_center);
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = state_at_cell.rho;
        state.mom_r(radial, theta, phi) = state_at_cell.mom_r;
        state.mom_theta(radial, theta, phi) = state_at_cell.mom_theta;
        state.mom_phi(radial, theta, phi) = state_at_cell.mom_phi;
        state.e_fluid_total(radial, theta, phi) = state_at_cell.e_fluid_total;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(
                dec3d::state::ElectronPressureFromChiE(state_at_cell.chi_e));

        before.rho(radial, theta, phi) = state.rho(radial, theta, phi);
        before.mom_r(radial, theta, phi) = state.mom_r(radial, theta, phi);
        before.mom_theta(radial, theta, phi) = state.mom_theta(radial, theta, phi);
        before.mom_phi(radial, theta, phi) = state.mom_phi(radial, theta, phi);
        before.e_fluid_total(radial, theta, phi) = state.e_fluid_total(radial, theta, phi);
        before.e_electron(radial, theta, phi) = state.e_electron(radial, theta, phi);
      }
    }
  }
}

[[nodiscard]] dec3d::hydro::HydroConservativeState NohInflowGhostState(
    double rho0,
    double pressure0,
    double radial_velocity0,
    double electron_energy_fraction) {
  const double e_internal = pressure0 / (dec3d::state::HydroIdealGasGamma() - 1.0);
  return dec3d::hydro::HydroConservativeState{
      rho0,
      rho0 * radial_velocity0,
      0.0,
      0.0,
      e_internal + 0.5 * rho0 * radial_velocity0 * radial_velocity0,
      dec3d::state::ChiEFromElectronPressure(electron_energy_fraction * pressure0)};
}

[[nodiscard]] double NohExactpackUpstreamDensity(
    double rho0,
    double radial_velocity0,
    double time_s,
    double radius) noexcept {
  if (!(rho0 > 0.0) || !(radius > 0.0) || !(time_s >= 0.0)) {
    return rho0;
  }
  const double speed = std::abs(radial_velocity0);
  const double compression = 1.0 + (speed * time_s / radius);
  return rho0 * compression * compression;
}

[[nodiscard]] dec3d::hydro::HydroConservativeState NohExactpackUpstreamInflowGhostState(
    double rho0,
    double pressure0,
    double radial_velocity0,
    double electron_energy_fraction,
    double time_s,
    double radius) {
  const double rho =
      NohExactpackUpstreamDensity(rho0, radial_velocity0, time_s, radius);
  const double e_internal = pressure0 / (dec3d::state::HydroIdealGasGamma() - 1.0);
  return dec3d::hydro::HydroConservativeState{
      rho,
      rho * radial_velocity0,
      0.0,
      0.0,
      e_internal + 0.5 * rho * radial_velocity0 * radial_velocity0,
      dec3d::state::ChiEFromElectronPressure(electron_energy_fraction * pressure0)};
}

void ApplyOuterInflowGhostIfNeeded(
    dec3d::hydro::RadialGhostOverride& ghost_override,
    bool is_outer_rank,
    std::size_t theta_cells,
    std::size_t phi_cells,
    const dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    double time_s,
    double rho0,
    double pressure0,
    double radial_velocity0,
    double electron_energy_fraction) {
  if (!is_outer_rank) {
    return;
  }

  ghost_override.has_outer_neighbor = true;
  ghost_override.outer_ghost_reuses_boundary_partition = true;
  ghost_override.outer_ghost_states = dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
      ghost_override.ghost_layers,
      theta_cells,
      phi_cells);
  const bool has_outer_spacing = local_geometry.radial_faces.size() >= 2u;
  const double outer_face_radius =
      has_outer_spacing ? local_geometry.radial_faces.back() : 1.0;
  const double outer_cell_width =
      has_outer_spacing
          ? local_geometry.radial_faces.back() -
                local_geometry.radial_faces[local_geometry.radial_faces.size() - 2u]
          : 0.0;
  for (std::size_t ghost = 0; ghost < ghost_override.ghost_layers; ++ghost) {
    const double ghost_radius =
        outer_face_radius +
        (static_cast<double>(ghost) + 0.5) * std::max(outer_cell_width, 0.0);
    const auto inflow_state = NohExactpackUpstreamInflowGhostState(
        rho0,
        pressure0,
        radial_velocity0,
        electron_energy_fraction,
        time_s,
        ghost_radius);
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        ghost_override.outer_ghost_states(ghost, theta, phi) = inflow_state;
      }
    }
  }
  ghost_override.report_line += "; outer_boundary=exactpack_noh_upstream_inflow";
}

void ApplyOuterNoh2InflowGhostIfNeeded(
    dec3d::hydro::RadialGhostOverride& ghost_override,
    bool is_outer_rank,
    std::size_t theta_cells,
    std::size_t phi_cells,
    const dec3d::mesh::SphericalGeometryMetadata& local_geometry,
    double time_s,
    double rho0,
    double specific_internal_energy0,
    double electron_energy_fraction) {
  if (!is_outer_rank) {
    return;
  }

  ghost_override.has_outer_neighbor = true;
  ghost_override.outer_ghost_reuses_boundary_partition = true;
  ghost_override.outer_ghost_states = dec3d::core::Array3D<dec3d::hydro::HydroConservativeState>(
      ghost_override.ghost_layers,
      theta_cells,
      phi_cells);
  const bool has_outer_spacing = local_geometry.radial_faces.size() >= 2u;
  const double outer_face_radius =
      has_outer_spacing ? local_geometry.radial_faces.back() : (1.0 - time_s);
  const double outer_cell_width =
      has_outer_spacing
          ? local_geometry.radial_faces.back() -
                local_geometry.radial_faces[local_geometry.radial_faces.size() - 2u]
          : 0.0;
  for (std::size_t ghost = 0; ghost < ghost_override.ghost_layers; ++ghost) {
    const double ghost_radius =
        outer_face_radius +
        (static_cast<double>(ghost) + 0.5) * std::max(outer_cell_width, 0.0);
    const auto inflow_state = Noh2ExactState(
        rho0,
        specific_internal_energy0,
        electron_energy_fraction,
        time_s,
        ghost_radius);
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        ghost_override.outer_ghost_states(ghost, theta, phi) = inflow_state;
      }
    }
  }
  ghost_override.report_line += "; outer_boundary=exactpack_noh2_inflow";
}

[[nodiscard]] bool WriteBudgetResidual(
    const std::filesystem::path& path,
    const dec3d::hydro::HydroBudgetResidualSummary& budget) {
  std::ofstream output(path);
  if (!output.is_open()) {
    return false;
  }
  output << std::setprecision(17);
  output << "old_mass=" << budget.old_mass << '\n';
  output << "flux_mass_delta=" << budget.flux_mass_delta << '\n';
  output << "source_mass_delta=" << budget.source_mass_delta << '\n';
  output << "new_mass=" << budget.new_mass << '\n';
  output << "mass_residual=" << budget.mass_residual << '\n';
  output << "old_mom_r=" << budget.old_mom_r << '\n';
  output << "flux_mom_r_delta=" << budget.flux_mom_r_delta << '\n';
  output << "source_mom_r_delta=" << budget.source_mom_r_delta << '\n';
  output << "new_mom_r=" << budget.new_mom_r << '\n';
  output << "mom_r_residual=" << budget.mom_r_residual << '\n';
  output << "old_e_fluid_total=" << budget.old_e_fluid_total << '\n';
  output << "flux_e_fluid_total_delta=" << budget.flux_e_fluid_total_delta << '\n';
  output << "source_e_fluid_total_delta=" << budget.source_e_fluid_total_delta << '\n';
  output << "new_e_fluid_total=" << budget.new_e_fluid_total << '\n';
  output << "e_fluid_total_residual=" << budget.e_fluid_total_residual << '\n';
  return true;
}

[[nodiscard]] bool WriteWallHeatingDiagnostics(
    const std::filesystem::path& path,
    const NohRadialProfileSample& final_profile) {
  std::ofstream output(path);
  if (!output.is_open() || final_profile.rho_avg.empty()) {
    return false;
  }

  constexpr std::size_t kOriginCells = 4u;
  constexpr double kExactPlateauDensity = 64.0;
  const std::size_t origin_cells =
      std::min(kOriginCells, final_profile.rho_avg.size());

  double origin_peak_density = final_profile.rho_avg.front();
  std::size_t origin_peak_index = 0u;
  for (std::size_t radial = 0; radial < origin_cells; ++radial) {
    if (final_profile.rho_avg[radial] > origin_peak_density) {
      origin_peak_density = final_profile.rho_avg[radial];
      origin_peak_index = radial;
    }
  }

  double plateau_sum = 0.0;
  std::size_t plateau_count = 0u;
  const double plateau_outer_radius = 0.8 * final_profile.exact_shock_radius;
  for (std::size_t radial = origin_cells; radial < final_profile.rho_avg.size(); ++radial) {
    if (final_profile.radial_centers[radial] >= plateau_outer_radius) {
      break;
    }
    plateau_sum += final_profile.rho_avg[radial];
    ++plateau_count;
  }

  const double plateau_mean =
      plateau_count > 0u
          ? plateau_sum / static_cast<double>(plateau_count)
          : std::numeric_limits<double>::quiet_NaN();
  const double first_cell_to_plateau_ratio =
      plateau_mean > 0.0 ? final_profile.rho_avg.front() / plateau_mean
                         : std::numeric_limits<double>::quiet_NaN();
  const double origin_peak_to_plateau_ratio =
      plateau_mean > 0.0 ? origin_peak_density / plateau_mean
                         : std::numeric_limits<double>::quiet_NaN();
  const double plateau_relative_error =
      std::isfinite(plateau_mean)
          ? std::abs(plateau_mean - kExactPlateauDensity) / kExactPlateauDensity
          : std::numeric_limits<double>::quiet_NaN();

  output << std::setprecision(17);
  output << "# kind=noh_wall_heating_diagnostics\n";
  output << "time_s=" << final_profile.time_s << '\n';
  output << "exact_shock_radius=" << final_profile.exact_shock_radius << '\n';
  output << "origin_window_cells=" << origin_cells << '\n';
  output << "origin_peak_radial_index=" << origin_peak_index << '\n';
  output << "origin_peak_radius=" << final_profile.radial_centers[origin_peak_index] << '\n';
  output << "origin_peak_density=" << origin_peak_density << '\n';
  output << "first_cell_density=" << final_profile.rho_avg.front() << '\n';
  output << "plateau_excludes_origin_cells=" << origin_cells << '\n';
  output << "plateau_outer_radius_fraction_of_exact_shock=0.8\n";
  output << "plateau_sample_count=" << plateau_count << '\n';
  output << "plateau_mean_density=" << plateau_mean << '\n';
  output << "exact_plateau_density=" << kExactPlateauDensity << '\n';
  output << "plateau_relative_error=" << plateau_relative_error << '\n';
  output << "first_cell_to_plateau_density_ratio=" << first_cell_to_plateau_ratio << '\n';
  output << "origin_peak_to_plateau_density_ratio=" << origin_peak_to_plateau_ratio << '\n';
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  int mpi_initialized = 0;

  try {
    MPI_Init(&argc, &argv);
    MPI_Initialized(&mpi_initialized);

    int rank = 0;
    int rank_count = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_count);

    constexpr int kExpectedRanks = DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT;
    if (rank_count != kExpectedRanks) {
      if (rank == 0) {
        std::cerr << "expected " << kExpectedRanks << " MPI ranks\n";
      }
      MPI_Finalize();
      return 1;
    }

    constexpr std::size_t kRadialCells = DEC3D_NOH_MPI24_RADIAL_CELLS;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr std::size_t kGhostLayers =
        DEC3D_NOH_MPI24_USE_PPM != 0 ? 3u : 1u;
    constexpr bool kUseAle = DEC3D_NOH_MPI24_USE_ALE != 0;
    constexpr bool kUseNoh2 = DEC3D_NOH_MPI24_USE_NOH2 != 0;
    constexpr double kRho0 = 1.0;
    constexpr double kPressure0 = 1.0e-6;
    constexpr double kRadialVelocity0 = -1.0;
    constexpr double kNoh2SpecificInternalEnergy0 =
        DEC3D_NOH_MPI24_NOH2_SPECIFIC_INTERNAL_ENERGY0;
    constexpr double kElectronEnergyFraction = 0.5;
    constexpr double kFinalTime = DEC3D_NOH_MPI24_FINAL_TIME;
    constexpr std::size_t kMaxSteps = DEC3D_NOH_MPI24_MAX_STEPS;
    constexpr bool kWriteStageAngularDebug =
        DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG != 0;
    constexpr std::size_t kStageDebugSteps = DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS;
    constexpr double kShockTolerance = 0.25;
    constexpr double kNoh2ProfileTolerance = 0.20;
    constexpr double kShellSymmetryTolerance = 1.0e-8;
    constexpr double kTangentialMomentumTolerance = 1.0e-10;

    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
    if (!geometry.is_valid()) {
      if (rank == 0) {
        std::cerr << "invalid geometry\n";
      }
      MPI_Finalize();
      return 1;
    }

    auto global_before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto global_seed = dec3d::state::CanonicalState::Create(global_before.layout);
    if (kUseNoh2) {
      SeedNoh2SphericalFlow(
          global_seed,
          global_before,
          geometry,
          kRho0,
          kNoh2SpecificInternalEnergy0,
          kElectronEnergyFraction);
    } else {
      SeedNohSphericalInflow(
          global_seed,
          global_before,
          kRho0,
          kPressure0,
          kRadialVelocity0,
          kElectronEnergyFraction);
    }

    const auto decomposition =
        dec3d::mesh::BuildRadialOwnership(kRadialCells, static_cast<std::size_t>(rank_count));
    if (!decomposition.is_valid()) {
      if (rank == 0) {
        std::cerr << decomposition.failure_reason << '\n';
      }
      MPI_Finalize();
      return 1;
    }

    const auto& local_slice = decomposition.slices[static_cast<std::size_t>(rank)];
    auto local_geometry =
        dec3d::testsupport::BuildLocalGeometrySlice(geometry, local_slice, kThetaCells, kPhiCells);
    auto local_state = dec3d::testsupport::BuildLocalStateSlice(global_seed, local_slice);
    auto current_global_radial_faces = geometry.radial_faces;

    dec3d::hydro::StaticGridHydroOptions options{};
    options.apply_geometric_source = DEC3D_NOH_MPI24_APPLY_GEOMETRIC_SOURCE != 0;
    options.apply_radial_sweep = true;
    options.apply_theta_sweep = DEC3D_NOH_MPI24_APPLY_THETA_SWEEP != 0;
    options.apply_phi_sweep = DEC3D_NOH_MPI24_APPLY_PHI_SWEEP != 0;
    options.use_ppm_reconstruction = DEC3D_NOH_MPI24_USE_PPM != 0;
    options.reconstruction_ghost_layers = kGhostLayers;
    options.apply_radial_ale_flux_correction = kUseAle;
    options.request_radial_ale_proposal = kUseAle;
    options.radial_ale_global_face_begin_index = local_slice.begin_index;
    options.use_macro_zoning = DEC3D_NOH_MPI24_USE_MACRO_ZONING != 0;
    options.macro_zoning_coarse_factor = 0.5;
    options.macro_zoning_use_radial_ppm = DEC3D_NOH_MPI24_MACRO_RADIAL_PPM != 0;
    options.macro_zoning_use_theta_ppm = DEC3D_NOH_MPI24_MACRO_THETA_PPM != 0;
    options.macro_zoning_use_phi_ppm = DEC3D_NOH_MPI24_MACRO_PHI_PPM != 0;
    options.debug_angular_stage_diagnostics = kWriteStageAngularDebug;
    options.debug_macro_radial_ppm_face_smoothness =
        DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS != 0;
    options.debug_macro_radial_ppm_face =
        DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE;

    std::vector<NohRadialProfileSample> profile_history;
    dec3d::hydro::HydroBudgetResidualSummary reduced_budget{};
    std::ofstream progress_output;
    std::ofstream shell_symmetry_output;
    std::ofstream stage_angular_debug_output;
    std::ofstream moving_mesh_output;
    const std::filesystem::path output_dir =
        DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL;
    if (rank == 0) {
      std::filesystem::remove_all(output_dir);
      std::filesystem::create_directories(output_dir);
      progress_output.open(output_dir / "dec3d.out", std::ios::out | std::ios::trunc);
      progress_output << "# kind=run_progress\n";
      progress_output << "# columns=step time_s status\n";
      progress_output.flush();
      shell_symmetry_output.open(output_dir / "shell_symmetry_vs_time.txt", std::ios::out | std::ios::trunc);
      shell_symmetry_output << "# kind=noh_shell_symmetry_history\n";
      shell_symmetry_output << "# columns=step time_s max_rho_angular_relative_spread "
                            << "max_spread_radial_index radial_center min_rho_theta min_rho_phi "
                            << "max_rho_theta max_rho_phi tangential_to_radial_momentum_ratio\n";
      shell_symmetry_output.flush();
      if (kUseAle) {
        moving_mesh_output.open(output_dir / "moving_mesh_diagnostics.txt", std::ios::out | std::ios::trunc);
        moving_mesh_output << "# kind=noh_moving_mesh_diagnostics\n";
        moving_mesh_output << "# columns=step time_s dt_s proposal_complete "
                           << "global_max_face_speed global_min_outer_face_speed "
                           << "global_max_outer_face_speed\n";
        moving_mesh_output.flush();
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (kWriteStageAngularDebug) {
      std::ostringstream rank_file;
      rank_file << "stage_angular_debug_rank"
                << std::setw(3) << std::setfill('0') << rank << ".txt";
      stage_angular_debug_output.open(output_dir / rank_file.str(), std::ios::out | std::ios::trunc);
      stage_angular_debug_output << "# kind=noh_stage_angular_debug\n";
      stage_angular_debug_output << "# columns=step time_s dt_s rank message\n";
      stage_angular_debug_output.flush();
    }

    double time_s = 0.0;
    std::size_t step = 0u;

    auto abort_with_rank = [&](const std::string& message) -> int {
      std::cerr << "rank " << rank << ": " << message << '\n';
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    };

    auto record_profile_snapshot = [&](std::size_t step_index, double sample_time) -> bool {
      auto local_rho = BuildLocalShellAverage(local_state.rho);
      auto local_radial_velocity = BuildLocalRadialVelocityAverage(local_state);
      auto local_velocity = BuildLocalVelocityMagnitudeAverage(local_state);
      auto local_mom_r = BuildLocalShellAverage(local_state.mom_r);
      auto local_e_total = BuildLocalShellAverage(local_state.e_fluid_total);
      auto local_symmetry = BuildLocalShellSymmetryProfile(local_state);

      std::vector<double> global_rho;
      std::vector<double> global_radial_velocity;
      std::vector<double> global_velocity;
      std::vector<double> global_mom_r;
      std::vector<double> global_e_total;
      std::vector<double> global_shell_spread;
      std::vector<int> global_min_theta;
      std::vector<int> global_min_phi;
      std::vector<int> global_max_theta;
      std::vector<int> global_max_phi;
      if (rank == 0) {
        global_rho.resize(kRadialCells);
        global_radial_velocity.resize(kRadialCells);
        global_velocity.resize(kRadialCells);
        global_mom_r.resize(kRadialCells);
        global_e_total.resize(kRadialCells);
        global_shell_spread.resize(kRadialCells);
        global_min_theta.resize(kRadialCells);
        global_min_phi.resize(kRadialCells);
        global_max_theta.resize(kRadialCells);
        global_max_phi.resize(kRadialCells);
      }

      GatherRadialAverageToRoot(local_rho, decomposition, rank, &global_rho);
      GatherRadialAverageToRoot(local_radial_velocity, decomposition, rank, &global_radial_velocity);
      GatherRadialAverageToRoot(local_velocity, decomposition, rank, &global_velocity);
      GatherRadialAverageToRoot(local_mom_r, decomposition, rank, &global_mom_r);
      GatherRadialAverageToRoot(local_e_total, decomposition, rank, &global_e_total);
      GatherRadialAverageToRoot(local_symmetry.spread, decomposition, rank, &global_shell_spread);
      GatherRadialIntToRoot(local_symmetry.min_theta, decomposition, rank, &global_min_theta);
      GatherRadialIntToRoot(local_symmetry.min_phi, decomposition, rank, &global_min_phi);
      GatherRadialIntToRoot(local_symmetry.max_theta, decomposition, rank, &global_max_theta);
      GatherRadialIntToRoot(local_symmetry.max_phi, decomposition, rank, &global_max_phi);
      double global_tangential_sum = 0.0;
      double global_radial_sum = 0.0;
      MPI_Reduce(
          &local_symmetry.tangential_sum,
          &global_tangential_sum,
          1,
          MPI_DOUBLE,
          MPI_SUM,
          0,
          MPI_COMM_WORLD);
      MPI_Reduce(
          &local_symmetry.radial_sum,
          &global_radial_sum,
          1,
          MPI_DOUBLE,
          MPI_SUM,
          0,
          MPI_COMM_WORLD);

      if (rank == 0) {
        std::size_t max_spread_radial = 0u;
        for (std::size_t radial = 1u; radial < global_shell_spread.size(); ++radial) {
          if (global_shell_spread[radial] > global_shell_spread[max_spread_radial]) {
            max_spread_radial = radial;
          }
        }
        const double max_spread_radial_center =
            0.5 * (current_global_radial_faces[max_spread_radial] +
                   current_global_radial_faces[max_spread_radial + 1u]);
        const double tangential_ratio =
            global_radial_sum > 0.0 ? global_tangential_sum / global_radial_sum : 0.0;
        shell_symmetry_output << step_index << ' '
                              << std::setprecision(17) << sample_time << ' '
                              << global_shell_spread[max_spread_radial] << ' '
                              << max_spread_radial << ' '
                              << max_spread_radial_center << ' '
                              << global_min_theta[max_spread_radial] << ' '
                              << global_min_phi[max_spread_radial] << ' '
                              << global_max_theta[max_spread_radial] << ' '
                              << global_max_phi[max_spread_radial] << ' '
                              << tangential_ratio << '\n';
        shell_symmetry_output.flush();

        NohRadialProfileSample sample;
        sample.step_index = step_index;
        sample.time_s = sample_time;
        sample.exact_shock_radius = kUseNoh2 ? 0.0 : sample_time / 3.0;
        sample.numerical_shock_radius =
            (!kUseNoh2 && sample_time > 0.0)
                ? EstimateNohShockRadiusFromShellAverage(
                      global_rho,
                      current_global_radial_faces,
                      sample.exact_shock_radius)
                : 0.0;
        sample.radial_centers.reserve(kRadialCells);
        for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
          sample.radial_centers.push_back(
              0.5 * (current_global_radial_faces[radial] +
                     current_global_radial_faces[radial + 1u]));
        }
        sample.rho_avg = std::move(global_rho);
        sample.radial_velocity_avg = std::move(global_radial_velocity);
        sample.velocity_magnitude_avg = std::move(global_velocity);
        sample.mom_r_avg = std::move(global_mom_r);
        sample.e_fluid_total_avg = std::move(global_e_total);

        std::ostringstream file_name;
        file_name << "radial_profile_step"
                  << std::setw(6) << std::setfill('0') << step_index
                  << ".txt";
        std::ostringstream label;
        label << "step_" << step_index;
        if (!WriteNohRadialProfileSampleFile(output_dir / file_name.str(), label.str().c_str(), sample)) {
          return false;
        }
        profile_history.push_back(std::move(sample));
      }
      return true;
    };

    if (!record_profile_snapshot(0u, 0.0)) {
      if (rank == 0) {
        std::cerr << "failed to record initial profile\n";
      }
      MPI_Finalize();
      return 1;
    }

    while (time_s + 1.0e-16 < kFinalTime && step < kMaxSteps) {
      dec3d::hydro::HydroOperator dt_hydro;
      dt_hydro.SetStaticGridOptions(options);
      const dec3d::core::StageContext probe_context{
          time_s,
          1.0e-12,
          static_cast<std::uint64_t>(step + 1u),
          dec3d::core::PhaseId::p1,
          "p1-v0.1",
          DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL,
          "ownership:rank" + std::to_string(rank),
          DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL};
      if (!dt_hydro.bind(probe_context, local_geometry, local_state)) {
        if (rank == 0) {
          std::cerr << "bind failed for dt estimate\n";
        }
        MPI_Finalize();
        return 1;
      }

      const auto local_dt = dt_hydro.estimate_dt();
      double step_dt = local_dt.hard_cap_dt;
      MPI_Allreduce(MPI_IN_PLACE, &step_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      step_dt = std::min(step_dt, kFinalTime - time_s);

      const auto rho_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.rho, kGhostLayers);
      const auto mom_r_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_r, kGhostLayers);
      const auto mom_theta_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_theta, kGhostLayers);
      const auto mom_phi_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.mom_phi, kGhostLayers);
      const auto e_total_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_fluid_total, kGhostLayers);
      const auto e_electron_exchange = dec3d::mesh::ExchangeRadialHaloField(local_state.e_electron, kGhostLayers);

      auto local_hydro_view = dec3d::state::BuildHydroWorkView(local_state);
      if (!local_hydro_view.is_complete()) {
        return abort_with_rank(local_hydro_view.failure_reason);
      }

      auto local_override = dec3d::testsupport::BuildRadialGhostOverride(
          rho_exchange,
          mom_r_exchange,
          mom_theta_exchange,
          mom_phi_exchange,
          e_total_exchange,
          e_electron_exchange);
      if (kUseNoh2) {
        ApplyOuterNoh2InflowGhostIfNeeded(
            local_override,
            rank == rank_count - 1,
            kThetaCells,
            kPhiCells,
            local_geometry,
            time_s,
            kRho0,
            kNoh2SpecificInternalEnergy0,
            kElectronEnergyFraction);
      } else {
        ApplyOuterInflowGhostIfNeeded(
            local_override,
            rank == rank_count - 1,
            kThetaCells,
            kPhiCells,
            local_geometry,
            time_s,
            kRho0,
            kPressure0,
            kRadialVelocity0,
            kElectronEnergyFraction);
      }

      dec3d::core::MeshUpdateProposal local_ale_proposal;
      const dec3d::core::MeshUpdateProposal* local_ale_proposal_ptr = nullptr;
      if (kUseAle) {
        const auto global_ale_result = dec3d::mesh::BuildRadialAleMpiGlobalProposal(
            local_state.rho,
            local_state.mom_r,
            decomposition,
            current_global_radial_faces,
            step_dt,
            1.0,
            MPI_COMM_WORLD);
        if (!global_ale_result.is_complete()) {
          if (rank == 0) {
            std::cerr << "global ALE proposal failed: "
                      << (global_ale_result.failure_reason.empty()
                              ? global_ale_result.report_line
                              : global_ale_result.failure_reason)
                      << '\n';
          }
          MPI_Finalize();
          return 1;
        }
        local_ale_proposal = global_ale_result.proposal;
        const auto local_ale_window = dec3d::mesh::BuildRadialAleLocalProposalWindow(
            local_ale_proposal,
            local_slice.begin_index,
            local_geometry.radial_faces.size());
        const int local_proposal_complete =
            local_ale_window.is_complete() ? 1 : 0;
        int global_proposal_complete = 0;
        MPI_Allreduce(
            &local_proposal_complete,
            &global_proposal_complete,
            1,
            MPI_INT,
            MPI_MIN,
            MPI_COMM_WORLD);
        if (global_proposal_complete == 0) {
          if (rank == 0) {
            std::cerr << "radial ALE global proposal local window is incomplete\n";
          }
          MPI_Finalize();
          return 1;
        }

        double local_max_face_speed = 0.0;
        for (const double face_speed : local_ale_window.radial_face_velocities) {
          local_max_face_speed = std::max(local_max_face_speed, std::abs(face_speed));
        }
        const bool owns_global_outer_face = local_slice.end_index == kRadialCells;
        const double local_outer_face_speed =
            !owns_global_outer_face || local_ale_window.radial_face_velocities.empty()
                ? 0.0
                : local_ale_window.radial_face_velocities.back();
        double global_max_face_speed = 0.0;
        double global_min_outer_face_speed = 0.0;
        double global_max_outer_face_speed = 0.0;
        MPI_Reduce(
            &local_max_face_speed,
            &global_max_face_speed,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0,
            MPI_COMM_WORLD);
        MPI_Reduce(
            &local_outer_face_speed,
            &global_min_outer_face_speed,
            1,
            MPI_DOUBLE,
            MPI_MIN,
            0,
            MPI_COMM_WORLD);
        MPI_Reduce(
            &local_outer_face_speed,
            &global_max_outer_face_speed,
            1,
            MPI_DOUBLE,
            MPI_MAX,
            0,
            MPI_COMM_WORLD);
        if (rank == 0) {
          moving_mesh_output << (step + 1u) << ' '
                             << std::setprecision(17) << time_s << ' '
                             << step_dt << ' '
                             << global_proposal_complete << ' '
                             << global_max_face_speed << ' '
                             << global_min_outer_face_speed << ' '
                             << global_max_outer_face_speed << '\n';
          moving_mesh_output.flush();
        }
        local_ale_proposal_ptr = &local_ale_proposal;
      }

      dec3d::hydro::StaticGridHydroResult local_result;
      if (kUseAle && options.use_macro_zoning) {
        const auto macro_ale_result = dec3d::hydro::AdvanceMacroAleHllcMpiStep(
            local_hydro_view,
            local_geometry,
            decomposition,
            local_ale_proposal,
            step_dt,
            local_override,
            options,
            MPI_COMM_WORLD);
        local_result = macro_ale_result.staged_hydro_result;
        if (!macro_ale_result.success) {
          return abort_with_rank(
              macro_ale_result.failure_reason.empty()
                  ? "macro ALE direct HLLC MPI step failed"
                  : macro_ale_result.failure_reason);
        }
        current_global_radial_faces = local_ale_proposal.proposed_radial_faces;
      } else {
        local_result = dec3d::hydro::AdvanceStaticGridHydro(
            local_hydro_view,
            local_geometry,
            step_dt,
            local_override,
            options,
            local_ale_proposal_ptr);
      }
      if (kWriteStageAngularDebug && step < kStageDebugSteps) {
        for (const auto& entry : local_result.diagnostics.entries) {
          if (entry.code == "p1.hydro.debug.angular.stage" ||
              entry.code == "p1.hydro.debug.ale.commit_density_balance" ||
              entry.code == "p1.hydro.debug.radial_ale.face_flux") {
            stage_angular_debug_output << (step + 1u) << ' '
                                       << std::setprecision(17) << time_s << ' '
                                       << step_dt << ' '
                                       << rank << ' '
                                       << entry.message << '\n';
          }
        }
        stage_angular_debug_output.flush();
      }
      if (!local_result.success) {
        return abort_with_rank(local_result.failure_reason);
      }

      const auto writeback = dec3d::state::CommitHydroWriteback(
          local_state,
          local_hydro_view,
          dec3d::state::BuildHydroAuthorizedWriteMask());
      if (!writeback.success) {
        return abort_with_rank(writeback.report_line.empty()
                                   ? writeback.failure_reason
                                   : writeback.report_line);
      }

      dec3d::testsupport::ReduceBudgetSummaryToRoot(local_result.budget, rank, &reduced_budget);

      time_s += step_dt;
      if (rank == 0) {
        progress_output << (step + 1u) << ' ' << std::setprecision(17) << time_s << " pending\n";
        progress_output.flush();
      }

      if (((step + 1u) % 20u) == 0u || time_s + 1.0e-16 >= kFinalTime) {
        if (!record_profile_snapshot(step + 1u, time_s)) {
          if (rank == 0) {
            std::cerr << "failed to record profile snapshot\n";
          }
          MPI_Finalize();
          return 1;
        }
      }

      ++step;
    }

    if (time_s + 1.0e-16 < kFinalTime) {
      if (rank == 0) {
        std::cerr << "Noh MPI benchmark hit max_steps before final_time\n";
      }
      MPI_Finalize();
      return 1;
    }

    auto gathered_after = dec3d::state::CanonicalState::Create(global_before.layout);
    dec3d::testsupport::GatherFieldToRoot(local_state.rho, decomposition, rank, gathered_after.rho);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_r, decomposition, rank, gathered_after.mom_r);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_theta, decomposition, rank, gathered_after.mom_theta);
    dec3d::testsupport::GatherFieldToRoot(local_state.mom_phi, decomposition, rank, gathered_after.mom_phi);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_fluid_total, decomposition, rank, gathered_after.e_fluid_total);
    dec3d::testsupport::GatherFieldToRoot(local_state.e_electron, decomposition, rank, gathered_after.e_electron);

    if (rank == 0) {
      if (!profile_history.empty()) {
        if (!WriteNohRadialProfileSampleFile(output_dir / "radial_profile_t0.txt", "t0", profile_history.front()) ||
            !WriteNohRadialProfileSampleFile(output_dir / "radial_profile_t1.txt", "t1", profile_history.back())) {
          std::cerr << "failed to write Noh endpoint profiles\n";
          MPI_Finalize();
          return 1;
        }
      }
      const auto& final_profile = profile_history.back();
      const bool wrote_profile_diagnostics =
          kUseNoh2
              ? WriteNoh2ProfileErrorHistory(
                    output_dir / "noh2_profile_error_vs_time.txt",
                    profile_history,
                    kRho0,
                    kNoh2SpecificInternalEnergy0)
              : (WriteNohShockHistory(output_dir / "shock_radius_vs_time.txt", profile_history) &&
                 WriteWallHeatingDiagnostics(output_dir / "wall_heating_diagnostics.txt", final_profile));
      if (!wrote_profile_diagnostics ||
          !WriteBudgetResidual(output_dir / "budget_residual.txt", reduced_budget)) {
        std::cerr << "failed to write Noh/Noh2 diagnostics\n";
        MPI_Finalize();
        return 1;
      }

      const auto symmetry = MeasureShellSymmetry(gathered_after);
      const auto noh2_profile_error =
          kUseNoh2
              ? MeasureNoh2ProfileError(
                    final_profile,
                    kRho0,
                    kNoh2SpecificInternalEnergy0)
              : Noh2ProfileErrorMetrics{};
      const double shock_relative_error =
          final_profile.exact_shock_radius > 0.0
              ? std::abs(final_profile.numerical_shock_radius - final_profile.exact_shock_radius) /
                    final_profile.exact_shock_radius
              : 0.0;
      const bool shock_radius_within_tolerance =
          kUseNoh2 ? true : shock_relative_error <= kShockTolerance;
      const bool noh2_profile_within_tolerance =
          !kUseNoh2 ||
          (noh2_profile_error.rho_relative_l1 <= kNoh2ProfileTolerance &&
           noh2_profile_error.radial_velocity_relative_l1 <= kNoh2ProfileTolerance &&
           noh2_profile_error.e_fluid_total_relative_l1 <= kNoh2ProfileTolerance);
      const bool shell_symmetry_preserved =
          symmetry.max_rho_angular_relative_spread <= kShellSymmetryTolerance;
      const bool tangential_momentum_quiet =
          symmetry.tangential_to_radial_momentum_ratio <= kTangentialMomentumTolerance;
      const bool noh_success =
          shock_radius_within_tolerance &&
          noh2_profile_within_tolerance &&
          shell_symmetry_preserved &&
          tangential_momentum_quiet;

      std::ostringstream report;
      report << std::setprecision(6)
             << "noh_success=" << (noh_success ? "true" : "false")
             << "; shock_radius_within_tolerance=" << (shock_radius_within_tolerance ? "true" : "false")
             << "; noh2_profile_within_tolerance=" << (noh2_profile_within_tolerance ? "true" : "false")
             << "; shell_symmetry_preserved=" << (shell_symmetry_preserved ? "true" : "false")
             << "; tangential_momentum_quiet=" << (tangential_momentum_quiet ? "true" : "false")
             << "; final_time_s=" << time_s
             << "; numerical_shock_radius=" << final_profile.numerical_shock_radius
             << "; exact_shock_radius=" << final_profile.exact_shock_radius
             << "; shock_radius_relative_error=" << shock_relative_error
             << "; noh2_rho_relative_l1=" << noh2_profile_error.rho_relative_l1
             << "; noh2_radial_velocity_relative_l1="
             << noh2_profile_error.radial_velocity_relative_l1
             << "; noh2_e_fluid_total_relative_l1="
             << noh2_profile_error.e_fluid_total_relative_l1
             << "; max_rho_angular_relative_spread=" << symmetry.max_rho_angular_relative_spread
             << "; max_spread_radial_index=" << symmetry.max_spread_radial_index
             << "; min_rho_theta_phi=" << symmetry.min_rho_theta_index << ',' << symmetry.min_rho_phi_index
             << "; max_rho_theta_phi=" << symmetry.max_rho_theta_index << ',' << symmetry.max_rho_phi_index
             << "; tangential_to_radial_momentum_ratio=" << symmetry.tangential_to_radial_momentum_ratio;

      std::ofstream summary(output_dir / "summary.txt");
      summary << std::setprecision(17);
      summary << "case=" << DEC3D_NOH_MPI24_CASE_NAME_LITERAL << '\n';
      summary << "case_model=" << (kUseNoh2 ? "noh2" : "noh") << '\n';
      summary << "ale_enabled=" << (kUseAle ? "true" : "false") << '\n';
      summary << "success=" << (noh_success ? "true" : "false") << '\n';
      summary << "shock_radius_within_tolerance=" << (shock_radius_within_tolerance ? "true" : "false") << '\n';
      summary << "noh2_profile_within_tolerance=" << (noh2_profile_within_tolerance ? "true" : "false") << '\n';
      summary << "shell_symmetry_preserved=" << (shell_symmetry_preserved ? "true" : "false") << '\n';
      summary << "tangential_momentum_quiet=" << (tangential_momentum_quiet ? "true" : "false") << '\n';
      summary << "final_time_s=" << time_s << '\n';
      summary << "numerical_shock_radius=" << final_profile.numerical_shock_radius << '\n';
      summary << "exact_shock_radius=" << final_profile.exact_shock_radius << '\n';
      summary << "shock_radius_relative_error=" << shock_relative_error << '\n';
      summary << "noh2_rho_relative_l1=" << noh2_profile_error.rho_relative_l1 << '\n';
      summary << "noh2_radial_velocity_relative_l1="
              << noh2_profile_error.radial_velocity_relative_l1 << '\n';
      summary << "noh2_e_fluid_total_relative_l1="
              << noh2_profile_error.e_fluid_total_relative_l1 << '\n';
      summary << "noh2_rho_relative_linf=" << noh2_profile_error.rho_relative_linf << '\n';
      summary << "noh2_radial_velocity_relative_linf="
              << noh2_profile_error.radial_velocity_relative_linf << '\n';
      summary << "noh2_e_fluid_total_relative_linf="
              << noh2_profile_error.e_fluid_total_relative_linf << '\n';
      summary << "max_rho_angular_relative_spread=" << symmetry.max_rho_angular_relative_spread << '\n';
      summary << "max_spread_radial_index=" << symmetry.max_spread_radial_index << '\n';
      summary << "min_rho_theta_phi=" << symmetry.min_rho_theta_index << ',' << symmetry.min_rho_phi_index << '\n';
      summary << "max_rho_theta_phi=" << symmetry.max_rho_theta_index << ',' << symmetry.max_rho_phi_index << '\n';
      summary << "tangential_to_radial_momentum_ratio=" << symmetry.tangential_to_radial_momentum_ratio << '\n';
      summary << "report_line=" << report.str() << '\n';

      std::ofstream manifest(output_dir / "case_manifest.txt");
      manifest << "# kind=case_manifest\n";
      manifest << "case=" << DEC3D_NOH_MPI24_CASE_NAME_LITERAL << '\n';
      manifest << "summary=summary.txt\n";
      manifest << "run_progress=dec3d.out\n";
      if (kUseNoh2) {
        manifest << "noh2_profile_error_history=noh2_profile_error_vs_time.txt\n";
      } else {
        manifest << "shock_radius_history=shock_radius_vs_time.txt\n";
      }
      manifest << "shell_symmetry_history=shell_symmetry_vs_time.txt\n";
      if (kUseAle) {
        manifest << "moving_mesh_diagnostics=moving_mesh_diagnostics.txt\n";
      }
      manifest << "radial_profile_series=radial_profile_step*.txt\n";
      manifest << "radial_profile_t0=radial_profile_t0.txt\n";
      manifest << "radial_profile_t1=radial_profile_t1.txt\n";
      manifest << "budget_residual=budget_residual.txt\n";
      if (!kUseNoh2) {
        manifest << "wall_heating_diagnostics=wall_heating_diagnostics.txt\n";
      }

      std::cout << report.str() << '\n';
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    if (mpi_initialized != 0) {
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    std::cerr << error.what() << '\n';
    return 1;
  }
}
