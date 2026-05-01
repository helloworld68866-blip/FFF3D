#include "mesh/macro_zoning/macro_zoning.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::mesh {
namespace {

[[nodiscard]] bool IsPositivePowerOfTwo(std::size_t value) noexcept {
  return value > 0 && std::has_single_bit(value);
}

[[nodiscard]] std::size_t LinearIndex(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  return ((radial * theta_cells) + theta) * phi_cells + phi;
}

[[nodiscard]] std::size_t SmallestPowerOfTwoAtLeast(double target_ratio, std::size_t max_factor) noexcept {
  std::size_t factor = 1;
  while (factor < max_factor && static_cast<double>(factor) + 1.0e-12 < target_ratio) {
    factor *= 2;
  }
  return factor;
}

[[nodiscard]] bool GeometryHasCellCounts(const SphericalGeometryMetadata& geometry) noexcept {
  return geometry.radial_faces.size() >= 2 &&
         geometry.theta_faces.size() >= 2 &&
         geometry.phi_faces.size() >= 2 &&
         geometry.cell_volumes.size() ==
             (geometry.radial_faces.size() - 1) *
                 (geometry.theta_faces.size() - 1) *
                 (geometry.phi_faces.size() - 1);
}

}  // namespace

bool MacroZoneThetaBand::is_valid(
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) const noexcept {
  return radial_index < radial_cells &&
         theta_begin < theta_end &&
         theta_end <= theta_cells &&
         theta_factor > 0 &&
         IsPositivePowerOfTwo(theta_factor) &&
         (theta_end - theta_begin) == theta_factor &&
         phi_factor > 0 &&
         phi_factor <= phi_cells &&
         IsPositivePowerOfTwo(phi_factor) &&
         delta_s_theta > 0.0 &&
         min_delta_s_phi > 0.0;
}

bool MacroZoneCell::is_valid(
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) const noexcept {
  return initialized &&
         radial_index < radial_cells &&
         theta_begin < theta_end &&
         theta_end <= theta_cells &&
         phi_begin < phi_end &&
         phi_end <= phi_cells &&
         theta_factor > 0 &&
         phi_factor > 0 &&
         (theta_end - theta_begin) == theta_factor &&
         (phi_end - phi_begin) == phi_factor &&
         IsPositivePowerOfTwo(theta_factor) &&
         IsPositivePowerOfTwo(phi_factor) &&
         total_volume > 0.0;
}

bool MacroZoneMap::is_complete() const noexcept {
  if (!valid ||
      fine_radial_cells == 0 ||
      fine_theta_cells == 0 ||
      fine_phi_cells == 0 ||
      fine_to_coarse.empty() ||
      fine_to_coarse.extent_r() != fine_radial_cells ||
      fine_to_coarse.extent_theta() != fine_theta_cells ||
      fine_to_coarse.extent_phi() != fine_phi_cells ||
      coarse_cells.empty() ||
      report_line.empty()) {
    return false;
  }

  for (const auto& band : theta_bands) {
    if (!band.is_valid(fine_radial_cells, fine_theta_cells, fine_phi_cells)) {
      return false;
    }
  }

  for (const auto& cell : coarse_cells) {
    if (!cell.is_valid(fine_radial_cells, fine_theta_cells, fine_phi_cells)) {
      return false;
    }
  }

  for (std::size_t radial = 0; radial < fine_radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < fine_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < fine_phi_cells; ++phi) {
        if (fine_to_coarse(radial, theta, phi) >= coarse_cells.size()) {
          return false;
        }
      }
    }
  }

  return true;
}

bool FineHydroPackageView::is_complete() const noexcept {
  if (rho == nullptr ||
      mom_r == nullptr ||
      mom_theta == nullptr ||
      mom_phi == nullptr ||
      e_fluid_total == nullptr ||
      chi_e == nullptr ||
      rho->empty()) {
    return false;
  }

  return rho->extent_r() == mom_r->extent_r() &&
         rho->extent_r() == mom_theta->extent_r() &&
         rho->extent_r() == mom_phi->extent_r() &&
         rho->extent_r() == e_fluid_total->extent_r() &&
         rho->extent_r() == chi_e->extent_r() &&
         rho->extent_theta() == mom_r->extent_theta() &&
         rho->extent_theta() == mom_theta->extent_theta() &&
         rho->extent_theta() == mom_phi->extent_theta() &&
         rho->extent_theta() == e_fluid_total->extent_theta() &&
         rho->extent_theta() == chi_e->extent_theta() &&
         rho->extent_phi() == mom_r->extent_phi() &&
         rho->extent_phi() == mom_theta->extent_phi() &&
         rho->extent_phi() == mom_phi->extent_phi() &&
         rho->extent_phi() == e_fluid_total->extent_phi() &&
         rho->extent_phi() == chi_e->extent_phi();
}

bool MacroZonedHydroPackage::is_complete() const noexcept {
  return valid &&
         short_lived_work_view &&
         !authoritative_backed &&
         map.is_complete() &&
         rho.size() == map.coarse_cells.size() &&
         mom_r.size() == map.coarse_cells.size() &&
         mom_theta.size() == map.coarse_cells.size() &&
         mom_phi.size() == map.coarse_cells.size() &&
         e_fluid_total.size() == map.coarse_cells.size() &&
         chi_e.size() == map.coarse_cells.size() &&
         report_line.size() > 0;
}

bool FineHydroPackageBuffers::is_complete() const noexcept {
  return !rho.empty() &&
         rho.extent_r() == mom_r.extent_r() &&
         rho.extent_r() == mom_theta.extent_r() &&
         rho.extent_r() == mom_phi.extent_r() &&
         rho.extent_r() == e_fluid_total.extent_r() &&
         rho.extent_r() == chi_e.extent_r() &&
         rho.extent_theta() == mom_r.extent_theta() &&
         rho.extent_theta() == mom_theta.extent_theta() &&
         rho.extent_theta() == mom_phi.extent_theta() &&
         rho.extent_theta() == e_fluid_total.extent_theta() &&
         rho.extent_theta() == chi_e.extent_theta() &&
         rho.extent_phi() == mom_r.extent_phi() &&
         rho.extent_phi() == mom_theta.extent_phi() &&
         rho.extent_phi() == mom_phi.extent_phi() &&
         rho.extent_phi() == e_fluid_total.extent_phi() &&
         rho.extent_phi() == chi_e.extent_phi();
}

MacroZoneMap DetectMacroZones(
    const SphericalGeometryMetadata& geometry,
    const MacroZoningDetectionOptions& options) noexcept {
  MacroZoneMap map;

  if (!geometry.is_valid() || !GeometryHasCellCounts(geometry) || options.coarse_factor <= 0.0) {
    map.failure_reason = "macro-zoning requires valid spherical geometry and positive coarse factor";
    return map;
  }

  map.fine_radial_cells = geometry.radial_faces.size() - 1;
  map.fine_theta_cells = geometry.theta_faces.size() - 1;
  map.fine_phi_cells = geometry.phi_faces.size() - 1;

  if (!IsPositivePowerOfTwo(map.fine_theta_cells) || !IsPositivePowerOfTwo(map.fine_phi_cells)) {
    map.failure_reason = "macro-zoning requires power-of-two theta/phi resolution for binary recombination";
    return map;
  }

  map.fine_to_coarse = dec3d::core::Array3D<std::size_t>(
      map.fine_radial_cells,
      map.fine_theta_cells,
      map.fine_phi_cells,
      InvalidMacroZoneIndex());

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    const double radial_inner = geometry.radial_faces[radial];
    const double radial_outer = geometry.radial_faces[radial + 1];
    const double radial_center = 0.5 * (radial_inner + radial_outer);
    const double delta_r = radial_outer - radial_inner;
    const double delta_s_threshold = options.coarse_factor * delta_r;
    const double delta_theta = geometry.theta_faces[1] - geometry.theta_faces[0];
    const double theta_arc = radial_center * delta_theta;
    const std::size_t theta_factor = SmallestPowerOfTwoAtLeast(
        delta_s_threshold / std::max(theta_arc, 1.0e-16),
        map.fine_theta_cells);

    if (theta_arc * static_cast<double>(theta_factor) + 1.0e-12 < delta_s_threshold) {
      // Near the origin the full theta extent can still be smaller than the
      // radial threshold.  The thesis macro-zone operation can only recombine
      // existing cells, so saturate at the largest legal angular zone instead
      // of failing the whole hydro step.
      ++map.theta_factor_saturated_count;
    }

    for (std::size_t theta_begin = 0; theta_begin < map.fine_theta_cells; theta_begin += theta_factor) {
      const std::size_t theta_end = theta_begin + theta_factor;

      double min_delta_s_phi = std::numeric_limits<double>::max();
      for (std::size_t theta = theta_begin; theta < theta_end; ++theta) {
        const double theta_center = 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1]);
        const double delta_phi = geometry.phi_faces[1] - geometry.phi_faces[0];
        const double delta_s_phi = radial_center * std::sin(theta_center) * delta_phi;
        min_delta_s_phi = std::min(min_delta_s_phi, delta_s_phi);
      }

      const std::size_t phi_factor = SmallestPowerOfTwoAtLeast(
          delta_s_threshold / std::max(min_delta_s_phi, 1.0e-16),
          map.fine_phi_cells);

      if (min_delta_s_phi * static_cast<double>(phi_factor) + 1.0e-12 < delta_s_threshold) {
        // At the polar axis, even an all-phi macro-zone may not reach dr/2 in
        // azimuthal arc length.  Full-phi saturation is the only available
        // topological coarsening there; treating this as a failure prevents
        // legitimate full-sphere runs from entering the hydro path.
        ++map.phi_factor_saturated_count;
      }

      map.theta_bands.push_back(MacroZoneThetaBand{
          radial,
          theta_begin,
          theta_end,
          theta_factor,
          phi_factor,
          theta_arc,
          min_delta_s_phi});

      for (std::size_t phi_begin = 0; phi_begin < map.fine_phi_cells; phi_begin += phi_factor) {
        const std::size_t phi_end = phi_begin + phi_factor;
        MacroZoneCell cell;
        cell.initialized = true;
        cell.radial_index = radial;
        cell.theta_begin = theta_begin;
        cell.theta_end = theta_end;
        cell.phi_begin = phi_begin;
        cell.phi_end = phi_end;
        cell.theta_factor = theta_factor;
        cell.phi_factor = phi_factor;

        for (std::size_t theta = theta_begin; theta < theta_end; ++theta) {
          for (std::size_t phi = phi_begin; phi < phi_end; ++phi) {
            const auto linear_index = LinearIndex(
                radial,
                theta,
                phi,
                map.fine_theta_cells,
                map.fine_phi_cells);
            cell.total_volume += geometry.cell_volumes.at(linear_index);
          }
        }

        const std::size_t coarse_index = map.coarse_cells.size();
        map.coarse_cells.push_back(cell);

        for (std::size_t theta = theta_begin; theta < theta_end; ++theta) {
          for (std::size_t phi = phi_begin; phi < phi_end; ++phi) {
            map.fine_to_coarse(radial, theta, phi) = coarse_index;
          }
        }
      }
    }
  }

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    std::size_t last_phi_factor = 0;
    bool have_band = false;
    for (const auto& band : map.theta_bands) {
      if (band.radial_index != radial) {
        continue;
      }
      if (have_band && band.phi_factor != last_phi_factor) {
        map.phi_factor_varies_by_theta_band = true;
      }
      have_band = true;
      last_phi_factor = band.phi_factor;
    }
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "macro_zoning_detected=true"
         << " radial_cells=" << map.fine_radial_cells
         << " theta_cells=" << map.fine_theta_cells
         << " phi_cells=" << map.fine_phi_cells
         << " coarse_cells=" << map.coarse_cells.size()
         << " theta_bands=" << map.theta_bands.size()
         << " phi_factor_varies_by_theta_band=" << (map.phi_factor_varies_by_theta_band ? "true" : "false")
         << " theta_factor_saturated_count=" << map.theta_factor_saturated_count
         << " phi_factor_saturated_count=" << map.phi_factor_saturated_count;
  map.report_line = report.str();
  map.valid = true;
  return map;
}

MacroZonedHydroPackage RestrictFineHydroPackage(
    const FineHydroPackageView& fine_view,
    const SphericalGeometryMetadata& geometry,
    const MacroZoneMap& map) noexcept {
  MacroZonedHydroPackage coarse_package;
  coarse_package.map = map;

  if (!fine_view.is_complete()) {
    coarse_package.failure_reason = "macro-zoning restrict requires a complete fine hydro package";
    return coarse_package;
  }

  if (!geometry.is_valid() || !GeometryHasCellCounts(geometry) || !map.is_complete()) {
    coarse_package.failure_reason = "macro-zoning restrict requires valid geometry and coarse map";
    return coarse_package;
  }

  if (fine_view.rho->extent_r() != map.fine_radial_cells ||
      fine_view.rho->extent_theta() != map.fine_theta_cells ||
      fine_view.rho->extent_phi() != map.fine_phi_cells) {
    coarse_package.failure_reason = "macro-zoning restrict requires fine package extents to match coarse map";
    return coarse_package;
  }

  coarse_package.rho.assign(map.coarse_cells.size(), 0.0);
  coarse_package.mom_r.assign(map.coarse_cells.size(), 0.0);
  coarse_package.mom_theta.assign(map.coarse_cells.size(), 0.0);
  coarse_package.mom_phi.assign(map.coarse_cells.size(), 0.0);
  coarse_package.e_fluid_total.assign(map.coarse_cells.size(), 0.0);
  coarse_package.chi_e.assign(map.coarse_cells.size(), 0.0);
  std::vector<bool> has_reference(map.coarse_cells.size(), false);
  std::vector<double> reference_rho(map.coarse_cells.size(), 0.0);
  std::vector<double> reference_mom_r(map.coarse_cells.size(), 0.0);
  std::vector<double> reference_mom_theta(map.coarse_cells.size(), 0.0);
  std::vector<double> reference_mom_phi(map.coarse_cells.size(), 0.0);
  std::vector<double> reference_e_fluid_total(map.coarse_cells.size(), 0.0);
  std::vector<double> reference_chi_e(map.coarse_cells.size(), 0.0);

  for (std::size_t radial = 0; radial < map.fine_radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < map.fine_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < map.fine_phi_cells; ++phi) {
        const std::size_t coarse_index = map.fine_to_coarse(radial, theta, phi);
        const double cell_volume = geometry.cell_volumes.at(
            LinearIndex(radial, theta, phi, map.fine_theta_cells, map.fine_phi_cells));
        const double rho = (*fine_view.rho)(radial, theta, phi);
        const double mom_r = (*fine_view.mom_r)(radial, theta, phi);
        const double mom_theta = (*fine_view.mom_theta)(radial, theta, phi);
        const double mom_phi = (*fine_view.mom_phi)(radial, theta, phi);
        const double e_fluid_total = (*fine_view.e_fluid_total)(radial, theta, phi);
        const double chi_e = (*fine_view.chi_e)(radial, theta, phi);
        if (!has_reference[coarse_index]) {
          has_reference[coarse_index] = true;
          reference_rho[coarse_index] = rho;
          reference_mom_r[coarse_index] = mom_r;
          reference_mom_theta[coarse_index] = mom_theta;
          reference_mom_phi[coarse_index] = mom_phi;
          reference_e_fluid_total[coarse_index] = e_fluid_total;
          reference_chi_e[coarse_index] = chi_e;
        }
        coarse_package.rho[coarse_index] +=
            (rho - reference_rho[coarse_index]) * cell_volume;
        coarse_package.mom_r[coarse_index] +=
            (mom_r - reference_mom_r[coarse_index]) * cell_volume;
        coarse_package.mom_theta[coarse_index] +=
            (mom_theta - reference_mom_theta[coarse_index]) * cell_volume;
        coarse_package.mom_phi[coarse_index] +=
            (mom_phi - reference_mom_phi[coarse_index]) * cell_volume;
        coarse_package.e_fluid_total[coarse_index] +=
            (e_fluid_total - reference_e_fluid_total[coarse_index]) *
            cell_volume;
        coarse_package.chi_e[coarse_index] +=
            (chi_e - reference_chi_e[coarse_index]) * cell_volume;
      }
    }
  }

  for (std::size_t coarse_index = 0; coarse_index < map.coarse_cells.size(); ++coarse_index) {
    const double coarse_volume = map.coarse_cells[coarse_index].total_volume;
    coarse_package.rho[coarse_index] =
        reference_rho[coarse_index] +
        coarse_package.rho[coarse_index] / coarse_volume;
    coarse_package.mom_r[coarse_index] =
        reference_mom_r[coarse_index] +
        coarse_package.mom_r[coarse_index] / coarse_volume;
    coarse_package.mom_theta[coarse_index] =
        reference_mom_theta[coarse_index] +
        coarse_package.mom_theta[coarse_index] / coarse_volume;
    coarse_package.mom_phi[coarse_index] =
        reference_mom_phi[coarse_index] +
        coarse_package.mom_phi[coarse_index] / coarse_volume;
    coarse_package.e_fluid_total[coarse_index] =
        reference_e_fluid_total[coarse_index] +
        coarse_package.e_fluid_total[coarse_index] / coarse_volume;
    coarse_package.chi_e[coarse_index] =
        reference_chi_e[coarse_index] +
        coarse_package.chi_e[coarse_index] / coarse_volume;
  }

  std::ostringstream report;
  report << std::setprecision(17)
         << "macro_zoning_restrict=true"
         << " coarse_cells=" << map.coarse_cells.size()
         << " short_lived_work_view=true"
         << " authoritative_backed=false";
  coarse_package.report_line = report.str();
  coarse_package.valid = true;
  return coarse_package;
}

FineHydroPackageBuffers ProlongMacroZoneHydroPackage(
    const MacroZonedHydroPackage& coarse_package) noexcept {
  FineHydroPackageBuffers fine_buffers;

  if (!coarse_package.is_complete()) {
    fine_buffers.failure_reason = "macro-zoning prolong requires a complete coarse hydro package";
    return fine_buffers;
  }

  fine_buffers.rho = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);
  fine_buffers.mom_r = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);
  fine_buffers.mom_theta = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);
  fine_buffers.mom_phi = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);
  fine_buffers.e_fluid_total = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);
  fine_buffers.chi_e = dec3d::core::Array3D<double>(
      coarse_package.map.fine_radial_cells,
      coarse_package.map.fine_theta_cells,
      coarse_package.map.fine_phi_cells,
      0.0);

  for (std::size_t radial = 0; radial < coarse_package.map.fine_radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < coarse_package.map.fine_theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < coarse_package.map.fine_phi_cells; ++phi) {
        const std::size_t coarse_index = coarse_package.map.fine_to_coarse(radial, theta, phi);
        fine_buffers.rho(radial, theta, phi) = coarse_package.rho[coarse_index];
        fine_buffers.mom_r(radial, theta, phi) = coarse_package.mom_r[coarse_index];
        fine_buffers.mom_theta(radial, theta, phi) = coarse_package.mom_theta[coarse_index];
        fine_buffers.mom_phi(radial, theta, phi) = coarse_package.mom_phi[coarse_index];
        fine_buffers.e_fluid_total(radial, theta, phi) = coarse_package.e_fluid_total[coarse_index];
        fine_buffers.chi_e(radial, theta, phi) = coarse_package.chi_e[coarse_index];
      }
    }
  }

  return fine_buffers;
}

MacroZoneAuthoritativeBoundaryResult RejectMacroZoneAuthoritativeCommit(
    const MacroZonedHydroPackage& coarse_package) noexcept {
  MacroZoneAuthoritativeBoundaryResult result;
  result.failure_reason =
      coarse_package.is_complete()
          ? "macro-zoned coarse package is a short-lived hydro work view and cannot become authoritative truth"
          : "macro-zoned coarse package is incomplete and cannot become authoritative truth";
  return result;
}

}  // namespace dec3d::mesh
