#include "hydro/driver/radial_overlap_remap.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>
#include <vector>

namespace dec3d::hydro {
namespace {

constexpr const char* kRemapOrder = "first_order_proposal_mapped_overlap";
constexpr double kCoverageTolerance = 1.0e-10;

[[nodiscard]] std::size_t LinearIndex(
    std::size_t radial,
    std::size_t theta,
    std::size_t phi,
    std::size_t theta_cells,
    std::size_t phi_cells) noexcept {
  return (radial * theta_cells + theta) * phi_cells + phi;
}

[[nodiscard]] bool SameFaces(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) noexcept {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    if (lhs[index] != rhs[index]) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] bool StrictlyIncreasing(
    const std::vector<double>& values) noexcept {
  if (values.size() < 2u) {
    return false;
  }
  for (std::size_t index = 1u; index < values.size(); ++index) {
    if (!(values[index] > values[index - 1u]) ||
        !std::isfinite(values[index])) {
      return false;
    }
  }
  return std::isfinite(values.front());
}

[[nodiscard]] double WedgeFactor(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double theta_factor =
      std::cos(geometry.theta_faces[theta]) -
      std::cos(geometry.theta_faces[theta + 1u]);
  const double phi_factor =
      geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return theta_factor * phi_factor / 3.0;
}

[[nodiscard]] double ShellWedgeVolume(
    double r_inner,
    double r_outer,
    double wedge_factor) noexcept {
  return (r_outer * r_outer * r_outer - r_inner * r_inner * r_inner) *
         wedge_factor;
}

[[nodiscard]] std::vector<double> ScaleRadiationBundle(
    const std::vector<double>& values,
    double factor) {
  std::vector<double> scaled;
  scaled.reserve(values.size());
  for (const double value : values) {
    scaled.push_back(value * factor);
  }
  return scaled;
}

[[nodiscard]] std::vector<double> AddRadiationBundles(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs) {
  if (lhs.size() != rhs.size()) {
    return {};
  }

  std::vector<double> sum(lhs.size(), 0.0);
  for (std::size_t index = 0; index < lhs.size(); ++index) {
    sum[index] = lhs[index] + rhs[index];
  }
  return sum;
}

[[nodiscard]] HydroConservativeState ScaleState(
    const HydroConservativeState& state,
    double factor) noexcept {
  return HydroConservativeState{
      state.rho * factor,
      state.mom_r * factor,
      state.mom_theta * factor,
      state.mom_phi * factor,
      state.e_fluid_total * factor,
      state.chi_e * factor,
      state.alpha_chi * factor,
      ScaleRadiationBundle(state.radiation_chi, factor)};
}

[[nodiscard]] HydroConservativeState AddState(
    const HydroConservativeState& lhs,
    const HydroConservativeState& rhs) noexcept {
  return HydroConservativeState{
      lhs.rho + rhs.rho,
      lhs.mom_r + rhs.mom_r,
      lhs.mom_theta + rhs.mom_theta,
      lhs.mom_phi + rhs.mom_phi,
      lhs.e_fluid_total + rhs.e_fluid_total,
      lhs.chi_e + rhs.chi_e,
      lhs.alpha_chi + rhs.alpha_chi,
      AddRadiationBundles(lhs.radiation_chi, rhs.radiation_chi)};
}

void AccumulateExtensiveTotals(
    const HydroConservativeState& state,
    double volume,
    HydroConservativeState& totals) noexcept {
  totals = AddState(totals, ScaleState(state, volume));
}

[[nodiscard]] double MaxDensityAngularSpread(
    const dec3d::core::Array3D<HydroConservativeState>& cells,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  double max_spread = 0.0;
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    double rho_min = std::numeric_limits<double>::infinity();
    double rho_max = -std::numeric_limits<double>::infinity();
    double rho_volume_sum = 0.0;
    double volume_sum = 0.0;
    for (std::size_t theta = 0; theta < cells.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < cells.extent_phi(); ++phi) {
        const std::size_t index = LinearIndex(
            radial,
            theta,
            phi,
            cells.extent_theta(),
            cells.extent_phi());
        const double volume = geometry.cell_volumes[index];
        const double rho = cells(radial, theta, phi).rho;
        rho_min = std::min(rho_min, rho);
        rho_max = std::max(rho_max, rho);
        rho_volume_sum += rho * volume;
        volume_sum += volume;
      }
    }
    if (volume_sum > 0.0) {
      const double mean = rho_volume_sum / volume_sum;
      const double spread =
          std::abs(mean) > 0.0 ? (rho_max - rho_min) / std::abs(mean) : 0.0;
      max_spread = std::max(max_spread, spread);
    }
  }
  return max_spread;
}

[[nodiscard]] double MaxTangentialMomentumRatio(
    const dec3d::core::Array3D<HydroConservativeState>& cells) noexcept {
  double max_ratio = 0.0;
  for (std::size_t radial = 0; radial < cells.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < cells.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < cells.extent_phi(); ++phi) {
        const auto& state = cells(radial, theta, phi);
        const double tangential =
            std::hypot(state.mom_theta, state.mom_phi);
        const double denominator =
            std::max(std::abs(state.mom_r), state.rho);
        if (denominator > 0.0) {
          max_ratio = std::max(max_ratio, tangential / denominator);
        }
      }
    }
  }
  return max_ratio;
}

[[nodiscard]] std::string BuildReportLine(
    const RadialOverlapRemapDiagnostics& diagnostics) {
  std::ostringstream report;
  report << std::setprecision(17)
         << "remap_order=" << diagnostics.remap_order
         << "; mass_residual=" << diagnostics.mass_residual
         << "; mom_r_residual=" << diagnostics.mom_r_residual
         << "; mom_theta_residual=" << diagnostics.mom_theta_residual
         << "; mom_phi_residual=" << diagnostics.mom_phi_residual
         << "; e_fluid_total_residual="
         << diagnostics.e_fluid_total_residual
         << "; chi_e_residual=" << diagnostics.chi_e_residual
         << "; max_volume_coverage_error="
         << diagnostics.max_volume_coverage_error
         << "; max_density_angular_spread_before_remap="
         << diagnostics.max_density_angular_spread_before_remap
         << "; max_density_angular_spread_after_remap="
         << diagnostics.max_density_angular_spread_after_remap
         << "; max_tangential_momentum_ratio_after_remap="
         << diagnostics.max_tangential_momentum_ratio_after_remap;
  return report.str();
}

}  // namespace

bool RadialOverlapRemapDiagnostics::is_complete() const noexcept {
  return executed &&
         remap_order == kRemapOrder &&
         std::isfinite(mass_residual) &&
         std::isfinite(mom_r_residual) &&
         std::isfinite(mom_theta_residual) &&
         std::isfinite(mom_phi_residual) &&
         std::isfinite(e_fluid_total_residual) &&
         std::isfinite(chi_e_residual) &&
         std::isfinite(max_volume_coverage_error) &&
         std::isfinite(max_density_angular_spread_before_remap) &&
         std::isfinite(max_density_angular_spread_after_remap) &&
         std::isfinite(max_tangential_momentum_ratio_after_remap) &&
         !report_line.empty();
}

bool RadialOverlapRemapResult::is_complete() const noexcept {
  return success &&
         remapped_cells.size() > 0u &&
         diagnostics.is_complete() &&
         failure_reason.empty();
}

RadialOverlapRemapResult RemapHydroStateRadiallyConservative(
    const dec3d::core::Array3D<HydroConservativeState>& old_cell_averages,
    const dec3d::mesh::SphericalGeometryMetadata& old_geometry,
    const dec3d::mesh::SphericalGeometryMetadata& new_geometry,
    const std::vector<double>& proposal_mapped_source_faces) noexcept {
  RadialOverlapRemapResult result;
  result.diagnostics.executed = true;
  result.diagnostics.remap_order = kRemapOrder;

  const std::size_t radial_cells = old_cell_averages.extent_r();
  const std::size_t theta_cells = old_cell_averages.extent_theta();
  const std::size_t phi_cells = old_cell_averages.extent_phi();

  const auto fail = [&result](std::string reason) {
    result.failure_reason = std::move(reason);
    result.diagnostics.report_line = BuildReportLine(result.diagnostics);
    return result;
  };

  if (old_cell_averages.empty() ||
      !old_geometry.is_valid() ||
      !new_geometry.is_valid()) {
    return fail("radial overlap remap requires complete state and geometry");
  }
  if (old_geometry.radial_faces.size() != radial_cells + 1u ||
      new_geometry.radial_faces.size() != radial_cells + 1u ||
      old_geometry.theta_faces.size() != theta_cells + 1u ||
      old_geometry.phi_faces.size() != phi_cells + 1u ||
      new_geometry.theta_faces.size() != theta_cells + 1u ||
      new_geometry.phi_faces.size() != phi_cells + 1u) {
    return fail("radial overlap remap geometry shape does not match state");
  }
  if (!SameFaces(old_geometry.theta_faces, new_geometry.theta_faces) ||
      !SameFaces(old_geometry.phi_faces, new_geometry.phi_faces)) {
    return fail("radial overlap remap requires unchanged angular geometry");
  }
  if (proposal_mapped_source_faces.size() != old_geometry.radial_faces.size() ||
      !StrictlyIncreasing(proposal_mapped_source_faces)) {
    return fail("radial overlap remap requires monotone proposal-mapped faces");
  }

  result.diagnostics.max_density_angular_spread_before_remap =
      MaxDensityAngularSpread(old_cell_averages, old_geometry);

  dec3d::core::Array3D<HydroConservativeState> remapped_extensive(
      radial_cells,
      theta_cells,
      phi_cells,
      HydroConservativeState{});
  result.remapped_cells = dec3d::core::Array3D<HydroConservativeState>(
      radial_cells,
      theta_cells,
      phi_cells,
      HydroConservativeState{});

  HydroConservativeState old_totals;
  HydroConservativeState new_totals;
  double max_coverage_error = 0.0;

  for (std::size_t theta = 0; theta < theta_cells; ++theta) {
    for (std::size_t phi = 0; phi < phi_cells; ++phi) {
      const double wedge = WedgeFactor(new_geometry, theta, phi);
      if (!(wedge > 0.0) || !std::isfinite(wedge)) {
        return fail("radial overlap remap encountered invalid angular wedge");
      }

      for (std::size_t old_radial = 0; old_radial < radial_cells; ++old_radial) {
        const std::size_t old_index = LinearIndex(
            old_radial,
            theta,
            phi,
            theta_cells,
            phi_cells);
        const double old_volume = old_geometry.cell_volumes[old_index];
        if (!(old_volume > 0.0) || !std::isfinite(old_volume)) {
          return fail("radial overlap remap encountered invalid old volume");
        }
        const auto old_extensive =
            ScaleState(old_cell_averages(old_radial, theta, phi), old_volume);
        old_totals = AddState(old_totals, old_extensive);

        const double source_inner =
            proposal_mapped_source_faces[old_radial];
        const double source_outer =
            proposal_mapped_source_faces[old_radial + 1u];
        const double source_volume =
            ShellWedgeVolume(source_inner, source_outer, wedge);
        if (!(source_volume > 0.0) || !std::isfinite(source_volume)) {
          return fail("radial overlap remap encountered invalid source volume");
        }

        double covered_volume = 0.0;
        for (std::size_t new_radial = 0; new_radial < radial_cells; ++new_radial) {
          const double overlap_inner =
              std::max(source_inner, new_geometry.radial_faces[new_radial]);
          const double overlap_outer =
              std::min(source_outer, new_geometry.radial_faces[new_radial + 1u]);
          if (overlap_outer <= overlap_inner) {
            continue;
          }
          const double overlap_volume =
              ShellWedgeVolume(overlap_inner, overlap_outer, wedge);
          if (!(overlap_volume > 0.0) || !std::isfinite(overlap_volume)) {
            return fail("radial overlap remap encountered invalid overlap volume");
          }
          covered_volume += overlap_volume;
          remapped_extensive(new_radial, theta, phi) = AddState(
              remapped_extensive(new_radial, theta, phi),
              ScaleState(old_extensive, overlap_volume / source_volume));
        }

        const double coverage_error =
            std::abs(covered_volume - source_volume) / source_volume;
        max_coverage_error = std::max(max_coverage_error, coverage_error);
      }
    }
  }

  result.diagnostics.max_volume_coverage_error = max_coverage_error;
  if (max_coverage_error > kCoverageTolerance) {
    return fail("radial overlap remap coverage is incomplete");
  }

  for (std::size_t radial = 0; radial < radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        const std::size_t index = LinearIndex(
            radial,
            theta,
            phi,
            theta_cells,
            phi_cells);
        const double new_volume = new_geometry.cell_volumes[index];
        if (!(new_volume > 0.0) || !std::isfinite(new_volume)) {
          return fail("radial overlap remap encountered invalid new volume");
        }
        const auto new_state =
            ScaleState(remapped_extensive(radial, theta, phi), 1.0 / new_volume);
        if (!new_state.is_finite() ||
            !RecoverPrimitiveState(new_state).is_physical()) {
          return fail("radial overlap remap produced a non-physical state");
        }
        result.remapped_cells(radial, theta, phi) = new_state;
        AccumulateExtensiveTotals(new_state, new_volume, new_totals);
      }
    }
  }

  result.diagnostics.mass_residual = new_totals.rho - old_totals.rho;
  result.diagnostics.mom_r_residual = new_totals.mom_r - old_totals.mom_r;
  result.diagnostics.mom_theta_residual =
      new_totals.mom_theta - old_totals.mom_theta;
  result.diagnostics.mom_phi_residual =
      new_totals.mom_phi - old_totals.mom_phi;
  result.diagnostics.e_fluid_total_residual =
      new_totals.e_fluid_total - old_totals.e_fluid_total;
  result.diagnostics.chi_e_residual =
      new_totals.chi_e - old_totals.chi_e;
  result.diagnostics.max_density_angular_spread_after_remap =
      MaxDensityAngularSpread(result.remapped_cells, new_geometry);
  result.diagnostics.max_tangential_momentum_ratio_after_remap =
      MaxTangentialMomentumRatio(result.remapped_cells);
  result.diagnostics.report_line = BuildReportLine(result.diagnostics);
  result.success = true;
  return result;
}

}  // namespace dec3d::hydro
