#include "hydro/source/geometric_source.hpp"

#include "hydro/riemann/hllc_solver.hpp"

#include <cmath>
#include <sstream>
#include <vector>

namespace dec3d::hydro {

namespace {

struct GridShape {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
};

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] HydroConservativeState LoadCellState(
    const dec3d::state::HydroStateView& hydro_view,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return {
      (*hydro_view.rho)(radial, theta, phi),
      (*hydro_view.mom_r)(radial, theta, phi),
      (*hydro_view.mom_theta)(radial, theta, phi),
      (*hydro_view.mom_phi)(radial, theta, phi),
      (*hydro_view.e_fluid_total)(radial, theta, phi),
      hydro_view.chi_e(radial, theta, phi),
      hydro_view.operator_local_alpha_chi
          ? hydro_view.alpha_chi(radial, theta, phi)
          : 0.0,
      [&]() {
        std::vector<double> bundle;
        bundle.reserve(hydro_view.radiation_chi.size());
        for (const auto& group_chi : hydro_view.radiation_chi) {
          bundle.push_back(group_chi(radial, theta, phi));
        }
        return bundle;
      }()};
}

[[nodiscard]] std::size_t LinearIndex(
    const GridShape& shape,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return ((radial * shape.theta_cells) + theta) * shape.phi_cells + phi;
}

[[nodiscard]] double CellCenterRadius(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

[[nodiscard]] double CellCenterTheta(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const GridShape& shape,
    std::size_t radial,
    std::size_t theta,
    std::size_t phi) noexcept {
  return geometry.cell_volumes[LinearIndex(shape, radial, theta, phi)];
}

[[nodiscard]] double RadialFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face,
    std::size_t theta,
    std::size_t phi) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  const double theta_lower = geometry.theta_faces[theta];
  const double theta_upper = geometry.theta_faces[theta + 1u];
  const double phi_lower = geometry.phi_faces[phi];
  const double phi_upper = geometry.phi_faces[phi + 1u];
  return radius * radius *
         (std::cos(theta_lower) - std::cos(theta_upper)) *
         (phi_upper - phi_lower);
}

[[nodiscard]] double RadialShellVolumeFactor(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  return ((radial_outer * radial_outer * radial_outer) -
          (radial_inner * radial_inner * radial_inner)) / 3.0;
}

[[nodiscard]] double RadialFaceAreaPerSolidAngle(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial_face) noexcept {
  const double radius = geometry.radial_faces[radial_face];
  return radius * radius;
}

[[nodiscard]] double ThetaFaceArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial,
    std::size_t theta_face,
    std::size_t phi) noexcept {
  const double radial_inner = geometry.radial_faces[radial];
  const double radial_outer = geometry.radial_faces[radial + 1u];
  const double phi_lower = geometry.phi_faces[phi];
  const double phi_upper = geometry.phi_faces[phi + 1u];
  return 0.5 *
         (radial_outer * radial_outer - radial_inner * radial_inner) *
         std::sin(geometry.theta_faces[theta_face]) *
         (phi_upper - phi_lower);
}

}  // namespace

bool GeometricSourceSnapshot::is_complete() const noexcept {
  return radial_cells > 0u &&
         theta_cells > 0u &&
         phi_cells > 0u &&
         states.size() == radial_cells * theta_cells * phi_cells;
}

bool GeometricSourceStepResult::is_complete() const noexcept {
  return !report_line.empty() && diagnostics.has_entries();
}

bool GeometricSourceTermsResult::is_complete(
    std::size_t expected_cell_count) const noexcept {
  return !report_line.empty() &&
         source_terms.size() == expected_cell_count &&
         updated_cell_count == expected_cell_count;
}

GeometricSourceSnapshot CaptureGeometricSourceSnapshot(
    const dec3d::state::HydroStateView& hydro_view) noexcept {
  GeometricSourceSnapshot snapshot;
  if (!hydro_view.is_complete()) {
    return snapshot;
  }

  snapshot.radial_cells = hydro_view.rho->extent_r();
  snapshot.theta_cells = hydro_view.rho->extent_theta();
  snapshot.phi_cells = hydro_view.rho->extent_phi();
  snapshot.states.resize(snapshot.radial_cells * snapshot.theta_cells * snapshot.phi_cells);

  const GridShape shape{
      snapshot.radial_cells,
      snapshot.theta_cells,
      snapshot.phi_cells};
  for (std::size_t radial = 0; radial < snapshot.radial_cells; ++radial) {
    for (std::size_t theta = 0; theta < snapshot.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < snapshot.phi_cells; ++phi) {
        snapshot.states[LinearIndex(shape, radial, theta, phi)] =
            LoadCellState(hydro_view, radial, theta, phi);
      }
    }
  }

  return snapshot;
}

GeometricSourceTermsResult ComputeGeometricSourceTerms(
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) noexcept {
  GeometricSourceTermsResult result;
  if (!source_snapshot.is_complete()) {
    result.failure_reason = "geometric source snapshot is incomplete";
  } else if (!geometry.is_valid()) {
    result.failure_reason = "spherical geometry metadata is invalid for geometric source step";
  } else {
    const GridShape shape{
        source_snapshot.radial_cells,
        source_snapshot.theta_cells,
        source_snapshot.phi_cells};
    result.source_terms.assign(source_snapshot.states.size(), HydroConservativeState{});

    auto load_snapshot =
        [&source_snapshot, &shape](std::size_t radial, std::size_t theta, std::size_t phi) noexcept {
          return source_snapshot.states[LinearIndex(shape, radial, theta, phi)];
        };

    for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
      if (!result.failure_reason.empty()) {
        break;
      }
      const double radius = CellCenterRadius(geometry, radial);
      if (!(radius > 0.0) || !std::isfinite(radius)) {
        result.failure_reason = "cell-center radius must remain positive for geometric source step";
        break;
      }

      for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
        const double theta_center = CellCenterTheta(geometry, theta);
        const double sine_theta = std::sin(theta_center);
        if (!(std::abs(sine_theta) > 1.0e-12) || !std::isfinite(theta_center)) {
          result.failure_reason = "cell-center theta must avoid singular geometric source evaluation";
          break;
        }

        const double cot_theta = std::cos(theta_center) / sine_theta;
        for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
          const auto conservative = load_snapshot(radial, theta, phi);
          const auto primitive = RecoverPrimitiveState(conservative);
          if (!primitive.is_physical()) {
            result.failure_reason = "geometric source step encountered a non-physical hydro state";
            break;
          }

          const double volume = CellVolume(geometry, shape, radial, theta, phi);
          const double radial_volume_factor = RadialShellVolumeFactor(geometry, radial);
          if (!(radial_volume_factor > 0.0) || !std::isfinite(radial_volume_factor)) {
            result.failure_reason =
                "radial shell volume factor must remain positive for geometric source step";
            break;
          }
          const double radial_pressure_source =
              primitive.pressure *
              (RadialFaceAreaPerSolidAngle(geometry, radial + 1u) -
               RadialFaceAreaPerSolidAngle(geometry, radial)) /
              radial_volume_factor;
          const double theta_pressure_source =
              primitive.pressure *
              (ThetaFaceArea(geometry, radial, theta + 1u, phi) -
               ThetaFaceArea(geometry, radial, theta, phi)) /
              volume;
          const double source_mom_r =
              radial_pressure_source +
              ((conservative.mom_theta * conservative.mom_theta +
                conservative.mom_phi * conservative.mom_phi) /
               (conservative.rho * radius));
          const double source_mom_theta =
              (-(conservative.mom_theta * conservative.mom_r) /
               (conservative.rho * radius)) +
              ((cot_theta / radius) *
               ((conservative.mom_phi * conservative.mom_phi) / conservative.rho)) +
              theta_pressure_source;
          const double source_mom_phi =
              (-(conservative.mom_r * conservative.mom_phi) /
               (conservative.rho * radius)) -
              ((cot_theta / radius) *
               ((conservative.mom_theta * conservative.mom_phi) / conservative.rho));

          auto& source = result.source_terms[LinearIndex(shape, radial, theta, phi)];
          source.mom_r = source_mom_r;
          source.mom_theta = source_mom_theta;
          source.mom_phi = source_mom_phi;
          ++result.updated_cell_count;
        }

        if (!result.failure_reason.empty()) {
          break;
        }
      }

      if (!result.failure_reason.empty()) {
        break;
      }
    }

    result.success = result.failure_reason.empty() &&
                     result.source_terms.size() == source_snapshot.states.size();
  }

  std::ostringstream report;
  report << "geometric_source_terms_success=" << (result.success ? "true" : "false")
         << "; updated_cell_count=" << result.updated_cell_count;
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

GeometricSourceStepResult ApplyGeometricSourceStep(
    dec3d::state::HydroStateView& hydro_view,
    const GeometricSourceSnapshot& source_snapshot,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    double dt_s) noexcept {
  GeometricSourceStepResult result;

  if (!hydro_view.is_complete()) {
    result.failure_reason = "hydro view is incomplete for geometric source step";
  } else if (!(dt_s > 0.0) || !std::isfinite(dt_s)) {
    result.failure_reason = "time step must be finite and positive for geometric source step";
  } else {
    const auto source_terms = ComputeGeometricSourceTerms(source_snapshot, geometry);
    if (!source_terms.success) {
      result.failure_reason = source_terms.failure_reason;
    } else if (!source_terms.is_complete(source_snapshot.states.size())) {
      result.failure_reason = "geometric source terms are incomplete";
    } else {
      const GridShape shape{
          source_snapshot.radial_cells,
          source_snapshot.theta_cells,
          source_snapshot.phi_cells};

      for (std::size_t radial = 0; radial < shape.radial_cells; ++radial) {
        if (!result.failure_reason.empty()) {
          break;
        }
        for (std::size_t theta = 0; theta < shape.theta_cells; ++theta) {
          for (std::size_t phi = 0; phi < shape.phi_cells; ++phi) {
            const auto source =
                source_terms.source_terms[LinearIndex(shape, radial, theta, phi)];
            HydroConservativeState updated = LoadCellState(hydro_view, radial, theta, phi);
            updated.mom_r += dt_s * source.mom_r;
            updated.mom_theta += dt_s * source.mom_theta;
            updated.mom_phi += dt_s * source.mom_phi;
            if (!updated.is_finite() || !RecoverPrimitiveState(updated).is_physical()) {
              result.failure_reason = "geometric source step produced a non-physical hydro state";
              break;
            }

            (*hydro_view.mom_r)(radial, theta, phi) = updated.mom_r;
            (*hydro_view.mom_theta)(radial, theta, phi) = updated.mom_theta;
            (*hydro_view.mom_phi)(radial, theta, phi) = updated.mom_phi;
            ++result.updated_cell_count;
          }
          if (!result.failure_reason.empty()) {
            break;
          }
        }
      }
    }
  }

  result.success = result.failure_reason.empty();
  if (result.success) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.source.geometric_step",
        "hydro geometric source step updated momentum without modifying rho, E_fluid_total, or chi_e");
  } else {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.source.geometric_step.failed",
        result.failure_reason.empty() ? "hydro geometric source step failed" : result.failure_reason);
  }

  std::ostringstream report;
  report << "geometric_source_step_success=" << (result.success ? "true" : "false")
         << "; updated_cell_count=" << result.updated_cell_count;
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

}  // namespace dec3d::hydro
