#include "core/array/array3d.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "hydro/riemann/hllc_solver.hpp"
#include "mesh/macro_zoning/macro_zoning.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

struct SymmetryMetrics {
  double max_rho_angular_relative_spread{0.0};
  double max_pressure_angular_relative_spread{0.0};
  double max_abs_mom_theta{0.0};
  double max_abs_mom_phi{0.0};
  std::size_t mom_theta_radial{0};
  std::size_t mom_theta_theta{0};
  std::size_t mom_theta_phi{0};
};

[[nodiscard]] std::string FormatMetrics(const char* label, const SymmetryMetrics& metrics) {
  std::ostringstream report;
  report << std::setprecision(17)
         << label
         << " rho_spread=" << metrics.max_rho_angular_relative_spread
         << " pressure_spread=" << metrics.max_pressure_angular_relative_spread
         << " max_abs_mom_theta=" << metrics.max_abs_mom_theta
         << " mom_theta_cell="
         << metrics.mom_theta_radial << ","
         << metrics.mom_theta_theta << ","
         << metrics.mom_theta_phi
         << " max_abs_mom_phi=" << metrics.max_abs_mom_phi;
  return report.str();
}

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& code) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] bool HasDiagnosticMessageFragment(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const std::string& code,
    const std::string& fragment) {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code && entry.message.find(fragment) != std::string::npos) {
      return true;
    }
  }
  return false;
}

[[nodiscard]] SymmetryMetrics MeasureSymmetry(
    const dec3d::core::Array3D<double>& rho,
    const dec3d::core::Array3D<double>& mom_r,
    const dec3d::core::Array3D<double>& mom_theta,
    const dec3d::core::Array3D<double>& mom_phi,
    const dec3d::core::Array3D<double>& e_fluid_total) {
  SymmetryMetrics metrics;
  for (std::size_t radial = 0; radial < rho.extent_r(); ++radial) {
    double min_rho = std::numeric_limits<double>::infinity();
    double max_rho = -std::numeric_limits<double>::infinity();
    double rho_sum = 0.0;
    double min_pressure = std::numeric_limits<double>::infinity();
    double max_pressure = -std::numeric_limits<double>::infinity();
    double pressure_sum = 0.0;
    for (std::size_t theta = 0; theta < rho.extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < rho.extent_phi(); ++phi) {
        const double density = rho(radial, theta, phi);
        min_rho = std::min(min_rho, density);
        max_rho = std::max(max_rho, density);
        rho_sum += density;
        const auto primitive = dec3d::hydro::RecoverPrimitiveState(
            dec3d::hydro::HydroConservativeState{
                density,
                mom_r(radial, theta, phi),
                mom_theta(radial, theta, phi),
                mom_phi(radial, theta, phi),
                e_fluid_total(radial, theta, phi),
                0.0});
        const double pressure = primitive.pressure;
        min_pressure = std::min(min_pressure, pressure);
        max_pressure = std::max(max_pressure, pressure);
        pressure_sum += pressure;
        const double abs_mom_theta = std::abs(mom_theta(radial, theta, phi));
        if (abs_mom_theta > metrics.max_abs_mom_theta) {
          metrics.max_abs_mom_theta = abs_mom_theta;
          metrics.mom_theta_radial = radial;
          metrics.mom_theta_theta = theta;
          metrics.mom_theta_phi = phi;
        }
        metrics.max_abs_mom_phi =
            std::max(metrics.max_abs_mom_phi, std::abs(mom_phi(radial, theta, phi)));
      }
    }
    const double angular_count = static_cast<double>(rho.extent_theta() * rho.extent_phi());
    const double mean_rho = rho_sum / angular_count;
    if (std::abs(mean_rho) > 0.0) {
      metrics.max_rho_angular_relative_spread = std::max(
          metrics.max_rho_angular_relative_spread,
          (max_rho - min_rho) / std::abs(mean_rho));
    }
    const double mean_pressure = pressure_sum / angular_count;
    if (std::abs(mean_pressure) > 0.0) {
      metrics.max_pressure_angular_relative_spread = std::max(
          metrics.max_pressure_angular_relative_spread,
          (max_pressure - min_pressure) / std::abs(mean_pressure));
    }
  }
  return metrics;
}

[[nodiscard]] dec3d::state::CanonicalState SeedSphericalPulse(
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  using dec3d::state::CanonicalState;
  using dec3d::state::CanonicalStateLayout;
  using dec3d::state::ElectronEnergyDensityFromPressure;

  const std::size_t radial_cells = geometry.radial_faces.size() - 1u;
  const std::size_t theta_cells = geometry.theta_faces.size() - 1u;
  const std::size_t phi_cells = geometry.phi_faces.size() - 1u;
  auto state = CanonicalState::Create(CanonicalStateLayout{radial_cells, theta_cells, phi_cells, 0u});

  for (std::size_t radial = 0; radial < radial_cells; ++radial) {
    const double radius = 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
    const double normalized = radius / 0.12;
    const double pressure = 1.0e-4 + 50.0 * std::exp(-(normalized * normalized));
    const double e_total = pressure / (dec3d::state::HydroIdealGasGamma() - 1.0);
    const double e_electron = ElectronEnergyDensityFromPressure(0.5 * pressure);
    for (std::size_t theta = 0; theta < theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < phi_cells; ++phi) {
        state.rho(radial, theta, phi) = 1.0;
        state.mom_r(radial, theta, phi) = 0.0;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) = e_total;
        state.e_electron(radial, theta, phi) = e_electron;
      }
    }
  }

  return state;
}

[[nodiscard]] double EstimateStepDt(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::hydro::StaticGridHydroOptions& options) {
  dec3d::hydro::HydroOperator hydro;
  hydro.SetStaticGridOptions(options);
  const dec3d::core::StageContext context{
      0.0,
      1.0e-12,
      1u,
      dec3d::core::PhaseId::p1,
      "p1-v0.1",
      "mesh:p1.spherical-one-step",
      "ownership:single-rank",
      "diagnostics:p1.spherical-one-step"};
  DEC3D_CHECK(hydro.bind(context, geometry, state));
  const auto advice = hydro.estimate_dt();
  DEC3D_CHECK(advice.is_complete());
  DEC3D_CHECK(advice.hard_cap_dt > 0.0 && std::isfinite(advice.hard_cap_dt));
  return 0.8 * advice.hard_cap_dt;
}

[[nodiscard]] dec3d::hydro::StaticGridHydroOptions BasePpmOptions() {
  dec3d::hydro::StaticGridHydroOptions options;
  options.apply_geometric_source = true;
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = true;
  options.use_ppm_reconstruction = true;
  options.reconstruction_ghost_layers = 3u;
  options.use_macro_zoning = false;
  options.macro_zoning_coarse_factor = 0.5;
  return options;
}

[[nodiscard]] SymmetryMetrics RunHydroVariant(
    const char* label,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::hydro::StaticGridHydroOptions& options) {
  auto state = SeedSphericalPulse(geometry);
  const double dt_s = EstimateStepDt(state, geometry, options);
  auto view = dec3d::state::BuildHydroWorkView(state);
  DEC3D_CHECK(view.is_complete());
  const auto result = dec3d::hydro::AdvanceStaticGridHydro(view, geometry, dt_s, options);
  if (!result.is_complete()) {
    std::ostringstream message;
    message << label << " hydro update failed: " << result.failure_reason;
    throw std::runtime_error(message.str());
  }
  if (options.apply_geometric_source) {
    constexpr const char* kGeometricSourceDiagnostic = "p1.hydro.source.geometric_step";
    if (!HasDiagnosticCode(result.diagnostics, kGeometricSourceDiagnostic)) {
      std::ostringstream message;
      message << label << " did not emit required geometric source diagnostic";
      throw std::runtime_error(message.str());
    }
    const bool expects_macro_coarse_source =
        options.use_macro_zoning && (options.apply_theta_sweep || options.apply_phi_sweep);
    if (expects_macro_coarse_source &&
        !HasDiagnosticMessageFragment(
            result.diagnostics,
            kGeometricSourceDiagnostic,
            "macro-zoned coarse-consistent source")) {
      std::ostringstream message;
      message << label << " did not execute the macro-zoned coarse-consistent source path";
      throw std::runtime_error(message.str());
    }
  }
  const auto metrics = MeasureSymmetry(
      *view.rho,
      *view.mom_r,
      *view.mom_theta,
      *view.mom_phi,
      *view.e_fluid_total);
  std::cout << FormatMetrics(label, metrics) << '\n';
  return metrics;
}

[[nodiscard]] SymmetryMetrics RunMacroProjectionOnly(
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  using dec3d::mesh::DetectMacroZones;
  using dec3d::mesh::FineHydroPackageView;
  using dec3d::mesh::ProlongMacroZoneHydroPackage;
  using dec3d::mesh::RestrictFineHydroPackage;

  auto state = SeedSphericalPulse(geometry);
  auto view = dec3d::state::BuildHydroWorkView(state);
  DEC3D_CHECK(view.is_complete());
  const auto map = DetectMacroZones(
      geometry,
      dec3d::mesh::MacroZoningDetectionOptions{0.5});
  DEC3D_CHECK(map.is_complete());
  const FineHydroPackageView fine_view{
      view.rho,
      view.mom_r,
      view.mom_theta,
      view.mom_phi,
      view.e_fluid_total,
      &view.chi_e};
  DEC3D_CHECK(fine_view.is_complete());
  const auto coarse = RestrictFineHydroPackage(fine_view, geometry, map);
  DEC3D_CHECK(coarse.is_complete());
  const auto prolonged = ProlongMacroZoneHydroPackage(coarse);
  DEC3D_CHECK(prolonged.is_complete());
  const auto metrics = MeasureSymmetry(
      prolonged.rho,
      prolonged.mom_r,
      prolonged.mom_theta,
      prolonged.mom_phi,
      prolonged.e_fluid_total);
  std::cout << FormatMetrics("macro_projection_only", metrics) << '\n';
  return metrics;
}

void RequireSphericalPreservation(const char* label, const SymmetryMetrics& metrics) {
  constexpr double kRhoSpreadTolerance = 1.0e-12;
  constexpr double kPressureSpreadTolerance = 1.0e-12;
  constexpr double kTangentialMomentumTolerance = 1.0e-12;
  if (metrics.max_rho_angular_relative_spread > kRhoSpreadTolerance ||
      metrics.max_pressure_angular_relative_spread > kPressureSpreadTolerance ||
      metrics.max_abs_mom_theta > kTangentialMomentumTolerance ||
      metrics.max_abs_mom_phi > kTangentialMomentumTolerance) {
    throw std::runtime_error(FormatMetrics(label, metrics));
  }
}

}  // namespace

int main() {
  try {
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;

    constexpr std::size_t kRadialCells = 64u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
    DEC3D_CHECK(geometry.is_valid());

    const auto projection_only = RunMacroProjectionOnly(geometry);

    auto fine_ppm_options = BasePpmOptions();
    const auto fine_ppm = RunHydroVariant("fine_ppm", geometry, fine_ppm_options);

    auto macro_angular_options = BasePpmOptions();
    macro_angular_options.use_macro_zoning = true;
    macro_angular_options.apply_radial_sweep = false;
    const auto macro_angular = RunHydroVariant(
        "macro_angular_only_no_radial_macro",
        geometry,
        macro_angular_options);

    auto macro_angular_no_source_options = macro_angular_options;
    macro_angular_no_source_options.apply_geometric_source = false;
    const auto macro_angular_no_source = RunHydroVariant(
        "macro_angular_only_no_source",
        geometry,
        macro_angular_no_source_options);

    auto macro_theta_only_options = BasePpmOptions();
    macro_theta_only_options.use_macro_zoning = true;
    macro_theta_only_options.apply_radial_sweep = false;
    macro_theta_only_options.apply_phi_sweep = false;
    const auto macro_theta_only = RunHydroVariant(
        "macro_theta_only",
        geometry,
        macro_theta_only_options);

    auto macro_theta_only_no_source_options = macro_theta_only_options;
    macro_theta_only_no_source_options.apply_geometric_source = false;
    const auto macro_theta_only_no_source = RunHydroVariant(
        "macro_theta_only_no_source",
        geometry,
        macro_theta_only_no_source_options);

    auto macro_phi_only_options = BasePpmOptions();
    macro_phi_only_options.use_macro_zoning = true;
    macro_phi_only_options.apply_radial_sweep = false;
    macro_phi_only_options.apply_theta_sweep = false;
    const auto macro_phi_only = RunHydroVariant(
        "macro_phi_only",
        geometry,
        macro_phi_only_options);

    auto macro_phi_only_no_source_options = macro_phi_only_options;
    macro_phi_only_no_source_options.apply_geometric_source = false;
    const auto macro_phi_only_no_source = RunHydroVariant(
        "macro_phi_only_no_source",
        geometry,
        macro_phi_only_no_source_options);

    auto macro_angular_first_order_options = macro_angular_options;
    macro_angular_first_order_options.use_ppm_reconstruction = false;
    macro_angular_first_order_options.reconstruction_ghost_layers = 1u;
    const auto macro_angular_first_order = RunHydroVariant(
        "macro_angular_first_order",
        geometry,
        macro_angular_first_order_options);

    auto macro_radial_only_options = BasePpmOptions();
    macro_radial_only_options.use_macro_zoning = true;
    macro_radial_only_options.apply_geometric_source = false;
    macro_radial_only_options.apply_theta_sweep = false;
    macro_radial_only_options.apply_phi_sweep = false;
    const auto macro_radial_only = RunHydroVariant(
        "macro_radial_only_no_angular",
        geometry,
        macro_radial_only_options);

    auto macro_radial_angular_options = BasePpmOptions();
    macro_radial_angular_options.use_macro_zoning = true;
    const auto macro_radial_angular = RunHydroVariant(
        "macro_radial_plus_angular",
        geometry,
        macro_radial_angular_options);

    RequireSphericalPreservation("macro_projection_only", projection_only);
    RequireSphericalPreservation("fine_ppm", fine_ppm);
    RequireSphericalPreservation("macro_angular_only_no_radial_macro", macro_angular);
    (void)macro_angular_no_source;
    RequireSphericalPreservation("macro_theta_only", macro_theta_only);
    (void)macro_theta_only_no_source;
    (void)macro_phi_only;
    (void)macro_phi_only_no_source;
    RequireSphericalPreservation("macro_angular_first_order", macro_angular_first_order);
    RequireSphericalPreservation("macro_radial_only_no_angular", macro_radial_only);
    RequireSphericalPreservation("macro_radial_plus_angular", macro_radial_angular);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
