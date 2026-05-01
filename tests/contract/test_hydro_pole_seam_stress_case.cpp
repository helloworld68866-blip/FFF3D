#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

namespace {

[[nodiscard]] bool HasDiagnosticCode(
    const dec3d::core::DiagnosticsPayload& diagnostics,
    const char* code) noexcept {
  for (const auto& entry : diagnostics.entries) {
    if (entry.code == code) {
      return true;
    }
  }

  return false;
}

[[nodiscard]] std::string ReadTextFile(const std::filesystem::path& path) {
  std::ifstream input(path);
  std::string content;
  std::string line;
  while (std::getline(input, line)) {
    content += line;
    content.push_back('\n');
  }

  return content;
}

void SeedPoleSeamStressCase(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  constexpr double kBaseRho = 1.0;
  constexpr double kBasePressure = 1.0;
  constexpr double kBaseElectronPressure = 0.4;
  constexpr double kAmplitude = 3.0e-2;
  constexpr double kSigma = 0.28;
  constexpr double kPolarOffset = 0.25;

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      const double theta_center =
          0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
      const double north = std::exp(
          -std::pow(theta_center - kPolarOffset, 2.0) / (2.0 * kSigma * kSigma));
      const double south = std::exp(
          -std::pow(theta_center - (std::numbers::pi - kPolarOffset), 2.0) /
          (2.0 * kSigma * kSigma));

      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double phi_center =
            0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
        const double phi_mode = std::cos(phi_center) + 0.5 * std::sin(2.0 * phi_center);
        const double delta = kAmplitude * (north - south) * phi_mode;
        const double pressure = kBasePressure * (1.0 + 0.15 * delta);
        const double electron_pressure = kBaseElectronPressure * (1.0 + 0.10 * delta);
        const double v_theta = 2.0e-2 * delta;
        const double v_phi = -2.5e-2 * delta;
        const double kinetic = 0.5 * kBaseRho * (v_theta * v_theta + v_phi * v_phi);

        state.rho(radial, theta, phi) = kBaseRho;
        state.mom_r(radial, theta, phi) = 0.0;
        state.mom_theta(radial, theta, phi) = kBaseRho * v_theta;
        state.mom_phi(radial, theta, phi) = kBaseRho * v_phi;
        state.e_fluid_total(radial, theta, phi) =
            pressure / (dec3d::state::HydroIdealGasGamma() - 1.0) + kinetic;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);
      }
    }
  }
}

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::EvaluatePoleAdjacentSeamStress;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::WritePoleAdjacentSeamStressOutputs;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;

    constexpr std::size_t kRadialCells = 3u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    const std::filesystem::path output_dir =
        "F:\\dec3d\\analysis\\output\\case4_pole_seam_stress";
    std::filesystem::remove_all(output_dir);

    auto before = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    SeedPoleSeamStressCase(before, geometry);
    SeedPoleSeamStressCase(state, geometry);

    HydroOperator hydro;
    const StageContext context{
        0.0,
        5.0e-4,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.case4.pole-stress",
        "ownership:rank0",
        "diagnostics:p1.hydro.case4"};
    DEC3D_CHECK(hydro.bind(context, geometry, state));
    const auto stage_result = hydro.advance();
    DEC3D_CHECK(stage_result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.theta_pole_contract.executed"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.phi_periodic_contract.executed"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.ghost.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.ghost.direction.phi"));
    dec3d::hydro::HydroBudgetResidualSummary budget_summary;
    DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
        stage_result.diagnostics,
        &budget_summary));
    DEC3D_CHECK(budget_summary.is_complete());
    DEC3D_CHECK(std::abs(budget_summary.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(budget_summary.e_fluid_total_residual) < 1.0e-10);

    const auto stress = EvaluatePoleAdjacentSeamStress(
        before,
        state,
        geometry,
        1u,
        1.0e-4,
        2.0e-1);
    DEC3D_CHECK(stress.is_complete());
    DEC3D_CHECK(WritePoleAdjacentSeamStressOutputs(
        before,
        state,
        geometry,
        1u,
        output_dir,
        stress,
        &budget_summary));
    DEC3D_CHECK(stress.success);

    DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "north_polar_slice_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "north_polar_slice_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "south_polar_slice_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "south_polar_slice_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_shell_map_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_shell_map_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "seam_polar_diagnostics.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

    const auto summary = ReadTextFile(output_dir / "summary.txt");
    DEC3D_CHECK(summary.find("case=case4_pole_seam_stress") != std::string::npos);
    DEC3D_CHECK(summary.find("seam_growth_bounded=true") != std::string::npos);

    const auto manifest = ReadTextFile(output_dir / "case_manifest.txt");
    DEC3D_CHECK(manifest.find("seam_polar_diagnostics=seam_polar_diagnostics.txt") != std::string::npos);
    DEC3D_CHECK(manifest.find("budget_residual=budget_residual.txt") != std::string::npos);

    const auto diagnostics = ReadTextFile(output_dir / "seam_polar_diagnostics.txt");
    DEC3D_CHECK(diagnostics.find("max_seam_velocity_jump_after") != std::string::npos);

    const auto budget = ReadTextFile(output_dir / "budget_residual.txt");
    DEC3D_CHECK(budget.find("# kind=budget_residual") != std::string::npos);
    DEC3D_CHECK(budget.find("e_fluid_total_residual=") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
