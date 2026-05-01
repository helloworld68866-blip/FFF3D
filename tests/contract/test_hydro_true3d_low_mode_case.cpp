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

void SeedTrue3DLowModeCase(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  constexpr double kBaseRho = 1.0;
  constexpr double kBasePressure = 1.0;
  constexpr double kBaseElectronPressure = 0.4;
  constexpr double kAmplitude = 2.0e-2;

  for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
      const double theta_center =
          0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
      for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
        const double phi_center =
            0.5 * (geometry.phi_faces[phi] + geometry.phi_faces[phi + 1u]);
        const double low_mode = std::sin(theta_center) * std::cos(phi_center);
        const double v_r = kAmplitude * low_mode;
        const double kinetic = 0.5 * kBaseRho * v_r * v_r;

        state.rho(radial, theta, phi) = kBaseRho;
        state.mom_r(radial, theta, phi) = kBaseRho * v_r;
        state.mom_theta(radial, theta, phi) = 0.0;
        state.mom_phi(radial, theta, phi) = 0.0;
        state.e_fluid_total(radial, theta, phi) =
            kBasePressure / (dec3d::state::HydroIdealGasGamma() - 1.0) + kinetic;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(kBaseElectronPressure);
      }
    }
  }
}

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::EvaluateTrue3DLowModeSanity;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::hydro::WriteTrue3DLowModeOutputs;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroAuthorizedWriteMask;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::CommitHydroWriteback;

    constexpr std::size_t kRadialCells = 4u;
    constexpr std::size_t kThetaCells = 6u;
    constexpr std::size_t kPhiCells = 8u;
    const std::filesystem::path output_dir =
        "F:\\dec3d\\analysis\\output\\case3_true3d_low_mode";
    std::filesystem::remove_all(output_dir);

    auto before = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto full_state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto radial_only_state = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    SeedTrue3DLowModeCase(before, geometry);
    SeedTrue3DLowModeCase(full_state, geometry);
    SeedTrue3DLowModeCase(radial_only_state, geometry);

    HydroOperator hydro;
    const StageContext context{
        0.0,
        5.0e-4,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.case3.low-mode",
        "ownership:rank0",
        "diagnostics:p1.hydro.case3"};
    DEC3D_CHECK(hydro.bind(context, geometry, full_state));
    const auto stage_result = hydro.advance();
    DEC3D_CHECK(stage_result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.direction.theta"));
    DEC3D_CHECK(HasDiagnosticCode(stage_result.diagnostics, "p1.hydro.direction.phi"));
    dec3d::hydro::HydroBudgetResidualSummary budget_summary;
    DEC3D_CHECK(dec3d::hydro::TryExtractHydroBudgetResidualSummary(
        stage_result.diagnostics,
        &budget_summary));
    DEC3D_CHECK(budget_summary.is_complete());
    DEC3D_CHECK(std::abs(budget_summary.mass_residual) < 1.0e-10);
    DEC3D_CHECK(std::abs(budget_summary.e_fluid_total_residual) < 1.0e-10);

    auto radial_only_view = BuildHydroWorkView(radial_only_state);
    DEC3D_CHECK(radial_only_view.is_complete());
    const auto radial_only_result = AdvanceStaticGridHydro(
        radial_only_view,
        geometry,
        context.dt_s,
        StaticGridHydroOptions{true, true, false, false});
    DEC3D_CHECK(radial_only_result.is_complete());
    const auto radial_only_writeback = CommitHydroWriteback(
        radial_only_state,
        radial_only_view,
        BuildHydroAuthorizedWriteMask());
    DEC3D_CHECK(radial_only_writeback.is_complete());
    DEC3D_CHECK(radial_only_writeback.success);

    const auto sanity = EvaluateTrue3DLowModeSanity(
        before,
        full_state,
        radial_only_state,
        geometry,
        2u,
        1.0e-8,
        1.0e-8);
    DEC3D_CHECK(sanity.is_complete());
    DEC3D_CHECK(WriteTrue3DLowModeOutputs(
        before,
        full_state,
        radial_only_state,
        geometry,
        2u,
        output_dir,
        sanity,
        &budget_summary));
    DEC3D_CHECK(sanity.success);

    DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_shell_map_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_shell_map_t1_full.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "hydro_only_te_shell_map_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "hydro_only_te_shell_map_t1_full.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_shell_map_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_shell_map_t1_full.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "mode_projection.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "radial_only_comparison.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

    const auto summary = ReadTextFile(output_dir / "summary.txt");
    DEC3D_CHECK(summary.find("case=case3_true3d_low_mode") != std::string::npos);
    DEC3D_CHECK(summary.find("full_directional_differs_from_radial_only=true") != std::string::npos);

    const auto manifest = ReadTextFile(output_dir / "case_manifest.txt");
    DEC3D_CHECK(manifest.find("mode_projection=mode_projection.txt") != std::string::npos);
    DEC3D_CHECK(manifest.find("budget_residual=budget_residual.txt") != std::string::npos);

    const auto mode_projection = ReadTextFile(output_dir / "mode_projection.txt");
    DEC3D_CHECK(mode_projection.find("# kind=mode_projection") != std::string::npos);

    const auto budget = ReadTextFile(output_dir / "budget_residual.txt");
    DEC3D_CHECK(budget.find("# kind=budget_residual") != std::string::npos);
    DEC3D_CHECK(budget.find("mass_residual=") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
