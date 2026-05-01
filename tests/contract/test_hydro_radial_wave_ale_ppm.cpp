#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ale/radial_ale.hpp"
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

void SeedRadialWave(
    dec3d::state::CanonicalState& state,
    dec3d::state::CanonicalState& before) {
  const dec3d::hydro::HydroPrimitiveState inner_state{
      1.0,
      -0.15,
      0.0,
      0.0,
      1.0,
      std::pow(0.4, 3.0 / 5.0)};
  const dec3d::hydro::HydroPrimitiveState outer_state{
      0.125,
      -0.05,
      0.0,
      0.0,
      0.1,
      std::pow(0.04, 3.0 / 5.0)};

  const auto inner_conservative = dec3d::hydro::MakeConservativeState(inner_state);
  const auto outer_conservative = dec3d::hydro::MakeConservativeState(outer_state);

  const std::size_t radial_midpoint = state.layout.radial_cells / 2u;
  for (std::size_t radial = 0; radial < state.layout.radial_cells; ++radial) {
    const bool inside = radial < radial_midpoint;
    const auto& cell = inside ? inner_conservative : outer_conservative;
    const double electron_pressure = inside ? 0.4 : 0.04;
    for (std::size_t theta = 0; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0; phi < state.layout.phi_cells; ++phi) {
        state.rho(radial, theta, phi) = cell.rho;
        state.mom_r(radial, theta, phi) = cell.mom_r;
        state.mom_theta(radial, theta, phi) = cell.mom_theta;
        state.mom_phi(radial, theta, phi) = cell.mom_phi;
        state.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
        state.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);

        before.rho(radial, theta, phi) = cell.rho;
        before.mom_r(radial, theta, phi) = cell.mom_r;
        before.mom_theta(radial, theta, phi) = cell.mom_theta;
        before.mom_phi(radial, theta, phi) = cell.mom_phi;
        before.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
        before.e_electron(radial, theta, phi) =
            dec3d::state::ElectronEnergyDensityFromPressure(electron_pressure);
      }
    }
  }
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::AdvanceStaticGridHydro;
    using dec3d::hydro::EvaluateRadialAleOnOff;
    using dec3d::hydro::StaticGridHydroOptions;
    using dec3d::hydro::WriteRadialAleOnOffOutputs;
    using dec3d::mesh::BuildRadialAleMeshUpdateProposal;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::BuildHydroAuthorizedWriteMask;
    using dec3d::state::BuildHydroWorkView;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::CommitHydroWriteback;

    constexpr std::size_t kRadialCells = 16u;
    constexpr std::size_t kThetaCells = 3u;
    constexpr std::size_t kPhiCells = 4u;
    constexpr std::size_t kGhostLayers = 3u;
    constexpr double kDt = 1.0e-3;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    auto before = CanonicalState::Create(
        CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto ale_off_state = CanonicalState::Create(before.layout);
    auto ale_on_state = CanonicalState::Create(before.layout);
    SeedRadialWave(ale_off_state, before);
    SeedRadialWave(ale_on_state, before);

    StaticGridHydroOptions ale_off_options;
    ale_off_options.apply_theta_sweep = false;
    ale_off_options.apply_phi_sweep = false;
    ale_off_options.use_ppm_reconstruction = true;
    ale_off_options.reconstruction_ghost_layers = kGhostLayers;

    StaticGridHydroOptions ale_on_options = ale_off_options;
    ale_on_options.apply_radial_ale_flux_correction = true;

    auto ale_off_view = BuildHydroWorkView(ale_off_state);
    auto ale_on_view = BuildHydroWorkView(ale_on_state);
    DEC3D_CHECK(ale_off_view.is_complete());
    DEC3D_CHECK(ale_on_view.is_complete());

    const auto proposal = BuildRadialAleMeshUpdateProposal(
        ale_on_state.rho,
        ale_on_state.mom_r,
        geometry,
        kDt);
    DEC3D_CHECK(proposal.is_complete(geometry.radial_faces.size()));

    const auto ale_off_result = AdvanceStaticGridHydro(
        ale_off_view,
        geometry,
        kDt,
        ale_off_options);
    const auto ale_on_result = AdvanceStaticGridHydro(
        ale_on_view,
        geometry,
        kDt,
        dec3d::hydro::RadialGhostOverride{},
        ale_on_options,
        &proposal);

    DEC3D_CHECK(ale_off_result.is_complete());
    DEC3D_CHECK(ale_on_result.is_complete());
    DEC3D_CHECK(HasDiagnosticCode(
        ale_off_result.diagnostics,
        "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_off_result.diagnostics,
        "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_off_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_on_result.diagnostics,
        "p1.hydro.reconstruction.ppm.executed"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_on_result.diagnostics,
        "p1.hydro.reconstruction.ppm.ng3"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_on_result.diagnostics,
        "p1.hydro.reconstruction.ppm.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(
        ale_on_result.diagnostics,
        "p1.hydro.ale.radial_flux_correction.executed"));

    DEC3D_CHECK(CommitHydroWriteback(
        ale_off_state,
        ale_off_view,
        BuildHydroAuthorizedWriteMask()).success);
    DEC3D_CHECK(CommitHydroWriteback(
        ale_on_state,
        ale_on_view,
        BuildHydroAuthorizedWriteMask()).success);

    const auto summary = EvaluateRadialAleOnOff(
        before,
        ale_off_state,
        ale_on_state,
        geometry,
        proposal,
        1.0e-12,
        true,
        kGhostLayers);
    DEC3D_CHECK(summary.is_complete());
    DEC3D_CHECK(summary.success);

    const std::filesystem::path output_dir =
        "F:\\dec3d\\analysis\\output\\case2_radial_wave_ale_ppm";
    std::filesystem::remove_all(output_dir);
    DEC3D_CHECK(WriteRadialAleOnOffOutputs(
        before,
        ale_off_state,
        ale_on_state,
        geometry,
        proposal,
        output_dir,
        summary,
        &ale_off_result.budget,
        &ale_on_result.budget,
        "case2_radial_wave_ale_ppm"));

    DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "moving_mesh_diagnostics.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual_comparison.txt"));
    const auto summary_text = ReadTextFile(output_dir / "summary.txt");
    DEC3D_CHECK(summary_text.find("case=case2_radial_wave_ale_ppm") != std::string::npos);
    DEC3D_CHECK(summary_text.find("success=true") != std::string::npos);
    DEC3D_CHECK(summary_text.find("ale_changes_solution=true") != std::string::npos);
    DEC3D_CHECK(summary_text.find("ppm_executed=true") != std::string::npos);
    DEC3D_CHECK(summary_text.find("reconstruction_ghost_layers=3") != std::string::npos);
    DEC3D_CHECK(summary_text.find("reconstruction_mode=radial_only_ppm_ale") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
