#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_numerical_checks.hpp"
#include "hydro/driver/hydro_operator.hpp"
#include "hydro/riemann/hllc_solver.hpp"
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

}  // namespace

int main() {
  try {
    using dec3d::core::PhaseId;
    using dec3d::core::StageContext;
    using dec3d::hydro::EvaluateSphericalRadialWaveMpiParity;
    using dec3d::hydro::HydroOperator;
    using dec3d::hydro::HydroPrimitiveState;
    using dec3d::hydro::MakeConservativeState;
    using dec3d::hydro::WriteSphericalRadialWaveMpiParityOutputs;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;
    using dec3d::state::CanonicalState;
    using dec3d::state::CanonicalStateLayout;
    using dec3d::state::ElectronEnergyDensityFromPressure;

    constexpr std::size_t kRadialCells = 16u;
    constexpr std::size_t kThetaCells = 3u;
    constexpr std::size_t kPhiCells = 4u;
    constexpr std::size_t kRankCount = 3u;
    const std::filesystem::path output_dir =
        "F:\\dec3d\\analysis\\output\\case2_radial_wave_mpi_parity";
    std::filesystem::remove_all(output_dir);

    auto before = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    auto single_rank_after = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});

    const HydroPrimitiveState inner_state{
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        std::pow(0.4, 3.0 / 5.0)};
    const HydroPrimitiveState outer_state{
        0.125,
        0.0,
        0.0,
        0.0,
        0.1,
        std::pow(0.04, 3.0 / 5.0)};

    const auto inner_conservative = MakeConservativeState(inner_state);
    const auto outer_conservative = MakeConservativeState(outer_state);

    for (std::size_t radial = 0; radial < kRadialCells; ++radial) {
      const bool inside = radial < (kRadialCells / 2u);
      const auto& cell = inside ? inner_conservative : outer_conservative;
      const double electron_pressure = inside ? 0.4 : 0.04;
      for (std::size_t theta = 0; theta < kThetaCells; ++theta) {
        for (std::size_t phi = 0; phi < kPhiCells; ++phi) {
          before.rho(radial, theta, phi) = cell.rho;
          before.mom_r(radial, theta, phi) = cell.mom_r;
          before.mom_theta(radial, theta, phi) = cell.mom_theta;
          before.mom_phi(radial, theta, phi) = cell.mom_phi;
          before.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
          before.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(electron_pressure);

          single_rank_after.rho(radial, theta, phi) = cell.rho;
          single_rank_after.mom_r(radial, theta, phi) = cell.mom_r;
          single_rank_after.mom_theta(radial, theta, phi) = cell.mom_theta;
          single_rank_after.mom_phi(radial, theta, phi) = cell.mom_phi;
          single_rank_after.e_fluid_total(radial, theta, phi) = cell.e_fluid_total;
          single_rank_after.e_electron(radial, theta, phi) =
              ElectronEnergyDensityFromPressure(electron_pressure);
        }
      }
    }

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.5, 1.5});
    DEC3D_CHECK(geometry.is_valid());

    const StageContext context{
        0.0,
        1.0e-3,
        1,
        PhaseId::p1,
        "p1-v0.1",
        "mesh:p1.radial-wave-mpi-parity",
        "ownership:rank0",
        "diagnostics:p1.hydro.radial-wave-mpi-parity"};
    DEC3D_CHECK(context.is_complete());

    HydroOperator hydro;
    DEC3D_CHECK(hydro.bind(context, geometry, single_rank_after));
    const auto single_rank_stage_result = hydro.advance();
    DEC3D_CHECK(single_rank_stage_result.is_semantically_complete());
    DEC3D_CHECK(HasDiagnosticCode(single_rank_stage_result.diagnostics, "p1.hydro.direction.radial"));
    DEC3D_CHECK(HasDiagnosticCode(single_rank_stage_result.diagnostics, "p1.hydro.ghost.direction.radial"));

    auto multi_rank_after = CanonicalState::Create(CanonicalStateLayout{kRadialCells, kThetaCells, kPhiCells, 0u});
    const auto parity = EvaluateSphericalRadialWaveMpiParity(
        before,
        single_rank_after,
        geometry,
        context.dt_s,
        kRankCount,
        1.0e-10,
        1.0e-10,
        &multi_rank_after);
    DEC3D_CHECK(parity.is_complete());
    DEC3D_CHECK(parity.success);
    DEC3D_CHECK(parity.rank_decomposition_valid);
    DEC3D_CHECK(parity.parity_within_tolerance);
    DEC3D_CHECK(parity.seam_jumps_bounded);

    DEC3D_CHECK(WriteSphericalRadialWaveMpiParityOutputs(
        single_rank_after,
        multi_rank_after,
        geometry,
        output_dir,
        parity));

    DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "seam_diagnostics.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "single_rank_radial_profile_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "multi_rank_radial_profile_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_shell_map_t1_single.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_shell_map_t1_multi.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_shell_map_t1_single.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_shell_map_t1_multi.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));

    const auto summary = ReadTextFile(output_dir / "summary.txt");
    DEC3D_CHECK(summary.find("case=case2_radial_wave_mpi_parity") != std::string::npos);
    DEC3D_CHECK(summary.find("rank_count=3") != std::string::npos);
    DEC3D_CHECK(summary.find("success=true") != std::string::npos);

    const auto seam = ReadTextFile(output_dir / "seam_diagnostics.txt");
    DEC3D_CHECK(seam.find("# kind=seam_diagnostics") != std::string::npos);

    const auto manifest = ReadTextFile(output_dir / "case_manifest.txt");
    DEC3D_CHECK(manifest.find("single_rank_profile=single_rank_radial_profile_t1.txt") != std::string::npos);
    DEC3D_CHECK(manifest.find("multi_rank_profile=multi_rank_radial_profile_t1.txt") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
