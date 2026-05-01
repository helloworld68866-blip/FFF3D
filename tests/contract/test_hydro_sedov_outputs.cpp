#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"
#include "hydro_benchmark_case_support.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
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

}  // namespace

int main() {
  try {
    using dec3d::hydro::EvaluateSedovSphericalBlast;
    using dec3d::hydro::WriteSedovSphericalBlastOutputs;
    using dec3d::mesh::BuildSphericalGeometry;
    using dec3d::mesh::SphericalMeshDescriptor;

    constexpr std::size_t kRadialCells = 256u;
    constexpr std::size_t kThetaCells = 8u;
    constexpr std::size_t kPhiCells = 8u;
    constexpr double kAmbientDensity = 1.0;
    constexpr double kAmbientPressure = 1.0e-5;
    constexpr double kBlastEnergy = 1.0;
    constexpr double kBlastRadius = 2.0e-2;
    constexpr double kElectronEnergyFraction = 0.5;
    constexpr double kFinalTime = 2.0e-3;
    const std::filesystem::path output_dir = "F:\\dec3d\\analysis\\output\\case_sedov_spherical";
    const std::filesystem::path progress_path = output_dir / "dec3d.out";

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
    DEC3D_CHECK(geometry.is_valid());

    const auto run = dec3d::testsupport::RunSedovBaseline(
        geometry,
        kRadialCells,
        kThetaCells,
        kPhiCells,
        kAmbientDensity,
        kAmbientPressure,
        kBlastEnergy,
        kBlastRadius,
        kElectronEnergyFraction,
        kFinalTime,
        2048u,
        &progress_path);
    if (!run.success) {
      throw std::runtime_error(run.failure_reason);
    }

    const auto check = EvaluateSedovSphericalBlast(
        run.before,
        run.after,
        geometry,
        kBlastEnergy,
        kAmbientDensity,
        run.final_time_s,
        0.15,
        1.0e-6,
        1.0e-8);
    DEC3D_CHECK(check.is_complete());

    std::filesystem::remove_all(output_dir);
    DEC3D_CHECK(WriteSedovSphericalBlastOutputs(
        run.before,
        run.after,
        geometry,
        output_dir,
        check,
        run.shock_radius_history,
        run.profile_history,
        &run.final_budget,
        "case_sedov_spherical"));

    DEC3D_CHECK(std::filesystem::exists(output_dir / "summary.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "case_manifest.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "dec3d.out"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shock_radius_vs_time.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "radial_profile_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "radial_profile_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "radial_profile_step000000.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "radial_profile_step000020.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "rho_rt_map.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "hydro_only_te_rt_map.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "velocity_magnitude_rt_map.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "mom_r_rt_map.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "e_fluid_total_rt_map.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "budget_residual.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_rho_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_rho_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_te_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_te_t1.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_velocity_t0.txt"));
    DEC3D_CHECK(std::filesystem::exists(output_dir / "shell_map_velocity_t1.txt"));

    const auto summary_text = ReadTextFile(output_dir / "summary.txt");
    DEC3D_CHECK(summary_text.find("case=case_sedov_spherical") != std::string::npos);
    DEC3D_CHECK(summary_text.find("shock_radius_within_tolerance=true") != std::string::npos);

    const auto manifest_text = ReadTextFile(output_dir / "case_manifest.txt");
    DEC3D_CHECK(manifest_text.find("run_progress=dec3d.out") != std::string::npos);
    DEC3D_CHECK(manifest_text.find("shock_radius_history=shock_radius_vs_time.txt") != std::string::npos);
    DEC3D_CHECK(manifest_text.find("radial_profile_series=radial_profile_step*.txt") != std::string::npos);
    DEC3D_CHECK(manifest_text.find("rho_rt_map=rho_rt_map.txt") != std::string::npos);
    DEC3D_CHECK(manifest_text.find("budget_residual=budget_residual.txt") != std::string::npos);

    const auto progress_text = ReadTextFile(output_dir / "dec3d.out");
    DEC3D_CHECK(progress_text.find("# kind=run_progress") != std::string::npos);
    DEC3D_CHECK(progress_text.find("pending") != std::string::npos);

    const auto shock_text = ReadTextFile(output_dir / "shock_radius_vs_time.txt");
    DEC3D_CHECK(shock_text.find("# kind=shock_radius_history") != std::string::npos);
    DEC3D_CHECK(shock_text.find("reference_radius") != std::string::npos);

    const auto rt_map_text = ReadTextFile(output_dir / "rho_rt_map.txt");
    DEC3D_CHECK(rt_map_text.find("# kind=radial_time_map") != std::string::npos);
    DEC3D_CHECK(rt_map_text.find("# field=rho") != std::string::npos);
    DEC3D_CHECK(rt_map_text.find("# radial_centers=") != std::string::npos);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
