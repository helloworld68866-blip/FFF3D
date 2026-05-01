#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "test_assert.hpp"
#include "hydro_benchmark_case_support.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main() {
  try {
    using dec3d::hydro::EvaluateSedovSphericalBlast;
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
    std::filesystem::create_directories(output_dir);

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
    DEC3D_CHECK(check.shell_symmetry_preserved);
    DEC3D_CHECK(check.tangential_momentum_quiet);
    DEC3D_CHECK(check.max_rho_angular_relative_spread <= 1.0e-6);
    DEC3D_CHECK(check.tangential_to_radial_momentum_ratio <= 1.0e-8);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
