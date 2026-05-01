#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/hydro_state/hydro_view.hpp"
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
    DEC3D_CHECK(run.final_budget.is_complete());

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
    DEC3D_CHECK(check.success);
    DEC3D_CHECK(check.shock_radius_within_tolerance);
    DEC3D_CHECK(std::abs(check.shock_radius_relative_error) <= 0.15);
    DEC3D_CHECK(check.numerical_shock_radius > geometry.radial_faces[2]);

    const auto synthetic_geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{6u, 1u, 1u, 0.0, 0.12});
    DEC3D_CHECK(synthetic_geometry.is_valid());
    auto synthetic_before = dec3d::state::CanonicalState::Create(
        dec3d::state::CanonicalStateLayout{6u, 1u, 1u, 0u});
    auto synthetic_after = dec3d::state::CanonicalState::Create(synthetic_before.layout);
    const double rho_profile[] = {0.2, 1.0, 2.5, 4.0, 1.0, 1.0};
    for (std::size_t radial = 0; radial < 6u; ++radial) {
      synthetic_before.rho(radial, 0u, 0u) = 1.0;
      synthetic_before.e_fluid_total(radial, 0u, 0u) = 1.0;
      synthetic_before.e_electron(radial, 0u, 0u) = 0.5;
      synthetic_after.rho(radial, 0u, 0u) = rho_profile[radial];
      synthetic_after.e_fluid_total(radial, 0u, 0u) = 1.0;
      synthetic_after.e_electron(radial, 0u, 0u) = 0.5;
    }
    const auto synthetic_check = EvaluateSedovSphericalBlast(
        synthetic_before,
        synthetic_after,
        synthetic_geometry,
        kBlastEnergy,
        kAmbientDensity,
        kFinalTime,
        1.0,
        1.0e-12,
        1.0e-12);
    DEC3D_CHECK(synthetic_check.is_complete());
    DEC3D_CHECK(std::abs(synthetic_check.numerical_shock_radius - 0.09) <= 1.0e-12);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
