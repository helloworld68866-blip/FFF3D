#include "hydro/driver/hydro_numerical_checks.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "hydro_benchmark_case_support.hpp"

#include <filesystem>
#include <iostream>

int main() {
  try {
    using dec3d::hydro::EvaluateSedovSphericalBlast;
    using dec3d::hydro::StaticGridHydroOptions;
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
    constexpr std::size_t kMaxSteps = 1000000u;

    const auto geometry = BuildSphericalGeometry(
        SphericalMeshDescriptor{kRadialCells, kThetaCells, kPhiCells, 0.0, 1.0});
    if (!geometry.is_valid()) {
      std::cerr << "invalid geometry\n";
      return 1;
    }

    StaticGridHydroOptions macro_options{};
    macro_options.apply_geometric_source = true;
    macro_options.apply_radial_sweep = true;
    macro_options.apply_theta_sweep = true;
    macro_options.apply_phi_sweep = true;
    macro_options.use_ppm_reconstruction = false;
    macro_options.use_macro_zoning = true;
    macro_options.macro_zoning_coarse_factor = 0.5;

    const std::filesystem::path output_dir =
        "F:\\dec3d\\analysis\\output\\case_sedov_spherical_macro_trial";
    const std::filesystem::path progress_path = output_dir / "dec3d.out";
    std::filesystem::remove_all(output_dir);
    std::filesystem::create_directories(output_dir);

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
        kMaxSteps,
        &progress_path,
        &macro_options);

    if (!run.success) {
      std::cerr << run.failure_reason << '\n';
      return 1;
    }

    const auto summary = EvaluateSedovSphericalBlast(
        run.before,
        run.after,
        geometry,
        kBlastEnergy,
        kAmbientDensity,
        run.final_time_s,
        0.25,
        1.0e-6,
        1.0e-8);

    if (!summary.is_complete()) {
      std::cerr << "sedov macro summary incomplete\n";
      return 1;
    }

    if (!WriteSedovSphericalBlastOutputs(
            run.before,
            run.after,
            geometry,
            output_dir,
            summary,
            run.shock_radius_history,
            run.profile_history,
            &run.final_budget,
            "case_sedov_spherical_macro_trial")) {
      std::cerr << "failed to write sedov macro outputs\n";
      return 1;
    }

    std::cout << summary.report_line << '\n';
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
