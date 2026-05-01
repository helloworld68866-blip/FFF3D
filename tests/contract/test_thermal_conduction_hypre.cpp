#include "test_assert.hpp"
#include "thermal_conduction_test_support.hpp"
#include "transport/diffusion/hypre_diffusion_solver.hpp"
#include "transport/thermal/thermal_conduction.hpp"

#include <cstddef>
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    using dec3d::transport::ApplyConstantKappaThermalConduction;
    using dec3d::transport::ApplyThermalConduction;
    using dec3d::transport::GenericDiffusionHypreSolveOptions;
    using dec3d::transport::SolveGenericDiffusionHypre;
    using dec3d::transport::ThermalConductionBackend;
    using dec3d::transport::ThermalConductionKappaModel;
    using dec3d::transport::ThermalConductionOptions;
    using dec3d::transport::ThermalConductionSolveResult;
    using dec3d::test_support::BuildPatchGeometry;
    using dec3d::test_support::BuildThermalState;
    using dec3d::test_support::CheckNear;
    using dec3d::test_support::Layout;
    using dec3d::test_support::PatchBoundary;

    const auto layout = Layout(2, 1, 1);
    const auto geometry = BuildPatchGeometry(layout, 1.0, 3.0, 0.5, 1.2);
    auto serial_state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);
    auto hypre_state = BuildThermalState(layout, 5.0, 1.0, 4.0, 2.0);

    auto options = ThermalConductionOptions{};
    options.dt_s = 0.01;
    options.kappa_model = ThermalConductionKappaModel::constant_user_supplied;
    options.kappa_e_cm_inv_s = 0.2;
    options.kappa_i_cm_inv_s = 0.3;
    options.backend = ThermalConductionBackend::hypre_parcsr_gmres_boomeramg;
    options.boundary_policy = PatchBoundary();

    auto serial_options = options;
    serial_options.backend = ThermalConductionBackend::serial_dense_reference;
    const auto serial_result =
        ApplyConstantKappaThermalConduction(serial_state, geometry, serial_options);
    DEC3D_CHECK(serial_result.success);

    const auto hook = [](const dec3d::transport::GenericDiffusionAssemblyResult& assembly) {
      const auto hypre = SolveGenericDiffusionHypre(
          assembly,
          GenericDiffusionHypreSolveOptions{});
      ThermalConductionSolveResult out;
      out.success = hypre.success;
      out.scalar_new = hypre.scalar_new;
      out.report_line = hypre.report_line;
      out.failure_diagnostics = hypre.failure_diagnostics;
      out.failure_reason = hypre.failure_reason;
      out.backend_executed = hypre.backend;
      out.residual_l2 = hypre.residual_l2;
      out.residual_linf = hypre.residual_linf;
      out.max_abs_delta = hypre.max_abs_delta;
      return out;
    };

    const auto result = ApplyConstantKappaThermalConduction(hypre_state, geometry, options, hook);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.report_line.find("backend_requested=hypre_parcsr_gmres_boomeramg") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("backend_executed=hypre_parcsr_gmres_boomeramg") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("fallback_used=false") != std::string::npos);
    DEC3D_CHECK(result.electron_solve_report.find("backend=hypre_parcsr_gmres_boomeramg") != std::string::npos);
    DEC3D_CHECK(result.ion_solve_report.find("backend=hypre_parcsr_gmres_boomeramg") != std::string::npos);

    for (std::size_t r = 0; r < layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < layout.theta_cells; ++t) {
        for (std::size_t p = 0; p < layout.phi_cells; ++p) {
          CheckNear(
              hypre_state.e_electron(r, t, p),
              serial_state.e_electron(r, t, p),
              1.0e-12,
              "HYPRE/serial e_electron parity");
          CheckNear(
              hypre_state.e_fluid_total(r, t, p),
              serial_state.e_fluid_total(r, t, p),
              1.0e-12,
              "HYPRE/serial e_fluid_total parity");
        }
      }
    }

    {
      auto serial_variable_state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);
      auto hypre_variable_state = BuildThermalState(layout, 0.01, 0.005, 0.008, 0.006);

      auto variable_options = options;
      variable_options.dt_s = 1.0e-16;
      variable_options.kappa_model = ThermalConductionKappaModel::spitzer_no_degeneracy;
      variable_options.kappa_e_cm_inv_s = 0.0;
      variable_options.kappa_i_cm_inv_s = 0.0;

      auto serial_variable_options = variable_options;
      serial_variable_options.backend = ThermalConductionBackend::serial_dense_reference;
      serial_variable_options.serial_reference_options.residual_tolerance = 1.0e-2;
      const auto serial_variable_result =
          ApplyThermalConduction(serial_variable_state, geometry, serial_variable_options);
      DEC3D_CHECK(serial_variable_result.success);
      DEC3D_CHECK(serial_variable_result.report_line.find("diagnostic_id=p2.thermal_conduction.variable_kappa") != std::string::npos);

      const auto hypre_variable_result =
          ApplyThermalConduction(hypre_variable_state, geometry, variable_options, hook);
      DEC3D_CHECK(hypre_variable_result.success);
      DEC3D_CHECK(hypre_variable_result.report_line.find("backend_requested=hypre_parcsr_gmres_boomeramg") != std::string::npos);
      DEC3D_CHECK(hypre_variable_result.report_line.find("backend_executed=hypre_parcsr_gmres_boomeramg") != std::string::npos);
      DEC3D_CHECK(hypre_variable_result.report_line.find("fallback_used=false") != std::string::npos);

      for (std::size_t r = 0; r < layout.radial_cells; ++r) {
        for (std::size_t t = 0; t < layout.theta_cells; ++t) {
          for (std::size_t p = 0; p < layout.phi_cells; ++p) {
            CheckNear(
                hypre_variable_state.e_electron(r, t, p),
                serial_variable_state.e_electron(r, t, p),
                1.0e-12,
                "variable HYPRE/serial e_electron parity");
            CheckNear(
                hypre_variable_state.e_fluid_total(r, t, p),
                serial_variable_state.e_fluid_total(r, t, p),
                1.0e-12,
                "variable HYPRE/serial e_fluid_total parity");
          }
        }
      }
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
