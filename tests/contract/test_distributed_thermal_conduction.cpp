#include "test_assert.hpp"
#include "transport/thermal/distributed_thermal_conduction.hpp"

#include <mpi.h>

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    int size = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
      dec3d::test::Fail("rank count", __FILE__, __LINE__, "requires mpiexec -n 2");
    }

    auto problem = dec3d::transport::BuildDistributedThermalConductionSmokeProblem(
        MPI_COMM_WORLD,
        4,
        4,
        4);
    auto result = dec3d::transport::ApplyDistributedVariableKappaThermalConduction(
        problem,
        dec3d::transport::DistributedThermalConductionOptions{
            1.0e-12,
            dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy,
            dec3d::transport::ElectronFluxLimiterModel::disabled,
            0.0});

    if (!result.success) {
      std::cerr << result.failure_diagnostics << '\n';
      std::cerr << result.electron_solve_report << '\n';
      std::cerr << result.ion_solve_report << '\n';
      std::cerr << result.coefficient_report << '\n';
    }
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.distributed_writeback);
    DEC3D_CHECK(!result.gathered_writeback_used);
    DEC3D_CHECK(result.local_owned_write_count > 0);
    DEC3D_CHECK(result.global_owned_write_count > result.local_owned_write_count);
    DEC3D_CHECK(result.report_line.find("global_stage_ok=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("global_publish_ok=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("diagnostic_id=p2.production.thermal_stage") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("stage_id=T") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("stage_order=H,T,E") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("backend_requested=hypre_parcsr_gmres_boomeramg") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("thermal_kappa_model_requested=spitzer_no_degeneracy") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("electron_matrix_report_present=true") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("ion_matrix_report_present=true") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("global_volume_integrated_thermal_energy_residual=") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("updated_fields=e_electron,e_fluid_total") !=
                std::string::npos);
    DEC3D_CHECK_EQ(result.updated_fields_label, std::string{"e_electron,e_fluid_total"});
    DEC3D_CHECK_EQ(
        result.updated_fields,
        dec3d::transport::kDistributedThermalWritesElectronEnergy |
            dec3d::transport::kDistributedThermalWritesFluidTotalEnergy);
    DEC3D_CHECK(dec3d::transport::ValidateDistributedThermalConductionDiagnostics(result));

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    problem = dec3d::transport::BuildDistributedThermalConductionSmokeProblem(
        MPI_COMM_WORLD,
        4,
        4,
        4);
    if (rank == 1) {
      problem.local_state.e_electron(0, 0, 0) = -1.0;
    }
    result = dec3d::transport::ApplyDistributedVariableKappaThermalConduction(
        problem,
        dec3d::transport::DistributedThermalConductionOptions{
            1.0e-12,
            dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy,
            dec3d::transport::ElectronFluxLimiterModel::disabled,
            0.0});
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK_EQ(result.updated_fields, dec3d::transport::kDistributedThermalWritesNoFields);
    DEC3D_CHECK_EQ(result.updated_fields_label, std::string{"none"});
    DEC3D_CHECK(result.failure_diagnostics.find("global_stage_ok=false") != std::string::npos);
    DEC3D_CHECK(result.failure_diagnostics.find("canonical_state_mutated=false") !=
                std::string::npos);

    problem = dec3d::transport::BuildDistributedThermalConductionSmokeProblem(
        MPI_COMM_WORLD,
        4,
        4,
        4);
    const auto limited_result = dec3d::transport::ApplyDistributedVariableKappaThermalConduction(
        problem,
        dec3d::transport::DistributedThermalConductionOptions{
            1.0e-12,
            dec3d::transport::ThermalConductionKappaModel::spitzer_no_degeneracy,
            dec3d::transport::ElectronFluxLimiterModel::minmax_old_time_face_effective_kappa,
            1.0e-8});
    if (!limited_result.success) {
      std::cerr << limited_result.failure_diagnostics << '\n';
      std::cerr << limited_result.flux_limiter_report << '\n';
    }
    DEC3D_CHECK(limited_result.success);
    DEC3D_CHECK(limited_result.flux_limiter_report.find("electron_flux_limiter_enabled=true") !=
                std::string::npos);
    DEC3D_CHECK(limited_result.flux_limiter_report.find("electron_flux_limiter_model=minmax_old_time_face_effective_kappa") !=
                std::string::npos);
    DEC3D_CHECK(limited_result.flux_limiter_report.find("limited_face_count_global=") !=
                std::string::npos);
    DEC3D_CHECK(limited_result.report_line.find("electron_flux_limiter_enabled=true") !=
                std::string::npos);

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
