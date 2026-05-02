#include "app/dec3d_app.hpp"
#include "test_assert.hpp"

#include <mpi.h>

#include <filesystem>
#include <fstream>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

void WriteText(const std::filesystem::path& path, const std::string& text) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  out << text;
}

std::string ReadText(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  std::stringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

std::filesystem::path TempRoot() {
  return std::filesystem::temp_directory_path() / "dec3d_p5_runtime_distributed_all_stages";
}

void RequireTwoRanks() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    dec3d::test::Fail("rank count", __FILE__, __LINE__, "requires mpiexec -n 2");
  }
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    RequireTwoRanks();
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const auto root = TempRoot();
    const auto deck = root / "case.in";
    const auto profile = root / "case.pro";
    const auto output = root / "output";

    if (rank == 0) {
      std::filesystem::remove_all(root);
      std::filesystem::create_directories(root);
      WriteText(deck, R"ini(
[run] # run
case_name = p5_runtime_distributed_all_stages # name
phase = P5 # phase
stage_order = H,T,E,R,A # full stage chain
step_count = 1 # one real runtime step
dt_mode = hydro_relaxed # hydro CFL dt
cfl = 0.4 # CFL
output_dir = OUTPUT_ROOT_REPLACED # output dir
[mesh] # mesh
geometry = spherical # geometry
radial_cells = 8 # radial
theta_cells = 4 # theta
phi_cells = 4 # phi
radial_min_cm = 0.0 # inner radius
radial_max_cm = 1.0e-2 # outer radius
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = true # exercise MPI ALE hydro
macro_zoning = true # exercise PPM/macro hydro path
[initial_condition] # init
profile_interpolation = linear # interpolation
profile_radius_unit = um # unit
outside_profile_policy = hard_fail # coverage
[physics] # physics
enable_hydro = true # H
enable_thermal = true # T
enable_equilibration = true # E
enable_radiation = true # R
enable_alpha = true # A
[radiation] # radiation
group_mode = explicit_frequency_groups # groups
group_edges_eV = 1,10,100 # two groups within TOPS photon grid
opacity_provider = tops_dt_tabulated # provider
radiation_initialization = local_blackbody # initial U_g
radiation_initial_blackbody_scale = 1.0 # scale
boundary_model = thesis_marshak_vacuum # R boundary
flux_limiter = harmonic_eq_5_209 # limiter token
[thermal] # thermal
kappa_model = lee_more_with_degeneracy # thermal kappa
electron_flux_limiter = disabled # thermal limiter disabled for backend wiring contract
[alpha] # alpha
composition_model = equimolar_dt_from_p2_recovery # composition
reactivity_model = bosch_hale_dt # reactivity
alpha_initialization = zero # alpha starts from birth source
[boundaries] # boundaries
inner_radial = scalar_origin_remap_required # origin
outer_radial = neumann_zero_flux # outer for non-radiation diffusion
theta = scalar_pole_remap_required # poles
phi = periodic # phi
[output] # output
write_profiles_every = 0 # legacy profile output disabled
write_diagnostics_every = 0 # legacy diagnostics disabled
write_restart_every = 0 # legacy restart disabled
write_dec3d_out_every_steps = 0 # dec3d.out disabled in MPI contract
field_checkpoint_every_steps = 0 # field checkpoint disabled
field_checkpoint_interval_s = 0.0 # time trigger disabled
field_checkpoint_prefix = fields # field prefix
field_checkpoint_format = csv3d # field format
restart_checkpoint_every_steps = 0 # restart disabled
restart_checkpoint_interval_s = 0.0 # restart time disabled
restart_checkpoint_prefix = restart # restart prefix
restart_checkpoint_format = dec3d_restart_text # restart format
)ini");

      auto deck_text = ReadText(deck);
      const auto marker = std::string("OUTPUT_ROOT_REPLACED");
      deck_text.replace(deck_text.find(marker), marker.size(), output.generic_string());
      WriteText(deck, deck_text);

      WriteText(profile, R"pro(
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale
0.0 10.0 1.0 1.0 0.0 0.0 0.0 0.0 1.0
100.0 10.0 1.0 1.0 0.0 0.0 0.0 0.0 1.0
)pro");
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::vector<std::string> args{
        "dec3d.exe",
        "-input",
        deck.string(),
        "-profile",
        profile.string()};
    const auto result = dec3d::app::RunDec3DCommandLine(args);
    if (result.exit_code != 0) {
      dec3d::test::Fail("p5 distributed runtime success",
                        __FILE__,
                        __LINE__,
                        result.failure_diagnostics);
    }

    DEC3D_CHECK(result.report_line.find("distributed_mpi_runtime=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("runtime_stage_backend=mpi_hypre_callable_runtime") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("hydro_runtime_mode=distributed_mpi_hydro") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("h_stage_backend=mpi_macro_ale_hllc") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("hydro_halo_exchange_used=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("hydro_scalar_bundle_halo=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("t_stage_backend=hypre_parcsr_gmres_boomeramg") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("r_stage_backend=hypre_parcsr_gmres_boomeramg") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("a_stage_backend=hypre_parcsr_gmres_boomeramg") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("rank0_gather_solve_used=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("local_slab_allgather_after_distributed_stages=false") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("runtime_scalar_output_uses_mpi_reduce=false") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("root_gather_for_checkpoint_outputs=false") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("gather_metadata_precomputed=true") !=
                std::string::npos);

    const auto axisym_root = root / "axisymmetric_2d";
    const auto axisym_deck = axisym_root / "axisymmetric_2d_noale_allstages_smoke.in";
    const auto axisym_profile = axisym_root / "axisymmetric_2d_noale_smoke.pro";
    const auto axisym_output = axisym_root / "output";
    if (rank == 0) {
      std::filesystem::create_directories(axisym_root);
      std::filesystem::copy_file("F:/dec3d/cases/axisymmetric_2d_noale_allstages_smoke.in",
                                 axisym_deck,
                                 std::filesystem::copy_options::overwrite_existing);
      std::filesystem::copy_file("F:/dec3d/cases/axisymmetric_2d_noale_smoke.pro",
                                 axisym_profile,
                                 std::filesystem::copy_options::overwrite_existing);
      auto axisym_deck_text = ReadText(axisym_deck);
      const auto old_output =
          std::string("F:/dec3d/analysis/output/axisymmetric_2d_noale_allstages_smoke");
      axisym_deck_text.replace(axisym_deck_text.find(old_output),
                               old_output.size(),
                               axisym_output.generic_string());
      WriteText(axisym_deck, axisym_deck_text);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::vector<std::string> axisym_args{
        "dec3d.exe",
        "-input",
        axisym_deck.string(),
        "-profile",
        axisym_profile.string()};
    const auto axisym_result = dec3d::app::RunDec3DCommandLine(axisym_args);
    if (axisym_result.exit_code != 0) {
      dec3d::test::Fail("axisymmetric distributed runtime success",
                        __FILE__,
                        __LINE__,
                        axisym_result.failure_diagnostics);
    }
    DEC3D_CHECK(axisym_result.report_line.find("runtime_steps_executed=2") !=
                std::string::npos);
    DEC3D_CHECK(axisym_result.report_line.find("mesh_dimensionality=axisymmetric_2d") !=
                std::string::npos);
    DEC3D_CHECK(axisym_result.report_line.find("phi_sweep_executed=false") !=
                std::string::npos);
    DEC3D_CHECK(axisym_result.report_line.find("axisymmetric_invariant_ok=true") !=
                std::string::npos);
    DEC3D_CHECK(axisym_result.report_line.find("global_phi_coupling_count=0") !=
                std::string::npos);
    DEC3D_CHECK(axisym_result.report_line.find("global_duplicate_column_row_count=0") !=
                std::string::npos);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
      std::filesystem::remove_all(root);
    }
    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
