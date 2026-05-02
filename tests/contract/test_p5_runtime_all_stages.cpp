#include "app/dec3d_app.hpp"
#include "io/input_deck.hpp"
#include "test_assert.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
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
  auto root = std::filesystem::temp_directory_path() / "dec3d_p5_runtime_all_stages";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  return root;
}

}  // namespace

int main() {
  const auto root = TempRoot();
  const auto deck = root / "case.in";
  const auto profile = root / "case.pro";
  const auto output = root / "output";

  WriteText(deck, R"ini(
[run] # run
case_name = p5_runtime_all_stages # name
phase = P5 # phase
stage_order = H,T,E,R,A # full stage chain
step_count = 3 # physical target stops before this upper step cap
target_time_s = 1.0e-20 # clamp final dt to hit this physical stop time
dt_mode = hydro_relaxed # hydro CFL dt
cfl = 0.4 # CFL
output_dir = OUTPUT_ROOT_REPLACED # output dir
[mesh] # mesh
geometry = spherical # geometry
radial_cells = 8 # radial
theta_cells = 2 # theta
phi_cells = 2 # phi
radial_min_cm = 0.0 # inner radius
radial_max_cm = 1.0e-2 # outer radius
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = false # keep geometry fixed in this contract
macro_zoning = false # keep H path small
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
electron_flux_limiter = minmax_old_time_face_effective_kappa # thermal limiter
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
write_profiles_every = 1 # legacy
write_diagnostics_every = 1 # legacy
write_restart_every = 0 # legacy restart disabled
write_dec3d_out_every_steps = 1 # every runtime step
field_checkpoint_every_steps = 1 # every runtime step
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

  const std::vector<std::string> argv{
      "dec3d.exe",
      "-input",
      deck.string(),
      "-profile",
      profile.string()};
  const auto result = dec3d::app::RunDec3DCommandLine(argv);
  assert(result.exit_code == 0);
  assert(result.report_line.find("diagnostic_id=p5.io.runtime_entry") != std::string::npos);
  assert(result.report_line.find("runtime_time_loop_executed=true") != std::string::npos);
  assert(result.report_line.find("runtime_stage_order_executed=H,T,E,R,A") != std::string::npos);
  assert(result.report_line.find("runtime_steps_executed=1") != std::string::npos);
  assert(result.report_line.find("target_time_s=") != std::string::npos);
  assert(result.report_line.find("target_time_reached=true") != std::string::npos);
  assert(result.report_line.find("h_stage_executed=true") != std::string::npos);
  assert(result.report_line.find("t_stage_executed=true") != std::string::npos);
  assert(result.report_line.find("e_stage_executed=true") != std::string::npos);
  assert(result.report_line.find("r_stage_executed=true") != std::string::npos);
  assert(result.report_line.find("a_stage_executed=true") != std::string::npos);
  const std::vector<std::string> hydro_timing_keys{
      "h_hydro_snapshot_wall_s=",
      "h_hydro_scratch_wall_s=",
      "h_hydro_radial_sweep_wall_s=",
      "h_hydro_macro_detect_wall_s=",
      "h_hydro_macro_restrict_wall_s=",
      "h_hydro_macro_update_wall_s=",
      "h_hydro_macro_radial_update_wall_s=",
      "h_hydro_macro_theta_update_wall_s=",
      "h_hydro_macro_phi_update_wall_s=",
      "h_hydro_macro_state_update_wall_s=",
      "h_hydro_macro_prolong_wall_s=",
      "h_hydro_theta_sweep_wall_s=",
      "h_hydro_phi_sweep_wall_s=",
      "h_hydro_commit_wall_s=",
      "h_hydro_source_wall_s=",
      "h_hydro_budget_wall_s=",
      "h_hydro_diagnostics_wall_s="};
  for (const auto& key : hydro_timing_keys) {
    DEC3D_CHECK(result.report_line.find(key) != std::string::npos);
  }

  assert(std::filesystem::exists(output / "dec3d.out"));
  assert(std::filesystem::exists(output / "fields_000000.snap"));
  assert(std::filesystem::exists(output / "fields_000001.snap"));
  assert(!std::filesystem::exists(output / "fields_000002.snap"));
  assert(!std::filesystem::exists(output / "restart_000000.snap"));
  assert(!std::filesystem::exists(output / "restart_000001.snap"));

  const auto log = ReadText(output / "dec3d.out");
  assert(log.find("initialized") != std::string::npos);
  assert(log.find("advanced") != std::string::npos);

  const auto noh_root = root / "noh";
  const auto noh_deck = noh_root / "case_noh.in";
  const auto noh_profile = noh_root / "case_noh.pro";
  const auto noh_output = noh_root / "output";
  WriteText(noh_deck, R"ini(
[run] # run
case_name = p5_noh_exact_inflow_h_only # name
phase = P5 # phase
stage_order = H # H-only runtime validation
step_count = 1 # one hydro step
dt_mode = hydro_relaxed # hydro CFL dt
cfl = 0.5 # Noh CFL
output_dir = NOH_OUTPUT_ROOT_REPLACED # output dir
[mesh] # mesh
geometry = spherical # geometry
radial_cells = 32 # radial
theta_cells = 2 # theta
phi_cells = 2 # phi
radial_min_cm = 0.0 # inner radius
radial_max_cm = 1.0 # outer radius
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = false # fixed mesh
macro_zoning = true # macro path
[initial_condition] # init
profile_interpolation = linear # interpolation
profile_radius_unit = um # unit
outside_profile_policy = hard_fail # coverage
[physics] # physics
enable_hydro = true # H
enable_thermal = false # no T
enable_equilibration = false # no E
enable_radiation = false # no R
enable_alpha = false # no A
[radiation] # radiation placeholder
group_mode = explicit_frequency_groups # groups
group_edges_eV = 1,10 # one placeholder group
opacity_provider = tops_dt_tabulated # provider placeholder
radiation_initialization = zero # no radiation energy
radiation_initial_blackbody_scale = 1.0 # scale
boundary_model = thesis_marshak_vacuum # placeholder
flux_limiter = harmonic_eq_5_209 # placeholder
[thermal] # thermal placeholder
kappa_model = lee_more_with_degeneracy # placeholder
electron_flux_limiter = minmax_old_time_face_effective_kappa # placeholder
[alpha] # alpha placeholder
composition_model = equimolar_dt_from_p2_recovery # placeholder
reactivity_model = bosch_hale_dt # placeholder
alpha_initialization = zero # no alpha
[boundaries] # boundaries
inner_radial = scalar_origin_remap_required # origin
outer_radial = noh_exact_inflow # exact spherical Noh inflow ghost
theta = scalar_pole_remap_required # poles
phi = periodic # periodic phi
[output] # output
write_profiles_every = 0 # legacy disabled
write_diagnostics_every = 0 # legacy disabled
write_restart_every = 0 # legacy disabled
write_dec3d_out_every_steps = 1 # runtime log
field_checkpoint_every_steps = 0 # no field snaps
field_checkpoint_interval_s = 0.0 # no time snaps
field_checkpoint_prefix = fields # field prefix
field_checkpoint_format = csv3d # field format
restart_checkpoint_every_steps = 0 # no restart
restart_checkpoint_interval_s = 0.0 # no restart time trigger
restart_checkpoint_prefix = restart # restart prefix
restart_checkpoint_format = dec3d_restart_text # restart format
)ini");

  auto noh_deck_text = ReadText(noh_deck);
  const auto noh_marker = std::string("NOH_OUTPUT_ROOT_REPLACED");
  noh_deck_text.replace(noh_deck_text.find(noh_marker),
                        noh_marker.size(),
                        noh_output.generic_string());
  WriteText(noh_deck, noh_deck_text);

  WriteText(noh_profile, R"pro(
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale
0.0 1.0 1.0e-9 1.0e-9 -1.0 0.0 0.0 0.0 1.0
1000000.0 1.0 1.0e-9 1.0e-9 -1.0 0.0 0.0 0.0 1.0
)pro");

  const auto noh_deck_result = dec3d::io::LoadInputDeck(noh_deck);
  DEC3D_CHECK(noh_deck_result.success);
  DEC3D_CHECK(noh_deck_result.config.boundaries.outer_radial == "noh_exact_inflow");
  DEC3D_CHECK(noh_deck_result.report_line.find("outer_radial=noh_exact_inflow") !=
              std::string::npos);
  DEC3D_CHECK(noh_deck_result.report_line.find("noh_exact_inflow_deck_enabled=true") !=
              std::string::npos);
  DEC3D_CHECK(noh_deck_result.report_line.find(
                  "noh_exact_inflow_chi_e_time_level=current_time_density_scaled") !=
              std::string::npos);
  DEC3D_CHECK(noh_deck_result.report_line.find(
                  "noh_exact_inflow_pressure_time_level=current_time_adiabatic") !=
              std::string::npos);

  const std::vector<std::string> noh_argv{
      "dec3d.exe",
      "-input",
      noh_deck.string(),
      "-profile",
      noh_profile.string()};
  const auto noh_result = dec3d::app::RunDec3DCommandLine(noh_argv);
  DEC3D_CHECK(noh_result.exit_code == 0);
  for (const auto& key : hydro_timing_keys) {
    DEC3D_CHECK(noh_result.report_line.find(key) != std::string::npos);
  }

  std::filesystem::remove_all(root);
  return 0;
}
