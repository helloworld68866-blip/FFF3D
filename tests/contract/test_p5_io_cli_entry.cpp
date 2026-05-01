#include "app/dec3d_app.hpp"

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
  auto root = std::filesystem::temp_directory_path() / "dec3d_p5_io_cli_contract";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  return root;
}

}  // namespace

int main() {
  const auto root = TempRoot();
  const auto deck = root / "case.in";
  const auto profile = root / "case.pro";

  WriteText(deck, R"ini(
# case.in: CLI initialization test
[run] # run
case_name = p5_io_cli # name
phase = P5 # phase
stage_order = H,T,E,R,A # order
step_count = 0 # no advance
dt_mode = hydro_relaxed # dt
cfl = 0.4 # CFL
output_dir = OUTPUT_ROOT_REPLACED # output
[mesh] # mesh
geometry = spherical # geometry
radial_cells = 2 # r
theta_cells = 2 # theta
phi_cells = 4 # phi
radial_min_cm = 0.0 # r min
radial_max_cm = 1.0e-2 # r max
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = true # moving
macro_zoning = true # macro
[initial_condition] # initial condition
profile_interpolation = linear # interpolation
profile_radius_unit = um # unit
outside_profile_policy = hard_fail # outside
[physics] # physics
enable_hydro = true # H
enable_thermal = true # T
enable_equilibration = true # E
enable_radiation = true # R
enable_alpha = true # A
[radiation] # radiation
group_mode = explicit_frequency_groups # groups
group_edges_eV = 1,10,100 # edges
opacity_provider = tops_dt_tabulated # opacity
radiation_initialization = local_blackbody # blackbody
radiation_initial_blackbody_scale = 1.0 # scale
boundary_model = thesis_marshak_vacuum # boundary
flux_limiter = harmonic_eq_5_209 # limiter
[thermal] # thermal
kappa_model = lee_more_with_degeneracy # kappa
electron_flux_limiter = minmax_old_time_face_effective_kappa # limiter
[alpha] # alpha
composition_model = equimolar_dt_from_p2_recovery # composition
reactivity_model = bosch_hale_dt # reactivity
alpha_initialization = zero # alpha
[boundaries] # boundaries
inner_radial = scalar_origin_remap_required # inner
outer_radial = neumann_zero_flux # outer
theta = scalar_pole_remap_required # theta
phi = periodic # phi
[output] # output
write_profiles_every = 1 # profiles
write_diagnostics_every = 1 # diagnostics
write_restart_every = 10 # restart
write_dec3d_out_every_steps = 1 # dec3d.out
field_checkpoint_every_steps = 1 # field checkpoint
field_checkpoint_interval_s = 0.0 # disable time field checkpoint
field_checkpoint_prefix = fields # field prefix
field_checkpoint_format = csv3d # field format
restart_checkpoint_every_steps = 1 # restart checkpoint
restart_checkpoint_interval_s = 0.0 # disable time restart checkpoint
restart_checkpoint_prefix = restart # restart prefix
restart_checkpoint_format = dec3d_restart_text # restart format
)ini");

  auto deck_text = ReadText(deck);
  const auto output_dir = (root / "output").generic_string();
  const auto marker = std::string("OUTPUT_ROOT_REPLACED");
  deck_text.replace(deck_text.find(marker), marker.size(), output_dir);
  WriteText(deck, deck_text);

  WriteText(profile, R"pro(
# case.pro: CLI profile
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s radiation_scale # header
0.0 10.0 5.0 5.0 0.0 1.0 # inner
100.0 20.0 1.0 1.0 0.0 1.0 # outer
)pro");

  const std::vector<std::string> argv{"dec3d.exe", "-input", deck.string(), "-profile", profile.string()};
  const auto result = dec3d::app::RunDec3DCommandLine(argv);
  assert(result.exit_code == 0);
  assert(result.report_line.find("diagnostic_id=p5.io.runtime_entry") != std::string::npos);
  assert(result.report_line.find("profile_file_source=command_line.profile") != std::string::npos);
  assert(result.report_line.find("canonical_state_initialized=true") != std::string::npos);
  assert(result.report_line.find("runtime_output_report_present=true") != std::string::npos);
  assert(std::filesystem::exists(root / "output" / "dec3d.out"));
  assert(std::filesystem::exists(root / "output" / "fields_000000.snap"));
  assert(std::filesystem::exists(root / "output" / "restart_000000.snap"));

  const std::vector<std::string> missing_profile{"dec3d.exe", "-input", deck.string()};
  const auto bad = dec3d::app::RunDec3DCommandLine(missing_profile);
  assert(bad.exit_code != 0);
  assert(bad.failure_diagnostics.find("missing -profile") != std::string::npos);

  std::filesystem::remove_all(root);
  return 0;
}
