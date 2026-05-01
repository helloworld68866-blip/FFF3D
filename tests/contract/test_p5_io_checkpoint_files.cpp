#include "initialization/profile_initializer.hpp"
#include "io/input_deck.hpp"
#include "io/radial_profile.hpp"
#include "io/runtime_output.hpp"
#include "test_assert.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

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

}  // namespace

int main() {
  const auto root = std::filesystem::temp_directory_path() / "dec3d_p5_io_checkpoint_files";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);

  const auto deck_path = root / "case.in";
  const auto profile_path = root / "case.pro";

  WriteText(deck_path, R"deck(
[run]
case_name = output_contract
phase = P5
stage_order = H,T,E,R,A
step_count = 0
dt_mode = hydro_relaxed
cfl = 0.4
output_dir = OUTPUT_ROOT_REPLACED
[mesh]
geometry = spherical
radial_cells = 2
theta_cells = 2
phi_cells = 4
radial_min_cm = 0.0
radial_max_cm = 1.0e-2
theta_min = 0.0
theta_max = 3.141592653589793
phi_min = 0.0
phi_max = 6.283185307179586
moving_mesh = true
macro_zoning = true
[initial_condition]
profile_interpolation = linear
profile_radius_unit = um
outside_profile_policy = hard_fail
[physics]
enable_hydro = true
enable_thermal = true
enable_equilibration = true
enable_radiation = true
enable_alpha = true
[radiation]
group_mode = explicit_frequency_groups
group_edges_eV = 1,10,100
opacity_provider = tops_dt_tabulated
radiation_initialization = zero
radiation_initial_blackbody_scale = 1.0
boundary_model = thesis_marshak_vacuum
flux_limiter = harmonic_eq_5_209
[thermal]
kappa_model = lee_more_with_degeneracy
electron_flux_limiter = minmax_old_time_face_effective_kappa
[alpha]
composition_model = equimolar_dt_from_p2_recovery
reactivity_model = bosch_hale_dt
alpha_initialization = profile
[boundaries]
inner_radial = scalar_origin_remap_required
outer_radial = neumann_zero_flux
theta = scalar_pole_remap_required
phi = periodic
[output]
write_profiles_every = 1
write_diagnostics_every = 1
write_restart_every = 10
write_dec3d_out_every_steps = 1
field_checkpoint_every_steps = 1
field_checkpoint_interval_s = 0.0
field_checkpoint_prefix = fields
field_checkpoint_format = csv3d
restart_checkpoint_every_steps = 1
restart_checkpoint_interval_s = 0.0
restart_checkpoint_prefix = restart
restart_checkpoint_format = dec3d_restart_text
history_profile_interval_s = 1.0e-12
history_profile_file = output_contract.his
)deck");

  auto deck_text = ReadText(deck_path);
  const auto output_dir = (root / "output").generic_string();
  const auto marker = std::string("OUTPUT_ROOT_REPLACED");
  deck_text.replace(deck_text.find(marker), marker.size(), output_dir);
  WriteText(deck_path, deck_text);

  WriteText(profile_path, R"profile(
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale
0.0 10.0 5.0 4.0 -1.0e7 0.0 0.0 0.0 1.0
100.0 20.0 3.0 2.0 -5.0e6 0.0 0.0 1.0e9 0.5
)profile");

  const auto deck = dec3d::io::LoadInputDeck(deck_path);
  assert(deck.success);
  const auto profile = dec3d::io::LoadRadialProfile(profile_path, "um");
  assert(profile.success);
  const auto init =
      dec3d::initialization::InitializeFromRadialProfile(deck.config, profile.profile, profile_path);
  assert(init.success);

  const auto write = dec3d::io::WriteInitialRuntimeOutputs(
      deck.config,
      init.state,
      init.geometry,
      init.group_layout,
      profile_path);
  DEC3D_CHECK(write.success);
  DEC3D_CHECK(write.dec3d_out_written);
  DEC3D_CHECK(write.field_checkpoint_written);
  DEC3D_CHECK(write.restart_checkpoint_written);
  DEC3D_CHECK(write.history_profile_written);
  DEC3D_CHECK(write.report_line.find("field_checkpoint_write_s=") != std::string::npos);
  DEC3D_CHECK(write.report_line.find("restart_checkpoint_write_s=") != std::string::npos);
  DEC3D_CHECK(write.report_line.find("history_profile_written=true") != std::string::npos);

  const auto dec3d_out = ReadText(root / "output" / "dec3d.out");
  assert(dec3d_out.find("step time_s dt_s wall_elapsed_s rho_min rho_max") != std::string::npos);
  assert(dec3d_out.find("field_checkpoint_write_s restart_checkpoint_write_s") !=
         std::string::npos);
  assert(dec3d_out.find("0 0") != std::string::npos);

  const auto field = ReadText(root / "output" / "fields_000000.snap");
  assert(field.find("diagnostic_id=p5.io.field_checkpoint") != std::string::npos);
  assert(field.find("restart_compatible=false") != std::string::npos);
  assert(field.find("field,group,radial,theta,phi,value") != std::string::npos);
  assert(field.find("Te_keV,-1,") != std::string::npos);
  assert(field.find("vr_cm_s,-1,") != std::string::npos);

  const auto restart = ReadText(root / "output" / "restart_000000.snap");
  assert(restart.find("diagnostic_id=p5.io.restart_checkpoint") != std::string::npos);
  assert(restart.find("restart_compatible=true") != std::string::npos);
  assert(restart.find("field,group,radial,theta,phi,value") != std::string::npos);
  assert(restart.find("rho,-1,") != std::string::npos);
  assert(restart.find("mom_r,-1,") != std::string::npos);
  assert(restart.find("e_electron,-1,") != std::string::npos);
  assert(restart.find("radiation_groups,0,") != std::string::npos);
  assert(restart.find("alpha_state,-1,") != std::string::npos);

  const auto history = ReadText(root / "output" / "output_contract.his");
  DEC3D_CHECK(history.find("# diagnostic_id=p5.io.history_profile") != std::string::npos);
  DEC3D_CHECK(history.find("step,time_s,r_index,r_cm,r_um,rho_g_cm3,Te_keV,Ti_keV,vr_cm_s,vtheta_cm_s,vphi_cm_s") != std::string::npos);
  DEC3D_CHECK(history.find("0,0,0,") != std::string::npos);
  DEC3D_CHECK(history.find("0,0,1,") != std::string::npos);

  auto interval_only_config = deck.config;
  interval_only_config.run.output_dir = (root / "output_interval_only").generic_string();
  interval_only_config.output.field_checkpoint_every_steps = 0;
  interval_only_config.output.field_checkpoint_interval_s = 1.0e-11;
  interval_only_config.output.restart_checkpoint_every_steps = 0;
  interval_only_config.output.restart_checkpoint_interval_s = 0.0;
  interval_only_config.output.history_profile_interval_s = 0.0;
  interval_only_config.output.history_profile_file = "";
  const auto interval_only_write = dec3d::io::WriteInitialRuntimeOutputs(
      interval_only_config,
      init.state,
      init.geometry,
      init.group_layout,
      profile_path);
  assert(interval_only_write.success);
  assert(interval_only_write.field_checkpoint_written);
  assert(!interval_only_write.restart_checkpoint_written);
  const auto interval_only_field =
      ReadText(root / "output_interval_only" / "fields_000000.snap");
  assert(interval_only_field.find("# step=0") != std::string::npos);
  assert(interval_only_field.find("# time_s=0") != std::string::npos);

  std::filesystem::remove_all(root);
  return 0;
}
