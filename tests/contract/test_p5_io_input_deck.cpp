#include "io/input_deck.hpp"
#include "test_assert.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

void WriteText(const std::filesystem::path& path, const std::string& text) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  out << text;
}

std::filesystem::path TempRoot() {
  auto root = std::filesystem::temp_directory_path() / "dec3d_p5_io_input_deck_contract";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  return root;
}

std::string ValidDeckText() {
  return R"ini(
# smoke.in: P5-IO input deck
[run] # run controls
case_name = p5_io_smoke # case name
phase = P5 # phase
stage_order = H,T,E,R,A # stage order
step_count = 0 # initialization-only step count
target_time_s = 1.0e-10 # optional physical stop time
dt_mode = hydro_relaxed # dt mode
cfl = 0.4 # CFL number
output_dir = analysis/output/p5_io_smoke # output directory

[mesh] # mesh controls
geometry = spherical # geometry
radial_cells = 4 # radial cell count
theta_cells = 2 # theta cell count
phi_cells = 4 # phi cell count
radial_min_cm = 0.0 # minimum radius
radial_max_cm = 1.0e-2 # maximum radius
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = true # moving mesh
macro_zoning = true # macro zoning

[initial_condition] # initial condition controls
profile_interpolation = linear # interpolation
profile_radius_unit = um # radius unit
outside_profile_policy = hard_fail # outside profile policy

[physics] # physics switches
enable_hydro = true # H
enable_thermal = true # T
enable_equilibration = true # E
enable_radiation = true # R
enable_alpha = true # A

[radiation] # radiation controls
group_mode = explicit_frequency_groups # group mode
group_edges_eV = 1,10,100 # group edges
opacity_provider = tops_dt_tabulated # opacity provider
radiation_initialization = local_blackbody # radiation initialization
radiation_initial_blackbody_scale = 1.0 # blackbody scale
boundary_model = thesis_marshak_vacuum # boundary model
flux_limiter = harmonic_eq_5_209 # flux limiter

[thermal] # thermal controls
kappa_model = lee_more_with_degeneracy # kappa model
electron_flux_limiter = minmax_old_time_face_effective_kappa # electron limiter

[alpha] # alpha controls
composition_model = equimolar_dt_from_p2_recovery # composition model
reactivity_model = bosch_hale_dt # reactivity model
alpha_initialization = zero # alpha initialization

[boundaries] # boundary controls
inner_radial = scalar_origin_remap_required # inner radial boundary
outer_radial = neumann_zero_flux # outer radial boundary
theta = scalar_pole_remap_required # theta boundary
phi = periodic # phi boundary

[output] # output controls
write_profiles_every = 1 # profile cadence
write_diagnostics_every = 1 # diagnostics cadence
write_restart_every = 10 # restart cadence
write_dec3d_out_every_steps = 1 # dec3d.out cadence
field_checkpoint_every_steps = 10 # field checkpoint step cadence
field_checkpoint_interval_s = 0.0 # disabled field checkpoint time cadence
field_checkpoint_prefix = fields # field checkpoint prefix
field_checkpoint_format = csv3d # field checkpoint format
restart_checkpoint_every_steps = 100 # restart checkpoint step cadence
restart_checkpoint_interval_s = 0.0 # disabled restart checkpoint time cadence
restart_checkpoint_prefix = restart # restart checkpoint prefix
restart_checkpoint_format = dec3d_restart_text # restart checkpoint format
history_profile_interval_s = 1.0e-10 # angular-average history profile time cadence
history_profile_file = p5_io_smoke.his # angular-average history profile file
)ini";
}

}  // namespace

int main() {
  const auto root = TempRoot();
  const auto deck_path = root / "smoke.in";
  const auto profile_path = root / "smoke.pro";
  WriteText(deck_path, ValidDeckText());
  WriteText(profile_path, "r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s\n0 1 1 1 0\n1 1 1 1 0\n");

  {
    const std::vector<std::string> argv{
        "dec3d.exe", "-input", deck_path.string(), "-profile", profile_path.string()};
    const auto args = dec3d::io::ParseP5IORuntimeArguments(argv);
    assert(args.success);
    assert(args.input_deck_path == deck_path);
    assert(args.profile_path == profile_path);
    assert(args.report_line.find("profile_file_source=command_line.profile") != std::string::npos);
    assert(args.report_line.find("implicit_profile_used=false") != std::string::npos);
  }

  {
    const auto deck = dec3d::io::LoadInputDeck(deck_path);
    DEC3D_CHECK(deck.success);
    DEC3D_CHECK(dec3d::io::ValidateInputDeckDiagnostics(deck));
    DEC3D_CHECK(deck.config.run.case_name == "p5_io_smoke");
    DEC3D_CHECK(deck.config.run.target_time_s == 1.0e-10);
    DEC3D_CHECK(deck.config.run.stage_order.size() == 5u);
    DEC3D_CHECK(deck.config.mesh.radial_cells == 4u);
    DEC3D_CHECK(deck.config.radiation.group_edges_eV.size() == 3u);
    DEC3D_CHECK(deck.config.initial_condition.profile_interpolation == "linear");
    DEC3D_CHECK(deck.config.output.write_dec3d_out_every_steps == 1);
    DEC3D_CHECK(deck.config.output.field_checkpoint_every_steps == 10);
    DEC3D_CHECK(deck.config.output.field_checkpoint_interval_s == 0.0);
    DEC3D_CHECK(deck.config.output.field_checkpoint_prefix == "fields");
    DEC3D_CHECK(deck.config.output.field_checkpoint_format == "csv3d");
    DEC3D_CHECK(deck.config.output.restart_checkpoint_every_steps == 100);
    DEC3D_CHECK(deck.config.output.restart_checkpoint_interval_s == 0.0);
    DEC3D_CHECK(deck.config.output.restart_checkpoint_prefix == "restart");
    DEC3D_CHECK(deck.config.output.restart_checkpoint_format == "dec3d_restart_text");
    DEC3D_CHECK(deck.config.output.history_profile_interval_s == 1.0e-10);
    DEC3D_CHECK(deck.config.output.history_profile_file == "p5_io_smoke.his");
    DEC3D_CHECK(deck.report_line.find("deck_profile_file_present=false") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("field_checkpoint_step_trigger_enabled=true") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("field_checkpoint_time_trigger_enabled=false") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("restart_checkpoint_step_trigger_enabled=true") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("restart_checkpoint_time_trigger_enabled=false") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("history_profile_time_trigger_enabled=true") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("history_profile_file=p5_io_smoke.his") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("target_time_trigger_enabled=true") != std::string::npos);
    DEC3D_CHECK(deck.report_line.find("target_time_s=1e-10") != std::string::npos);
  }

  {
    const std::vector<std::string> argv{"dec3d.exe", "-input", deck_path.string()};
    const auto missing_profile = dec3d::io::ParseP5IORuntimeArguments(argv);
    assert(!missing_profile.success);
    assert(missing_profile.failure_reason.find("missing -profile") != std::string::npos);
  }

  {
    const auto bad_path = root / "bad_profile_in_deck.in";
    WriteText(bad_path, ValidDeckText() + "\nprofile_file = profiles/hidden.pro # forbidden\n");
    const auto bad_deck = dec3d::io::LoadInputDeck(bad_path);
    assert(!bad_deck.success);
    assert(bad_deck.failure_reason.find("profile_file") != std::string::npos);
  }

  {
    const auto bad_order_path = root / "bad_order.in";
    std::string bad = ValidDeckText();
    const auto pos = bad.find("stage_order = H,T,E,R,A");
    bad.replace(pos, std::string("stage_order = H,T,E,R,A").size(), "stage_order = H,T,Q");
    WriteText(bad_order_path, bad);
    const auto bad_deck = dec3d::io::LoadInputDeck(bad_order_path);
    assert(!bad_deck.success);
    assert(bad_deck.failure_reason.find("stage_order") != std::string::npos);
  }

  std::filesystem::remove_all(root);
  return 0;
}
