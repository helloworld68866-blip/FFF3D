#include "initialization/profile_initializer.hpp"
#include "io/input_deck.hpp"
#include "io/radial_profile.hpp"
#include "test_assert.hpp"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

void WriteText(const std::filesystem::path& path, const std::string& text) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  out << text;
}

std::filesystem::path TempRoot() {
  auto root = std::filesystem::temp_directory_path() / "dec3d_p5_io_initializer_contract";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  return root;
}

std::string DeckText(const std::string& radiation_init, const std::string& alpha_init) {
  return R"ini(
# init.in: initialization test input
[run] # run
case_name = p5_io_init # case name
phase = P5 # phase
stage_order = H,T,E,R,A # order
step_count = 0 # no advance
dt_mode = hydro_relaxed # dt
cfl = 0.4 # CFL
output_dir = analysis/output/p5_io_init # output

[mesh] # mesh
geometry = spherical # spherical
radial_cells = 2 # radial
theta_cells = 2 # theta
phi_cells = 4 # phi
radial_min_cm = 0.0 # r min
radial_max_cm = 1.0e-2 # r max
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = true # moving mesh
macro_zoning = true # macro

[initial_condition] # initial condition
profile_interpolation = linear # interpolation
profile_radius_unit = um # radius unit
outside_profile_policy = hard_fail # outside policy

[physics] # physics
enable_hydro = true # H
enable_thermal = true # T
enable_equilibration = true # E
enable_radiation = true # R
enable_alpha = true # A

[radiation] # radiation
group_mode = explicit_frequency_groups # explicit groups
group_edges_eV = 1,10,100 # eV edges
opacity_provider = tops_dt_tabulated # provider
)ini" + std::string("radiation_initialization = ") + radiation_init + R"ini( # radiation init
radiation_initial_blackbody_scale = 2.0 # scale
boundary_model = thesis_marshak_vacuum # boundary
flux_limiter = harmonic_eq_5_209 # limiter

[thermal] # thermal
kappa_model = lee_more_with_degeneracy # model
electron_flux_limiter = minmax_old_time_face_effective_kappa # limiter

[alpha] # alpha
composition_model = equimolar_dt_from_p2_recovery # composition
reactivity_model = bosch_hale_dt # reactivity
)ini" + std::string("alpha_initialization = ") + alpha_init + R"ini( # alpha init

[boundaries] # boundaries
inner_radial = scalar_origin_remap_required # origin
outer_radial = neumann_zero_flux # outer
theta = scalar_pole_remap_required # pole
phi = periodic # periodic

[output] # output
write_profiles_every = 1 # profiles
write_diagnostics_every = 1 # diagnostics
write_restart_every = 10 # restart
write_dec3d_out_every_steps = 1 # dec3d.out
field_checkpoint_every_steps = 0 # disable field checkpoint in initializer test
field_checkpoint_interval_s = 0.0 # disable field checkpoint time trigger
field_checkpoint_prefix = fields # field prefix
field_checkpoint_format = csv3d # field format
restart_checkpoint_every_steps = 0 # disable restart checkpoint in initializer test
restart_checkpoint_interval_s = 0.0 # disable restart checkpoint time trigger
restart_checkpoint_prefix = restart # restart prefix
restart_checkpoint_format = dec3d_restart_text # restart format
)ini";
}

}  // namespace

int main() {
  const auto root = TempRoot();
  const auto deck_path = root / "init.in";
  const auto profile_path = root / "init.pro";
  WriteText(deck_path, DeckText("scaled_local_blackbody", "profile"));
  WriteText(profile_path, R"pro(
# init.pro: initialization profile
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale # header
0.0 10.0 5.0 4.0 -1.0e7 1.0e5 2.0e5 1.0e8 0.5 # inner
100.0 20.0 1.0 1.0 0.0 0.0 0.0 2.0e8 0.25 # outer
)pro");

  const auto deck = dec3d::io::LoadInputDeck(deck_path);
  assert(deck.success);
  const auto profile =
      dec3d::io::LoadRadialProfile(profile_path, deck.config.initial_condition.profile_radius_unit);
  assert(profile.success);

  const auto init =
      dec3d::initialization::InitializeFromRadialProfile(deck.config, profile.profile, profile_path);
  assert(init.success);
  assert(dec3d::initialization::ValidateProfileInitializationDiagnostics(init));
  assert(init.state.layout.radial_cells == 2u);
  assert(init.state.layout.theta_cells == 2u);
  assert(init.state.layout.phi_cells == 4u);
  assert(init.state.layout.radiation_group_count == 2u);
  assert(init.state.radiation_groups.size() == 2u);
  assert(init.state.last_authoritative_write_mask != 0u);
  assert(init.report_line.find("canonical_state_initialized=true") != std::string::npos);
  assert(init.report_line.find("radiation_blackbody_provider_executed=true") != std::string::npos);
  assert(init.report_line.find("group_edges_unit=eV") != std::string::npos);
  assert(init.report_line.find("group_count=2") != std::string::npos);
  assert(init.report_line.find("double_kB_guard_passed=true") != std::string::npos);
  assert(init.report_line.find("alpha_initialization=profile") != std::string::npos);

  for (std::size_t r = 0; r < init.state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < init.state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < init.state.layout.phi_cells; ++p) {
        assert(init.state.rho(r, t, p) > 0.0);
        assert(init.state.e_electron(r, t, p) > 0.0);
        assert(init.state.e_fluid_total(r, t, p) > init.state.e_electron(r, t, p));
        assert(init.state.alpha_state.storage(r, t, p) >= 0.0);
        for (const auto& group : init.state.radiation_groups) {
          assert(group(r, t, p) >= 0.0);
        }
      }
    }
  }

  WriteText(deck_path, DeckText("zero", "zero"));
  const auto zero_deck = dec3d::io::LoadInputDeck(deck_path);
  assert(zero_deck.success);
  const auto zero_init = dec3d::initialization::InitializeFromRadialProfile(
      zero_deck.config, profile.profile, profile_path);
  assert(zero_init.success);
  assert(zero_init.report_line.find("radiation_initialization=zero") != std::string::npos);
  assert(zero_init.report_line.find("radiation_blackbody_provider_executed=false") !=
         std::string::npos);
  assert(zero_init.report_line.find("alpha_initialization=zero") != std::string::npos);
  for (const auto& group : zero_init.state.radiation_groups) {
    assert(group(0u, 0u, 0u) == 0.0);
  }
  assert(zero_init.state.alpha_state.storage(0u, 0u, 0u) == 0.0);

  const std::string p2_velocity_deck = R"ini(
# init_p2.in: initialization P2 perturbation test input
[run] # run
case_name = p5_io_init_p2 # case name
phase = P5 # phase
stage_order = H,T,E,R,A # order
step_count = 0 # no advance
dt_mode = hydro_relaxed # dt
cfl = 0.4 # CFL
output_dir = analysis/output/p5_io_init_p2 # output

[mesh] # mesh
geometry = spherical # spherical
radial_cells = 2 # radial
theta_cells = 4 # theta
phi_cells = 4 # phi
radial_min_cm = 0.0 # r min
radial_max_cm = 1.0e-2 # r max
theta_min = 0.0 # theta min
theta_max = 3.141592653589793 # theta max
phi_min = 0.0 # phi min
phi_max = 6.283185307179586 # phi max
moving_mesh = false # no ALE
macro_zoning = true # macro

[initial_condition] # initial condition
profile_interpolation = linear # interpolation
profile_radius_unit = um # radius unit
outside_profile_policy = hard_fail # outside policy

[perturbation] # P2 single-mode perturbation
enabled = true # enable perturbation
type = single_mode_radial_velocity # thesis-style radial velocity perturbation
ell = 2 # P2 mode
m = 0 # axisymmetric mode
amplitude = 0.2 # fractional amplitude
r0_cm = 5.0e-3 # inner shell surface radius
target = radial_velocity_cm_s # perturb radial velocity only

[physics] # physics
enable_hydro = true # H
enable_thermal = true # T
enable_equilibration = true # E
enable_radiation = true # R
enable_alpha = true # A

[radiation] # radiation
group_mode = explicit_frequency_groups # explicit groups
group_edges_eV = 1,10,100 # eV edges
opacity_provider = tops_dt_tabulated # provider
radiation_initialization = zero # radiation init
radiation_initial_blackbody_scale = 1.0 # scale
boundary_model = thesis_marshak_vacuum # boundary
flux_limiter = harmonic_eq_5_209 # limiter

[thermal] # thermal
kappa_model = lee_more_with_degeneracy # model
electron_flux_limiter = minmax_old_time_face_effective_kappa # limiter

[alpha] # alpha
composition_model = equimolar_dt_from_p2_recovery # composition
reactivity_model = bosch_hale_dt # reactivity
alpha_initialization = zero # alpha init

[boundaries] # boundaries
inner_radial = scalar_origin_remap_required # origin
outer_radial = neumann_zero_flux # outer
theta = scalar_pole_remap_required # pole
phi = periodic # periodic

[output] # output
write_profiles_every = 1 # profiles
write_diagnostics_every = 1 # diagnostics
write_restart_every = 10 # restart
write_dec3d_out_every_steps = 1 # dec3d.out
field_checkpoint_every_steps = 0 # disable field checkpoint in initializer test
field_checkpoint_interval_s = 0.0 # disable field checkpoint time trigger
field_checkpoint_prefix = fields # field prefix
field_checkpoint_format = csv3d # field format
restart_checkpoint_every_steps = 0 # disable restart checkpoint in initializer test
restart_checkpoint_interval_s = 0.0 # disable restart checkpoint time trigger
restart_checkpoint_prefix = restart # restart prefix
restart_checkpoint_format = dec3d_restart_text # restart format
)ini";
  WriteText(deck_path, p2_velocity_deck);
  WriteText(profile_path, R"pro(
# init_p2.pro: monotone radial profile for P2 perturbation
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale # header
0.0 10.0 5.0 4.0 -1.0e7 0.0 0.0 0.0 1.0 # inner
100.0 20.0 1.0 1.0 0.0 0.0 0.0 0.0 1.0 # outer
)pro");
  const auto p2_deck = dec3d::io::LoadInputDeck(deck_path);
  DEC3D_CHECK(p2_deck.success);
  const auto p2_profile =
      dec3d::io::LoadRadialProfile(profile_path, p2_deck.config.initial_condition.profile_radius_unit);
  DEC3D_CHECK(p2_profile.success);
  const auto p2_init = dec3d::initialization::InitializeFromRadialProfile(
      p2_deck.config, p2_profile.profile, profile_path);
  DEC3D_CHECK(p2_init.success);
  DEC3D_CHECK(p2_init.report_line.find("initial_perturbation_enabled=true") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find(
                  "initial_perturbation_type=single_mode_radial_velocity") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find("initial_perturbation_l=2") != std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find("initial_perturbation_m=0") != std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find(
                  "initial_perturbation_target=radial_velocity_cm_s") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find("initial_perturbation_r0_cm=0.005") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find("density_angular_perturbation=false") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.report_line.find("temperature_angular_perturbation=false") !=
              std::string::npos);
  DEC3D_CHECK(p2_init.state.rho(1u, 0u, 0u) == p2_init.state.rho(1u, 1u, 0u));
  DEC3D_CHECK(p2_init.state.e_electron(1u, 0u, 0u) ==
              p2_init.state.e_electron(1u, 1u, 0u));
  DEC3D_CHECK(p2_init.state.mom_r(0u, 0u, 0u) > p2_init.state.mom_r(0u, 1u, 0u));

  {
    auto config = p2_deck.config;
    config.mesh.dimensionality = dec3d::io::MeshDimensionality::axisymmetric_2d;
    config.mesh.phi_cells = 1u;
    config.mesh.moving_mesh = false;
    const auto result = dec3d::initialization::InitializeFromRadialProfile(
        config, p2_profile.profile, profile_path);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.report_line.find("initial_perturbation_normalization=raw_Pl") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("axisymmetric_mom_phi_zero=true") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("initial_max_abs_mom_phi=0") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("initial_max_abs_v_phi=0") !=
                std::string::npos);
    DEC3D_CHECK(result.state.layout.phi_cells == 1u);
    for (std::size_t r = 0; r < result.state.layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < result.state.layout.theta_cells; ++t) {
        DEC3D_CHECK(result.state.mom_phi(r, t, 0u) == 0.0);
      }
    }
  }

  {
    auto config = p2_deck.config;
    config.mesh.dimensionality = dec3d::io::MeshDimensionality::axisymmetric_2d;
    config.mesh.phi_cells = 1u;
    config.mesh.moving_mesh = false;
    auto profile_with_vp = p2_profile.profile;
    for (auto& row : profile_with_vp.rows) {
      row.vp_cm_s = 1.0e5;
    }
    const auto result = dec3d::initialization::InitializeFromRadialProfile(
        config, profile_with_vp, profile_path);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("axisymmetric_2d requires zero vp_cm_s") !=
                std::string::npos);
  }

  auto old_coordinate_deck = p2_velocity_deck;
  const std::string velocity_type = "type = single_mode_radial_velocity";
  const auto velocity_type_pos = old_coordinate_deck.find(velocity_type);
  DEC3D_CHECK(velocity_type_pos != std::string::npos);
  old_coordinate_deck.replace(
      velocity_type_pos,
      velocity_type.size(),
      "type = single_mode_radial_coordinate");
  WriteText(deck_path, old_coordinate_deck);
  const auto old_coordinate_load = dec3d::io::LoadInputDeck(deck_path);
  DEC3D_CHECK(!old_coordinate_load.success);

  std::filesystem::remove_all(root);
  return 0;
}
