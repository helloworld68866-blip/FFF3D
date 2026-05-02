#pragma once

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace dec3d::io {

struct RuntimeArgumentResult {
  bool success{false};
  std::filesystem::path input_deck_path;
  std::filesystem::path profile_path;
  std::filesystem::path restart_path;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

struct RunConfig {
  std::string case_name;
  std::string phase;
  std::vector<char> stage_order;
  int step_count{0};
  double target_time_s{0.0};
  std::string dt_mode;
  double cfl{0.0};
  std::string output_dir;
};

struct MeshConfig {
  std::string geometry;
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
  double radial_min_cm{0.0};
  double radial_max_cm{0.0};
  double theta_min{0.0};
  double theta_max{0.0};
  double phi_min{0.0};
  double phi_max{0.0};
  bool moving_mesh{false};
  bool macro_zoning{false};
};

struct InitialConditionConfig {
  std::string profile_interpolation;
  std::string profile_radius_unit;
  std::string outside_profile_policy;
};

struct PerturbationConfig {
  bool enabled{false};
  std::string type;
  int ell{0};
  int m{0};
  double amplitude{0.0};
  double r0_cm{0.0};
  std::string target;
};

struct PhysicsConfig {
  bool enable_hydro{false};
  bool enable_thermal{false};
  bool enable_equilibration{false};
  bool enable_radiation{false};
  bool enable_alpha{false};
};

struct RadiationConfig {
  std::string group_mode;
  std::vector<double> group_edges_eV;
  std::string opacity_provider;
  std::string radiation_initialization;
  double radiation_initial_blackbody_scale{1.0};
  std::string boundary_model;
  std::string flux_limiter;
};

struct ThermalConfig {
  std::string kappa_model;
  std::string electron_flux_limiter;
};

struct AlphaConfig {
  std::string composition_model;
  std::string reactivity_model;
  std::string alpha_initialization;
};

struct BoundaryConfig {
  std::string inner_radial;
  std::string outer_radial;
  std::string theta;
  std::string phi;
};

struct OutputConfig {
  int write_profiles_every{0};
  int write_diagnostics_every{0};
  int write_restart_every{0};
  int write_dec3d_out_every_steps{1};
  int field_checkpoint_every_steps{0};
  double field_checkpoint_interval_s{0.0};
  std::string field_checkpoint_prefix{"fields"};
  std::string field_checkpoint_format{"csv3d"};
  int restart_checkpoint_every_steps{0};
  double restart_checkpoint_interval_s{0.0};
  std::string restart_checkpoint_prefix{"restart"};
  std::string restart_checkpoint_format{"dec3d_restart_text"};
  double history_profile_interval_s{0.0};
  std::string history_profile_file;
};

struct InputDeckConfig {
  RunConfig run;
  MeshConfig mesh;
  InitialConditionConfig initial_condition;
  PerturbationConfig perturbation;
  PhysicsConfig physics;
  RadiationConfig radiation;
  ThermalConfig thermal;
  AlphaConfig alpha;
  BoundaryConfig boundaries;
  OutputConfig output;
};

struct InputDeckLoadResult {
  bool success{false};
  InputDeckConfig config;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] RuntimeArgumentResult ParseP5IORuntimeArguments(
    const std::vector<std::string>& argv) noexcept;

[[nodiscard]] InputDeckLoadResult LoadInputDeck(
    const std::filesystem::path& input_deck_path) noexcept;

[[nodiscard]] bool ValidateInputDeckDiagnostics(
    const InputDeckLoadResult& result) noexcept;

}  // namespace dec3d::io
