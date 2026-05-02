#include "io/input_deck.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string_view>

namespace dec3d::io {
namespace {

[[nodiscard]] std::string Trim(std::string_view text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string_view::npos) {
    return {};
  }
  const auto last = text.find_last_not_of(" \t\r\n");
  return std::string(text.substr(first, last - first + 1u));
}

[[nodiscard]] std::string StripComment(std::string_view line) {
  const auto pos = line.find('#');
  return Trim(pos == std::string_view::npos ? line : line.substr(0u, pos));
}

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

template <typename T>
[[nodiscard]] bool ParseNumber(const std::string& text, T& out) {
  std::istringstream in(text);
  in >> out;
  return in && in.eof();
}

[[nodiscard]] bool ParseBool(const std::string& text, bool& out) {
  if (text == "true") {
    out = true;
    return true;
  }
  if (text == "false") {
    out = false;
    return true;
  }
  return false;
}

[[nodiscard]] std::vector<std::string> SplitComma(const std::string& text) {
  std::vector<std::string> parts;
  std::stringstream stream(text);
  std::string part;
  while (std::getline(stream, part, ',')) {
    const auto trimmed = Trim(part);
    if (!trimmed.empty()) {
      parts.push_back(trimmed);
    }
  }
  return parts;
}

[[nodiscard]] bool ValidStage(char c) noexcept {
  return c == 'H' || c == 'T' || c == 'E' || c == 'R' || c == 'A';
}

[[nodiscard]] bool StrictlyIncreasing(const std::vector<double>& values) noexcept {
  if (values.size() < 2u) {
    return false;
  }
  for (std::size_t i = 1u; i < values.size(); ++i) {
    if (!(values[i] > values[i - 1u])) {
      return false;
    }
  }
  return true;
}

[[nodiscard]] std::string BuildRuntimeReport(const std::filesystem::path& input,
                                             const std::filesystem::path& profile,
                                             const std::filesystem::path& restart) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "diagnostic_id=p5.io.runtime_arguments"
      << "; input_deck_path=" << input.string()
      << "; profile_file=" << profile.string()
      << "; profile_file_resolved=" << profile.string()
      << "; profile_file_source=command_line.profile"
      << "; deck_profile_file_present=false"
      << "; implicit_profile_used=false"
      << "; restart_file_present=" << (!restart.empty() ? "true" : "false");
  if (!restart.empty()) {
    out << "; restart_file=" << restart.string()
        << "; restart_file_source=command_line.restart";
  }
  return out.str();
}

RuntimeArgumentResult FailRuntime(std::string reason) {
  RuntimeArgumentResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = "diagnostic_id=p5.io.runtime_arguments.failure; failure_reason=" +
                               result.failure_reason;
  return result;
}

InputDeckLoadResult FailDeck(std::string reason) {
  InputDeckLoadResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = "diagnostic_id=p5.io.input_deck.failure; failure_reason=" +
                               result.failure_reason;
  return result;
}

[[nodiscard]] bool RequireString(const std::map<std::string, std::string>& values,
                                 const std::string& key,
                                 std::string& out,
                                 std::string& failure) {
  const auto it = values.find(key);
  if (it == values.end() || it->second.empty()) {
    failure = "missing " + key;
    return false;
  }
  out = it->second;
  return true;
}

[[nodiscard]] bool RequireInt(const std::map<std::string, std::string>& values,
                              const std::string& key,
                              int& out,
                              std::string& failure) {
  const auto it = values.find(key);
  if (it == values.end() || !ParseNumber(it->second, out)) {
    failure = "invalid " + key;
    return false;
  }
  return true;
}

[[nodiscard]] bool RequireSize(const std::map<std::string, std::string>& values,
                               const std::string& key,
                               std::size_t& out,
                               std::string& failure) {
  int parsed = 0;
  if (!RequireInt(values, key, parsed, failure) || parsed <= 0) {
    failure = "invalid " + key;
    return false;
  }
  out = static_cast<std::size_t>(parsed);
  return true;
}

[[nodiscard]] bool RequireDouble(const std::map<std::string, std::string>& values,
                                 const std::string& key,
                                 double& out,
                                 std::string& failure) {
  const auto it = values.find(key);
  if (it == values.end() || !ParseNumber(it->second, out) || !std::isfinite(out)) {
    failure = "invalid " + key;
    return false;
  }
  return true;
}

[[nodiscard]] bool OptionalDouble(const std::map<std::string, std::string>& values,
                                  const std::string& key,
                                  double default_value,
                                  double& out,
                                  std::string& failure) {
  const auto it = values.find(key);
  if (it == values.end()) {
    out = default_value;
    return true;
  }
  if (!ParseNumber(it->second, out) || !std::isfinite(out)) {
    failure = "invalid " + key;
    return false;
  }
  return true;
}

[[nodiscard]] void OptionalString(const std::map<std::string, std::string>& values,
                                  const std::string& key,
                                  std::string& out) {
  const auto it = values.find(key);
  if (it != values.end()) {
    out = it->second;
  }
}

[[nodiscard]] bool RequireBool(const std::map<std::string, std::string>& values,
                               const std::string& key,
                               bool& out,
                               std::string& failure) {
  const auto it = values.find(key);
  if (it == values.end() || !ParseBool(it->second, out)) {
    failure = "invalid " + key;
    return false;
  }
  return true;
}

[[nodiscard]] std::string BuildDeckReport(const InputDeckConfig& config,
                                          const std::filesystem::path& path) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "diagnostic_id=p5.io.input_deck"
      << "; input_deck_path=" << path.string()
      << "; case_name=" << config.run.case_name
      << "; phase=" << config.run.phase
      << "; target_time_trigger_enabled="
      << (config.run.target_time_s > 0.0 ? "true" : "false")
      << "; target_time_s=" << config.run.target_time_s
      << "; profile_file_source=command_line.profile"
      << "; deck_profile_file_present=false"
      << "; implicit_profile_used=false"
      << "; profile_interpolation=" << config.initial_condition.profile_interpolation
      << "; profile_radius_unit=" << config.initial_condition.profile_radius_unit
      << "; outside_profile_policy=" << config.initial_condition.outside_profile_policy
      << "; inner_radial=" << config.boundaries.inner_radial
      << "; outer_radial=" << config.boundaries.outer_radial
      << "; noh_exact_inflow_deck_enabled="
      << (config.boundaries.outer_radial == "noh_exact_inflow" ? "true" : "false");
  if (config.boundaries.outer_radial == "noh_exact_inflow") {
    out << "; noh_exact_inflow_pressure_time_level=current_time_adiabatic"
        << "; noh_exact_inflow_chi_e_time_level=current_time_density_scaled";
  }
  out
      << "; initial_perturbation_enabled="
      << (config.perturbation.enabled ? "true" : "false");
  if (config.perturbation.enabled) {
    out << "; initial_perturbation_type=" << config.perturbation.type
        << "; initial_perturbation_l=" << config.perturbation.ell
        << "; initial_perturbation_m=" << config.perturbation.m
        << "; initial_perturbation_amplitude=" << config.perturbation.amplitude
        << "; initial_perturbation_r0_cm=" << config.perturbation.r0_cm
        << "; initial_perturbation_target=" << config.perturbation.target;
  }
  out
      << "; comment_policy=hash_full_line_and_trailing"
      << "; generated_file_chinese_comments_present=true"
      << "; chinese_comments_runtime_contract=false"
      << "; group_edges_unit=eV"
      << "; group_count=" << (config.radiation.group_edges_eV.empty()
                                  ? 0u
                                  : config.radiation.group_edges_eV.size() - 1u)
      << "; write_dec3d_out_every_steps=" << config.output.write_dec3d_out_every_steps
      << "; field_checkpoint_step_trigger_enabled="
      << (config.output.field_checkpoint_every_steps > 0 ? "true" : "false")
      << "; field_checkpoint_time_trigger_enabled="
      << (config.output.field_checkpoint_interval_s > 0.0 ? "true" : "false")
      << "; restart_checkpoint_step_trigger_enabled="
      << (config.output.restart_checkpoint_every_steps > 0 ? "true" : "false")
      << "; restart_checkpoint_time_trigger_enabled="
      << (config.output.restart_checkpoint_interval_s > 0.0 ? "true" : "false")
      << "; history_profile_time_trigger_enabled="
      << (config.output.history_profile_interval_s > 0.0 ? "true" : "false")
      << "; history_profile_interval_s=" << config.output.history_profile_interval_s
      << "; history_profile_file="
      << (config.output.history_profile_file.empty()
              ? config.run.case_name + ".his"
              : config.output.history_profile_file)
      << "; deck_validation_success=true";
  return out.str();
}

}  // namespace

RuntimeArgumentResult ParseP5IORuntimeArguments(const std::vector<std::string>& argv) noexcept {
  std::filesystem::path input;
  std::filesystem::path profile;
  std::filesystem::path restart;
  for (std::size_t i = 1u; i < argv.size(); ++i) {
    if (argv[i] == "-input") {
      if (i + 1u >= argv.size() || argv[i + 1u].empty()) {
        return FailRuntime("missing -input value");
      }
      input = argv[++i];
      continue;
    }
    if (argv[i] == "-profile") {
      if (i + 1u >= argv.size() || argv[i + 1u].empty()) {
        return FailRuntime("missing -profile value");
      }
      profile = argv[++i];
      continue;
    }
    if (argv[i] == "-restart") {
      if (i + 1u >= argv.size() || argv[i + 1u].empty()) {
        return FailRuntime("missing -restart value");
      }
      restart = argv[++i];
      continue;
    }
    return FailRuntime("unknown argument " + argv[i]);
  }
  if (input.empty()) {
    return FailRuntime("missing -input");
  }
  if (profile.empty()) {
    return FailRuntime("missing -profile");
  }
  RuntimeArgumentResult result;
  result.success = true;
  result.input_deck_path = input;
  result.profile_path = profile;
  result.restart_path = restart;
  result.report_line = BuildRuntimeReport(input, profile, restart);
  return result;
}

InputDeckLoadResult LoadInputDeck(const std::filesystem::path& input_deck_path) noexcept {
  std::ifstream in(input_deck_path);
  if (!in) {
    return FailDeck("missing input deck");
  }

  std::map<std::string, std::map<std::string, std::string>> sections;
  std::string section;
  std::string raw;
  while (std::getline(in, raw)) {
    const auto line = StripComment(raw);
    if (line.empty()) {
      continue;
    }
    if (line.front() == '[' && line.back() == ']') {
      section = Trim(std::string_view(line).substr(1u, line.size() - 2u));
      if (section.empty()) {
        return FailDeck("empty section");
      }
      continue;
    }
    const auto eq = line.find('=');
    if (eq == std::string::npos || section.empty()) {
      return FailDeck("invalid deck line");
    }
    const auto key = Trim(std::string_view(line).substr(0u, eq));
    const auto value = Trim(std::string_view(line).substr(eq + 1u));
    if (key == "profile_file") {
      return FailDeck("profile_file is forbidden in .in; pass -profile at runtime");
    }
    auto& values = sections[section];
    if (values.contains(key)) {
      return FailDeck("duplicate key " + section + "." + key);
    }
    values[key] = value;
  }

  InputDeckLoadResult result;
  std::string failure;
  auto section_values = [&](const std::string& name) -> const std::map<std::string, std::string>* {
    const auto it = sections.find(name);
    return it == sections.end() ? nullptr : &it->second;
  };

  const auto* run = section_values("run");
  const auto* mesh = section_values("mesh");
  const auto* init = section_values("initial_condition");
  const auto* perturbation = section_values("perturbation");
  const auto* physics = section_values("physics");
  const auto* radiation = section_values("radiation");
  const auto* thermal = section_values("thermal");
  const auto* alpha = section_values("alpha");
  const auto* boundaries = section_values("boundaries");
  const auto* output = section_values("output");
  if (!run || !mesh || !init || !physics || !radiation || !thermal || !alpha || !boundaries ||
      !output) {
    return FailDeck("missing required sections");
  }

  auto& cfg = result.config;
  if (!RequireString(*run, "case_name", cfg.run.case_name, failure) ||
      !RequireString(*run, "phase", cfg.run.phase, failure) ||
      !RequireInt(*run, "step_count", cfg.run.step_count, failure) ||
      !OptionalDouble(*run, "target_time_s", 0.0, cfg.run.target_time_s, failure) ||
      !RequireString(*run, "dt_mode", cfg.run.dt_mode, failure) ||
      !RequireDouble(*run, "cfl", cfg.run.cfl, failure) ||
      !RequireString(*run, "output_dir", cfg.run.output_dir, failure)) {
    return FailDeck(failure);
  }
  const auto order_it = run->find("stage_order");
  if (order_it == run->end()) {
    return FailDeck("missing stage_order");
  }
  for (const auto& part : SplitComma(order_it->second)) {
    if (part.size() != 1u || !ValidStage(part.front())) {
      return FailDeck("invalid stage_order");
    }
    cfg.run.stage_order.push_back(part.front());
  }
  if (cfg.run.stage_order.empty() || cfg.run.phase != "P5") {
    return FailDeck("invalid stage_order or phase");
  }

  if (!RequireString(*mesh, "geometry", cfg.mesh.geometry, failure) ||
      !RequireSize(*mesh, "radial_cells", cfg.mesh.radial_cells, failure) ||
      !RequireSize(*mesh, "theta_cells", cfg.mesh.theta_cells, failure) ||
      !RequireSize(*mesh, "phi_cells", cfg.mesh.phi_cells, failure) ||
      !RequireDouble(*mesh, "radial_min_cm", cfg.mesh.radial_min_cm, failure) ||
      !RequireDouble(*mesh, "radial_max_cm", cfg.mesh.radial_max_cm, failure) ||
      !RequireDouble(*mesh, "theta_min", cfg.mesh.theta_min, failure) ||
      !RequireDouble(*mesh, "theta_max", cfg.mesh.theta_max, failure) ||
      !RequireDouble(*mesh, "phi_min", cfg.mesh.phi_min, failure) ||
      !RequireDouble(*mesh, "phi_max", cfg.mesh.phi_max, failure) ||
      !RequireBool(*mesh, "moving_mesh", cfg.mesh.moving_mesh, failure) ||
      !RequireBool(*mesh, "macro_zoning", cfg.mesh.macro_zoning, failure)) {
    return FailDeck(failure);
  }
  if (cfg.mesh.geometry != "spherical" || !(cfg.mesh.radial_max_cm > cfg.mesh.radial_min_cm) ||
      cfg.mesh.phi_cells % 2u != 0u) {
    return FailDeck("invalid mesh");
  }

  if (!RequireString(*init, "profile_interpolation", cfg.initial_condition.profile_interpolation,
                     failure) ||
      !RequireString(*init, "profile_radius_unit", cfg.initial_condition.profile_radius_unit,
                     failure) ||
      !RequireString(*init, "outside_profile_policy", cfg.initial_condition.outside_profile_policy,
                     failure)) {
    return FailDeck(failure);
  }
  if (cfg.initial_condition.profile_interpolation != "linear" ||
      cfg.initial_condition.profile_radius_unit != "um" ||
      cfg.initial_condition.outside_profile_policy != "hard_fail") {
    return FailDeck("unsupported initial_condition option");
  }
  if (perturbation) {
    if (!RequireBool(*perturbation, "enabled", cfg.perturbation.enabled, failure)) {
      return FailDeck(failure);
    }
    if (cfg.perturbation.enabled) {
      if (!RequireString(*perturbation, "type", cfg.perturbation.type, failure) ||
          !RequireInt(*perturbation, "ell", cfg.perturbation.ell, failure) ||
          !RequireInt(*perturbation, "m", cfg.perturbation.m, failure) ||
          !RequireDouble(*perturbation, "amplitude", cfg.perturbation.amplitude, failure) ||
          !RequireDouble(*perturbation, "r0_cm", cfg.perturbation.r0_cm, failure) ||
          !RequireString(*perturbation, "target", cfg.perturbation.target, failure)) {
        return FailDeck(failure);
      }
      if (cfg.perturbation.type != "single_mode_radial_velocity" ||
          cfg.perturbation.ell <= 0 || cfg.perturbation.m != 0 ||
          !std::isfinite(cfg.perturbation.amplitude) ||
          cfg.perturbation.amplitude < 0.0 ||
          !std::isfinite(cfg.perturbation.r0_cm) ||
          !(cfg.perturbation.r0_cm > 0.0) ||
          cfg.perturbation.target != "radial_velocity_cm_s") {
        return FailDeck("unsupported perturbation option");
      }
    }
  }

  if (!RequireBool(*physics, "enable_hydro", cfg.physics.enable_hydro, failure) ||
      !RequireBool(*physics, "enable_thermal", cfg.physics.enable_thermal, failure) ||
      !RequireBool(*physics, "enable_equilibration", cfg.physics.enable_equilibration, failure) ||
      !RequireBool(*physics, "enable_radiation", cfg.physics.enable_radiation, failure) ||
      !RequireBool(*physics, "enable_alpha", cfg.physics.enable_alpha, failure)) {
    return FailDeck(failure);
  }

  if (!RequireString(*radiation, "group_mode", cfg.radiation.group_mode, failure) ||
      !RequireString(*radiation, "opacity_provider", cfg.radiation.opacity_provider, failure) ||
      !RequireString(*radiation, "radiation_initialization",
                     cfg.radiation.radiation_initialization, failure) ||
      !RequireDouble(*radiation, "radiation_initial_blackbody_scale",
                     cfg.radiation.radiation_initial_blackbody_scale, failure) ||
      !RequireString(*radiation, "boundary_model", cfg.radiation.boundary_model, failure) ||
      !RequireString(*radiation, "flux_limiter", cfg.radiation.flux_limiter, failure)) {
    return FailDeck(failure);
  }
  const auto edges_it = radiation->find("group_edges_eV");
  if (edges_it == radiation->end()) {
    return FailDeck("missing group_edges_eV");
  }
  for (const auto& part : SplitComma(edges_it->second)) {
    double value = 0.0;
    if (!ParseNumber(part, value) || !std::isfinite(value) || value <= 0.0) {
      return FailDeck("invalid group_edges_eV");
    }
    cfg.radiation.group_edges_eV.push_back(value);
  }
  if (cfg.radiation.group_mode != "explicit_frequency_groups" ||
      !StrictlyIncreasing(cfg.radiation.group_edges_eV)) {
    return FailDeck("invalid group_edges_eV");
  }
  if (cfg.radiation.radiation_initialization != "local_blackbody" &&
      cfg.radiation.radiation_initialization != "scaled_local_blackbody" &&
      cfg.radiation.radiation_initialization != "zero") {
    return FailDeck("unsupported radiation_initialization");
  }

  if (!RequireString(*thermal, "kappa_model", cfg.thermal.kappa_model, failure) ||
      !RequireString(*thermal, "electron_flux_limiter", cfg.thermal.electron_flux_limiter,
                     failure) ||
      !RequireString(*alpha, "composition_model", cfg.alpha.composition_model, failure) ||
      !RequireString(*alpha, "reactivity_model", cfg.alpha.reactivity_model, failure) ||
      !RequireString(*alpha, "alpha_initialization", cfg.alpha.alpha_initialization, failure) ||
      !RequireString(*boundaries, "inner_radial", cfg.boundaries.inner_radial, failure) ||
      !RequireString(*boundaries, "outer_radial", cfg.boundaries.outer_radial, failure) ||
      !RequireString(*boundaries, "theta", cfg.boundaries.theta, failure) ||
      !RequireString(*boundaries, "phi", cfg.boundaries.phi, failure) ||
      !RequireInt(*output, "write_profiles_every", cfg.output.write_profiles_every, failure) ||
      !RequireInt(*output, "write_diagnostics_every", cfg.output.write_diagnostics_every,
                  failure) ||
      !RequireInt(*output, "write_restart_every", cfg.output.write_restart_every, failure) ||
      !RequireInt(*output, "write_dec3d_out_every_steps",
                  cfg.output.write_dec3d_out_every_steps, failure) ||
      !RequireInt(*output, "field_checkpoint_every_steps",
                  cfg.output.field_checkpoint_every_steps, failure) ||
      !RequireDouble(*output, "field_checkpoint_interval_s",
                     cfg.output.field_checkpoint_interval_s, failure) ||
      !RequireString(*output, "field_checkpoint_prefix", cfg.output.field_checkpoint_prefix,
                     failure) ||
      !RequireString(*output, "field_checkpoint_format", cfg.output.field_checkpoint_format,
                     failure) ||
      !RequireInt(*output, "restart_checkpoint_every_steps",
                  cfg.output.restart_checkpoint_every_steps, failure) ||
      !RequireDouble(*output, "restart_checkpoint_interval_s",
                     cfg.output.restart_checkpoint_interval_s, failure) ||
      !RequireString(*output, "restart_checkpoint_prefix", cfg.output.restart_checkpoint_prefix,
                     failure) ||
      !RequireString(*output, "restart_checkpoint_format", cfg.output.restart_checkpoint_format,
                     failure)) {
    return FailDeck(failure);
  }
  if (!OptionalDouble(*output, "history_profile_interval_s", 0.0,
                      cfg.output.history_profile_interval_s, failure)) {
    return FailDeck(failure);
  }
  OptionalString(*output, "history_profile_file", cfg.output.history_profile_file);
  if (cfg.alpha.alpha_initialization != "zero" && cfg.alpha.alpha_initialization != "profile") {
    return FailDeck("unsupported alpha_initialization");
  }
  if (cfg.output.write_profiles_every < 0 || cfg.output.write_diagnostics_every < 0 ||
      cfg.output.write_restart_every < 0 || cfg.output.write_dec3d_out_every_steps < 0 ||
      cfg.output.field_checkpoint_every_steps < 0 ||
      cfg.output.restart_checkpoint_every_steps < 0 ||
      cfg.output.field_checkpoint_interval_s < 0.0 ||
      cfg.output.restart_checkpoint_interval_s < 0.0 ||
      cfg.output.history_profile_interval_s < 0.0 ||
      cfg.run.target_time_s < 0.0) {
    return FailDeck("output intervals must be nonnegative");
  }
  if (cfg.output.field_checkpoint_format != "csv3d") {
    return FailDeck("unsupported field_checkpoint_format");
  }
  if (cfg.output.restart_checkpoint_format != "dec3d_restart_text") {
    return FailDeck("unsupported restart_checkpoint_format");
  }

  result.success = true;
  result.report_line = BuildDeckReport(cfg, input_deck_path);
  return result;
}

bool ValidateInputDeckDiagnostics(const InputDeckLoadResult& result) noexcept {
  if (!result.success) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p5.io.input_deck") &&
         Contains(line, "profile_file_source=command_line.profile") &&
         Contains(line, "deck_profile_file_present=false") &&
         Contains(line, "implicit_profile_used=false") &&
         Contains(line, "chinese_comments_runtime_contract=false") &&
         Contains(line, "deck_validation_success=true");
}

}  // namespace dec3d::io
