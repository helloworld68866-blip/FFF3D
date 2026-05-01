#include "io/runtime_output.hpp"

#include "physics/units/physical_constants.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace dec3d::io {
namespace {

[[nodiscard]] double SecondsSince(std::chrono::steady_clock::time_point start) noexcept {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

RuntimeOutputWriteResult FailOutput(std::string reason) {
  RuntimeOutputWriteResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics =
      "diagnostic_id=p5.io.runtime_output.failure; failure_reason=" +
      result.failure_reason;
  return result;
}

[[nodiscard]] bool CheckpointTriggerEnabled(
    const CheckpointTriggerConfig& config) noexcept {
  return config.every_steps > 0 || config.interval_s > 0.0;
}

[[nodiscard]] std::filesystem::path HistoryProfilePath(
    const InputDeckConfig& config) {
  const std::filesystem::path output_dir(config.run.output_dir);
  const std::string filename = config.output.history_profile_file.empty()
                                   ? config.run.case_name + ".his"
                                   : config.output.history_profile_file;
  return output_dir / filename;
}

[[nodiscard]] std::size_t LinearIndex(
    const dec3d::state::CanonicalStateLayout& layout,
    std::size_t r,
    std::size_t t,
    std::size_t p) noexcept {
  return (r * layout.theta_cells + t) * layout.phi_cells + p;
}

void WriteFieldRows(
    std::ofstream& out,
    const char* name,
    int group,
    const dec3d::core::Array3D<double>& field) {
  out << std::setprecision(17);
  for (std::size_t r = 0u; r < field.extent_r(); ++r) {
    for (std::size_t t = 0u; t < field.extent_theta(); ++t) {
      for (std::size_t p = 0u; p < field.extent_phi(); ++p) {
        out << name << ',' << group << ',' << r << ',' << t << ',' << p << ','
            << field(r, t, p) << '\n';
      }
    }
  }
}

void WriteDerivedRow(
    std::ofstream& out,
    const char* name,
    std::size_t r,
    std::size_t t,
    std::size_t p,
    double value) {
  out << name << ",-1," << r << ',' << t << ',' << p << ',' << value << '\n';
}

void WriteSnapshotMetadata(
    std::ofstream& out,
    const char* diagnostic_id,
    bool restart_compatible,
    const dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s) {
  out << "# diagnostic_id=" << diagnostic_id << '\n'
      << "# format_version=1\n"
      << "# restart_compatible=" << (restart_compatible ? "true" : "false") << '\n'
      << "# step=" << step << '\n'
      << "# time_s=" << time_s << '\n'
      << "# radial_cells=" << state.layout.radial_cells << '\n'
      << "# theta_cells=" << state.layout.theta_cells << '\n'
      << "# phi_cells=" << state.layout.phi_cells << '\n'
      << "# radiation_group_count=" << group_layout.group_count << '\n'
      << "# source_profile=" << profile_path.string() << '\n'
      << "field,group,radial,theta,phi,value\n";
}

struct ScalarExtrema {
  double minimum{std::numeric_limits<double>::infinity()};
  double maximum{-std::numeric_limits<double>::infinity()};
};

void Observe(ScalarExtrema& extrema, double value) noexcept {
  extrema.minimum = std::min(extrema.minimum, value);
  extrema.maximum = std::max(extrema.maximum, value);
}

[[nodiscard]] bool RecoverForRuntimeOutput(
    const dec3d::state::CanonicalState& state,
    dec3d::state::ThermodynamicRecoveryResult& recovery,
    std::string& failure_reason) {
  recovery = dec3d::state::RecoverThermodynamicState(state);
  if (!recovery.success || !dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovery)) {
    failure_reason = recovery.failure_reason.empty() ? "thermodynamic recovery failed"
                                                     : recovery.failure_reason;
    return false;
  }
  return true;
}

[[nodiscard]] RuntimeOutputExtrema ToPublicExtrema(
    const ScalarExtrema& rho,
    const ScalarExtrema& te,
    const ScalarExtrema& ti) noexcept {
  RuntimeOutputExtrema extrema;
  extrema.rho_min = rho.minimum;
  extrema.rho_max = rho.maximum;
  extrema.te_min_keV = te.minimum;
  extrema.te_max_keV = te.maximum;
  extrema.ti_min_keV = ti.minimum;
  extrema.ti_max_keV = ti.maximum;
  return extrema;
}

void ComputeExtremaFromRecovery(
    const dec3d::state::CanonicalState& state,
    const dec3d::state::ThermodynamicRecoveryResult& recovery,
    ScalarExtrema& rho,
    ScalarExtrema& te,
    ScalarExtrema& ti) noexcept {
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0u; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0u; p < state.layout.phi_cells; ++p) {
        Observe(rho, state.rho(r, t, p));
        Observe(te, recovery.cells(r, t, p).t_e_keV);
        Observe(ti, recovery.cells(r, t, p).t_i_keV);
      }
    }
  }
}

[[nodiscard]] bool WriteFieldCheckpointFile(
    const std::filesystem::path& path,
    const dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s,
    const dec3d::state::ThermodynamicRecoveryResult& recovery,
    std::string& failure_reason) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    failure_reason = "failed to open field checkpoint";
    return false;
  }
  out << std::setprecision(17);
  WriteSnapshotMetadata(out, "p5.io.field_checkpoint", false, state, group_layout, profile_path,
                        step, time_s);
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0u; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0u; p < state.layout.phi_cells; ++p) {
        const double density = state.rho(r, t, p);
        const double inv_rho = density > 0.0 ? 1.0 / density : 0.0;
        WriteDerivedRow(out, "rho", r, t, p, density);
        WriteDerivedRow(out, "Te_keV", r, t, p, recovery.cells(r, t, p).t_e_keV);
        WriteDerivedRow(out, "Ti_keV", r, t, p, recovery.cells(r, t, p).t_i_keV);
        WriteDerivedRow(out, "vr_cm_s", r, t, p, state.mom_r(r, t, p) * inv_rho);
        WriteDerivedRow(out, "vtheta_cm_s", r, t, p, state.mom_theta(r, t, p) * inv_rho);
        WriteDerivedRow(out, "vphi_cm_s", r, t, p, state.mom_phi(r, t, p) * inv_rho);
      }
    }
  }
  return true;
}

[[nodiscard]] bool WriteFieldCheckpointFile(
    const std::filesystem::path& path,
    const dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s,
    std::string& failure_reason) {
  dec3d::state::ThermodynamicRecoveryResult recovery;
  if (!RecoverForRuntimeOutput(state, recovery, failure_reason)) {
    return false;
  }
  return WriteFieldCheckpointFile(
      path, state, group_layout, profile_path, step, time_s, recovery, failure_reason);
}

[[nodiscard]] bool WriteHistoryProfileFile(
    const std::filesystem::path& path,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t step,
    double time_s,
    const dec3d::state::ThermodynamicRecoveryResult& recovery,
    bool truncate,
    std::string& failure_reason) {
  const auto& layout = state.layout;
  const std::size_t cell_count =
      layout.radial_cells * layout.theta_cells * layout.phi_cells;
  if (!geometry.valid ||
      geometry.radial_faces.size() != layout.radial_cells + 1u ||
      geometry.cell_volumes.size() != cell_count) {
    failure_reason = "invalid geometry for history profile";
    return false;
  }

  const bool write_header = truncate || !std::filesystem::exists(path);
  const auto mode = std::ios::binary | (truncate ? std::ios::trunc : std::ios::app);
  std::ofstream out(path, mode);
  if (!out) {
    failure_reason = "failed to open history profile";
    return false;
  }
  out << std::setprecision(17);
  if (write_header) {
    out << "# diagnostic_id=p5.io.history_profile\n"
        << "# format_version=1\n"
        << "# angular_average=volume_weighted\n"
        << "# columns=step,time_s,r_index,r_cm,r_um,rho_g_cm3,Te_keV,Ti_keV,"
           "vr_cm_s,vtheta_cm_s,vphi_cm_s\n"
        << "step,time_s,r_index,r_cm,r_um,rho_g_cm3,Te_keV,Ti_keV,"
           "vr_cm_s,vtheta_cm_s,vphi_cm_s\n";
  }

  for (std::size_t r = 0u; r < layout.radial_cells; ++r) {
    double volume_sum = 0.0;
    double rho_sum = 0.0;
    double te_sum = 0.0;
    double ti_sum = 0.0;
    double vr_sum = 0.0;
    double vt_sum = 0.0;
    double vp_sum = 0.0;
    for (std::size_t t = 0u; t < layout.theta_cells; ++t) {
      for (std::size_t p = 0u; p < layout.phi_cells; ++p) {
        const auto linear = LinearIndex(layout, r, t, p);
        const double volume = geometry.cell_volumes[linear];
        const double rho = state.rho(r, t, p);
        if (!(volume > 0.0) || !std::isfinite(volume) ||
            !(rho > 0.0) || !std::isfinite(rho)) {
          failure_reason = "invalid cell state for history profile";
          return false;
        }
        const double inv_rho = 1.0 / rho;
        volume_sum += volume;
        rho_sum += rho * volume;
        te_sum += recovery.cells(r, t, p).t_e_keV * volume;
        ti_sum += recovery.cells(r, t, p).t_i_keV * volume;
        vr_sum += state.mom_r(r, t, p) * inv_rho * volume;
        vt_sum += state.mom_theta(r, t, p) * inv_rho * volume;
        vp_sum += state.mom_phi(r, t, p) * inv_rho * volume;
      }
    }
    if (!(volume_sum > 0.0) || !std::isfinite(volume_sum)) {
      failure_reason = "invalid radial shell volume for history profile";
      return false;
    }
    const double inv_volume = 1.0 / volume_sum;
    const double r_cm = 0.5 * (geometry.radial_faces[r] + geometry.radial_faces[r + 1u]);
    out << step << ',' << time_s << ',' << r << ','
        << r_cm << ',' << r_cm * 1.0e4 << ','
        << rho_sum * inv_volume << ','
        << te_sum * inv_volume << ','
        << ti_sum * inv_volume << ','
        << vr_sum * inv_volume << ','
        << vt_sum * inv_volume << ','
        << vp_sum * inv_volume << '\n';
  }
  return true;
}

[[nodiscard]] bool WriteRestartCheckpointFile(
    const std::filesystem::path& path,
    const dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s,
    std::string& failure_reason) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    failure_reason = "failed to open restart checkpoint";
    return false;
  }
  out << std::setprecision(17);
  WriteSnapshotMetadata(out, "p5.io.restart_checkpoint", true, state, group_layout, profile_path,
                        step, time_s);
  WriteFieldRows(out, "rho", -1, state.rho);
  WriteFieldRows(out, "mom_r", -1, state.mom_r);
  WriteFieldRows(out, "mom_theta", -1, state.mom_theta);
  WriteFieldRows(out, "mom_phi", -1, state.mom_phi);
  WriteFieldRows(out, "e_electron", -1, state.e_electron);
  WriteFieldRows(out, "e_fluid_total", -1, state.e_fluid_total);
  for (std::size_t g = 0u; g < state.radiation_groups.size(); ++g) {
    WriteFieldRows(out, "radiation_groups", static_cast<int>(g), state.radiation_groups[g]);
  }
  WriteFieldRows(out, "alpha_state", -1, state.alpha_state.storage);
  return true;
}

}  // namespace

bool ShouldWriteCheckpoint(
    const CheckpointTriggerConfig& config,
    std::size_t step,
    double time_s,
    double last_written_time_s) noexcept {
  const bool step_enabled = config.every_steps > 0;
  const bool time_enabled = config.interval_s > 0.0;
  if (!step_enabled && !time_enabled) {
    return false;
  }
  if (step_enabled && step % static_cast<std::size_t>(config.every_steps) == 0u) {
    return true;
  }
  if (time_enabled && time_s + 1.0e-30 >= last_written_time_s + config.interval_s) {
    return true;
  }
  return false;
}

std::string FormatStepSnapshotName(const std::string& prefix, std::size_t step) {
  std::ostringstream name;
  name << prefix << '_' << std::setw(6) << std::setfill('0') << step << ".snap";
  return name.str();
}

CheckpointTriggerConfig FieldCheckpointTrigger(const OutputConfig& output) noexcept {
  return CheckpointTriggerConfig{
      output.field_checkpoint_every_steps,
      output.field_checkpoint_interval_s};
}

CheckpointTriggerConfig RestartCheckpointTrigger(const OutputConfig& output) noexcept {
  return CheckpointTriggerConfig{
      output.restart_checkpoint_every_steps,
      output.restart_checkpoint_interval_s};
}

CheckpointTriggerConfig HistoryProfileTrigger(const OutputConfig& output) noexcept {
  return CheckpointTriggerConfig{0, output.history_profile_interval_s};
}

RuntimeOutputExtremaResult ComputeRuntimeOutputExtrema(
    const dec3d::state::CanonicalState& state) noexcept {
  RuntimeOutputExtremaResult result;
  try {
    ScalarExtrema rho;
    ScalarExtrema te;
    ScalarExtrema ti;
    std::string failure;
    dec3d::state::ThermodynamicRecoveryResult recovery;
    if (!RecoverForRuntimeOutput(state, recovery, failure)) {
      result.failure_reason = failure;
      result.failure_diagnostics =
          "diagnostic_id=p5.io.runtime_output_extrema.failure; failure_reason=" +
          result.failure_reason;
      return result;
    }
    ComputeExtremaFromRecovery(state, recovery, rho, te, ti);
    result.success = true;
    result.extrema = ToPublicExtrema(rho, te, ti);
    std::ostringstream report;
    report << std::setprecision(17)
           << "diagnostic_id=p5.io.runtime_output_extrema"
           << "; rho_min=" << result.extrema.rho_min
           << "; rho_max=" << result.extrema.rho_max
           << "; Te_min_keV=" << result.extrema.te_min_keV
           << "; Te_max_keV=" << result.extrema.te_max_keV
           << "; Ti_min_keV=" << result.extrema.ti_min_keV
           << "; Ti_max_keV=" << result.extrema.ti_max_keV;
    result.report_line = report.str();
    return result;
  } catch (const std::exception& e) {
    result.failure_reason = e.what();
  } catch (...) {
    result.failure_reason = "unknown runtime output extrema failure";
  }
  result.failure_diagnostics =
      "diagnostic_id=p5.io.runtime_output_extrema.failure; failure_reason=" +
      result.failure_reason;
  return result;
}

RuntimeOutputWriteResult WriteRuntimeStepDec3DOut(
    const InputDeckConfig& config,
    std::size_t step,
    double time_s,
    double dt_s,
    double wall_elapsed_s,
    const std::string& status,
    const RuntimeOutputExtrema& extrema) noexcept {
  try {
    RuntimeOutputWriteResult result;
    result.success = true;
    if (config.output.write_dec3d_out_every_steps <= 0 ||
        step % static_cast<std::size_t>(config.output.write_dec3d_out_every_steps) != 0u) {
      result.report_line =
          "diagnostic_id=p5.io.runtime_step_output; dec3d_out_written=false";
      return result;
    }

    const std::filesystem::path output_dir(config.run.output_dir);
    std::filesystem::create_directories(output_dir);
    const auto output_start = std::chrono::steady_clock::now();
    const auto write_start = std::chrono::steady_clock::now();
    std::ofstream log(output_dir / "dec3d.out", std::ios::binary | std::ios::app);
    if (!log) {
      return FailOutput("failed to open dec3d.out");
    }
    log << std::setprecision(17);
    log << step << ' ' << time_s << ' ' << dt_s << ' ' << wall_elapsed_s << ' '
        << extrema.rho_min << ' ' << extrema.rho_max << ' '
        << extrema.te_min_keV << ' ' << extrema.te_max_keV << ' '
        << extrema.ti_min_keV << ' ' << extrema.ti_max_keV << ' '
        << 0.0 << ' ' << 0.0 << ' ' << status << '\n';
    result.dec3d_out_written = true;
    result.dec3d_out_write_s = SecondsSince(write_start);
    result.runtime_output_total_s = SecondsSince(output_start);
    std::ostringstream report;
    report << std::setprecision(17)
           << "diagnostic_id=p5.io.runtime_step_output"
           << "; step=" << step
           << "; time_s=" << time_s
           << "; dec3d_out_written=true"
           << "; field_checkpoint_written=false"
           << "; restart_checkpoint_written=false"
           << "; history_profile_written=false"
           << "; dec3d_out_write_s=" << result.dec3d_out_write_s
           << "; field_checkpoint_write_s=0"
           << "; restart_checkpoint_write_s=0"
           << "; history_profile_write_s=0"
           << "; runtime_output_total_s=" << result.runtime_output_total_s;
    result.report_line = report.str();
    return result;
  } catch (const std::exception& e) {
    return FailOutput(e.what());
  } catch (...) {
    return FailOutput("unknown runtime step dec3d.out failure");
  }
}

RuntimeOutputWriteResult WriteInitialRuntimeOutputs(
    const InputDeckConfig& config,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path) noexcept {
  try {
    const std::filesystem::path output_dir(config.run.output_dir);
    const auto output_start = std::chrono::steady_clock::now();
    std::filesystem::create_directories(output_dir);

    ScalarExtrema rho;
    ScalarExtrema te;
    ScalarExtrema ti;
    std::string failure;
    dec3d::state::ThermodynamicRecoveryResult recovery;
    if (!RecoverForRuntimeOutput(state, recovery, failure)) {
      return FailOutput(failure);
    }
    ComputeExtremaFromRecovery(state, recovery, rho, te, ti);
    const auto extrema = ToPublicExtrema(rho, te, ti);

    RuntimeOutputWriteResult result;
    result.success = true;

    if (config.output.write_dec3d_out_every_steps > 0) {
      const auto write_start = std::chrono::steady_clock::now();
      std::ofstream log(output_dir / "dec3d.out", std::ios::binary | std::ios::trunc);
      if (!log) {
        return FailOutput("failed to open dec3d.out");
      }
      log << std::setprecision(17);
      log << "step time_s dt_s wall_elapsed_s rho_min rho_max "
          << "Te_min_keV Te_max_keV Ti_min_keV Ti_max_keV "
          << "field_checkpoint_write_s restart_checkpoint_write_s status\n";
      log << 0 << ' ' << 0.0 << ' ' << 0.0 << ' ' << 0.0 << ' '
          << extrema.rho_min << ' ' << extrema.rho_max << ' '
          << extrema.te_min_keV << ' ' << extrema.te_max_keV << ' '
          << extrema.ti_min_keV << ' ' << extrema.ti_max_keV << ' '
          << 0.0 << ' ' << 0.0 << " initialized\n";
      result.dec3d_out_written = true;
      result.dec3d_out_write_s = SecondsSince(write_start);
    }

    if (CheckpointTriggerEnabled(FieldCheckpointTrigger(config.output))) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteFieldCheckpointFile(
              output_dir / FormatStepSnapshotName(config.output.field_checkpoint_prefix, 0u),
              state, group_layout, profile_path, 0u, 0.0, recovery, failure)) {
        return FailOutput(failure);
      }
      result.field_checkpoint_written = true;
      result.field_checkpoint_write_s = SecondsSince(write_start);
    }

    if (CheckpointTriggerEnabled(RestartCheckpointTrigger(config.output))) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteRestartCheckpointFile(
              output_dir / FormatStepSnapshotName(config.output.restart_checkpoint_prefix, 0u),
              state, group_layout, profile_path, 0u, 0.0, failure)) {
        return FailOutput(failure);
      }
      result.restart_checkpoint_written = true;
      result.restart_checkpoint_write_s = SecondsSince(write_start);
    }

    if (CheckpointTriggerEnabled(HistoryProfileTrigger(config.output))) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteHistoryProfileFile(
              HistoryProfilePath(config), state, geometry, 0u, 0.0, recovery,
              true, failure)) {
        return FailOutput(failure);
      }
      result.history_profile_written = true;
      result.history_profile_write_s = SecondsSince(write_start);
    }
    result.runtime_output_total_s = SecondsSince(output_start);

    std::ostringstream report;
    report << std::setprecision(17)
           << "diagnostic_id=p5.io.runtime_output"
           << "; dec3d_out_written=" << (result.dec3d_out_written ? "true" : "false")
           << "; field_checkpoint_written="
           << (result.field_checkpoint_written ? "true" : "false")
           << "; restart_checkpoint_written="
           << (result.restart_checkpoint_written ? "true" : "false")
           << "; history_profile_written="
           << (result.history_profile_written ? "true" : "false")
           << "; field_checkpoint_trigger_steps_enabled="
           << (config.output.field_checkpoint_every_steps > 0 ? "true" : "false")
           << "; field_checkpoint_trigger_time_enabled="
           << (config.output.field_checkpoint_interval_s > 0.0 ? "true" : "false")
           << "; restart_checkpoint_trigger_steps_enabled="
           << (config.output.restart_checkpoint_every_steps > 0 ? "true" : "false")
           << "; restart_checkpoint_trigger_time_enabled="
           << (config.output.restart_checkpoint_interval_s > 0.0 ? "true" : "false")
           << "; history_profile_time_trigger_enabled="
           << (config.output.history_profile_interval_s > 0.0 ? "true" : "false")
           << "; dec3d_out_write_s=" << result.dec3d_out_write_s
           << "; field_checkpoint_write_s=" << result.field_checkpoint_write_s
           << "; restart_checkpoint_write_s=" << result.restart_checkpoint_write_s
           << "; history_profile_write_s=" << result.history_profile_write_s
           << "; runtime_output_total_s=" << result.runtime_output_total_s;
    result.report_line = report.str();
    return result;
  } catch (const std::exception& e) {
    return FailOutput(e.what());
  } catch (...) {
    return FailOutput("unknown runtime output failure");
  }
}

RuntimeOutputWriteResult WriteRuntimeStepOutputs(
    const InputDeckConfig& config,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s,
    double dt_s,
    double wall_elapsed_s,
    const std::string& status,
    RuntimeOutputStepState& output_state) noexcept {
  try {
    const std::filesystem::path output_dir(config.run.output_dir);
    const auto output_start = std::chrono::steady_clock::now();
    std::filesystem::create_directories(output_dir);
    RuntimeOutputWriteResult result;
    result.success = true;

    ScalarExtrema rho;
    ScalarExtrema te;
    ScalarExtrema ti;
    std::string failure;
    dec3d::state::ThermodynamicRecoveryResult recovery;
    if (!RecoverForRuntimeOutput(state, recovery, failure)) {
      return FailOutput(failure);
    }
    ComputeExtremaFromRecovery(state, recovery, rho, te, ti);
    const auto extrema = ToPublicExtrema(rho, te, ti);

    if (ShouldWriteCheckpoint(FieldCheckpointTrigger(config.output), step, time_s,
                              output_state.last_field_checkpoint_time_s)) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteFieldCheckpointFile(
              output_dir / FormatStepSnapshotName(config.output.field_checkpoint_prefix, step),
              state, group_layout, profile_path, step, time_s, recovery, failure)) {
        return FailOutput(failure);
      }
      output_state.last_field_checkpoint_time_s = time_s;
      result.field_checkpoint_written = true;
      result.field_checkpoint_write_s = SecondsSince(write_start);
    }

    if (ShouldWriteCheckpoint(RestartCheckpointTrigger(config.output), step, time_s,
                              output_state.last_restart_checkpoint_time_s)) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteRestartCheckpointFile(
              output_dir / FormatStepSnapshotName(config.output.restart_checkpoint_prefix, step),
              state, group_layout, profile_path, step, time_s, failure)) {
        return FailOutput(failure);
      }
      output_state.last_restart_checkpoint_time_s = time_s;
      result.restart_checkpoint_written = true;
      result.restart_checkpoint_write_s = SecondsSince(write_start);
    }

    if (ShouldWriteCheckpoint(HistoryProfileTrigger(config.output), step, time_s,
                              output_state.last_history_profile_time_s)) {
      const auto write_start = std::chrono::steady_clock::now();
      if (!WriteHistoryProfileFile(
              HistoryProfilePath(config), state, geometry, step, time_s, recovery,
              false, failure)) {
        return FailOutput(failure);
      }
      output_state.last_history_profile_time_s = time_s;
      result.history_profile_written = true;
      result.history_profile_write_s = SecondsSince(write_start);
    }

    if (config.output.write_dec3d_out_every_steps > 0 &&
        step % static_cast<std::size_t>(config.output.write_dec3d_out_every_steps) == 0u) {
      const auto write_start = std::chrono::steady_clock::now();
      std::ofstream log(output_dir / "dec3d.out", std::ios::binary | std::ios::app);
      if (!log) {
        return FailOutput("failed to open dec3d.out");
      }
      log << std::setprecision(17);
      log << step << ' ' << time_s << ' ' << dt_s << ' ' << wall_elapsed_s << ' '
          << extrema.rho_min << ' ' << extrema.rho_max << ' '
          << extrema.te_min_keV << ' ' << extrema.te_max_keV << ' '
          << extrema.ti_min_keV << ' ' << extrema.ti_max_keV << ' '
          << result.field_checkpoint_write_s << ' '
          << result.restart_checkpoint_write_s << ' ' << status << '\n';
      result.dec3d_out_written = true;
      result.dec3d_out_write_s = SecondsSince(write_start);
    }
    result.runtime_output_total_s = SecondsSince(output_start);

    std::ostringstream report;
    report << std::setprecision(17)
           << "diagnostic_id=p5.io.runtime_step_output"
           << "; step=" << step
           << "; time_s=" << time_s
           << "; dec3d_out_written=" << (result.dec3d_out_written ? "true" : "false")
           << "; field_checkpoint_written="
           << (result.field_checkpoint_written ? "true" : "false")
           << "; restart_checkpoint_written="
           << (result.restart_checkpoint_written ? "true" : "false")
           << "; history_profile_written="
           << (result.history_profile_written ? "true" : "false")
           << "; dec3d_out_write_s=" << result.dec3d_out_write_s
           << "; field_checkpoint_write_s=" << result.field_checkpoint_write_s
           << "; restart_checkpoint_write_s=" << result.restart_checkpoint_write_s
           << "; history_profile_write_s=" << result.history_profile_write_s
           << "; runtime_output_total_s=" << result.runtime_output_total_s;
    result.report_line = report.str();
    return result;
  } catch (const std::exception& e) {
    return FailOutput(e.what());
  } catch (...) {
    return FailOutput("unknown runtime step output failure");
  }
}

}  // namespace dec3d::io
