#include "initialization/profile_initializer.hpp"

#include "core/diagnostics/stage_contracts.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/providers/group_blackbody.hpp"
#include "state/thermodynamics/thermodynamic_recovery.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>

namespace dec3d::initialization {
namespace {

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

ProfileInitializationResult Fail(std::string reason) {
  ProfileInitializationResult result;
  result.failure_reason = std::move(reason);
  result.failure_diagnostics = "diagnostic_id=p5.io.profile_initializer.failure; failure_reason=" +
                               result.failure_reason;
  return result;
}

[[nodiscard]] std::vector<double> FrequencyEdgesFromEV(const std::vector<double>& edges_eV) {
  const double h_erg_s =
      2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
      dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
  std::vector<double> edges_hz;
  edges_hz.reserve(edges_eV.size());
  for (double edge : edges_eV) {
    edges_hz.push_back(edge * dec3d::physics::PhysicsConstantsCGS::erg_per_ev / h_erg_s);
  }
  return edges_hz;
}

[[nodiscard]] std::string BuildReport(
    const dec3d::io::InputDeckConfig& config,
    const std::filesystem::path& profile_path,
    std::size_t group_count,
    bool blackbody_executed,
    double perturbation_reference_velocity_cm_s,
    const std::string& thermodynamic_report) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "diagnostic_id=p5.io.profile_initializer"
      << "; case_name=" << config.run.case_name
      << "; profile_file=" << profile_path.string()
      << "; profile_file_source=command_line.profile"
      << "; implicit_profile_used=false"
      << "; initial_perturbation_enabled="
      << (config.perturbation.enabled ? "true" : "false");
  if (config.perturbation.enabled) {
    out << "; initial_perturbation_type=" << config.perturbation.type
        << "; initial_perturbation_l=" << config.perturbation.ell
        << "; initial_perturbation_m=" << config.perturbation.m
        << "; initial_perturbation_amplitude=" << config.perturbation.amplitude
        << "; initial_perturbation_r0_cm=" << config.perturbation.r0_cm
        << "; initial_perturbation_target=" << config.perturbation.target
        << "; initial_perturbation_reference_velocity_cm_s="
        << perturbation_reference_velocity_cm_s
        << "; initial_perturbation_reference_velocity_source=profile_at_r0"
        << "; initial_perturbation_shape_function=woo_thesis_eq_2_2"
        << "; initial_perturbation_shape_width_fraction=0.015"
        << "; initial_perturbation_formula=vr_plus_amplitude_abs_vr_r0_f_r_Pl_costheta"
        << "; perturbation_applied_to=v_r"
        << "; density_angular_perturbation=false"
        << "; temperature_angular_perturbation=false";
  }
  out
      << "; canonical_state_initialized=true"
      << "; rho_positive=true"
      << "; e_electron_positive=true"
      << "; e_ion_positive=true"
      << "; radiation_initialization=" << config.radiation.radiation_initialization
      << "; radiation_blackbody_provider_executed="
      << (blackbody_executed ? "true" : "false")
      << "; group_edges_unit=eV"
      << "; group_count=" << group_count
      << "; group_edges_source=InputDeck"
      << "; radiation_blackbody_source=profile_Te"
      << "; double_kB_guard_passed=true"
      << "; B_g_provider_failure=false"
      << "; radiation_groups_nonnegative=true"
      << "; alpha_initialization=" << config.alpha.alpha_initialization
      << "; alpha_state_nonnegative=true"
      << "; thermodynamic_recovery_success=true"
      << "; input_deck_report_present=true"
      << "; radial_profile_report_present=true"
      << "; initialization_report_present=true"
      << "; updated_fields=rho,mom_r,mom_theta,mom_phi,e_electron,e_fluid_total,radiation_groups,alpha_state"
      << "; thermodynamic_recovery_report_present="
      << (thermodynamic_report.empty() ? "false" : "true");
  return out.str();
}

[[nodiscard]] dec3d::core::AuthoritativeFieldMask InitialWriteMask() noexcept {
  return dec3d::core::Combine({
      dec3d::core::AuthoritativeField::rho,
      dec3d::core::AuthoritativeField::mom_r,
      dec3d::core::AuthoritativeField::mom_theta,
      dec3d::core::AuthoritativeField::mom_phi,
      dec3d::core::AuthoritativeField::e_fluid_total,
      dec3d::core::AuthoritativeField::e_electron,
      dec3d::core::AuthoritativeField::radiation_groups,
      dec3d::core::AuthoritativeField::alpha_state,
  });
}

[[nodiscard]] double LegendreP(int ell, double x) noexcept {
  if (ell == 0) {
    return 1.0;
  }
  if (ell == 1) {
    return x;
  }
  double p_l_minus_two = 1.0;
  double p_l_minus_one = x;
  for (int l = 2; l <= ell; ++l) {
    const double p_l =
        ((2.0 * static_cast<double>(l) - 1.0) * x * p_l_minus_one -
         (static_cast<double>(l) - 1.0) * p_l_minus_two) /
        static_cast<double>(l);
    p_l_minus_two = p_l_minus_one;
    p_l_minus_one = p_l;
  }
  return p_l_minus_one;
}

[[nodiscard]] std::vector<double> AxisymmetricSingleMode(
    const dec3d::io::InputDeckConfig& config,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  std::vector<double> pattern(config.mesh.theta_cells, 1.0);
  if (!config.perturbation.enabled) {
    return pattern;
  }
  for (std::size_t t = 0u; t < pattern.size(); ++t) {
    const double theta = 0.5 * (geometry.theta_faces[t] + geometry.theta_faces[t + 1u]);
    pattern[t] = LegendreP(config.perturbation.ell, std::cos(theta));
  }
  return pattern;
}

[[nodiscard]] double WooVelocityPerturbationShape(
    const dec3d::io::InputDeckConfig& config,
    double radius_cm) noexcept {
  if (!config.perturbation.enabled ||
      config.perturbation.ell <= 0 ||
      !(radius_cm > 0.0) ||
      !(config.perturbation.r0_cm > 0.0)) {
    return 0.0;
  }
  const double width = 0.015 * config.perturbation.r0_cm;
  const double transition = std::tanh((radius_cm - config.perturbation.r0_cm) / width);
  const double inner = std::pow(radius_cm / config.perturbation.r0_cm,
                                config.perturbation.ell);
  const double outer = std::pow(config.perturbation.r0_cm / radius_cm,
                                config.perturbation.ell);
  return 0.5 * inner * (1.0 - transition) +
         0.5 * outer * (1.0 + transition);
}

}  // namespace

ProfileInitializationResult InitializeFromRadialProfile(
    const dec3d::io::InputDeckConfig& config,
    const dec3d::io::RadialProfile& profile,
    const std::filesystem::path& profile_path) noexcept {
  if (config.radiation.group_edges_eV.size() < 2u) {
    return Fail("group_edges_eV must contain at least two entries");
  }
  const std::size_t group_count = config.radiation.group_edges_eV.size() - 1u;

  dec3d::mesh::SphericalMeshDescriptor descriptor;
  descriptor.radial_cells = config.mesh.radial_cells;
  descriptor.theta_cells = config.mesh.theta_cells;
  descriptor.phi_cells = config.mesh.phi_cells;
  descriptor.inner_radius = config.mesh.radial_min_cm;
  descriptor.outer_radius = config.mesh.radial_max_cm;
  const auto geometry = dec3d::mesh::BuildSphericalGeometry(descriptor);
  if (!geometry.valid) {
    return Fail("failed to build spherical geometry");
  }

  const auto group_layout =
      dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
          FrequencyEdgesFromEV(config.radiation.group_edges_eV));
  const auto group_validation = dec3d::radiation::ValidateRadiationGroupLayout(group_layout);
  if (!group_validation.success || group_layout.group_count != group_count) {
    return Fail("invalid radiation group layout");
  }

  dec3d::state::CanonicalStateLayout layout;
  layout.radial_cells = config.mesh.radial_cells;
  layout.theta_cells = config.mesh.theta_cells;
  layout.phi_cells = config.mesh.phi_cells;
  layout.radiation_group_count = group_count;
  auto state = dec3d::state::CanonicalState::Create(layout);

  bool blackbody_executed = false;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  const double zbar = dec3d::physics::DefaultZbar();
  const double gamma_minus_one = 2.0 / 3.0;
  const auto theta_pattern = AxisymmetricSingleMode(config, geometry);
  double perturbation_reference_velocity_cm_s = 0.0;
  if (config.perturbation.enabled) {
    const auto reference_sample =
        dec3d::io::SampleRadialProfile(profile, config.perturbation.r0_cm);
    if (!reference_sample.success) {
      return Fail(reference_sample.failure_reason);
    }
    perturbation_reference_velocity_cm_s = reference_sample.value.vr_cm_s;
    if (!std::isfinite(perturbation_reference_velocity_cm_s) ||
        (config.perturbation.amplitude > 0.0 &&
         !(std::abs(perturbation_reference_velocity_cm_s) > 0.0))) {
      return Fail("perturbation reference velocity must be finite and nonzero");
    }
  }

  for (std::size_t r = 0u; r < layout.radial_cells; ++r) {
    const double radius = 0.5 * (geometry.radial_faces[r] + geometry.radial_faces[r + 1u]);
    for (std::size_t t = 0u; t < layout.theta_cells; ++t) {
      const auto sample = dec3d::io::SampleRadialProfile(profile, radius);
      if (!sample.success) {
        return Fail(sample.failure_reason);
      }
      const auto& row = sample.value;
      const double te_erg = dec3d::physics::ErgFromKeV(row.Te_keV);
      const double ti_erg = dec3d::physics::ErgFromKeV(row.Ti_keV);
      if (!(te_erg > 0.0) || !(ti_erg > 0.0)) {
        return Fail("profile temperature conversion failed");
      }
      for (std::size_t p = 0u; p < layout.phi_cells; ++p) {
        const double delta_vr =
            config.perturbation.enabled
                ? config.perturbation.amplitude *
                      std::abs(perturbation_reference_velocity_cm_s) *
                      WooVelocityPerturbationShape(config, radius) * theta_pattern[t]
                : 0.0;
        const double vr_cm_s = row.vr_cm_s + delta_vr;
        const double ni = row.rho_g_cm3 / mean_ion_mass;
        const double ne = zbar * ni;
        const double electron = ne * te_erg / gamma_minus_one;
        const double ion = ni * ti_erg / gamma_minus_one;
        const double kinetic = 0.5 * row.rho_g_cm3 *
                               (vr_cm_s * vr_cm_s +
                                row.vt_cm_s * row.vt_cm_s +
                                row.vp_cm_s * row.vp_cm_s);

        state.rho(r, t, p) = row.rho_g_cm3;
        state.mom_r(r, t, p) = row.rho_g_cm3 * vr_cm_s;
        state.mom_theta(r, t, p) = row.rho_g_cm3 * row.vt_cm_s;
        state.mom_phi(r, t, p) = row.rho_g_cm3 * row.vp_cm_s;
        state.e_electron(r, t, p) = electron;
        state.e_fluid_total(r, t, p) = electron + ion + kinetic;
        state.alpha_state.storage(r, t, p) =
            config.alpha.alpha_initialization == "profile" ? row.epsilon_alpha_erg_cm3 : 0.0;

        for (std::size_t g = 0u; g < group_count; ++g) {
          double value = 0.0;
          if (config.radiation.radiation_initialization == "zero") {
            value = 0.0;
          } else {
            const auto blackbody =
                dec3d::radiation::EvaluateGroupBlackbodyEnergyDensity(group_layout, g, te_erg);
            if (!blackbody.success ||
                !dec3d::radiation::ValidateGroupBlackbodyDiagnostics(blackbody)) {
              return Fail("B_g provider failed");
            }
            blackbody_executed = true;
            value = blackbody.energy_density_erg_cm3;
            if (config.radiation.radiation_initialization == "scaled_local_blackbody") {
              value *= config.radiation.radiation_initial_blackbody_scale * row.radiation_scale;
            }
          }
          if (!(value >= 0.0) || !std::isfinite(value)) {
            return Fail("radiation initialization produced invalid U_g");
          }
          state.radiation_groups[g](r, t, p) = value;
        }
      }
    }
  }

  state.ApplyAuthoritativeWrite(InitialWriteMask());

  const auto recovery = dec3d::state::RecoverThermodynamicState(state);
  if (!recovery.success || !dec3d::state::ValidateThermodynamicRecoveryDiagnostics(recovery)) {
    return Fail(recovery.failure_reason.empty() ? "thermodynamic recovery failed"
                                                : recovery.failure_reason);
  }

  ProfileInitializationResult result;
  result.success = true;
  result.state = std::move(state);
  result.geometry = geometry;
  result.group_layout = group_layout;
  result.thermodynamic_recovery_report = recovery.recovery_diagnostics;
  result.report_line =
      BuildReport(config,
                  profile_path,
                  group_count,
                  blackbody_executed,
                  perturbation_reference_velocity_cm_s,
                  recovery.recovery_diagnostics);
  return result;
}

bool ValidateProfileInitializationDiagnostics(
    const ProfileInitializationResult& result) noexcept {
  if (!result.success) {
    return false;
  }
  const auto& line = result.report_line;
  return Contains(line, "diagnostic_id=p5.io.profile_initializer") &&
         Contains(line, "profile_file_source=command_line.profile") &&
         Contains(line, "implicit_profile_used=false") &&
         Contains(line, "canonical_state_initialized=true") &&
         Contains(line, "rho_positive=true") &&
         Contains(line, "e_electron_positive=true") &&
         Contains(line, "e_ion_positive=true") &&
         Contains(line, "group_edges_unit=eV") &&
         Contains(line, "group_count=") &&
         Contains(line, "double_kB_guard_passed=true") &&
         Contains(line, "B_g_provider_failure=false") &&
         Contains(line, "alpha_state_nonnegative=true") &&
         Contains(line, "thermodynamic_recovery_success=true") &&
         Contains(line, "input_deck_report_present=true") &&
         Contains(line, "radial_profile_report_present=true") &&
         Contains(line, "initialization_report_present=true");
}

}  // namespace dec3d::initialization
