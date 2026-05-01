#include "state/hydro_state/hydro_view.hpp"

#include <cmath>
#include <sstream>

namespace dec3d::state {

namespace {

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] bool MatchesShape(
    const dec3d::core::Array3D<double>& lhs,
    const dec3d::core::Array3D<double>& rhs) noexcept {
  return lhs.extent_r() == rhs.extent_r() &&
         lhs.extent_theta() == rhs.extent_theta() &&
         lhs.extent_phi() == rhs.extent_phi();
}

[[nodiscard]] double RadiationPressureScalarFromEnergyDensity(double ug) noexcept {
  return ug < 0.0 ? std::nan("") : std::pow(ug / 3.0, 3.0 / 4.0);
}

[[nodiscard]] double RadiationEnergyDensityFromPressureScalar(double chi_rad) noexcept {
  return chi_rad < 0.0 ? std::nan("") : 3.0 * std::pow(chi_rad, 4.0 / 3.0);
}

[[nodiscard]] double AlphaPressureScalarFromEnergyDensity(double epsilon_alpha) noexcept {
  const double p_alpha = (2.0 / 3.0) * epsilon_alpha;
  return p_alpha < 0.0 || !std::isfinite(p_alpha)
             ? std::nan("")
             : std::pow(p_alpha, 3.0 / 5.0);
}

[[nodiscard]] double AlphaEnergyDensityFromPressureScalar(double chi_alpha) noexcept {
  return chi_alpha < 0.0 || !std::isfinite(chi_alpha)
             ? std::nan("")
             : 1.5 * std::pow(chi_alpha, 5.0 / 3.0);
}

void PopulateChiEFromElectronChannel(HydroStateView& view) noexcept {
  if (view.rho == nullptr || view.e_electron == nullptr) {
    view.failure_reason = "hydro view is missing rho or electron storage";
    view.report_line = "hydro_view_complete=false; failure_reason=missing_hydro_storage";
    return;
  }

  view.chi_e = dec3d::core::Array3D<double>(
      view.rho->extent_r(),
      view.rho->extent_theta(),
      view.rho->extent_phi(),
      0.0);

  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const double electron_pressure =
            ElectronPressureFromElectronEnergyDensity((*view.e_electron)(radial, theta, phi));
        if (electron_pressure < 0.0) {
          view.failure_reason = "electron pressure must remain non-negative for chi_e recovery";
          std::ostringstream report;
          report << "hydro_view_complete=false"
                 << "; failure_reason=" << view.failure_reason
                 << "; failing_cell=" << radial << "," << theta << "," << phi;
          view.report_line = report.str();
          return;
        }

        view.chi_e(radial, theta, phi) = ChiEFromElectronPressure(electron_pressure);
      }
    }
  }

  std::ostringstream report;
  report << "hydro_view_complete=true"
         << "; operator_local_chi_e=true"
         << "; transactional_work_copy=" << (view.transactional_work_copy ? "true" : "false")
         << "; cell_count=" << view.rho->size();
  view.report_line = report.str();
}

void PopulateAlphaChiFromAlphaState(
    HydroStateView& view,
    const dec3d::state::CanonicalState& state,
    HydroAlphaViewOptions options) noexcept {
  if (view.rho == nullptr) {
    view.failure_reason = "hydro view is missing rho storage for alpha scalar recovery";
    view.report_line = "hydro_view_complete=false; failure_reason=missing_hydro_storage";
    return;
  }
  if (!std::isfinite(options.alpha_energy_floor_erg_cm3) ||
      options.alpha_energy_floor_erg_cm3 < 0.0) {
    view.failure_reason = "alpha hydro energy floor must be finite and non-negative";
    view.report_line = "hydro_view_complete=false; failure_reason=invalid_alpha_energy_floor";
    return;
  }
  if (!state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state) ||
      !MatchesShape(*view.rho, state.alpha_state.storage)) {
    view.failure_reason = "alpha_state storage does not match hydro view shape";
    view.report_line = "hydro_view_complete=false; failure_reason=alpha_state_shape_mismatch";
    return;
  }

  view.alpha_chi = dec3d::core::Array3D<double>(
      view.rho->extent_r(),
      view.rho->extent_theta(),
      view.rho->extent_phi(),
      0.0);

  for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
    for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
      for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
        const double epsilon = state.alpha_state.storage(radial, theta, phi);
        if (!std::isfinite(epsilon) || epsilon < options.alpha_energy_floor_erg_cm3) {
          view.failure_reason = "alpha_state is non-finite or below alpha hydro energy floor";
          std::ostringstream report;
          report << "hydro_view_complete=false"
                 << "; failure_reason=" << view.failure_reason
                 << "; failing_cell=" << radial << "," << theta << "," << phi;
          view.report_line = report.str();
          return;
        }
        const double chi = AlphaPressureScalarFromEnergyDensity(epsilon);
        if (!(std::isfinite(chi) && chi >= 0.0)) {
          view.failure_reason = "alpha pressure scalar must remain finite and non-negative";
          std::ostringstream report;
          report << "hydro_view_complete=false"
                 << "; failure_reason=" << view.failure_reason
                 << "; failing_cell=" << radial << "," << theta << "," << phi;
          view.report_line = report.str();
          return;
        }
        view.alpha_chi(radial, theta, phi) = chi;
      }
    }
  }

  view.operator_local_alpha_chi = true;
}

void PopulateRadiationChiFromRadiationGroups(
    HydroStateView& view,
    const std::vector<dec3d::core::Array3D<double>>& radiation_groups) noexcept {
  if (view.rho == nullptr) {
    view.failure_reason = "hydro view is missing rho storage for radiation scalar recovery";
    view.report_line = "hydro_view_complete=false; failure_reason=missing_hydro_storage";
    return;
  }

  view.radiation_chi.clear();
  view.radiation_chi.reserve(radiation_groups.size());

  for (std::size_t group = 0; group < radiation_groups.size(); ++group) {
    if (!MatchesShape(*view.rho, radiation_groups[group])) {
      view.failure_reason = "radiation group storage does not match hydro view shape";
      std::ostringstream report;
      report << "hydro_view_complete=false"
             << "; failure_reason=" << view.failure_reason
             << "; failing_radiation_group=" << group;
      view.report_line = report.str();
      return;
    }

    dec3d::core::Array3D<double> group_chi(
        view.rho->extent_r(),
        view.rho->extent_theta(),
        view.rho->extent_phi(),
        0.0);
    for (std::size_t radial = 0; radial < view.rho->extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < view.rho->extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < view.rho->extent_phi(); ++phi) {
          const double ug = radiation_groups[group](radial, theta, phi);
          if (!(ug >= 0.0) || !std::isfinite(ug)) {
            view.failure_reason = "radiation energy density must remain finite and non-negative";
            std::ostringstream report;
            report << "hydro_view_complete=false"
                   << "; failure_reason=" << view.failure_reason
                   << "; failing_radiation_group=" << group
                   << "; failing_cell=" << radial << "," << theta << "," << phi;
            view.report_line = report.str();
            return;
          }

          group_chi(radial, theta, phi) = RadiationPressureScalarFromEnergyDensity(ug);
        }
      }
    }

    view.radiation_chi.push_back(std::move(group_chi));
  }
}

void AppendRadiationHydroTermsReport(HydroStateView& view) {
  if (!view.failure_reason.empty()) {
    return;
  }

  std::ostringstream report;
  report << view.report_line
         << "; operator_local_radiation_chi=true"
         << "; radiation_group_count=" << view.radiation_chi.size()
         << "; advected_radiation_scalar=P_g_power_3_over_4"
         << "; passive_Ug_advection=false";
  view.report_line = report.str();
}

}  // namespace

dec3d::core::AuthoritativeFieldMask BuildHydroAuthorizedWriteMask() noexcept {
  using dec3d::core::AuthoritativeField;
  using dec3d::core::Combine;

  return Combine({
      AuthoritativeField::rho,
      AuthoritativeField::mom_r,
      AuthoritativeField::mom_theta,
      AuthoritativeField::mom_phi,
      AuthoritativeField::e_fluid_total,
      AuthoritativeField::e_electron});
}

dec3d::core::AuthoritativeFieldMask BuildHydroElectronWritebackMask() noexcept {
  return dec3d::core::ToMask(dec3d::core::AuthoritativeField::e_electron);
}

double ElectronPressureFromElectronEnergyDensity(double electron_energy_density) noexcept {
  return (HydroIdealGasGamma() - 1.0) * electron_energy_density;
}

double ElectronEnergyDensityFromPressure(double electron_pressure) noexcept {
  return electron_pressure / (HydroIdealGasGamma() - 1.0);
}

double ChiEFromElectronPressure(double electron_pressure) noexcept {
  return electron_pressure < 0.0 ? std::nan("") : std::pow(electron_pressure, 3.0 / 5.0);
}

double ElectronPressureFromChiE(double chi_e) noexcept {
  return chi_e < 0.0 ? std::nan("") : std::pow(chi_e, 5.0 / 3.0);
}

bool HydroStateView::is_complete() const noexcept {
  if (rho == nullptr ||
      mom_r == nullptr ||
      mom_theta == nullptr ||
      mom_phi == nullptr ||
      e_fluid_total == nullptr ||
      e_electron == nullptr ||
      !operator_local_chi_e ||
      !operator_local_radiation_chi ||
      !MatchesShape(*rho, chi_e) ||
      (operator_local_alpha_chi && !MatchesShape(*rho, alpha_chi)) ||
      !MatchesShape(*rho, *mom_r) ||
      !MatchesShape(*rho, *mom_theta) ||
      !MatchesShape(*rho, *mom_phi) ||
      !MatchesShape(*rho, *e_fluid_total) ||
      !MatchesShape(*rho, *e_electron) ||
      report_line.empty()) {
    return false;
  }

  for (const auto& group_chi : radiation_chi) {
    if (!MatchesShape(*rho, group_chi)) {
      return false;
    }
  }

  return true;
}

bool HydroWritebackResult::is_complete() const noexcept {
  return !report_line.empty() && diagnostics.has_entries();
}

HydroStateView BuildHydroStateView(CanonicalState& state) noexcept {
  return BuildHydroStateView(state, HydroAlphaViewOptions{});
}

HydroStateView BuildHydroStateView(
    CanonicalState& state,
    HydroAlphaViewOptions alpha_options) noexcept {
  HydroStateView view;
  view.rho = &state.rho;
  view.mom_r = &state.mom_r;
  view.mom_theta = &state.mom_theta;
  view.mom_phi = &state.mom_phi;
  view.e_fluid_total = &state.e_fluid_total;
  view.e_electron = &state.e_electron;

  if (!(state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::rho) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_r) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_theta) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_phi) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_fluid_total) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_electron))) {
    view.failure_reason = "canonical state is missing hydro authoritative storage";
    view.report_line = "hydro_view_complete=false; failure_reason=missing_hydro_authoritative_storage";
    return view;
  }

  PopulateChiEFromElectronChannel(view);
  if (alpha_options.enabled && view.failure_reason.empty()) {
    PopulateAlphaChiFromAlphaState(view, state, alpha_options);
  }
  if (view.failure_reason.empty()) {
    PopulateRadiationChiFromRadiationGroups(view, state.radiation_groups);
  }
  AppendRadiationHydroTermsReport(view);
  return view;
}

HydroStateView BuildHydroWorkView(const CanonicalState& state) noexcept {
  return BuildHydroWorkView(state, HydroAlphaViewOptions{});
}

HydroStateView BuildHydroWorkView(
    const CanonicalState& state,
    HydroAlphaViewOptions alpha_options) noexcept {
  HydroStateView view;
  view.transactional_work_copy = true;

  if (!(state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::rho) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_r) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_theta) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_phi) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_fluid_total) &&
        state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_electron))) {
    view.failure_reason = "canonical state is missing hydro authoritative storage";
    view.report_line = "hydro_view_complete=false; failure_reason=missing_hydro_authoritative_storage";
    return view;
  }

  view.owned_rho = std::make_unique<dec3d::core::Array3D<double>>(state.rho);
  view.owned_mom_r = std::make_unique<dec3d::core::Array3D<double>>(state.mom_r);
  view.owned_mom_theta = std::make_unique<dec3d::core::Array3D<double>>(state.mom_theta);
  view.owned_mom_phi = std::make_unique<dec3d::core::Array3D<double>>(state.mom_phi);
  view.owned_e_fluid_total = std::make_unique<dec3d::core::Array3D<double>>(state.e_fluid_total);
  view.owned_e_electron = std::make_unique<dec3d::core::Array3D<double>>(state.e_electron);
  view.rho = view.owned_rho.get();
  view.mom_r = view.owned_mom_r.get();
  view.mom_theta = view.owned_mom_theta.get();
  view.mom_phi = view.owned_mom_phi.get();
  view.e_fluid_total = view.owned_e_fluid_total.get();
  view.e_electron = view.owned_e_electron.get();

  PopulateChiEFromElectronChannel(view);
  if (alpha_options.enabled && view.failure_reason.empty()) {
    PopulateAlphaChiFromAlphaState(view, state, alpha_options);
  }
  if (view.failure_reason.empty()) {
    PopulateRadiationChiFromRadiationGroups(view, state.radiation_groups);
  }
  AppendRadiationHydroTermsReport(view);
  return view;
}

HydroWritebackResult CommitHydroWriteback(
    CanonicalState& state,
    const HydroStateView& hydro_view,
    dec3d::core::AuthoritativeFieldMask updated_fields) noexcept {
  HydroWritebackResult result;

  if (!hydro_view.is_complete()) {
    result.failure_reason = "hydro state view is incomplete";
  } else if (!MatchesShape(state.rho, hydro_view.chi_e)) {
    result.failure_reason = "hydro chi_e storage does not match canonical state shape";
  } else if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::alpha_state) &&
             (!hydro_view.operator_local_alpha_chi ||
              !MatchesShape(state.alpha_state.storage, hydro_view.alpha_chi))) {
    result.failure_reason = "hydro alpha scalar storage does not match canonical alpha_state";
  } else if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::radiation_groups) &&
             hydro_view.radiation_chi.size() != state.radiation_groups.size()) {
    result.failure_reason = "hydro radiation scalar bundle size does not match canonical radiation groups";
  } else {
    CanonicalState candidate = state;

    if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::rho)) {
      candidate.rho = *hydro_view.rho;
    }
    if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_r)) {
      candidate.mom_r = *hydro_view.mom_r;
    }
    if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_theta)) {
      candidate.mom_theta = *hydro_view.mom_theta;
    }
    if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_phi)) {
      candidate.mom_phi = *hydro_view.mom_phi;
    }
    if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::e_fluid_total)) {
      candidate.e_fluid_total = *hydro_view.e_fluid_total;
    }

    for (std::size_t radial = 0; radial < state.rho.extent_r(); ++radial) {
      for (std::size_t theta = 0; theta < state.rho.extent_theta(); ++theta) {
        for (std::size_t phi = 0; phi < state.rho.extent_phi(); ++phi) {
          const double chi_e = hydro_view.chi_e(radial, theta, phi);
          if (chi_e < 0.0) {
            result.failure_reason = "chi_e must remain non-negative during hydro writeback";
            std::ostringstream report;
            report << "hydro_writeback_success=false"
                   << "; failure_reason=" << result.failure_reason
                   << "; failing_cell=" << radial << "," << theta << "," << phi;
            result.report_line = report.str();
            AppendDiagnostic(
                result.diagnostics,
                "p1.hydro.writeback.failed",
                result.failure_reason);
            return result;
          }

          if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::e_electron)) {
            candidate.e_electron(radial, theta, phi) =
                ElectronEnergyDensityFromPressure(ElectronPressureFromChiE(chi_e));
          }

          if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::alpha_state)) {
            const double chi_alpha = hydro_view.alpha_chi(radial, theta, phi);
            const double epsilon_alpha = AlphaEnergyDensityFromPressureScalar(chi_alpha);
            if (!(std::isfinite(epsilon_alpha) && epsilon_alpha >= 0.0)) {
              result.failure_reason =
                  "alpha pressure scalar must recover finite non-negative alpha_state during hydro writeback";
              std::ostringstream report;
              report << "hydro_writeback_success=false"
                     << "; failure_reason=" << result.failure_reason
                     << "; failing_cell=" << radial << "," << theta << "," << phi;
              result.report_line = report.str();
              AppendDiagnostic(
                  result.diagnostics,
                  "p1.hydro.writeback.failed",
                  result.failure_reason);
              return result;
            }
            candidate.alpha_state.storage(radial, theta, phi) = epsilon_alpha;
          }

          if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::radiation_groups)) {
            for (std::size_t group = 0; group < hydro_view.radiation_chi.size(); ++group) {
              const double chi_rad = hydro_view.radiation_chi[group](radial, theta, phi);
              if (chi_rad < 0.0 || !std::isfinite(chi_rad)) {
                result.failure_reason =
                    "radiation pressure scalar must remain finite and non-negative during hydro writeback";
                std::ostringstream report;
                report << "hydro_writeback_success=false"
                       << "; failure_reason=" << result.failure_reason
                       << "; failing_radiation_group=" << group
                       << "; failing_cell=" << radial << "," << theta << "," << phi;
                result.report_line = report.str();
                AppendDiagnostic(
                    result.diagnostics,
                    "p1.hydro.writeback.failed",
                    result.failure_reason);
                return result;
              }

              candidate.radiation_groups[group](radial, theta, phi) =
                  RadiationEnergyDensityFromPressureScalar(chi_rad);
            }
          }
        }
      }
    }

    const auto feasibility = ValidateThermodynamicFeasibility(
        candidate,
        ThermodynamicFeasibilityThresholds{0.0, 0.0});
    if (!feasibility.success) {
      result.failure_reason = feasibility.failure_reason.empty()
                                  ? "hydro writeback failed thermodynamic feasibility validation"
                                  : feasibility.failure_reason;
    } else {
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::rho)) {
        state.rho = std::move(candidate.rho);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_r)) {
        state.mom_r = std::move(candidate.mom_r);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_theta)) {
        state.mom_theta = std::move(candidate.mom_theta);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::mom_phi)) {
        state.mom_phi = std::move(candidate.mom_phi);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::e_fluid_total)) {
        state.e_fluid_total = std::move(candidate.e_fluid_total);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::e_electron)) {
        state.e_electron = std::move(candidate.e_electron);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::alpha_state)) {
        state.alpha_state = std::move(candidate.alpha_state);
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::radiation_groups)) {
        state.radiation_groups = std::move(candidate.radiation_groups);
      }

      result.success = true;
      result.electron_channel_updated =
          dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::e_electron);
      result.updated_fields = updated_fields;
      state.ApplyAuthoritativeWrite(result.updated_fields);
      AppendDiagnostic(
          result.diagnostics,
          "p1.hydro.writeback.electron_channel",
          "hydro electron channel updated through chi_e -> Pe -> E_electron");
      AppendDiagnostic(
          result.diagnostics,
          "p1.hydro.writeback.thermodynamic_feasibility",
          "hydro writeback preserved canonical thermodynamic feasibility");
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::radiation_groups)) {
        AppendDiagnostic(
            result.diagnostics,
            "p3.radiation.hydro_terms.writeback",
            "hydro writeback recovered radiation_groups from advected P_g^(3/4) scalar bundle");
      }
      if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::alpha_state)) {
        AppendDiagnostic(
            result.diagnostics,
            "p4.alpha.hydro_terms.writeback",
            "hydro writeback recovered alpha_state from advected P_alpha^(3/5) scalar");
      }
    }
  }

  if (!result.success && result.diagnostics.entries.empty()) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.writeback.failed",
        result.failure_reason.empty() ? "hydro writeback failed" : result.failure_reason);
  }

  std::ostringstream report;
  report << "hydro_writeback_success=" << (result.success ? "true" : "false")
         << "; electron_channel_updated=" << (result.electron_channel_updated ? "true" : "false")
         << "; updated_fields=" << result.updated_fields
         << "; radiation_group_count=" << hydro_view.radiation_chi.size()
         << "; advected_radiation_scalar=P_g_power_3_over_4"
         << "; passive_Ug_advection=false";
  if (dec3d::core::MaskContains(updated_fields, dec3d::core::AuthoritativeField::alpha_state)) {
    report << "; alpha_hydro_terms_enabled=true"
           << "; advected_alpha_scalar=P_alpha_power_3_over_5"
           << "; passive_epsilon_alpha_advection=false"
           << "; h_stage_updated_fields_contains_alpha_state=true"
           << "; alpha_state_epoch=post_H_committed"
           << "; alpha_hydro_terms_share_hydro_scalar_bundle=true";
  }
  if (!result.failure_reason.empty()) {
    report << "; failure_reason=" << result.failure_reason;
  }
  result.report_line = report.str();

  return result;
}

}  // namespace dec3d::state
