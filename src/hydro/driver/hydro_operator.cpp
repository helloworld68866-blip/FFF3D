#include "hydro/driver/hydro_operator.hpp"

#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/ale/radial_ale.hpp"

#include <cmath>
#include <limits>
#include <sstream>

namespace dec3d::hydro {

namespace {

constexpr const char* kHydroImplementationId = "p1.hydro.operator.hllc_static_grid";
constexpr const char* kHydroImplementationIdPpm = "p1.hydro.operator.ppm_static_grid";
constexpr double kHydroExplicitCflPrefactor = 0.5;

struct ExplicitDtSummary {
  bool success{false};
  double hard_cap_dt{0.0};
  std::string reason;
  std::string evidence;
};

void AppendDiagnostic(
    dec3d::core::DiagnosticsPayload& diagnostics,
    std::string code,
    std::string message) {
  diagnostics.entries.push_back({std::move(code), std::move(message)});
}

[[nodiscard]] double CellCenteredRadius(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t radial) noexcept {
  return 0.5 * (geometry.radial_faces[radial] + geometry.radial_faces[radial + 1u]);
}

[[nodiscard]] double CellCenteredTheta(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta) noexcept {
  return 0.5 * (geometry.theta_faces[theta] + geometry.theta_faces[theta + 1u]);
}

[[nodiscard]] std::string DtEvidenceForOptions(const StaticGridHydroOptions& options) {
  if (options.apply_radial_sweep && options.apply_theta_sweep && options.apply_phi_sweep) {
    return "p1.hydro.dt.explicit_cfl.full_direction";
  }
  if (options.apply_radial_sweep && !options.apply_theta_sweep && !options.apply_phi_sweep) {
    return "p1.hydro.dt.explicit_cfl.radial_only";
  }

  std::ostringstream evidence;
  evidence << "p1.hydro.dt.explicit_cfl";
  if (options.apply_radial_sweep) {
    evidence << ".radial";
  }
  if (options.apply_theta_sweep) {
    evidence << ".theta";
  }
  if (options.apply_phi_sweep) {
    evidence << ".phi";
  }
  return evidence.str();
}

[[nodiscard]] ExplicitDtSummary ComputeExplicitDtAdvice(
    const dec3d::state::HydroStateView& hydro_view,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const StaticGridHydroOptions& options) noexcept {
  ExplicitDtSummary summary;
  if (!hydro_view.is_complete()) {
    summary.reason = "hydro estimate_dt requires a complete hydro view";
    summary.evidence = "p1.hydro.dt.explicit_cfl.view_incomplete";
    return summary;
  }
  if (!geometry.is_valid()) {
    summary.reason = "hydro estimate_dt requires valid spherical geometry";
    summary.evidence = "p1.hydro.dt.explicit_cfl.geometry_invalid";
    return summary;
  }
  if (!(options.apply_radial_sweep || options.apply_theta_sweep || options.apply_phi_sweep)) {
    summary.reason = "hydro estimate_dt requires at least one active directional sweep";
    summary.evidence = "p1.hydro.dt.explicit_cfl.no_active_sweeps";
    return summary;
  }

  double limiting_dt = std::numeric_limits<double>::infinity();
  const char* limiting_direction = "none";
  std::size_t limiting_radial = 0u;
  std::size_t limiting_theta = 0u;
  std::size_t limiting_phi = 0u;
  double limiting_width = 0.0;
  double limiting_signal_speed = 0.0;

  for (std::size_t radial = 0; radial < hydro_view.rho->extent_r(); ++radial) {
    const double dr = geometry.radial_faces[radial + 1u] - geometry.radial_faces[radial];
    const double radius = CellCenteredRadius(geometry, radial);
    for (std::size_t theta = 0; theta < hydro_view.rho->extent_theta(); ++theta) {
      const double theta_center = CellCenteredTheta(geometry, theta);
      const double dtheta =
          radius * (geometry.theta_faces[theta + 1u] - geometry.theta_faces[theta]);
      const double dphi_base =
          radius * std::sin(theta_center);

      for (std::size_t phi = 0; phi < hydro_view.rho->extent_phi(); ++phi) {
        const double dphi =
            dphi_base * (geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi]);
        const HydroConservativeState conservative{
            (*hydro_view.rho)(radial, theta, phi),
            (*hydro_view.mom_r)(radial, theta, phi),
            (*hydro_view.mom_theta)(radial, theta, phi),
            (*hydro_view.mom_phi)(radial, theta, phi),
            (*hydro_view.e_fluid_total)(radial, theta, phi),
            hydro_view.chi_e(radial, theta, phi)};
        const auto primitive = RecoverPrimitiveState(conservative);
        if (!primitive.is_physical()) {
          summary.reason = "hydro estimate_dt encountered a non-physical cell state";
          summary.evidence = "p1.hydro.dt.explicit_cfl.non_physical_state";
          return summary;
        }

        const double sound_speed =
            std::sqrt(dec3d::state::HydroIdealGasGamma() * primitive.pressure / primitive.rho);
        const double radial_signal = std::abs(primitive.v_r) + sound_speed;
        const double theta_signal = std::abs(primitive.v_theta) + sound_speed;
        const double phi_signal = std::abs(primitive.v_phi) + sound_speed;

        const auto update_limit = [&](bool enabled,
                                      double width,
                                      double signal_speed,
                                      const char* direction_name) {
          if (!enabled || !(width > 0.0) || !(signal_speed > 0.0) ||
              !std::isfinite(width) || !std::isfinite(signal_speed)) {
            return;
          }

          const double candidate_dt = kHydroExplicitCflPrefactor * width / signal_speed;
          if (candidate_dt < limiting_dt) {
            limiting_dt = candidate_dt;
            limiting_direction = direction_name;
            limiting_radial = radial;
            limiting_theta = theta;
            limiting_phi = phi;
            limiting_width = width;
            limiting_signal_speed = signal_speed;
          }
        };

        update_limit(options.apply_radial_sweep, dr, radial_signal, "radial");
        update_limit(options.apply_theta_sweep, dtheta, theta_signal, "theta");
        update_limit(options.apply_phi_sweep, dphi, phi_signal, "phi");
      }
    }
  }

  if (!(limiting_dt > 0.0) || !std::isfinite(limiting_dt)) {
    summary.reason = "hydro estimate_dt could not derive a finite explicit CFL hard cap";
    summary.evidence = "p1.hydro.dt.explicit_cfl.unresolved";
    return summary;
  }

  std::ostringstream reason;
  reason << "hydro explicit CFL hard cap from active spherical "
         << (options.apply_radial_sweep && options.apply_theta_sweep && options.apply_phi_sweep
                 ? "full-direction"
                 : "directional")
         << " geometry; limiting_direction=" << limiting_direction
         << "; limiting_cell=" << limiting_radial << "," << limiting_theta << "," << limiting_phi
         << "; effective_width=" << limiting_width
         << "; signal_speed=" << limiting_signal_speed
         << "; cfl_prefactor=" << kHydroExplicitCflPrefactor;

  summary.success = true;
  summary.hard_cap_dt = limiting_dt;
  summary.reason = reason.str();
  summary.evidence = DtEvidenceForOptions(options);
  return summary;
}

[[nodiscard]] dec3d::core::StageResult FailedHydroStageResult(
    std::string failure_reason,
    dec3d::core::DiagnosticsPayload diagnostics,
    bool entered_stage,
    std::uint64_t touched_cell_count,
    const char* implementation_id) {
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.stage.failed",
      failure_reason.empty() ? "hydro stage failed" : failure_reason);

  dec3d::core::ExecutionEvidence execution_evidence;
  execution_evidence.entered_stage = entered_stage;
  execution_evidence.implementation_id = implementation_id == nullptr
                                             ? kHydroImplementationId
                                             : implementation_id;
  execution_evidence.touched_cell_count = touched_cell_count;

  return dec3d::core::StageResult::Failed(
      std::move(failure_reason),
      std::move(diagnostics),
      std::move(execution_evidence));
}

}  // namespace

bool HydroOperator::bind(
    const dec3d::core::StageContext& context,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    dec3d::state::CanonicalState& state) noexcept {
  bound_ = false;
  geometry_ = nullptr;
  state_ = nullptr;
  context_ = {};

  if (!context.is_complete() ||
      context.phase_id != dec3d::core::PhaseId::p1 ||
      !geometry.is_valid()) {
    return false;
  }

  context_ = context;
  geometry_ = &geometry;
  state_ = &state;
  bound_ = has_allocated_state();
  return bound_;
}

dec3d::core::DtAdvice HydroOperator::estimate_dt() const noexcept {
  dec3d::core::DtAdvice advice;
  if (!bound_ || geometry_ == nullptr || state_ == nullptr) {
    advice.hard_cap_dt = std::numeric_limits<double>::infinity();
    advice.soft_advice_dt = std::numeric_limits<double>::infinity();
    advice.reason = "hydro scaffold is not bound";
    advice.evidence = "p1.hydro.dt.scaffold.unbound";
    return advice;
  }

  const auto hydro_view = dec3d::state::BuildHydroStateView(*state_);
  const auto explicit_dt = ComputeExplicitDtAdvice(
      hydro_view,
      *geometry_,
      static_grid_options_);
  if (!explicit_dt.success) {
    advice.hard_cap_dt = std::numeric_limits<double>::infinity();
    advice.soft_advice_dt = std::numeric_limits<double>::infinity();
    advice.reason = explicit_dt.reason.empty() ? "hydro estimate_dt failed" : explicit_dt.reason;
    advice.evidence = explicit_dt.evidence.empty() ? "p1.hydro.dt.explicit_cfl.failed"
                                                   : explicit_dt.evidence;
    return advice;
  }

  advice.hard_cap_dt = explicit_dt.hard_cap_dt;
  advice.soft_advice_dt = explicit_dt.hard_cap_dt;
  advice.reason = explicit_dt.reason;
  advice.evidence = explicit_dt.evidence;
  return advice;
}

dec3d::core::StageResult HydroOperator::advance() noexcept {
  if (!bound_ || geometry_ == nullptr || state_ == nullptr) {
    dec3d::core::DiagnosticsPayload diagnostics;
    return FailedHydroStageResult(
        "hydro operator is not bound",
        std::move(diagnostics),
        false,
        0u,
        static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                    : kHydroImplementationId);
  }

  dec3d::state::HydroAlphaViewOptions alpha_view_options;
  alpha_view_options.enabled = static_grid_options_.enable_alpha_hydro_terms;
  auto hydro_view = dec3d::state::BuildHydroWorkView(*state_, alpha_view_options);
  if (!hydro_view.is_complete()) {
    dec3d::core::DiagnosticsPayload diagnostics;
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.view_build_failed",
        hydro_view.failure_reason.empty() ? "hydro work view construction failed"
                                          : hydro_view.failure_reason);
    return FailedHydroStageResult(
        hydro_view.failure_reason.empty() ? "hydro work view construction failed"
                                          : hydro_view.failure_reason,
        std::move(diagnostics),
        true,
        state_->rho.size(),
        static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                    : kHydroImplementationId);
  }

  const std::size_t radial_cells = hydro_view.rho->extent_r();
  if (!(radial_cells > 0u) || !(context_.dt_s > 0.0)) {
    dec3d::core::DiagnosticsPayload diagnostics;
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.precondition_failed",
        "hydro stage requires positive radial extent and positive dt");
    return FailedHydroStageResult(
        "hydro stage requires positive radial extent and positive dt",
        std::move(diagnostics),
        true,
        state_->rho.size(),
        static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                    : kHydroImplementationId);
  }

  std::optional<dec3d::core::MeshUpdateProposal> radial_ale_proposal;
  if (static_grid_options_.apply_radial_ale_flux_correction ||
      static_grid_options_.request_radial_ale_proposal) {
    radial_ale_proposal = dec3d::mesh::BuildRadialAleMeshUpdateProposal(
        state_->rho,
        state_->mom_r,
        *geometry_,
        context_.dt_s);
    if (!radial_ale_proposal->is_complete(geometry_->radial_faces.size())) {
      dec3d::core::DiagnosticsPayload diagnostics;
      AppendDiagnostic(
          diagnostics,
          "p1.hydro.ale.proposal.failed",
          radial_ale_proposal->summary.empty()
              ? "radial ALE proposal construction failed"
              : radial_ale_proposal->summary);
      return FailedHydroStageResult(
          radial_ale_proposal->summary.empty()
              ? "radial ALE proposal construction failed"
              : radial_ale_proposal->summary,
          std::move(diagnostics),
          true,
          state_->rho.size(),
          static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                      : kHydroImplementationId);
    }
  }

  const auto static_grid = AdvanceStaticGridHydro(
      hydro_view,
      *geometry_,
      context_.dt_s,
      static_grid_options_,
      radial_ale_proposal.has_value() ? &*radial_ale_proposal : nullptr);
  if (!static_grid.success) {
    dec3d::core::DiagnosticsPayload diagnostics = static_grid.diagnostics;
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.static_grid_failed",
        static_grid.failure_reason.empty() ? "static-grid hydro failed"
                                           : static_grid.failure_reason);
    return FailedHydroStageResult(
        static_grid.failure_reason.empty() ? "static-grid hydro failed"
                                           : static_grid.failure_reason,
        std::move(diagnostics),
        true,
        state_->rho.size(),
        static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                    : kHydroImplementationId);
  }

  auto updated_fields = dec3d::state::BuildHydroAuthorizedWriteMask();
  if (static_grid_options_.enable_radiation_hydro_terms) {
    updated_fields |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::radiation_groups);
  }
  if (static_grid_options_.enable_alpha_hydro_terms) {
    updated_fields |= dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state);
  }

  const auto writeback = dec3d::state::CommitHydroWriteback(
      *state_,
      hydro_view,
      updated_fields);
  if (!writeback.success) {
    dec3d::core::DiagnosticsPayload diagnostics = writeback.diagnostics;
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.writeback_failed",
        writeback.failure_reason.empty() ? "hydro writeback failed"
                                         : writeback.failure_reason);
    return FailedHydroStageResult(
        writeback.failure_reason.empty() ? "hydro writeback failed"
                                         : writeback.failure_reason,
        std::move(diagnostics),
        true,
        state_->rho.size(),
        static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                    : kHydroImplementationId);
  }

  dec3d::core::DiagnosticsPayload diagnostics;
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.stage.executed",
      "hydro scaffold executed through the runtime-owned hydro path");
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.stage.view_writeback",
      "hydro scaffold executed through hydro-local chi_e recovery and authorized electron writeback");
  AppendDiagnostic(
      diagnostics,
      "p1.hydro.stage.hllc_static_grid",
      "hydro stage executed static-grid HLLC flux assembly");
  if (static_grid_options_.use_ppm_reconstruction) {
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.ppm_reconstruction",
        "hydro stage executed direction-local characteristic traced-interface PPM reconstruction");
  }
  if (static_grid_options_.use_macro_zoning) {
    AppendDiagnostic(
        diagnostics,
        "p1.hydro.stage.macro_zoning",
        "hydro stage executed single-rank macro-zoned coarse/fine angular coupling");
  }
  if (static_grid_options_.enable_alpha_hydro_terms) {
    AppendDiagnostic(
        diagnostics,
        "p4.alpha.hydro_terms.enabled",
        "hydro stage transported alpha pressure scalar P_alpha^(3/5)");
    AppendDiagnostic(
        diagnostics,
        "p4.alpha.hydro_terms.write_set",
        "h_stage_updated_fields_contains_alpha_state=true; alpha_state_epoch=post_H_committed; alpha_hydro_terms_share_hydro_scalar_bundle=true");
  }
  diagnostics.entries.insert(
      diagnostics.entries.end(),
      static_grid.diagnostics.entries.begin(),
      static_grid.diagnostics.entries.end());
  diagnostics.entries.insert(
      diagnostics.entries.end(),
      writeback.diagnostics.entries.begin(),
      writeback.diagnostics.entries.end());

  dec3d::core::ExecutionEvidence execution_evidence;
  execution_evidence.entered_stage = true;
  execution_evidence.implementation_id =
      static_grid_options_.use_ppm_reconstruction ? kHydroImplementationIdPpm
                                                  : kHydroImplementationId;
  execution_evidence.touched_cell_count = state_->rho.size();

  auto result = dec3d::core::StageResult::Successful(
      writeback.updated_fields,
      diagnostics,
      execution_evidence);

  if (static_grid_options_.request_radial_ale_proposal && radial_ale_proposal.has_value()) {
    AppendDiagnostic(
        result.diagnostics,
        "p1.hydro.ale.proposal.generated",
        radial_ale_proposal->summary);
    result.mesh_update_proposal = *radial_ale_proposal;
  }

  return result;
}

bool HydroOperator::has_allocated_state() const noexcept {
  return state_ != nullptr &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::rho) &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_r) &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_theta) &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::mom_phi) &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_fluid_total) &&
         state_->HasAuthoritativeStorage(dec3d::core::AuthoritativeField::e_electron);
}

}  // namespace dec3d::hydro
