#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/hydro_boundary_ghosts.hpp"
#include "hydro/riemann/hllc_solver.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::hydro {

struct DirectionalPrimitiveState {
  double rho{0.0};
  double v_n{0.0};
  double v_t1{0.0};
  double v_t2{0.0};
  double pressure{0.0};
  double chi_e{0.0};
  double alpha_chi{0.0};
  std::vector<double> radiation_chi;

  [[nodiscard]] bool is_physical() const noexcept;
};

struct PpmSmoothnessProbe {
  int component{-1};
  bool sign_reversal{false};
  double q_im1{0.0};
  double q_i{0.0};
  double q_ip1{0.0};
  double delta_left{0.0};
  double delta_right{0.0};
  double primary_slope{0.0};
  double reversed_slope{0.0};
  double tolerance{0.0};
};

struct PpmInterfaceState {
  DirectionalPrimitiveState left;
  DirectionalPrimitiveState right;
  bool downgraded_to_first_order{false};
  DirectionalPrimitiveState left_edge_before_trace;
  DirectionalPrimitiveState right_edge_before_trace;
  bool edge_state_diagnostics_available{false};
  bool left_profile_troubled{false};
  bool right_profile_troubled{false};
  bool left_trace_failed{false};
  bool right_trace_failed{false};
  int left_troubled_component{-1};
  int right_troubled_component{-1};
  double left_trouble_q_im1{0.0};
  double left_trouble_q_i{0.0};
  double left_trouble_q_ip1{0.0};
  double left_trouble_total_variation{0.0};
  double left_trouble_end_to_end{0.0};
  double right_trouble_q_im1{0.0};
  double right_trouble_q_i{0.0};
  double right_trouble_q_ip1{0.0};
  double right_trouble_total_variation{0.0};
  double right_trouble_end_to_end{0.0};
  PpmSmoothnessProbe left_component1_smoothness;
  PpmSmoothnessProbe right_component1_smoothness;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct PpmLineReconstructionResult {
  bool success{false};
  HydroDirection direction{HydroDirection::radial};
  std::size_t ghost_layers_consumed{0};
  std::size_t interface_count{0};
  std::size_t downgraded_interface_count{0};
  std::vector<PpmInterfaceState> interfaces;
  dec3d::core::DiagnosticsPayload diagnostics;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete(std::size_t expected_interface_count) const noexcept;
};

[[nodiscard]] DirectionalPrimitiveState ToDirectionalPrimitive(
    const HydroConservativeState& canonical_state,
    HydroDirection direction) noexcept;

[[nodiscard]] HydroConservativeState MakeDirectionalConservativeState(
    const DirectionalPrimitiveState& directional_primitive) noexcept;

[[nodiscard]] PpmLineReconstructionResult ReconstructPpmLine(
    const std::vector<HydroConservativeState>& ghosted_line,
    HydroDirection direction,
    std::size_t interior_cells,
    std::size_t ghost_layers,
    double dt_s,
    const std::vector<double>& effective_cell_widths,
    const std::vector<double>* moving_face_speeds = nullptr) noexcept;

void AccumulateShellWidePpmFallbackMask(
    const PpmLineReconstructionResult& line_reconstruction,
    std::vector<bool>& shell_face_mask) noexcept;

[[nodiscard]] std::vector<bool> BuildShellWidePpmFallbackMask(
    const std::vector<PpmLineReconstructionResult>& angular_line_reconstructions,
    std::size_t interface_count) noexcept;

}  // namespace dec3d::hydro
