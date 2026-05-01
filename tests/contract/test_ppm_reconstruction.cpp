#include "hydro/reconstruction/ppm_reconstruction.hpp"
#include "test_assert.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

[[nodiscard]] dec3d::hydro::HydroConservativeState MakeState(
    double rho,
    double v_r,
    double pressure,
    double chi_e) noexcept {
  return dec3d::hydro::MakeConservativeState(
      dec3d::hydro::HydroPrimitiveState{
          rho,
          v_r,
          0.0,
          0.0,
          pressure,
          chi_e});
}

}  // namespace

int main() {
  try {
    using dec3d::hydro::HydroDirection;
    using dec3d::hydro::ReconstructPpmLine;

    const std::vector<dec3d::hydro::HydroConservativeState> monotone_line{
        MakeState(0.80, 0.00, 0.90, 0.50),
        MakeState(0.90, 0.00, 0.95, 0.55),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.10, 0.05, 1.05, 0.65),
        MakeState(1.20, 0.10, 1.10, 0.70),
        MakeState(1.30, 0.15, 1.15, 0.75),
        MakeState(1.40, 0.20, 1.20, 0.80),
        MakeState(1.50, 0.25, 1.25, 0.85)
    };

    const auto insufficient = ReconstructPpmLine(
        monotone_line,
        HydroDirection::radial,
        2u,
        2u,
        1.0e-4,
        std::vector<double>(monotone_line.size(), 0.1));
    DEC3D_CHECK(!insufficient.success);
    DEC3D_CHECK(!insufficient.failure_reason.empty());

    const std::vector<double> effective_widths(monotone_line.size(), 0.1);
    const auto base = ReconstructPpmLine(
        monotone_line,
        HydroDirection::radial,
        2u,
        3u,
        1.0e-4,
        effective_widths);
    DEC3D_CHECK(base.is_complete(3u));
    DEC3D_CHECK_EQ(base.ghost_layers_consumed, 3u);
    bool saw_distinct_traced_states = false;
    for (const auto& interface_state : base.interfaces) {
      if (!interface_state.downgraded_to_first_order &&
          (std::abs(interface_state.left.rho - interface_state.right.rho) > 1.0e-8 ||
           std::abs(interface_state.left.v_n - interface_state.right.v_n) > 1.0e-8 ||
           std::abs(interface_state.left.pressure - interface_state.right.pressure) > 1.0e-8 ||
           std::abs(interface_state.left.chi_e - interface_state.right.chi_e) > 1.0e-8)) {
        saw_distinct_traced_states = true;
        break;
      }
    }
    DEC3D_CHECK(saw_distinct_traced_states);

    const auto traced_with_small_dt = ReconstructPpmLine(
        monotone_line,
        HydroDirection::radial,
        2u,
        3u,
        1.0e-5,
        effective_widths);
    const auto traced_with_large_dt = ReconstructPpmLine(
        monotone_line,
        HydroDirection::radial,
        2u,
        3u,
        4.0e-4,
        effective_widths);
    DEC3D_CHECK(traced_with_small_dt.is_complete(3u));
    DEC3D_CHECK(traced_with_large_dt.is_complete(3u));
    double dt_sensitive_shift = 0.0;
    for (std::size_t interface_index = 0; interface_index < traced_with_small_dt.interfaces.size(); ++interface_index) {
      const auto& small_dt_state = traced_with_small_dt.interfaces[interface_index];
      const auto& large_dt_state = traced_with_large_dt.interfaces[interface_index];
      dt_sensitive_shift += std::abs(small_dt_state.left.rho - large_dt_state.left.rho);
      dt_sensitive_shift += std::abs(small_dt_state.left.v_n - large_dt_state.left.v_n);
      dt_sensitive_shift += std::abs(small_dt_state.left.pressure - large_dt_state.left.pressure);
      dt_sensitive_shift += std::abs(small_dt_state.right.rho - large_dt_state.right.rho);
      dt_sensitive_shift += std::abs(small_dt_state.right.v_n - large_dt_state.right.v_n);
      dt_sensitive_shift += std::abs(small_dt_state.right.pressure - large_dt_state.right.pressure);
    }
    DEC3D_CHECK(dt_sensitive_shift > 1.0e-10);

    const std::vector<dec3d::hydro::HydroConservativeState> cold_pressure_only_line{
        MakeState(1.0, -1.0, 1.0e-6 - 7.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 - 5.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 - 3.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 - 1.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 + 1.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 + 3.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 + 5.0e-16, 0.50),
        MakeState(1.0, -1.0, 1.0e-6 + 7.0e-16, 0.50)
    };
    const auto cold_pressure_only = ReconstructPpmLine(
        cold_pressure_only_line,
        HydroDirection::radial,
        2u,
        3u,
        1.0e-6,
        std::vector<double>(cold_pressure_only_line.size(), 0.1));
    DEC3D_CHECK(cold_pressure_only.is_complete(3u));
    for (const auto& interface_state : cold_pressure_only.interfaces) {
      DEC3D_CHECK(interface_state.edge_state_diagnostics_available);
      DEC3D_CHECK(std::abs(interface_state.left_edge_before_trace.rho - 1.0) <= 1.0e-15);
      DEC3D_CHECK(std::abs(interface_state.right_edge_before_trace.rho - 1.0) <= 1.0e-15);
      DEC3D_CHECK(std::abs(interface_state.left_edge_before_trace.v_n + 1.0) <= 1.0e-15);
      DEC3D_CHECK(std::abs(interface_state.right_edge_before_trace.v_n + 1.0) <= 1.0e-15);
      DEC3D_CHECK(std::abs(interface_state.left_edge_before_trace.chi_e - 0.50) <= 1.0e-15);
      DEC3D_CHECK(std::abs(interface_state.right_edge_before_trace.chi_e - 0.50) <= 1.0e-15);
    }

    auto perturbed_line = monotone_line;
    perturbed_line[5] = MakeState(1.28, 0.11, 1.11, 0.74);
    const auto perturbed = ReconstructPpmLine(
        perturbed_line,
        HydroDirection::radial,
        2u,
        3u,
        1.0e-4,
        effective_widths);
    DEC3D_CHECK(perturbed.is_complete(3u));

    double total_interface_shift = 0.0;
    for (std::size_t interface_index = 0; interface_index < base.interfaces.size(); ++interface_index) {
      total_interface_shift +=
          std::abs(base.interfaces[interface_index].left.rho -
                   perturbed.interfaces[interface_index].left.rho) +
          std::abs(base.interfaces[interface_index].left.v_n -
                   perturbed.interfaces[interface_index].left.v_n) +
          std::abs(base.interfaces[interface_index].right.rho -
                   perturbed.interfaces[interface_index].right.rho) +
          std::abs(base.interfaces[interface_index].right.v_n -
                   perturbed.interfaces[interface_index].right.v_n);
      if (base.interfaces[interface_index].downgraded_to_first_order !=
          perturbed.interfaces[interface_index].downgraded_to_first_order) {
        total_interface_shift += 1.0;
      }
    }
    DEC3D_CHECK(total_interface_shift > 1.0e-8);

    const std::vector<dec3d::hydro::HydroConservativeState> near_monotone_outer_roundoff_line{
        MakeState(1.0401393842326945, -0.9999999330784480, 1.0679200000000000e-6, 0.50),
        MakeState(1.0319977031220682, -0.9999999466065290, 1.0540423762694930e-6, 0.50),
        MakeState(1.0238560220114403, -0.9999999601346100, 1.0401830912905772e-6, 0.50),
        MakeState(1.0158100556769198, -0.9999999735452161, 1.0265593098542545e-6, 0.50),
        MakeState(1.0078585291303450, -0.9999999868334754, 1.0131665281557882e-6, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50)
    };
    const auto near_monotone_outer_roundoff = ReconstructPpmLine(
        near_monotone_outer_roundoff_line,
        HydroDirection::radial,
        2u,
        3u,
        3.0127404516800261e-4,
        std::vector<double>(near_monotone_outer_roundoff_line.size(), 1.0 / 256.0));
    DEC3D_CHECK(near_monotone_outer_roundoff.is_complete(3u));
    DEC3D_CHECK(!near_monotone_outer_roundoff.interfaces[1].downgraded_to_first_order);
    DEC3D_CHECK(!near_monotone_outer_roundoff.interfaces[1].left_profile_troubled);

    const std::vector<dec3d::hydro::HydroConservativeState> near_monotone_outer_step59_roundoff_line{
        MakeState(1.0401393842326945, -0.9999999330784480, 1.0679200000000000e-6, 0.50),
        MakeState(1.0319956705031565, -0.9999999466570508, 1.0540386533956318e-6, 0.50),
        MakeState(1.0238557543663105, -0.9999999601421921, 1.0401825989436733e-6, 0.50),
        MakeState(1.0158100341982284, -0.9999999735459044, 1.0265592701822850e-6, 0.50),
        MakeState(1.0078585283354005, -0.9999999868335034, 1.0131665269715504e-6, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50),
        MakeState(1.0000000000000000, -1.0000000000000000, 9.9999999999174837e-7, 0.50)
    };
    const auto near_monotone_outer_step59_roundoff = ReconstructPpmLine(
        near_monotone_outer_step59_roundoff_line,
        HydroDirection::radial,
        2u,
        3u,
        2.0112527252721320e-4,
        std::vector<double>(near_monotone_outer_step59_roundoff_line.size(), 1.0 / 256.0));
    DEC3D_CHECK(near_monotone_outer_step59_roundoff.is_complete(3u));
    DEC3D_CHECK(!near_monotone_outer_step59_roundoff.interfaces[1].downgraded_to_first_order);
    DEC3D_CHECK(!near_monotone_outer_step59_roundoff.interfaces[1].left_profile_troubled);

    const double dr = 1.0 / 256.0;
    std::vector<dec3d::hydro::HydroConservativeState> homologous_compression_line;
    homologous_compression_line.reserve(8u);
    for (int cell = -3; cell < 5; ++cell) {
      const double signed_radius = (static_cast<double>(cell) + 0.5) * dr;
      homologous_compression_line.push_back(
          MakeState(1.0, -signed_radius, 2.0 / 3.0, 0.50));
    }
    const std::vector<double> homologous_widths(
        homologous_compression_line.size(),
        dr);
    const std::vector<double> moving_face_speeds{
        -0.0,
        -dr,
        -2.0 * dr,
    };
    const auto moving_trace = ReconstructPpmLine(
        homologous_compression_line,
        HydroDirection::radial,
        2u,
        3u,
        1.0e-4,
        homologous_widths,
        &moving_face_speeds);
    DEC3D_CHECK(moving_trace.is_complete(3u));
    const auto left_face1 =
        MakeDirectionalConservativeState(moving_trace.interfaces[1].left);
    const auto right_face1 =
        MakeDirectionalConservativeState(moving_trace.interfaces[1].right);
    const auto moving_flux_face1 =
        dec3d::hydro::SolveMovingInterfaceHllcFlux(
            left_face1,
            right_face1,
            moving_face_speeds[1]);
    DEC3D_CHECK(moving_flux_face1.is_complete());
    DEC3D_CHECK(std::abs(moving_flux_face1.flux.rho) <= 1.0e-10);

    const std::vector<dec3d::hydro::HydroConservativeState> oscillatory_line{
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 0.02, 0.10),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60),
        MakeState(1.00, 0.00, 1.00, 0.60)
    };
    const auto troubled = ReconstructPpmLine(
        oscillatory_line,
        HydroDirection::radial,
        3u,
        3u,
        1.0e-4,
        std::vector<double>(oscillatory_line.size(), 0.1));
    DEC3D_CHECK(troubled.is_complete(4u));
    DEC3D_CHECK(troubled.downgraded_interface_count > 0u);
    bool saw_downgraded = false;
    for (const auto& interface_state : troubled.interfaces) {
      if (interface_state.downgraded_to_first_order) {
        saw_downgraded = true;
        DEC3D_CHECK(interface_state.left.is_physical());
        DEC3D_CHECK(interface_state.right.is_physical());
      }
    }
    DEC3D_CHECK(saw_downgraded);

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
