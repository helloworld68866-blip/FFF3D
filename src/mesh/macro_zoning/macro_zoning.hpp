#pragma once

#include "core/array/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::mesh {

[[nodiscard]] constexpr std::size_t InvalidMacroZoneIndex() noexcept {
  return static_cast<std::size_t>(-1);
}

struct MacroZoningDetectionOptions {
  double coarse_factor{0.5};
};

struct MacroZoneThetaBand {
  std::size_t radial_index{0};
  std::size_t theta_begin{0};
  std::size_t theta_end{0};
  std::size_t theta_factor{0};
  std::size_t phi_factor{0};
  double delta_s_theta{0.0};
  double min_delta_s_phi{0.0};

  [[nodiscard]] bool is_valid(std::size_t radial_cells, std::size_t theta_cells, std::size_t phi_cells) const noexcept;
};

struct MacroZoneCell {
  bool initialized{false};
  std::size_t radial_index{0};
  std::size_t theta_begin{0};
  std::size_t theta_end{0};
  std::size_t phi_begin{0};
  std::size_t phi_end{0};
  std::size_t theta_factor{0};
  std::size_t phi_factor{0};
  double total_volume{0.0};

  [[nodiscard]] bool is_valid(std::size_t radial_cells, std::size_t theta_cells, std::size_t phi_cells) const noexcept;
};

struct MacroZoneMap {
  bool valid{false};
  std::size_t fine_radial_cells{0};
  std::size_t fine_theta_cells{0};
  std::size_t fine_phi_cells{0};
  dec3d::core::Array3D<std::size_t> fine_to_coarse;
  std::vector<MacroZoneThetaBand> theta_bands;
  std::vector<MacroZoneCell> coarse_cells;
  bool phi_factor_varies_by_theta_band{false};
  std::size_t theta_factor_saturated_count{0u};
  std::size_t phi_factor_saturated_count{0u};
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct FineHydroPackageView {
  const dec3d::core::Array3D<double>* rho{nullptr};
  const dec3d::core::Array3D<double>* mom_r{nullptr};
  const dec3d::core::Array3D<double>* mom_theta{nullptr};
  const dec3d::core::Array3D<double>* mom_phi{nullptr};
  const dec3d::core::Array3D<double>* e_fluid_total{nullptr};
  const dec3d::core::Array3D<double>* chi_e{nullptr};

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MacroZonedHydroPackage {
  bool valid{false};
  bool short_lived_work_view{true};
  bool authoritative_backed{false};
  MacroZoneMap map;
  std::vector<double> rho;
  std::vector<double> mom_r;
  std::vector<double> mom_theta;
  std::vector<double> mom_phi;
  std::vector<double> e_fluid_total;
  std::vector<double> chi_e;
  std::string failure_reason;
  std::string report_line;

  [[nodiscard]] bool is_complete() const noexcept;
  [[nodiscard]] bool is_authoritative_truth() const noexcept { return authoritative_backed; }
};

struct FineHydroPackageBuffers {
  dec3d::core::Array3D<double> rho;
  dec3d::core::Array3D<double> mom_r;
  dec3d::core::Array3D<double> mom_theta;
  dec3d::core::Array3D<double> mom_phi;
  dec3d::core::Array3D<double> e_fluid_total;
  dec3d::core::Array3D<double> chi_e;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

struct MacroZoneAuthoritativeBoundaryResult {
  bool success{false};
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept { return !success && !failure_reason.empty(); }
};

[[nodiscard]] MacroZoneMap DetectMacroZones(
    const SphericalGeometryMetadata& geometry,
    const MacroZoningDetectionOptions& options = {}) noexcept;

[[nodiscard]] MacroZonedHydroPackage RestrictFineHydroPackage(
    const FineHydroPackageView& fine_view,
    const SphericalGeometryMetadata& geometry,
    const MacroZoneMap& map) noexcept;

[[nodiscard]] FineHydroPackageBuffers ProlongMacroZoneHydroPackage(
    const MacroZonedHydroPackage& coarse_package) noexcept;

[[nodiscard]] MacroZoneAuthoritativeBoundaryResult RejectMacroZoneAuthoritativeCommit(
    const MacroZonedHydroPackage& coarse_package) noexcept;

}  // namespace dec3d::mesh
