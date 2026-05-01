#pragma once

#include "core/diagnostics/stage_contracts.hpp"
#include "hydro/driver/static_grid_hydro.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"

namespace dec3d::hydro {

class HydroOperator {
 public:
  [[nodiscard]] bool bind(
      const dec3d::core::StageContext& context,
      const dec3d::mesh::SphericalGeometryMetadata& geometry,
      dec3d::state::CanonicalState& state) noexcept;

  [[nodiscard]] bool is_bound() const noexcept { return bound_; }

  [[nodiscard]] dec3d::core::DtAdvice estimate_dt() const noexcept;

  [[nodiscard]] dec3d::core::StageResult advance() noexcept;

  void SetStaticGridOptions(const StaticGridHydroOptions& options) noexcept {
    static_grid_options_ = options;
  }

 private:
  [[nodiscard]] bool has_allocated_state() const noexcept;

  bool bound_{false};
  dec3d::core::StageContext context_{};
  const dec3d::mesh::SphericalGeometryMetadata* geometry_{nullptr};
  dec3d::state::CanonicalState* state_{nullptr};
  StaticGridHydroOptions static_grid_options_{};
};

}  // namespace dec3d::hydro
