#pragma once

#include "core/array/array3d.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cstddef>
#include <string>

namespace dec3d::state {

struct ThermodynamicRecoveryOptions {
  double electron_energy_floor{0.0};
  double ion_energy_floor{0.0};
  double zbar{1.0};
  double mean_ion_mass_g{0.0};
};

struct ThermodynamicCellState {
  double rho_g_per_cm3{0.0};
  double kinetic_energy_density_erg_per_cm3{0.0};
  double e_electron_erg_per_cm3{0.0};
  double e_ion_erg_per_cm3{0.0};
  double p_e_erg_per_cm3{0.0};
  double p_i_erg_per_cm3{0.0};
  double n_e_cm3{0.0};
  double n_i_cm3{0.0};
  double t_e_erg_per_particle{0.0};
  double t_i_erg_per_particle{0.0};
  double t_e_keV{0.0};
  double t_i_keV{0.0};
  double mean_ion_mass_g{0.0};
  double zbar{0.0};
  double energy_partition_residual{0.0};
};

struct ThermodynamicRecoveryResult {
  bool success{false};
  std::size_t cell_count{0};
  dec3d::core::Array3D<ThermodynamicCellState> cells;
  std::string unit_diagnostics;
  std::string recovery_diagnostics;
  std::string failure_diagnostics;
  std::string failure_reason;

  [[nodiscard]] bool is_complete() const noexcept;
};

[[nodiscard]] ThermodynamicRecoveryResult RecoverThermodynamicState(
    const CanonicalState& state,
    ThermodynamicRecoveryOptions options = {}) noexcept;

[[nodiscard]] bool ValidateThermodynamicRecoveryDiagnostics(
    const ThermodynamicRecoveryResult& result) noexcept;

}  // namespace dec3d::state
