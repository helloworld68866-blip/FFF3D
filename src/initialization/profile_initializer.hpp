#pragma once

#include "io/input_deck.hpp"
#include "io/radial_profile.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <filesystem>
#include <string>

namespace dec3d::initialization {

struct ProfileInitializationResult {
  bool success{false};
  dec3d::state::CanonicalState state;
  dec3d::mesh::SphericalGeometryMetadata geometry;
  dec3d::radiation::RadiationGroupLayout group_layout;
  std::string report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
  std::string thermodynamic_recovery_report;
};

[[nodiscard]] ProfileInitializationResult InitializeFromRadialProfile(
    const dec3d::io::InputDeckConfig& config,
    const dec3d::io::RadialProfile& profile,
    const std::filesystem::path& profile_path) noexcept;

[[nodiscard]] bool ValidateProfileInitializationDiagnostics(
    const ProfileInitializationResult& result) noexcept;

}  // namespace dec3d::initialization
