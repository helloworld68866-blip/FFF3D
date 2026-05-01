#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace dec3d::mesh {

struct SphericalMeshDescriptor {
  std::size_t radial_cells{0};
  std::size_t theta_cells{0};
  std::size_t phi_cells{0};
  double inner_radius{0.0};
  double outer_radius{1.0};

  [[nodiscard]] bool is_valid() const noexcept {
    return radial_cells > 0 && theta_cells > 0 && phi_cells > 0 &&
           inner_radius >= 0.0 && outer_radius > inner_radius;
  }
};

struct SphericalGeometryMetadata {
  bool valid{false};
  std::vector<double> radial_faces;
  std::vector<double> theta_faces;
  std::vector<double> phi_faces;
  std::vector<double> cell_volumes;
  double global_volume{0.0};

  [[nodiscard]] bool is_valid() const noexcept { return valid; }
};

struct MeshBuildEvidence {
  bool built_by_new_code{false};
  std::string implementation_id;

  [[nodiscard]] bool is_present() const noexcept {
    return built_by_new_code && !implementation_id.empty();
  }
};

[[nodiscard]] SphericalGeometryMetadata BuildSphericalGeometry(
    const SphericalMeshDescriptor& descriptor) noexcept;

[[nodiscard]] MeshBuildEvidence BuildMeshSubstrateEvidence(
    const SphericalMeshDescriptor& descriptor) noexcept;

}  // namespace dec3d::mesh
