#include "hydro/reconstruction/ppm_reconstruction.hpp"
#include "test_assert.hpp"

#include <vector>

int main() {
  using dec3d::hydro::BuildShellWidePpmFallbackMask;
  using dec3d::hydro::PpmLineReconstructionResult;

  std::vector<PpmLineReconstructionResult> angular_lines(4);
  for (auto& line : angular_lines) {
    line.success = true;
    line.interface_count = 3u;
    line.interfaces.resize(3u);
  }

  angular_lines[2].interfaces[1].downgraded_to_first_order = true;

  const auto mask = BuildShellWidePpmFallbackMask(angular_lines, 3u);

  DEC3D_CHECK_EQ(mask.size(), 3u);
  DEC3D_CHECK(!mask[0]);
  DEC3D_CHECK(mask[1]);
  DEC3D_CHECK(!mask[2]);

  return 0;
}
