#pragma once

#include <string>
#include <vector>

namespace dec3d::app {

struct Dec3DAppResult {
  int exit_code{1};
  std::string report_line;
  std::string output_report_line;
  std::string failure_reason;
  std::string failure_diagnostics;
};

[[nodiscard]] Dec3DAppResult RunDec3DCommandLine(
    const std::vector<std::string>& argv) noexcept;

}  // namespace dec3d::app
