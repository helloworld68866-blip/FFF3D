#include "app/dec3d_app.hpp"

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#ifdef DEC3D_ENABLE_HYPRE
#include <mpi.h>
#endif

int main(int argc, char** argv) {
#ifdef DEC3D_ENABLE_HYPRE
  int mpi_initialized = 0;
  MPI_Initialized(&mpi_initialized);
  const bool initialized_here = mpi_initialized == 0;
  if (initialized_here) {
    MPI_Init(&argc, &argv);
  }
#endif
  std::vector<std::string> args;
  args.reserve(static_cast<std::size_t>(argc));
  for (int i = 0; i < argc; ++i) {
    args.emplace_back(argv[i]);
  }
  const auto result = dec3d::app::RunDec3DCommandLine(args);
  if (result.exit_code == 0) {
    std::cout << result.report_line << '\n';
  } else {
    std::cerr << result.failure_diagnostics << '\n';
  }
#ifdef DEC3D_ENABLE_HYPRE
  int mpi_finalized = 0;
  MPI_Finalized(&mpi_finalized);
  if (initialized_here && mpi_finalized == 0) {
    MPI_Finalize();
  }
#endif
  return result.exit_code;
}
