#define DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT 1
#define DEC3D_NOH_MPI24_RADIAL_CELLS 256u
#define DEC3D_NOH_MPI24_USE_ALE 1
#define DEC3D_NOH_MPI24_FINAL_TIME 0.2
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_ale_nr256_p1e6_t02_mpi1"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_ale_nr256_p1e6_t02_mpi1"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-ale.nr256.t02.mpi1"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-ale.nr256.t02.mpi1"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
