#define DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT 1
#define DEC3D_NOH_MPI24_RADIAL_CELLS 256u
#define DEC3D_NOH_MPI24_FINAL_TIME 0.105
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_nr256_p1e6_t0105_mpi1"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_nr256_p1e6_t0105_mpi1"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm.nr256.t0105.mpi1"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm.nr256.t0105.mpi1"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
