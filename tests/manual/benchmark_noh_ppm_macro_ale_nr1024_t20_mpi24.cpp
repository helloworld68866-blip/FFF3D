#define DEC3D_NOH_MPI24_RADIAL_CELLS 1024u
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_USE_ALE 1
#define DEC3D_NOH_MPI24_FINAL_TIME 2.0
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_macro_ale_nr1024_p1e6_t20_mpi24"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_macro_ale_nr1024_p1e6_t20_mpi24"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-macro-ale.nr1024.t20.mpi24"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-macro-ale.nr1024.t20.mpi24"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
