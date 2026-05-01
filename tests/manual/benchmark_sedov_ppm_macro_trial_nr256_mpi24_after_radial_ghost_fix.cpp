#define DEC3D_SEDOV_MPI16_RADIAL_CELLS 256u
#define DEC3D_SEDOV_MPI16_USE_MACRO_ZONING 1
#define DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT 24
#define DEC3D_SEDOV_MPI16_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_sedov_spherical_ppm_macro_nr256_mpi24_after_radial_ghost_fix"
#define DEC3D_SEDOV_MPI16_CASE_NAME_LITERAL \
  "case_sedov_spherical_ppm_macro_nr256_mpi24_after_radial_ghost_fix"
#define DEC3D_SEDOV_MPI16_MESH_HANDLE_LITERAL \
  "mesh:p1.sedov.ppm-macro.nr256.mpi24.after-radial-ghost-fix"
#define DEC3D_SEDOV_MPI16_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.sedov.ppm-macro.nr256.mpi24.after-radial-ghost-fix"

#include "benchmark_sedov_ppm_trial_nr512_mpi16.cpp"
