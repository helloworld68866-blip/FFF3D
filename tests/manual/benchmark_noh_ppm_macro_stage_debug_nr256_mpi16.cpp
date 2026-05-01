#define DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT 16
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_USE_PPM 1
#define DEC3D_NOH_MPI24_FINAL_TIME 0.08
#define DEC3D_NOH_MPI24_MAX_STEPS 200u
#define DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG 1
#define DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS 200u
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_macro_stage_debug_nr256_p1e6_mpi16"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_macro_stage_debug_nr256_p1e6_mpi16"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-macro.stage-debug.nr256.mpi16"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-macro.stage-debug.nr256.mpi16"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
