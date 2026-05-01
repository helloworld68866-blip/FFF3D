#define DEC3D_NOH_MPI24_USE_NOH2 1
#define DEC3D_NOH_MPI24_RADIAL_CELLS 256u
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_USE_ALE 1
#define DEC3D_NOH_MPI24_FINAL_TIME 0.003
#define DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG 1
#define DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS 4u
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS 1
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE 1u
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh2_spherical_ppm_macro_ale_stage_debug_nr256_t003_mpi24"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh2_spherical_ppm_macro_ale_stage_debug_nr256_t003_mpi24"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh2.ppm-macro-ale-stage-debug.nr256.t003.mpi24"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh2.ppm-macro-ale-stage-debug.nr256.t003.mpi24"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
