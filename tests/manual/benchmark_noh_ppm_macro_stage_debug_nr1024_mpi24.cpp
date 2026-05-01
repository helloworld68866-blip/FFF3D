#define DEC3D_NOH_MPI24_RADIAL_CELLS 1024u
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_FINAL_TIME 0.025
#define DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG 1
#define DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS 220u
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS 1
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE 36u
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_macro_nr1024_p1e6_mpi24_stage_debug"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_macro_nr1024_p1e6_mpi24_stage_debug"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-macro.nr1024.mpi24.stage-debug"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-macro.nr1024.mpi24.stage-debug"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
