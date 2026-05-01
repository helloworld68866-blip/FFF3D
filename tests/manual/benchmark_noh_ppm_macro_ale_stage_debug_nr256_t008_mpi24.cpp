#define DEC3D_NOH_MPI24_RADIAL_CELLS 256u
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_USE_ALE 1
#define DEC3D_NOH_MPI24_APPLY_GEOMETRIC_SOURCE 1
#define DEC3D_NOH_MPI24_APPLY_THETA_SWEEP 1
#define DEC3D_NOH_MPI24_APPLY_PHI_SWEEP 1
#define DEC3D_NOH_MPI24_FINAL_TIME 0.1
#define DEC3D_NOH_MPI24_MAX_STEPS 1000000u
#define DEC3D_NOH_MPI24_WRITE_STAGE_ANGULAR_DEBUG 1
#define DEC3D_NOH_MPI24_STAGE_DEBUG_STEPS 220u
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE_SMOOTHNESS 0
#define DEC3D_NOH_MPI24_DEBUG_MACRO_RADIAL_PPM_FACE 0u
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_macro_ale_nr256_p1e6_t010_mpi24_stage_ops"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_macro_ale_nr256_p1e6_t010_mpi24_stage_ops"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-macro-ale.nr256.t010.stage-ops.mpi24"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-macro-ale.nr256.t010.stage-ops.mpi24"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
