#define DEC3D_NOH_MPI24_EXPECTED_RANK_COUNT 16
#define DEC3D_NOH_MPI24_USE_MACRO_ZONING 1
#define DEC3D_NOH_MPI24_USE_PPM 1
#define DEC3D_NOH_MPI24_MACRO_RADIAL_PPM 1
#define DEC3D_NOH_MPI24_MACRO_THETA_PPM 0
#define DEC3D_NOH_MPI24_MACRO_PHI_PPM 0
#define DEC3D_NOH_MPI24_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_noh_spherical_ppm_macro_radial_ppm_ang_first_order_nr256_p1e6_mpi16"
#define DEC3D_NOH_MPI24_CASE_NAME_LITERAL \
  "case_noh_spherical_ppm_macro_radial_ppm_ang_first_order_nr256_p1e6_mpi16"
#define DEC3D_NOH_MPI24_MESH_HANDLE_LITERAL \
  "mesh:p1.noh.ppm-macro.radial-ppm-ang-first-order.nr256.mpi16"
#define DEC3D_NOH_MPI24_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.noh.ppm-macro.radial-ppm-ang-first-order.nr256.mpi16"

#include "benchmark_noh_ppm_nr256_mpi24.cpp"
