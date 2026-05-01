#define DEC3D_SEDOV_MPI16_RADIAL_CELLS 1024u
#define DEC3D_SEDOV_MPI16_USE_MACRO_ZONING 1
#define DEC3D_SEDOV_MPI16_EXPECTED_RANK_COUNT 24
#define DEC3D_SEDOV_MPI16_AMBIENT_PRESSURE 1.0e-8
#define DEC3D_SEDOV_MPI16_BLAST_RADIUS 5.0e-3
#define DEC3D_SEDOV_MPI16_OUTPUT_DIR_LITERAL \
  "F:\\dec3d\\analysis\\output\\case_sedov_spherical_ppm_macro_nr1024_br0005_p1e8_mpi24"
#define DEC3D_SEDOV_MPI16_CASE_NAME_LITERAL \
  "case_sedov_spherical_ppm_macro_nr1024_br0005_p1e8_mpi24"
#define DEC3D_SEDOV_MPI16_MESH_HANDLE_LITERAL \
  "mesh:p1.sedov.ppm-macro.nr1024.br0005.p1e8.mpi24"
#define DEC3D_SEDOV_MPI16_DIAGNOSTICS_HANDLE_LITERAL \
  "diagnostics:p1.hydro.sedov.ppm-macro.nr1024.br0005.p1e8.mpi24"

#include "benchmark_sedov_ppm_trial_nr512_mpi16.cpp"
