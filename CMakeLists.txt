cmake_minimum_required(VERSION 3.24)

project(dec3d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(DEC3D_BUILD_TESTS "Build DEC3D unit, contract, and acceptance tests." ON)
option(DEC3D_ENABLE_HYPRE "Build the real HYPRE diffusion adapter." OFF)

find_package(MPI COMPONENTS CXX)

if(MPI_CXX_FOUND)
  message(STATUS "DEC3D MPI support enabled")
else()
  message(STATUS "DEC3D MPI support not found; MPI hydro tests will be skipped")
endif()

if(DEC3D_ENABLE_HYPRE)
  if(NOT MPI_CXX_FOUND)
    message(FATAL_ERROR "DEC3D_ENABLE_HYPRE requires MPI CXX support because HYPRE IJ/ParCSR uses MPI communicators")
  endif()
  find_package(HYPRE CONFIG REQUIRED)
  if(TARGET HYPRE::HYPRE)
    set(DEC3D_HYPRE_TARGET HYPRE::HYPRE)
  elseif(TARGET HYPRE::hypre)
    set(DEC3D_HYPRE_TARGET HYPRE::hypre)
  elseif(TARGET hypre)
    set(DEC3D_HYPRE_TARGET hypre)
  else()
    message(FATAL_ERROR "DEC3D_ENABLE_HYPRE requires a real imported HYPRE target")
  endif()
  message(STATUS "DEC3D HYPRE support enabled")
else()
  message(STATUS "DEC3D HYPRE support disabled")
endif()

function(dec3d_apply_warnings target_name)
  if(MSVC)
    target_compile_options(${target_name} PRIVATE /W4 /permissive- /EHsc)
  else()
    target_compile_options(${target_name} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
endfunction()

add_library(dec3d_core STATIC
  src/core/array/array3d.cpp
  src/core/diagnostics/stage_contracts.cpp
  src/core/logging/log_event.cpp
)
target_include_directories(dec3d_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
dec3d_apply_warnings(dec3d_core)

add_library(dec3d_io STATIC
  src/io/input_deck.cpp
  src/io/radial_profile.cpp
  src/io/runtime_output.cpp
)
target_include_directories(dec3d_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_io PUBLIC dec3d_core)
dec3d_apply_warnings(dec3d_io)

set(DEC3D_MESH_SOURCES
  src/mesh/ale/radial_ale.cpp
  src/mesh/boundary/spherical_scalar_remap.cpp
  src/mesh/ghost/ghost_topology.cpp
  src/mesh/macro_zoning/macro_zoning.cpp
  src/mesh/ownership/radial_ownership.cpp
  src/mesh/spherical/spherical_mesh.cpp
)

if(MPI_CXX_FOUND)
  list(APPEND DEC3D_MESH_SOURCES
    src/mesh/ale/radial_ale_mpi.cpp
    src/mesh/ghost/hydro_halo_exchange.cpp
  )
endif()

add_library(dec3d_mesh STATIC
  ${DEC3D_MESH_SOURCES}
)
target_include_directories(dec3d_mesh PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_mesh PUBLIC dec3d_core)
if(MPI_CXX_FOUND)
  target_link_libraries(dec3d_mesh PRIVATE MPI::MPI_CXX)
endif()
dec3d_apply_warnings(dec3d_mesh)

add_library(dec3d_state STATIC
  src/state/canonical_state/canonical_state.cpp
  src/state/canonical_state/state_anchor.cpp
  src/state/canonical_state/state_utilities.cpp
  src/state/hydro_state/hydro_view.cpp
  src/state/thermodynamics/electron_ion_equilibration.cpp
  src/state/thermodynamics/electron_ion_tau.cpp
  src/state/thermodynamics/thermodynamic_recovery.cpp
)
target_include_directories(dec3d_state PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_state PUBLIC dec3d_core)
dec3d_apply_warnings(dec3d_state)

target_link_libraries(dec3d_io PUBLIC dec3d_state)

add_library(dec3d_alpha STATIC
  src/alpha/alpha_contract.cpp
  src/alpha/alpha_coefficients.cpp
  src/alpha/alpha_hydro_terms.cpp
)
target_include_directories(dec3d_alpha PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_alpha PUBLIC dec3d_core dec3d_state)
dec3d_apply_warnings(dec3d_alpha)

add_library(dec3d_transport STATIC
  src/transport/diffusion/generic_diffusion.cpp
  src/transport/thermal/electron_flux_limiter.cpp
  src/transport/thermal/thermal_conduction.cpp
  src/transport/thermal/thermal_conductivity.cpp
)
target_include_directories(dec3d_transport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_transport PUBLIC dec3d_core dec3d_mesh)
dec3d_apply_warnings(dec3d_transport)

add_library(dec3d_alpha_transport STATIC
  src/alpha/alpha_operator.cpp
)
target_include_directories(dec3d_alpha_transport PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_alpha_transport PUBLIC dec3d_alpha dec3d_transport dec3d_mesh dec3d_state)
dec3d_apply_warnings(dec3d_alpha_transport)

add_library(dec3d_radiation STATIC
  src/radiation/groups/radiation_group_layout.cpp
  src/radiation/hydro_terms/radiation_hydro_terms.cpp
  src/radiation/providers/group_blackbody.cpp
  src/radiation/providers/tops_opacity_provider.cpp
  src/radiation/transport/one_group_gray_diffusion.cpp
  src/radiation/transport/multigroup_gray_radiation.cpp
  src/radiation/transport/radiation_flux_limiter.cpp
)
target_include_directories(dec3d_radiation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_radiation PUBLIC dec3d_core dec3d_mesh dec3d_state dec3d_transport)
if(MPI_CXX_FOUND)
  target_link_libraries(dec3d_radiation PUBLIC MPI::MPI_CXX)
endif()
dec3d_apply_warnings(dec3d_radiation)

add_library(dec3d_initialization STATIC
  src/initialization/profile_initializer.cpp
)
target_include_directories(dec3d_initialization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_initialization
  PUBLIC
    dec3d_core
    dec3d_io
    dec3d_mesh
    dec3d_state
    dec3d_radiation
    dec3d_alpha
)
dec3d_apply_warnings(dec3d_initialization)

add_library(dec3d_app STATIC
  src/app/dec3d_app.cpp
)
target_include_directories(dec3d_app PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_app PUBLIC dec3d_io dec3d_initialization)
dec3d_apply_warnings(dec3d_app)

add_executable(dec3d
  src/app/dec3d_main.cpp
)
target_link_libraries(dec3d PRIVATE dec3d_app)
dec3d_apply_warnings(dec3d)

add_library(dec3d_benchmarks STATIC
  src/benchmarks/p2_thermal_benchmark.cpp
  src/benchmarks/p3_radiation_benchmark.cpp
  src/benchmarks/p3_radiation_benchmark_runner.cpp
  src/benchmarks/p4_alpha_benchmark.cpp
  src/benchmarks/p4_alpha_benchmark_runner.cpp
  src/benchmarks/p5_integrated_benchmark.cpp
)
target_include_directories(dec3d_benchmarks PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_benchmarks PUBLIC dec3d_core dec3d_mesh dec3d_state dec3d_transport dec3d_radiation dec3d_alpha_transport dec3d_alpha)
dec3d_apply_warnings(dec3d_benchmarks)

if(DEC3D_ENABLE_HYPRE)
  add_library(dec3d_transport_hypre STATIC
    src/transport/diffusion/hypre_diffusion_solver.cpp
    src/transport/diffusion/hypre_distributed_diffusion_solver.cpp
    src/transport/thermal/distributed_thermal_conduction.cpp
  )
  target_include_directories(dec3d_transport_hypre PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(dec3d_transport_hypre
    PUBLIC
      dec3d_transport
      dec3d_state
      MPI::MPI_CXX
      ${DEC3D_HYPRE_TARGET}
  )
  dec3d_apply_warnings(dec3d_transport_hypre)

  target_sources(dec3d_alpha_transport PRIVATE
    src/alpha/distributed_alpha_operator.cpp
  )
  target_link_libraries(dec3d_alpha_transport PUBLIC
    dec3d_transport_hypre
    MPI::MPI_CXX
  )

  add_library(dec3d_benchmarks_hypre STATIC
    src/benchmarks/p2_thermal_benchmark_runner.cpp
    src/benchmarks/p4_alpha_benchmark_runner_hypre.cpp
    src/benchmarks/p5_integrated_benchmark_runner_hypre.cpp
  )
  target_include_directories(dec3d_benchmarks_hypre PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(dec3d_benchmarks_hypre
    PUBLIC
      dec3d_benchmarks
      dec3d_alpha_transport
      dec3d_alpha
      dec3d_transport_hypre
      dec3d_runtime
      dec3d_state
      MPI::MPI_CXX
      ${DEC3D_HYPRE_TARGET}
  )
  dec3d_apply_warnings(dec3d_benchmarks_hypre)

  add_library(dec3d_radiation_hypre STATIC
    src/radiation/transport/distributed_multigroup_radiation.cpp
  )
  target_include_directories(dec3d_radiation_hypre PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(dec3d_radiation_hypre
    PUBLIC
      dec3d_radiation
      dec3d_transport_hypre
      dec3d_state
      MPI::MPI_CXX
      ${DEC3D_HYPRE_TARGET}
  )
  dec3d_apply_warnings(dec3d_radiation_hypre)

  target_compile_definitions(dec3d_app PUBLIC DEC3D_ENABLE_HYPRE=1)
  target_link_libraries(dec3d_app PUBLIC
    dec3d_transport_hypre
    dec3d_radiation_hypre
    MPI::MPI_CXX
    ${DEC3D_HYPRE_TARGET}
  )
  target_compile_definitions(dec3d PRIVATE DEC3D_ENABLE_HYPRE=1)
  target_link_libraries(dec3d PRIVATE MPI::MPI_CXX ${DEC3D_HYPRE_TARGET})

endif()

add_library(dec3d_runtime STATIC
  src/runtime/checkpoint_contract.cpp
  src/runtime/p2_production_orchestrator.cpp
  src/runtime/p3_production_orchestrator.cpp
  src/runtime/p4_production_orchestrator.cpp
  src/runtime/p1_acceptance.cpp
  src/runtime/p0_acceptance.cpp
  src/runtime/runtime_scaffold.cpp
)
target_include_directories(dec3d_runtime PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_runtime PUBLIC dec3d_core dec3d_mesh dec3d_state)
dec3d_apply_warnings(dec3d_runtime)

add_library(dec3d_hydro STATIC
  src/hydro/driver/hydro_numerical_checks.cpp
  src/hydro/reconstruction/ppm_reconstruction.cpp
  src/hydro/riemann/hllc_solver.cpp
  src/hydro/source/geometric_source.cpp
  src/hydro/driver/hydro_boundary_ghosts.cpp
  src/hydro/driver/macro_ale_hllc.cpp
  src/hydro/driver/radial_boundary_contract.cpp
  src/hydro/driver/radial_overlap_remap.cpp
  src/hydro/driver/static_grid_hydro.cpp
  src/hydro/driver/theta_pole_boundary.cpp
  src/hydro/driver/hydro_operator.cpp
)
target_include_directories(dec3d_hydro PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(dec3d_hydro PUBLIC dec3d_core dec3d_mesh dec3d_state)
if(MPI_CXX_FOUND)
  target_sources(dec3d_hydro PRIVATE
    src/hydro/driver/macro_ale_hllc_mpi.cpp
    src/hydro/driver/macro_ale_mpi_compatibility.cpp
  )
  target_link_libraries(dec3d_hydro PUBLIC MPI::MPI_CXX)
endif()
dec3d_apply_warnings(dec3d_hydro)

target_link_libraries(dec3d_app
  PUBLIC
    dec3d_hydro
    dec3d_transport
    dec3d_radiation
    dec3d_alpha_transport
    dec3d_state
)

if(DEC3D_BUILD_TESTS)
  enable_testing()

  function(dec3d_add_test test_name)
    add_executable(${test_name} ${ARGN})
    target_include_directories(${test_name}
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(${test_name}
      PRIVATE
        dec3d_hydro
        dec3d_runtime
        dec3d_transport
        dec3d_state
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(${test_name})
    add_test(NAME ${test_name} COMMAND ${test_name})
  endfunction()

  function(dec3d_add_mpi_test test_name rank_count)
    add_executable(${test_name} ${ARGN})
    target_include_directories(${test_name}
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(${test_name}
      PRIVATE
        MPI::MPI_CXX
        dec3d_hydro
        dec3d_runtime
        dec3d_state
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(${test_name})

    set(_mpi_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _mpi_command ${MPIEXEC_NUMPROC_FLAG} ${rank_count})
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _mpi_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _mpi_command $<TARGET_FILE:${test_name}>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _mpi_command ${MPIEXEC_POSTFLAGS})
    endif()

    add_test(NAME ${test_name} COMMAND ${_mpi_command})
  endfunction()

  dec3d_add_test(dec3d_unit_array3d
    tests/unit/test_array3d.cpp
  )

  dec3d_add_test(dec3d_contract_stage_result
    tests/contract/test_stage_result.cpp
  )

  dec3d_add_test(dec3d_contract_cache_invalidation
    tests/contract/test_cache_invalidation.cpp
  )

  dec3d_add_test(dec3d_contract_spherical_geometry
    tests/contract/test_spherical_geometry.cpp
  )

  dec3d_add_test(dec3d_contract_spherical_scalar_remap
    tests/contract/test_spherical_scalar_remap.cpp
  )

  dec3d_add_test(dec3d_contract_radial_ownership
    tests/contract/test_radial_ownership.cpp
  )

  dec3d_add_test(dec3d_contract_macro_zoning_detection
    tests/contract/test_macro_zoning_detection.cpp
  )

  dec3d_add_test(dec3d_contract_macro_zoning_conservative_transfer
    tests/contract/test_macro_zoning_conservative_transfer.cpp
  )

  dec3d_add_test(dec3d_contract_macro_zoning_not_authoritative
    tests/contract/test_macro_zoning_not_authoritative.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_zoning_single_rank
    tests/contract/test_hydro_macro_zoning_single_rank.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_zoning_ppm_single_rank
    tests/contract/test_hydro_macro_zoning_ppm_single_rank.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_zoning_transaction
    tests/contract/test_hydro_macro_zoning_transaction.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_zoning_budget
    tests/contract/test_hydro_macro_zoning_budget.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_zoning_spherical_one_step
    tests/contract/test_hydro_macro_zoning_spherical_one_step.cpp
  )

  dec3d_add_test(dec3d_contract_mesh_runtime_report
    tests/contract/test_mesh_runtime_report.cpp
  )

  dec3d_add_test(dec3d_contract_canonical_state_storage
    tests/contract/test_canonical_state_storage.cpp
  )

  dec3d_add_test(dec3d_contract_canonical_state_cached_commit
    tests/contract/test_canonical_state_cached_commit.cpp
  )

  dec3d_add_test(dec3d_contract_canonical_state_invalidation
    tests/contract/test_canonical_state_invalidation.cpp
  )

  dec3d_add_test(dec3d_contract_thermodynamic_feasibility
    tests/contract/test_thermodynamic_feasibility.cpp
  )

  dec3d_add_test(dec3d_contract_physics_units_constants
    tests/contract/test_physics_units_constants.cpp
  )

  dec3d_add_test(dec3d_contract_thermodynamic_recovery
    tests/contract/test_thermodynamic_recovery.cpp
  )

  dec3d_add_test(dec3d_contract_electron_ion_equilibration
    tests/contract/test_electron_ion_equilibration.cpp
  )

  dec3d_add_test(dec3d_contract_electron_ion_tau_provider
    tests/contract/test_electron_ion_tau_provider.cpp
  )

  dec3d_add_test(dec3d_contract_generic_implicit_diffusion
    tests/contract/test_generic_implicit_diffusion.cpp
  )

  dec3d_add_test(dec3d_contract_thermal_conduction_constant_kappa
    tests/contract/test_thermal_conduction_constant_kappa.cpp
  )

  dec3d_add_test(dec3d_contract_thermal_conductivity_coefficients
    tests/contract/test_thermal_conductivity_coefficients.cpp
  )

  dec3d_add_test(dec3d_contract_electron_flux_limiter
    tests/contract/test_electron_flux_limiter.cpp
  )

  dec3d_add_test(dec3d_contract_p2_production_orchestrator
    tests/contract/test_p2_production_orchestrator.cpp
  )

  dec3d_add_test(dec3d_contract_p3_production_orchestrator
    tests/contract/test_p3_production_orchestrator.cpp
  )

  dec3d_add_test(dec3d_contract_p4_production_orchestrator
    tests/contract/test_p4_production_orchestrator.cpp
  )
  target_link_libraries(dec3d_contract_p4_production_orchestrator
    PRIVATE
      dec3d_runtime
      dec3d_core
  )

  dec3d_add_test(dec3d_contract_p2_thermal_benchmark_contract
    tests/contract/test_p2_thermal_benchmark_contract.cpp
  )
  target_link_libraries(dec3d_contract_p2_thermal_benchmark_contract PRIVATE dec3d_benchmarks)

  dec3d_add_test(dec3d_contract_p3_radiation_benchmark_contract
    tests/contract/test_p3_radiation_benchmark_contract.cpp
  )
  target_link_libraries(dec3d_contract_p3_radiation_benchmark_contract PRIVATE dec3d_benchmarks dec3d_radiation)

  dec3d_add_test(dec3d_contract_p4_alpha_benchmark_contract
    tests/contract/test_p4_alpha_benchmark_contract.cpp
  )
  target_link_libraries(dec3d_contract_p4_alpha_benchmark_contract
    PRIVATE
      dec3d_benchmarks
      dec3d_alpha_transport
      dec3d_alpha
      dec3d_state
      dec3d_mesh
  )

  dec3d_add_test(dec3d_contract_p5_integrated_benchmark_contract
    tests/contract/test_p5_integrated_benchmark_contract.cpp
  )
  target_link_libraries(dec3d_contract_p5_integrated_benchmark_contract
    PRIVATE
      dec3d_benchmarks
  )

  dec3d_add_test(dec3d_contract_p5_io_input_deck
    tests/contract/test_p5_io_input_deck.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_input_deck PRIVATE dec3d_io)

  dec3d_add_test(dec3d_contract_p5_io_radial_profile
    tests/contract/test_p5_io_radial_profile.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_radial_profile PRIVATE dec3d_io)

  dec3d_add_test(dec3d_contract_p5_io_output_schedule
    tests/contract/test_p5_io_output_schedule.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_output_schedule PRIVATE dec3d_io dec3d_core)

  dec3d_add_test(dec3d_contract_p5_io_checkpoint_files
    tests/contract/test_p5_io_checkpoint_files.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_checkpoint_files
    PRIVATE
      dec3d_initialization
      dec3d_io
      dec3d_radiation
      dec3d_state
      dec3d_mesh
      dec3d_core
  )

  dec3d_add_test(dec3d_contract_p5_io_profile_initializer
    tests/contract/test_p5_io_profile_initializer.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_profile_initializer
    PRIVATE
      dec3d_initialization
      dec3d_io
      dec3d_radiation
      dec3d_state
  )

  dec3d_add_test(dec3d_contract_p5_io_cli_entry
    tests/contract/test_p5_io_cli_entry.cpp
  )
  target_link_libraries(dec3d_contract_p5_io_cli_entry PRIVATE dec3d_app)

  dec3d_add_test(dec3d_contract_p5_runtime_all_stages
    tests/contract/test_p5_runtime_all_stages.cpp
  )
  target_link_libraries(dec3d_contract_p5_runtime_all_stages PRIVATE dec3d_app)

  add_test(
    NAME dec3d_contract_p5_io_executable_smoke
    COMMAND dec3d
      -input ${CMAKE_CURRENT_SOURCE_DIR}/cases/p5_io_smoke.in
      -profile ${CMAKE_CURRENT_SOURCE_DIR}/cases/p5_io_smoke.pro
  )

  dec3d_add_test(dec3d_contract_radiation_group_layout
    tests/contract/test_radiation_group_layout.cpp
  )
  target_link_libraries(dec3d_contract_radiation_group_layout PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_radiation_one_group_gray_diffusion
    tests/contract/test_radiation_one_group_gray_diffusion.cpp
  )
  target_link_libraries(dec3d_contract_radiation_one_group_gray_diffusion PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_radiation_multigroup_gray_coupling
    tests/contract/test_radiation_multigroup_gray_coupling.cpp
  )
  target_link_libraries(dec3d_contract_radiation_multigroup_gray_coupling PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_radiation_tops_opacity_provider
    tests/contract/test_radiation_tops_opacity_provider.cpp
  )
  target_link_libraries(dec3d_contract_radiation_tops_opacity_provider PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_radiation_flux_limiter
    tests/contract/test_radiation_flux_limiter.cpp
  )
  target_link_libraries(dec3d_contract_radiation_flux_limiter PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_radiation_hydro_terms
    tests/contract/test_radiation_hydro_terms.cpp
  )
  target_link_libraries(dec3d_contract_radiation_hydro_terms PRIVATE dec3d_radiation dec3d_state)

  dec3d_add_test(dec3d_contract_alpha_contract
    tests/contract/test_alpha_contract.cpp
  )
  target_link_libraries(dec3d_contract_alpha_contract PRIVATE dec3d_alpha dec3d_state)

  dec3d_add_test(dec3d_contract_alpha_hydro_terms
    tests/contract/test_alpha_hydro_terms.cpp
  )
  target_link_libraries(dec3d_contract_alpha_hydro_terms PRIVATE
    dec3d_alpha
    dec3d_state
    dec3d_hydro
    dec3d_mesh
  )

  dec3d_add_test(dec3d_contract_alpha_one_group_operator
    tests/contract/test_alpha_one_group_operator.cpp
  )
  target_link_libraries(dec3d_contract_alpha_one_group_operator PRIVATE
    dec3d_alpha_transport
    dec3d_alpha
    dec3d_transport
    dec3d_state
    dec3d_mesh
  )

  dec3d_add_test(dec3d_contract_thermal_conduction_variable_kappa
    tests/contract/test_thermal_conduction_variable_kappa.cpp
  )

  if(DEC3D_ENABLE_HYPRE)
    add_executable(dec3d_contract_hypre_diffusion_solver
      tests/contract/test_hypre_diffusion_solver.cpp
    )
    target_include_directories(dec3d_contract_hypre_diffusion_solver
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(dec3d_contract_hypre_diffusion_solver
      PRIVATE
        dec3d_transport_hypre
        dec3d_transport
        MPI::MPI_CXX
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(dec3d_contract_hypre_diffusion_solver)
    add_test(
      NAME dec3d_contract_hypre_diffusion_solver
      COMMAND dec3d_contract_hypre_diffusion_solver
    )

    set(_hypre_mpi_reject_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _hypre_mpi_reject_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _hypre_mpi_reject_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _hypre_mpi_reject_command $<TARGET_FILE:dec3d_contract_hypre_diffusion_solver>)
    list(APPEND _hypre_mpi_reject_command --expect-world-reject)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _hypre_mpi_reject_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_hypre_diffusion_solver_world_reject
      COMMAND ${_hypre_mpi_reject_command}
    )

    add_executable(dec3d_contract_hypre_distributed_diffusion_solver
      tests/contract/test_hypre_distributed_diffusion_solver.cpp
    )
    target_include_directories(dec3d_contract_hypre_distributed_diffusion_solver
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_hypre_distributed_diffusion_solver
      PRIVATE
        dec3d_transport_hypre
        dec3d_transport
        MPI::MPI_CXX
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(dec3d_contract_hypre_distributed_diffusion_solver)

    set(_hypre_distributed_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _hypre_distributed_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _hypre_distributed_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _hypre_distributed_command $<TARGET_FILE:dec3d_contract_hypre_distributed_diffusion_solver>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _hypre_distributed_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_hypre_distributed_diffusion_solver
      COMMAND ${_hypre_distributed_command}
    )

    add_executable(dec3d_contract_thermal_conduction_hypre
      tests/contract/test_thermal_conduction_hypre.cpp
    )
    target_include_directories(dec3d_contract_thermal_conduction_hypre
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_thermal_conduction_hypre
      PRIVATE
        dec3d_transport_hypre
        dec3d_transport
        MPI::MPI_CXX
        dec3d_state
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(dec3d_contract_thermal_conduction_hypre)
    add_test(
      NAME dec3d_contract_thermal_conduction_hypre
      COMMAND dec3d_contract_thermal_conduction_hypre
    )

    add_executable(dec3d_contract_distributed_thermal_conduction
      tests/contract/test_distributed_thermal_conduction.cpp
    )
    target_include_directories(dec3d_contract_distributed_thermal_conduction
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_distributed_thermal_conduction
      PRIVATE
        dec3d_transport_hypre
        dec3d_transport
        MPI::MPI_CXX
        dec3d_state
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(dec3d_contract_distributed_thermal_conduction)
    set(_distributed_thermal_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _distributed_thermal_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _distributed_thermal_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _distributed_thermal_command $<TARGET_FILE:dec3d_contract_distributed_thermal_conduction>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _distributed_thermal_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_distributed_thermal_conduction
      COMMAND ${_distributed_thermal_command}
    )

    add_executable(dec3d_contract_radiation_distributed_provider_fed_multigroup
      tests/contract/test_radiation_distributed_provider_fed_multigroup.cpp
    )
    target_include_directories(dec3d_contract_radiation_distributed_provider_fed_multigroup
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_radiation_distributed_provider_fed_multigroup
      PRIVATE
        dec3d_radiation_hypre
        dec3d_radiation
        dec3d_transport_hypre
        dec3d_transport
        MPI::MPI_CXX
        dec3d_state
        dec3d_mesh
        dec3d_core
    )
    dec3d_apply_warnings(dec3d_contract_radiation_distributed_provider_fed_multigroup)
    set(_distributed_radiation_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _distributed_radiation_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _distributed_radiation_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _distributed_radiation_command $<TARGET_FILE:dec3d_contract_radiation_distributed_provider_fed_multigroup>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _distributed_radiation_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_radiation_distributed_provider_fed_multigroup
      COMMAND ${_distributed_radiation_command}
    )

    add_executable(dec3d_contract_alpha_distributed_one_group_operator
      tests/contract/test_alpha_distributed_one_group_operator.cpp
    )
    target_include_directories(dec3d_contract_alpha_distributed_one_group_operator
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_alpha_distributed_one_group_operator
      PRIVATE
        dec3d_alpha_transport
        dec3d_alpha
        dec3d_transport_hypre
        dec3d_transport
        dec3d_state
        dec3d_mesh
        MPI::MPI_CXX
    )
    dec3d_apply_warnings(dec3d_contract_alpha_distributed_one_group_operator)

    set(_alpha_distributed_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _alpha_distributed_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _alpha_distributed_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _alpha_distributed_command $<TARGET_FILE:dec3d_contract_alpha_distributed_one_group_operator>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _alpha_distributed_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_alpha_distributed_one_group_operator
      COMMAND ${_alpha_distributed_command}
    )

    add_executable(dec3d_contract_p4_alpha_benchmark_hypre
      tests/contract/test_p4_alpha_benchmark_hypre.cpp
    )
    target_include_directories(dec3d_contract_p4_alpha_benchmark_hypre
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
    )
    target_link_libraries(dec3d_contract_p4_alpha_benchmark_hypre
      PRIVATE
        dec3d_benchmarks_hypre
        dec3d_benchmarks
        dec3d_alpha_transport
        dec3d_alpha
        dec3d_transport_hypre
        dec3d_transport
        dec3d_state
        dec3d_mesh
        MPI::MPI_CXX
    )
    dec3d_apply_warnings(dec3d_contract_p4_alpha_benchmark_hypre)
    set(_p4_alpha_benchmark_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _p4_alpha_benchmark_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _p4_alpha_benchmark_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _p4_alpha_benchmark_command $<TARGET_FILE:dec3d_contract_p4_alpha_benchmark_hypre>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _p4_alpha_benchmark_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_p4_alpha_benchmark_hypre
      COMMAND ${_p4_alpha_benchmark_command}
    )

    add_executable(dec3d_contract_p5_integrated_benchmark_hypre
      tests/contract/test_p5_integrated_benchmark_hypre.cpp
    )
    target_include_directories(dec3d_contract_p5_integrated_benchmark_hypre
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(dec3d_contract_p5_integrated_benchmark_hypre
      PRIVATE
        dec3d_benchmarks_hypre
        dec3d_benchmarks
        dec3d_runtime
        dec3d_hydro
        dec3d_radiation_hypre
        dec3d_alpha_transport
        dec3d_transport_hypre
        dec3d_state
        dec3d_mesh
        dec3d_core
        MPI::MPI_CXX
        ${DEC3D_HYPRE_TARGET}
    )
    dec3d_apply_warnings(dec3d_contract_p5_integrated_benchmark_hypre)

    set(_p5_integrated_benchmark_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _p5_integrated_benchmark_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _p5_integrated_benchmark_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _p5_integrated_benchmark_command $<TARGET_FILE:dec3d_contract_p5_integrated_benchmark_hypre>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _p5_integrated_benchmark_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_p5_integrated_benchmark_hypre
      COMMAND ${_p5_integrated_benchmark_command}
    )

    add_executable(dec3d_contract_p5_runtime_distributed_all_stages
      tests/contract/test_p5_runtime_distributed_all_stages.cpp
    )
    target_include_directories(dec3d_contract_p5_runtime_distributed_all_stages
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(dec3d_contract_p5_runtime_distributed_all_stages
      PRIVATE
        dec3d_app
        dec3d_radiation_hypre
        dec3d_alpha_transport
        dec3d_transport_hypre
        dec3d_transport
        dec3d_state
        dec3d_mesh
        dec3d_core
        MPI::MPI_CXX
        ${DEC3D_HYPRE_TARGET}
    )
    dec3d_apply_warnings(dec3d_contract_p5_runtime_distributed_all_stages)

    set(_p5_runtime_distributed_command ${MPIEXEC_EXECUTABLE})
    if(MPIEXEC_NUMPROC_FLAG)
      list(APPEND _p5_runtime_distributed_command ${MPIEXEC_NUMPROC_FLAG} 2)
    endif()
    if(MPIEXEC_PREFLAGS)
      list(APPEND _p5_runtime_distributed_command ${MPIEXEC_PREFLAGS})
    endif()
    list(APPEND _p5_runtime_distributed_command $<TARGET_FILE:dec3d_contract_p5_runtime_distributed_all_stages>)
    if(MPIEXEC_POSTFLAGS)
      list(APPEND _p5_runtime_distributed_command ${MPIEXEC_POSTFLAGS})
    endif()
    add_test(
      NAME dec3d_contract_p5_runtime_distributed_all_stages
      COMMAND ${_p5_runtime_distributed_command}
    )

    add_executable(dec3d_contract_p2_thermal_benchmark_hypre
      tests/contract/test_p2_thermal_benchmark_hypre.cpp
    )
    target_include_directories(dec3d_contract_p2_thermal_benchmark_hypre
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/contract
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_support
    )
    target_link_libraries(dec3d_contract_p2_thermal_benchmark_hypre
      PRIVATE
        dec3d_benchmarks_hypre
        dec3d_benchmarks
        MPI::MPI_CXX
    )
    dec3d_apply_warnings(dec3d_contract_p2_thermal_benchmark_hypre)
    set(_p2_thermal_benchmark_command
      ${MPIEXEC_EXECUTABLE}
      ${MPIEXEC_NUMPROC_FLAG}
      2
      ${MPIEXEC_PREFLAGS}
    )
    list(APPEND _p2_thermal_benchmark_command $<TARGET_FILE:dec3d_contract_p2_thermal_benchmark_hypre>)
    list(APPEND _p2_thermal_benchmark_command ${MPIEXEC_POSTFLAGS})
    add_test(
      NAME dec3d_contract_p2_thermal_benchmark_hypre
      COMMAND ${_p2_thermal_benchmark_command}
    )

  endif()

  dec3d_add_test(dec3d_contract_canonical_state_layout
    tests/contract/test_canonical_state_layout.cpp
  )

  dec3d_add_test(dec3d_contract_canonical_state_utilities
    tests/contract/test_canonical_state_utilities.cpp
  )

  dec3d_add_test(dec3d_contract_state_restart_boundary
    tests/contract/test_state_restart_boundary.cpp
  )

  dec3d_add_test(dec3d_contract_runtime_stage_registry
    tests/contract/test_runtime_stage_registry.cpp
  )

  dec3d_add_test(dec3d_contract_runtime_required_stage_validation
    tests/contract/test_runtime_required_stage_validation.cpp
  )

  dec3d_add_test(dec3d_contract_runtime_stage_machine
    tests/contract/test_runtime_stage_machine.cpp
  )

  dec3d_add_test(dec3d_contract_runtime_stage_report_aggregation
    tests/contract/test_runtime_stage_report_aggregation.cpp
  )

  dec3d_add_test(dec3d_contract_runtime_stage_result_ingestion
    tests/contract/test_runtime_stage_result_ingestion.cpp
  )

  dec3d_add_test(dec3d_contract_checkpoint_contract
    tests/contract/test_checkpoint_contract.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_view
    tests/contract/test_hydro_view.cpp
  )

  dec3d_add_test(dec3d_contract_hllc_solver
    tests/contract/test_hllc_solver.cpp
  )
  target_link_libraries(dec3d_contract_hllc_solver PRIVATE dec3d_radiation)

  dec3d_add_test(dec3d_contract_hllc_moving_interface
    tests/contract/test_hllc_moving_interface.cpp
  )

  dec3d_add_test(dec3d_contract_ppm_reconstruction
    tests/contract/test_ppm_reconstruction.cpp
  )

  dec3d_add_test(dec3d_contract_ppm_shell_wide_fallback_mask
    tests/contract/test_ppm_shell_wide_fallback_mask.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_geometric_source
    tests/contract/test_hydro_geometric_source.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_theta_pole_boundary
    tests/contract/test_hydro_theta_pole_boundary.cpp
  )

  dec3d_add_test(dec3d_contract_radial_ale_runtime_commit
    tests/contract/test_radial_ale_runtime_commit.cpp
  )

  dec3d_add_test(dec3d_contract_radial_ale_global_proposal
    tests/contract/test_radial_ale_global_proposal.cpp
  )

  dec3d_add_test(dec3d_contract_radial_ale_flux
    tests/contract/test_radial_ale_flux.cpp
  )

  dec3d_add_test(dec3d_contract_radial_ale_extensive_update
    tests/contract/test_radial_ale_extensive_update.cpp
  )

  dec3d_add_test(dec3d_contract_macro_ale_radial_overlap_remap
    tests/contract/test_macro_ale_radial_overlap_remap.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_ale_single_rank
    tests/contract/test_hydro_macro_ale_single_rank.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_ale_staged_writeback
    tests/contract/test_hydro_macro_ale_staged_writeback.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_ale_direct_hllc_no_fallback
    tests/contract/test_hydro_macro_ale_direct_hllc_no_fallback.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_macro_ale_direct_hllc_single_rank
    tests/contract/test_hydro_macro_ale_direct_hllc_single_rank.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_radial_wave_ale_on_off
    tests/contract/test_hydro_radial_wave_ale_on_off.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_radial_wave_ale_ppm
    tests/contract/test_hydro_radial_wave_ale_ppm.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_radial_ppm_shell_wide_fallback
    tests/contract/test_hydro_radial_ppm_shell_wide_fallback.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_phi_periodic_ghost_fill
    tests/contract/test_hydro_phi_periodic_ghost_fill.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_theta_pole_ghost_fill_ng3
    tests/contract/test_hydro_theta_pole_ghost_fill_ng3.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_direction_ready_ghosts
    tests/contract/test_hydro_direction_ready_ghosts.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_spherical_radial_wave_sanity
    tests/contract/test_hydro_spherical_radial_wave_sanity.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_radial_wave_mpi_parity
    tests/contract/test_hydro_radial_wave_mpi_parity.cpp
  )

  if(MPI_CXX_FOUND)
    dec3d_add_mpi_test(dec3d_contract_hydro_halo_exchange 3
      tests/contract/test_hydro_halo_exchange.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_radial_ale_mpi_global_proposal 3
      tests/contract/test_radial_ale_mpi_global_proposal.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_radial_wave_mpi_real 3
      tests/contract/test_hydro_radial_wave_mpi_real.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_radial_wave_mpi_real_ppm 3
      tests/contract/test_hydro_radial_wave_mpi_real_ppm.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_true3d_low_mode_mpi_real 3
      tests/contract/test_hydro_true3d_low_mode_mpi_real.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_true3d_low_mode_mpi_real_ppm 3
      tests/contract/test_hydro_true3d_low_mode_mpi_real_ppm.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_true3d_low_mode_mpi_real_macro 3
      tests/contract/test_hydro_true3d_low_mode_mpi_real_macro.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_pole_seam_stress_mpi_real_macro 3
      tests/contract/test_hydro_pole_seam_stress_mpi_real_macro.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_pole_seam_stress_mpi_real_macro_ppm 3
      tests/contract/test_hydro_pole_seam_stress_mpi_real_macro_ppm.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_pole_seam_stress_mpi_real 3
      tests/contract/test_hydro_pole_seam_stress_mpi_real.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_macro_zoning_mpi_parity 3
      tests/contract/test_hydro_macro_zoning_mpi_parity.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_macro_zoning_mpi_real_ppm 3
      tests/contract/test_hydro_macro_zoning_mpi_real_ppm.cpp
    )

    dec3d_add_mpi_test(dec3d_contract_hydro_macro_zoning_mpi_coarse_ghost_bootstrap 3
      tests/contract/test_hydro_macro_zoning_mpi_coarse_ghost_bootstrap.cpp
    )

  dec3d_add_mpi_test(dec3d_contract_hydro_macro_zoning_mpi_transaction 3
    tests/contract/test_hydro_macro_zoning_mpi_transaction.cpp
  )

  dec3d_add_mpi_test(dec3d_contract_hydro_macro_ale_mpi_compatibility 3
    tests/contract/test_hydro_macro_ale_mpi_compatibility.cpp
  )

  dec3d_add_mpi_test(dec3d_contract_hydro_macro_ale_direct_hllc_mpi 3
    tests/contract/test_hydro_macro_ale_direct_hllc_mpi.cpp
  )

  dec3d_add_mpi_test(dec3d_contract_hydro_macro_ale_mpi_smoke 3
    tests/contract/test_hydro_macro_ale_mpi_smoke.cpp
  )
endif()

  dec3d_add_test(dec3d_contract_hydro_budget_residual
    tests/contract/test_hydro_budget_residual.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_dt_advice
    tests/contract/test_hydro_dt_advice.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_true3d_low_mode_case
    tests/contract/test_hydro_true3d_low_mode_case.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_pole_seam_stress_case
    tests/contract/test_hydro_pole_seam_stress_case.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_radial_boundary_contract
    tests/contract/test_hydro_radial_boundary_contract.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_origin_radial_ghost_remap
    tests/contract/test_hydro_origin_radial_ghost_remap.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_uniform_state_preservation
    tests/contract/test_hydro_uniform_state_preservation.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_uniform_state_flux_only
    tests/contract/test_hydro_uniform_state_flux_only.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_directional_distinguishability
    tests/contract/test_hydro_directional_distinguishability.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_pole_stability
    tests/contract/test_hydro_pole_stability.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_operator
    tests/contract/test_hydro_operator_contract.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_operator_ppm_reconstruction
    tests/contract/test_hydro_operator_ppm_reconstruction.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_operator_transaction
    tests/contract/test_hydro_operator_transaction.cpp
  )

  dec3d_add_test(dec3d_contract_hydro_directional_static_grid
    tests/contract/test_hydro_directional_static_grid.cpp
  )

  dec3d_add_test(dec3d_acceptance_phase_p0
    tests/acceptance/test_phase_p0_acceptance.cpp
  )

  dec3d_add_test(dec3d_acceptance_phase_p1
    tests/acceptance/test_phase_p1_acceptance.cpp
  )
endif()
