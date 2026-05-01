#include "alpha/distributed_alpha_operator.hpp"
#include "alpha/alpha_operator.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void RequireTwoRanks() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    dec3d::test::Fail("rank count", __FILE__, __LINE__, "requires mpiexec -n 2");
  }
}

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    dec3d::test::Fail(label,
                      __FILE__,
                      __LINE__,
                      std::to_string(actual) + " vs " + std::to_string(expected));
  }
}

dec3d::mesh::SphericalGeometryMetadata MakePatchGeometry(
    std::size_t radial_cells,
    std::size_t theta_cells,
    std::size_t phi_cells) {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces.resize(radial_cells + 1u);
  geometry.theta_faces.resize(theta_cells + 1u);
  geometry.phi_faces.resize(phi_cells + 1u);
  geometry.cell_volumes.resize(radial_cells * theta_cells * phi_cells, 0.0);

  for (std::size_t r = 0; r <= radial_cells; ++r) {
    geometry.radial_faces[r] =
        0.1 + 0.9 * static_cast<double>(r) / static_cast<double>(radial_cells);
  }
  for (std::size_t t = 0; t <= theta_cells; ++t) {
    geometry.theta_faces[t] =
        0.4 + 0.8 * static_cast<double>(t) / static_cast<double>(theta_cells);
  }
  for (std::size_t p = 0; p <= phi_cells; ++p) {
    geometry.phi_faces[p] =
        2.0 * kPi * static_cast<double>(p) / static_cast<double>(phi_cells);
  }

  std::size_t linear = 0u;
  geometry.global_volume = 0.0;
  for (std::size_t r = 0; r < radial_cells; ++r) {
    const double radial =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) /
        3.0;
    for (std::size_t t = 0; t < theta_cells; ++t) {
      const double polar =
          std::cos(geometry.theta_faces[t]) -
          std::cos(geometry.theta_faces[t + 1u]);
      for (std::size_t p = 0; p < phi_cells; ++p) {
        const double azimuth =
            geometry.phi_faces[p + 1u] - geometry.phi_faces[p];
        geometry.cell_volumes[linear] = radial * polar * azimuth;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

dec3d::transport::GenericDiffusionBoundaryPolicy PatchBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::interior_patch_no_origin,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

void FillAlphaLocalState(
    dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    bool invalid_for_provider = false) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    const std::size_t global_r = ownership.global_radial_begin + r;
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double rho = invalid_for_provider ? -1.0 : 50.0 + static_cast<double>(global_r);
        const double te = dec3d::physics::ErgFromKeV(2.0);
        const double ti = dec3d::physics::ErgFromKeV(3.0);
        const double ni = rho / dec3d::physics::DefaultMeanDTIonMassG();
        const double ne = ni;
        state.rho(r, t, p) = rho;
        state.mom_r(r, t, p) = 0.0;
        state.mom_theta(r, t, p) = 0.0;
        state.mom_phi(r, t, p) = 0.0;
        state.e_electron(r, t, p) = ne * te / gamma_minus_one;
        state.e_fluid_total(r, t, p) =
            state.e_electron(r, t, p) + ni * ti / gamma_minus_one;
        state.alpha_state.storage(r, t, p) =
            1.0e10 + 1.0e9 * static_cast<double>(global_r);
      }
    }
  }
}

struct AlphaProblemFixture {
  dec3d::state::CanonicalState state;
  dec3d::alpha::DistributedOneGroupAlphaTransportProblem problem;
};

void InitializeProblem(
    AlphaProblemFixture& fixture,
    double dt_s,
    bool rank_one_invalid_for_provider = false) {
  fixture.problem.communicator = MPI_COMM_WORLD;
  fixture.problem.ownership =
      dec3d::transport::BuildDistributedDiffusionRowOwnership(MPI_COMM_WORLD, 4, 2, 2);
  DEC3D_CHECK(fixture.problem.ownership.success);
  const std::size_t local_radial =
      fixture.problem.ownership.global_radial_end -
      fixture.problem.ownership.global_radial_begin;
  fixture.state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{local_radial, 2, 2, 0});
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  FillAlphaLocalState(
      fixture.state,
      fixture.problem.ownership,
      rank_one_invalid_for_provider && rank == 1);
  fixture.problem.local_state = &fixture.state;
  fixture.problem.global_geometry = MakePatchGeometry(4, 2, 2);
  fixture.problem.dt_s = dt_s;
}

dec3d::alpha::DistributedOneGroupAlphaTransportOptions MakeOptions() {
  dec3d::alpha::DistributedOneGroupAlphaTransportOptions options;
  options.provider_options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.provider_options.tau_model =
      dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.provider_options.reactivity_model =
      dec3d::alpha::AlphaReactivityModel::bosch_hale_dt;
  options.boundary_policy = PatchBoundary();
  options.solve_options.relative_tolerance = 1.0e-8;
  options.solve_options.max_iterations = 100;
  return options;
}

dec3d::state::CanonicalState MakeSerialState() {
  dec3d::transport::DistributedDiffusionRowOwnership ownership;
  ownership.global_radial_begin = 0;
  ownership.global_radial_end = 4;
  ownership.global_theta_cells = 2;
  ownership.global_phi_cells = 2;
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{4, 2, 2, 0});
  FillAlphaLocalState(state, ownership);
  return state;
}

dec3d::alpha::OneGroupAlphaTransportOptions MakeSerialOptions(double dt_s) {
  dec3d::alpha::OneGroupAlphaTransportOptions options;
  options.dt_s = dt_s;
  options.provider_options.composition_model =
      dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
  options.provider_options.tau_model =
      dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
  options.provider_options.reactivity_model =
      dec3d::alpha::AlphaReactivityModel::bosch_hale_dt;
  options.boundary_policy = PatchBoundary();
  options.serial_reference_options.residual_tolerance = 1.0e20;
  return options;
}

void CheckAlphaStateUnchanged(
    const dec3d::state::CanonicalState& state,
    const dec3d::core::Array3D<double>& before_alpha,
    const dec3d::core::Array3D<double>& before_electron,
    const dec3d::core::Array3D<double>& before_total) {
  DEC3D_CHECK_EQ(state.last_authoritative_write_mask, 0u);
  for (std::size_t r = 0; r < before_alpha.extent_r(); ++r) {
    for (std::size_t t = 0; t < before_alpha.extent_theta(); ++t) {
      for (std::size_t p = 0; p < before_alpha.extent_phi(); ++p) {
        DEC3D_CHECK_EQ(state.alpha_state.storage(r, t, p), before_alpha(r, t, p));
        DEC3D_CHECK_EQ(state.e_electron(r, t, p), before_electron(r, t, p));
        DEC3D_CHECK_EQ(state.e_fluid_total(r, t, p), before_total(r, t, p));
      }
    }
  }
}

void RequireDistributedSuccess(
    const dec3d::alpha::DistributedOneGroupAlphaTransportResult& result) {
  if (!result.success) {
    dec3d::test::Fail("distributed alpha success",
                      __FILE__,
                      __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    RequireTwoRanks();

    {
      AlphaProblemFixture fixture;
      InitializeProblem(fixture, 0.0);
      auto& problem = fixture.problem;
      const auto before = fixture.state.alpha_state.storage;
      const auto result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(problem, MakeOptions());
      DEC3D_CHECK(result.success);
      DEC3D_CHECK(dec3d::alpha::ValidateDistributedOneGroupAlphaDiagnostics(result));
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK(!result.metadata_written);
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::alpha::kDistributedAlphaWritesNoFields);
      DEC3D_CHECK_EQ(result.solve_count, std::size_t{0});
      DEC3D_CHECK(!result.distributed_alpha_solve_executed);
      DEC3D_CHECK(!result.hypre_solve_executed);
      DEC3D_CHECK(result.report_line.find("provider_skipped_for_no_change=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_group_solve_count") == std::string::npos);
      for (std::size_t r = 0; r < before.extent_r(); ++r) {
        for (std::size_t t = 0; t < before.extent_theta(); ++t) {
          for (std::size_t p = 0; p < before.extent_phi(); ++p) {
            DEC3D_CHECK_EQ(fixture.state.alpha_state.storage(r, t, p), before(r, t, p));
          }
        }
      }
    }

    {
      AlphaProblemFixture fixture;
      InitializeProblem(fixture, 1.0e-18);
      const auto before_alpha = fixture.state.alpha_state.storage;
      const auto before_electron = fixture.state.e_electron;
      const auto before_total = fixture.state.e_fluid_total;
      const auto result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(
              fixture.problem, MakeOptions());
      RequireDistributedSuccess(result);
      DEC3D_CHECK(dec3d::alpha::ValidateDistributedOneGroupAlphaDiagnostics(result));
      DEC3D_CHECK(result.canonical_state_mutated);
      DEC3D_CHECK(result.metadata_written);
      DEC3D_CHECK(result.local_stage_ok);
      DEC3D_CHECK(result.global_stage_ok);
      DEC3D_CHECK(result.local_publishable);
      DEC3D_CHECK(result.global_publish_ok);
      DEC3D_CHECK_EQ(result.solve_count, std::size_t{1});
      DEC3D_CHECK(result.distributed_alpha_solve_executed);
      DEC3D_CHECK(result.hypre_solve_executed);
      DEC3D_CHECK(result.seam_coefficient_halo_exchanged);
      DEC3D_CHECK(result.off_rank_face_conductance_uses_neighbor_Dalpha);
      DEC3D_CHECK(result.global_off_rank_column_count > 0u);
      DEC3D_CHECK(result.report_line.find("coefficient_scope=owned_slab") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("owned_slab_writeback_only=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("backend_executed=hypre_parcsr_gmres_boomeramg") !=
                  std::string::npos);
      DEC3D_CHECK(result.updated_fields ==
                  (dec3d::alpha::kDistributedAlphaWritesAlphaState |
                   dec3d::alpha::kDistributedAlphaWritesElectronEnergy |
                   dec3d::alpha::kDistributedAlphaWritesFluidTotalEnergy));
      DEC3D_CHECK(std::abs(result.alpha_equation_budget_residual_global) <=
                  1.0e-6 * std::max(1.0, std::abs(result.birth_source_total_global)));
      DEC3D_CHECK(std::abs(result.alpha_electron_exchange_residual_global) <=
                  1.0e-6 * std::max(1.0, std::abs(result.drag_deposition_total_global)));
      bool any_alpha_changed = false;
      bool any_electron_changed = false;
      for (std::size_t r = 0; r < before_alpha.extent_r(); ++r) {
        for (std::size_t t = 0; t < before_alpha.extent_theta(); ++t) {
          for (std::size_t p = 0; p < before_alpha.extent_phi(); ++p) {
            any_alpha_changed =
                any_alpha_changed ||
                fixture.state.alpha_state.storage(r, t, p) != before_alpha(r, t, p);
            any_electron_changed =
                any_electron_changed ||
                fixture.state.e_electron(r, t, p) != before_electron(r, t, p);
            DEC3D_CHECK(fixture.state.e_fluid_total(r, t, p) != before_total(r, t, p));
          }
        }
      }
      DEC3D_CHECK(any_alpha_changed);
      DEC3D_CHECK(any_electron_changed);
    }

    {
      constexpr double kDt = 1.0e-20;
      AlphaProblemFixture distributed;
      InitializeProblem(distributed, kDt);
      auto serial = MakeSerialState();

      auto serial_options = MakeSerialOptions(kDt);
      const auto serial_result =
          dec3d::alpha::ApplyOneGroupAlphaTransport(
              serial, distributed.problem.global_geometry, serial_options);
      if (!serial_result.success) {
        dec3d::test::Fail("serial alpha baseline success",
                          __FILE__,
                          __LINE__,
                          serial_result.failure_reason + " | " +
                              serial_result.failure_diagnostics);
      }

      auto distributed_options = MakeOptions();
      distributed_options.solve_options.relative_tolerance = 1.0e-10;
      distributed_options.solve_options.max_iterations = 200;
      const auto distributed_result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(
              distributed.problem, distributed_options);
      RequireDistributedSuccess(distributed_result);
      for (std::size_t r = 0; r < distributed.state.layout.radial_cells; ++r) {
        const std::size_t global_r =
            distributed.problem.ownership.global_radial_begin + r;
        for (std::size_t t = 0; t < distributed.state.layout.theta_cells; ++t) {
          for (std::size_t p = 0; p < distributed.state.layout.phi_cells; ++p) {
            CheckNear(distributed.state.alpha_state.storage(r, t, p),
                      serial.alpha_state.storage(global_r, t, p),
                      1.0e-3,
                      "distributed/serial alpha_state parity");
            CheckNear(distributed.state.e_electron(r, t, p),
                      serial.e_electron(global_r, t, p),
                      1.0e4,
                      "distributed/serial electron parity");
            CheckNear(distributed.state.e_fluid_total(r, t, p),
                      serial.e_fluid_total(global_r, t, p),
                      1.0e4,
                      "distributed/serial total parity");
          }
        }
      }
    }

    {
      AlphaProblemFixture fixture;
      InitializeProblem(fixture, 1.0e-18, true);
      const auto before_alpha = fixture.state.alpha_state.storage;
      const auto before_electron = fixture.state.e_electron;
      const auto before_total = fixture.state.e_fluid_total;
      const auto result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(
              fixture.problem, MakeOptions());
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!result.global_stage_ok);
      DEC3D_CHECK(!result.global_publish_ok);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::alpha::kDistributedAlphaWritesNoFields);
      DEC3D_CHECK(result.failure_diagnostics.find("first_failing_rank=1") !=
                  std::string::npos);
      CheckAlphaStateUnchanged(fixture.state, before_alpha, before_electron, before_total);
    }

    {
      AlphaProblemFixture fixture;
      InitializeProblem(fixture, 1.0e-18);
      const auto before_alpha = fixture.state.alpha_state.storage;
      const auto before_electron = fixture.state.e_electron;
      const auto before_total = fixture.state.e_fluid_total;
      auto options = MakeOptions();
      options.force_publish_preflight_failure_for_test = true;
      const auto result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(fixture.problem, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(result.global_stage_ok);
      DEC3D_CHECK(!result.global_publish_ok);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::alpha::kDistributedAlphaWritesNoFields);
      DEC3D_CHECK(result.failure_diagnostics.find("publish preflight") !=
                  std::string::npos);
      CheckAlphaStateUnchanged(fixture.state, before_alpha, before_electron, before_total);
    }

    {
      AlphaProblemFixture fixture;
      InitializeProblem(fixture, 1.0e-18);
      auto result =
          dec3d::alpha::ApplyDistributedOneGroupAlphaTransport(
              fixture.problem, MakeOptions());
      RequireDistributedSuccess(result);
      DEC3D_CHECK(dec3d::alpha::ValidateDistributedOneGroupAlphaDiagnostics(result));
      const std::string token = "off_rank_face_conductance_uses_neighbor_Dalpha=";
      const auto pos = result.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      result.report_line.erase(pos, token.size());
      DEC3D_CHECK(!dec3d::alpha::ValidateDistributedOneGroupAlphaDiagnostics(result));
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
}
