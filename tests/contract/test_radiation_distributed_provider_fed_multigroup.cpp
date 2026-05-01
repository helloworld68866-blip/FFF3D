#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/distributed_multigroup_radiation.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "test_assert.hpp"
#include "transport/diffusion/hypre_distributed_diffusion_solver.hpp"

#include <mpi.h>

#include <cmath>
#include <exception>
#include <iostream>
#include <string>

namespace {

constexpr const char* kTableRoot = "F:/dec3d/data/opacities/tops_dt_2026_04_27";
constexpr double kPi = 3.141592653589793238462643383279502884;

void RequireTwoRanks() {
  int size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 2) {
    dec3d::test::Fail("rank count", __FILE__, __LINE__, "requires mpiexec -n 2");
  }
}

double PlanckConstantErgS() {
  return 2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
         dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
}

dec3d::radiation::RadiationGroupLayout MakeTwoGroupLayout() {
  const double kev_to_hz =
      dec3d::physics::PhysicsConstantsCGS::erg_per_kev / PlanckConstantErgS();
  return dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
      {0.002 * kev_to_hz, 0.02 * kev_to_hz, 0.2 * kev_to_hz});
}

dec3d::radiation::RadiationGroupLayout MakeOneGroupLayout() {
  const double kev_to_hz =
      dec3d::physics::PhysicsConstantsCGS::erg_per_kev / PlanckConstantErgS();
  return dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
      {0.002 * kev_to_hz, 0.2 * kev_to_hz});
}

dec3d::radiation::TopsOpacityProviderOptions MakeProviderOptions() {
  dec3d::radiation::TopsOpacityProviderOptions options;
  options.table_root = kTableRoot;
  options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  options.density_clip_policy = dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return options;
}

dec3d::mesh::SphericalGeometryMetadata MakeFullSphereGeometry(
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
    geometry.radial_faces[r] = 0.1 + 0.9 * static_cast<double>(r) /
                                         static_cast<double>(radial_cells);
  }
  for (std::size_t t = 0; t <= theta_cells; ++t) {
    geometry.theta_faces[t] = 0.3 + 1.0 * static_cast<double>(t) /
                                        static_cast<double>(theta_cells);
  }
  for (std::size_t p = 0; p <= phi_cells; ++p) {
    geometry.phi_faces[p] = 2.0 * kPi * static_cast<double>(p) /
                            static_cast<double>(phi_cells);
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
        const double azimuth = geometry.phi_faces[p + 1u] - geometry.phi_faces[p];
        geometry.cell_volumes[linear] = radial * polar * azimuth;
        geometry.global_volume += geometry.cell_volumes[linear];
        ++linear;
      }
    }
  }
  return geometry;
}

dec3d::transport::GenericDiffusionBoundaryPolicy MakePatchBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::interior_patch_no_origin,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

void FillLocalState(
    dec3d::state::CanonicalState& state,
    const dec3d::transport::DistributedDiffusionRowOwnership& ownership,
    bool invalid_for_provider = false) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    const std::size_t global_r = ownership.global_radial_begin + r;
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double rho = invalid_for_provider ? 1.0e30 : 1.0 + 0.25 * static_cast<double>(global_r);
        const double te_kev = 0.5 + 0.05 * static_cast<double>(global_r);
        const double ti_kev = 0.4;
        const double te = dec3d::physics::ErgFromKeV(te_kev);
        const double ti = dec3d::physics::ErgFromKeV(ti_kev);
        const double ni = rho / dec3d::physics::DefaultMeanDTIonMassG();
        const double ne = ni;
        state.rho(r, t, p) = rho;
        state.mom_r(r, t, p) = 0.0;
        state.mom_theta(r, t, p) = 0.0;
        state.mom_phi(r, t, p) = 0.0;
        state.e_electron(r, t, p) = ne * te / gamma_minus_one;
        state.e_fluid_total(r, t, p) =
            state.e_electron(r, t, p) + ni * ti / gamma_minus_one;
        for (std::size_t g = 0; g < state.radiation_groups.size(); ++g) {
          state.radiation_groups[g](r, t, p) =
              1.0e8 + 1.0e7 * static_cast<double>(g) +
              1.0e6 * static_cast<double>(global_r);
        }
      }
    }
  }
}

dec3d::state::CanonicalState MakeSerialState(
    const dec3d::radiation::RadiationGroupLayout& group_layout) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{4, 2, 2, group_layout.group_count});
  dec3d::transport::DistributedDiffusionRowOwnership ownership;
  ownership.global_radial_begin = 0;
  ownership.global_radial_end = 4;
  ownership.global_theta_cells = 2;
  ownership.global_phi_cells = 2;
  FillLocalState(state, ownership);
  return state;
}

dec3d::radiation::DistributedMultigroupRadiationProblem MakeProblem(
    MPI_Comm communicator,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const dec3d::radiation::TopsOpacityTable& table,
    double dt_s,
    bool invalid_rank_one = false) {
  dec3d::radiation::DistributedMultigroupRadiationProblem problem;
  problem.ownership =
      dec3d::transport::BuildDistributedDiffusionRowOwnership(communicator, 4, 2, 2);
  const std::size_t local_radial =
      problem.ownership.global_radial_end - problem.ownership.global_radial_begin;
  problem.local_state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{local_radial, 2, 2, group_layout.group_count});
  int rank = 0;
  MPI_Comm_rank(communicator, &rank);
  FillLocalState(problem.local_state, problem.ownership, invalid_rank_one && rank == 1);
  problem.global_geometry = MakeFullSphereGeometry(4, 2, 2);
  problem.boundary_policy = MakePatchBoundary();
  problem.group_layout = group_layout;
  problem.opacity_table = &table;
  problem.dt_s = dt_s;
  return problem;
}

dec3d::radiation::DistributedMultigroupRadiationOptions MakeOptions() {
  dec3d::radiation::DistributedMultigroupRadiationOptions options;
  options.provider_options = MakeProviderOptions();
  options.boundary_model =
      dec3d::radiation::DistributedRadiationBoundaryModel::contract_zero_flux_or_scalar_remap;
  options.solve_options.relative_tolerance = 1.0e-8;
  options.solve_options.max_iterations = 100;
  return options;
}

void CheckSuccess(const dec3d::radiation::DistributedMultigroupRadiationResult& result) {
  if (!result.success) {
    dec3d::test::Fail("distributed provider-fed radiation success",
                      __FILE__,
                      __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
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

}  // namespace

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  try {
    RequireTwoRanks();
    const auto table_result = dec3d::radiation::LoadTopsOpacityTable(MakeProviderOptions());
    DEC3D_CHECK(table_result.success);
    const auto group_layout = MakeTwoGroupLayout();

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 0.0);
      const auto before = problem.local_state;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, MakeOptions());
      CheckSuccess(result);
      DEC3D_CHECK_EQ(result.per_group_solve_count, std::size_t{0});
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::radiation::kDistributedRadiationWritesNoFields);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK(result.global_stage_ok);
      DEC3D_CHECK(result.global_publish_ok);
      DEC3D_CHECK_EQ(result.delta_radiation_total_all_groups, 0.0);
      DEC3D_CHECK_EQ(result.delta_electron_total, 0.0);
      DEC3D_CHECK(problem.local_state.last_authoritative_write_mask == 0u);
      DEC3D_CHECK(problem.local_state.e_electron(0, 0, 0) == before.e_electron(0, 0, 0));
      DEC3D_CHECK(result.report_line.find("table_cache_scope=explicit_handle") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("table_handle_supplied=true") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_step_csv_io=false") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_lookup_full_table_scan=false") !=
                  std::string::npos);
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, MakeOptions());
      CheckSuccess(result);
      DEC3D_CHECK_EQ(result.per_group_solve_count, group_layout.group_count);
      DEC3D_CHECK(result.updated_fields ==
                  (dec3d::radiation::kDistributedRadiationWritesRadiationGroups |
                   dec3d::radiation::kDistributedRadiationWritesElectronEnergy |
                   dec3d::radiation::kDistributedRadiationWritesFluidTotalEnergy));
      DEC3D_CHECK(result.canonical_state_mutated);
      DEC3D_CHECK(result.seam_coefficient_halo_exchanged);
      DEC3D_CHECK(result.off_rank_face_conductance_uses_neighbor_Dbar);
      DEC3D_CHECK(result.global_off_rank_column_count > 0);
      DEC3D_CHECK(result.report_line.find("provider_fed=true") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("backend_executed=hypre_parcsr_gmres_boomeramg") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_group_provider_report_present=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("seam_coefficient_halo_exchanged=true") !=
                  std::string::npos);
      DEC3D_CHECK_EQ(result.delta_radiation_total_by_group.size(), group_layout.group_count);
      DEC3D_CHECK_EQ(result.source_gain_radiation_total_by_group.size(), group_layout.group_count);
      DEC3D_CHECK_EQ(result.boundary_leak_total_by_group.size(), group_layout.group_count);
      double delta_sum = 0.0;
      double source_sum = 0.0;
      double leak_sum = 0.0;
      for (std::size_t g = 0; g < group_layout.group_count; ++g) {
        delta_sum += result.delta_radiation_total_by_group[g];
        source_sum += result.source_gain_radiation_total_by_group[g];
        leak_sum += result.boundary_leak_total_by_group[g];
      }
      CheckNear(delta_sum, result.delta_radiation_total_all_groups, 1.0e-8, "delta by-group sum");
      CheckNear(source_sum, result.source_gain_radiation_total_all_groups, 1.0e-8, "source by-group sum");
      CheckNear(leak_sum, result.boundary_leak_total_all_groups, 1.0e-8, "leak by-group sum");
      DEC3D_CHECK(result.report_line.find("per_group_budget_count=2") != std::string::npos);
      DEC3D_CHECK(result.report_line.find("thermodynamic_recovery_report_present=true") !=
                  std::string::npos);
      DEC3D_CHECK(dec3d::radiation::ValidateDistributedMultigroupRadiationDiagnostics(result));
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      auto options = MakeOptions();
      options.radiation_flux_limiter.model =
          dec3d::radiation::RadiationFluxLimiterModel::minmax_eq_5_211;
      options.radiation_flux_limiter.U_floor_erg_per_cm3 = 1.0e-30;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, options);
      CheckSuccess(result);
      DEC3D_CHECK(result.radiation_flux_limiter_enabled);
      DEC3D_CHECK(result.radiation_flux_limiter_model == "minmax_eq_5_211");
      DEC3D_CHECK(result.report_line.find("radiation_flux_limiter_enabled=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("radiation_flux_limiter_model=minmax_eq_5_211") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("face_limiter_time_level=old_time_lagged") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("U_floor_is_numerical_regularization=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_group_limiter_report_present=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("outer_marshak_uses_face_effective_D=") !=
                  std::string::npos);
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 0.0);
      auto options = MakeOptions();
      options.radiation_flux_limiter.model =
          dec3d::radiation::RadiationFluxLimiterModel::harmonic_eq_5_209;
      options.radiation_flux_limiter.U_floor_erg_per_cm3 = 1.0e-30;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, options);
      CheckSuccess(result);
      DEC3D_CHECK_EQ(result.per_group_solve_count, std::size_t{0});
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::radiation::kDistributedRadiationWritesNoFields);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK(result.report_line.find("radiation_flux_limiter_enabled=true") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_group_limiter_report_present=false") !=
                  std::string::npos);
    }

    {
      auto baseline = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      auto disabled = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      auto baseline_options = MakeOptions();
      auto disabled_options = MakeOptions();
      disabled_options.radiation_flux_limiter.model =
          dec3d::radiation::RadiationFluxLimiterModel::disabled;
      const auto baseline_result =
          dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
              baseline, baseline_options);
      const auto disabled_result =
          dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
              disabled, disabled_options);
      CheckSuccess(baseline_result);
      CheckSuccess(disabled_result);
      for (std::size_t g = 0; g < group_layout.group_count; ++g) {
        for (std::size_t r = 0; r < baseline.local_state.layout.radial_cells; ++r) {
          for (std::size_t t = 0; t < baseline.local_state.layout.theta_cells; ++t) {
            for (std::size_t p = 0; p < baseline.local_state.layout.phi_cells; ++p) {
              CheckNear(disabled.local_state.radiation_groups[g](r, t, p),
                        baseline.local_state.radiation_groups[g](r, t, p),
                        0.0,
                        "disabled limiter radiation parity");
            }
          }
        }
      }
      DEC3D_CHECK(disabled_result.report_line.find(
                      "radiation_flux_limiter_model=disabled") != std::string::npos);
      DEC3D_CHECK(disabled_result.report_line.find(
                      "radiation_flux_limiter_enabled=false") != std::string::npos);
    }

    {
      const auto one_group_layout = MakeOneGroupLayout();
      auto distributed =
          MakeProblem(MPI_COMM_WORLD, one_group_layout, table_result.table, 1.0e-24);
      auto serial = MakeSerialState(one_group_layout);
      auto provider = dec3d::radiation::BuildRadiationCoefficientArrays(
          serial,
          distributed.global_geometry,
          one_group_layout,
          table_result.table,
          MakeProviderOptions());
      DEC3D_CHECK(provider.success);

      dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions serial_options;
      serial_options.dt_s = distributed.dt_s;
      serial_options.group_layout = one_group_layout;
      serial_options.coefficients = provider.coefficients;
      serial_options.boundary_policy = MakePatchBoundary();
      serial_options.radiation_boundary_model =
          dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;
      serial_options.serial_reference_options.row_limit = 128;
      serial_options.serial_reference_options.residual_tolerance = 1.0e20;
      const auto serial_result = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
          serial, distributed.global_geometry, serial_options);
      if (!serial_result.success) {
        dec3d::test::Fail("serial provider-fed baseline success",
                          __FILE__,
                          __LINE__,
                          serial_result.failure_reason + " | " +
                              serial_result.failure_diagnostics);
      }

      auto distributed_options = MakeOptions();
      distributed_options.boundary_model =
          dec3d::radiation::DistributedRadiationBoundaryModel::thesis_marshak_vacuum;
      const auto distributed_result =
          dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
              distributed, distributed_options);
      CheckSuccess(distributed_result);
      DEC3D_CHECK_EQ(distributed_result.group_count, std::size_t{1});
      DEC3D_CHECK_EQ(distributed_result.per_group_solve_count, std::size_t{1});
      for (std::size_t r = 0; r < distributed.local_state.layout.radial_cells; ++r) {
        const std::size_t global_r = distributed.ownership.global_radial_begin + r;
        for (std::size_t t = 0; t < distributed.local_state.layout.theta_cells; ++t) {
          for (std::size_t p = 0; p < distributed.local_state.layout.phi_cells; ++p) {
            CheckNear(distributed.local_state.radiation_groups[0](r, t, p),
                      serial.radiation_groups[0](global_r, t, p),
                      1.0e-2,
                      "distributed/serial one-group radiation parity");
            CheckNear(distributed.local_state.e_electron(r, t, p),
                      serial.e_electron(global_r, t, p),
                      1.0e2,
                      "distributed/serial one-group electron parity");
            CheckNear(distributed.local_state.e_fluid_total(r, t, p),
                      serial.e_fluid_total(global_r, t, p),
                      1.0e2,
                      "distributed/serial one-group total parity");
          }
        }
      }
    }

    {
      const auto one_group_layout = MakeOneGroupLayout();
      auto problem = MakeProblem(MPI_COMM_WORLD, one_group_layout, table_result.table, 1.0e-18);
      auto options = MakeOptions();
      options.radiation_flux_limiter.model =
          dec3d::radiation::RadiationFluxLimiterModel::harmonic_eq_5_209;
      options.radiation_flux_limiter.U_floor_erg_per_cm3 = 1.0e-30;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, options);
      CheckSuccess(result);
      DEC3D_CHECK_EQ(result.group_count, std::size_t{1});
      DEC3D_CHECK_EQ(result.per_group_solve_count, std::size_t{1});
      DEC3D_CHECK(result.report_line.find(
                      "radiation_flux_limiter_model=harmonic_eq_5_209") !=
                  std::string::npos);
      DEC3D_CHECK(result.report_line.find("per_group_limiter_report_present=true") !=
                  std::string::npos);
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      const auto before = problem.local_state;
      auto options = MakeOptions();
      options.ion_energy_floor_erg_per_cm3 = 1.0e30;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, options);
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!result.global_stage_ok);
      DEC3D_CHECK(!result.global_publish_ok);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::radiation::kDistributedRadiationWritesNoFields);
      DEC3D_CHECK(!result.thermodynamic_recovery_report.empty());
      DEC3D_CHECK(result.failure_diagnostics.find("thermodynamic") != std::string::npos);
      DEC3D_CHECK(problem.local_state.last_authoritative_write_mask == 0u);
      DEC3D_CHECK(problem.local_state.e_electron(0, 0, 0) == before.e_electron(0, 0, 0));
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18, true);
      const auto before = problem.local_state;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, MakeOptions());
      DEC3D_CHECK(!result.success);
      DEC3D_CHECK(!result.global_stage_ok);
      DEC3D_CHECK(!result.global_publish_ok);
      DEC3D_CHECK_EQ(result.updated_fields, dec3d::radiation::kDistributedRadiationWritesNoFields);
      DEC3D_CHECK(!result.canonical_state_mutated);
      DEC3D_CHECK(problem.local_state.last_authoritative_write_mask == 0u);
      DEC3D_CHECK(problem.local_state.e_electron(0, 0, 0) == before.e_electron(0, 0, 0));
      DEC3D_CHECK(result.failure_diagnostics.find("first_failing_rank=") != std::string::npos);
      DEC3D_CHECK(result.failure_diagnostics.find("first_failing_group_index=") !=
                  std::string::npos);
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      auto options = MakeOptions();
      options.boundary_model =
          dec3d::radiation::DistributedRadiationBoundaryModel::thesis_marshak_vacuum;
      const auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, options);
      CheckSuccess(result);
      DEC3D_CHECK(result.report_line.find("marshak_only_on_outer_global_rank=true") !=
                  std::string::npos);
    }

    {
      auto problem = MakeProblem(MPI_COMM_WORLD, group_layout, table_result.table, 1.0e-18);
      auto result = dec3d::radiation::ApplyDistributedProviderFedMultigroupRadiation(
          problem, MakeOptions());
      CheckSuccess(result);
      DEC3D_CHECK(dec3d::radiation::ValidateDistributedMultigroupRadiationDiagnostics(result));
      const std::string token = "off_rank_face_conductance_uses_neighbor_Dbar=";
      const auto pos = result.report_line.find(token);
      DEC3D_CHECK(pos != std::string::npos);
      result.report_line.erase(pos, token.size());
      DEC3D_CHECK(!dec3d::radiation::ValidateDistributedMultigroupRadiationDiagnostics(result));
    }

    MPI_Finalize();
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    MPI_Finalize();
    return 1;
  }
}
