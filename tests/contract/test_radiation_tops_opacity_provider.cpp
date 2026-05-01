#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/groups/radiation_group_layout.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

namespace {

constexpr const char* kTableRoot = "F:/dec3d/data/opacities/tops_dt_2026_04_27";

double PlanckConstantErgS() {
  return 2.0 * dec3d::physics::PhysicsConstantsCGS::pi *
         dec3d::physics::PhysicsConstantsCGS::hbar_erg_s;
}

dec3d::mesh::SphericalGeometryMetadata MakeTwoCellGeometry() {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces = {0.1, 0.2, 0.3};
  geometry.theta_faces = {0.5, 1.0};
  geometry.phi_faces = {0.0, 2.0 * dec3d::physics::PhysicsConstantsCGS::pi};
  geometry.cell_volumes.resize(2u);
  for (std::size_t r = 0; r < 2u; ++r) {
    const double radial =
        (std::pow(geometry.radial_faces[r + 1u], 3) -
         std::pow(geometry.radial_faces[r], 3)) /
        3.0;
    const double polar =
        std::cos(geometry.theta_faces[0]) - std::cos(geometry.theta_faces[1]);
    const double azimuthal = geometry.phi_faces[1] - geometry.phi_faces[0];
    geometry.cell_volumes[r] = radial * polar * azimuthal;
  }
  geometry.global_volume = geometry.cell_volumes[0] + geometry.cell_volumes[1];
  return geometry;
}

dec3d::transport::GenericDiffusionBoundaryPolicy MakePatchBoundary() {
  dec3d::transport::GenericDiffusionBoundaryPolicy policy;
  policy.inner_radial =
      dec3d::transport::DiffusionBoundaryKind::interior_patch_no_origin;
  policy.outer_radial = dec3d::transport::DiffusionBoundaryKind::neumann_zero_flux;
  policy.theta_lower =
      dec3d::transport::DiffusionBoundaryKind::interior_patch_no_pole;
  policy.theta_upper =
      dec3d::transport::DiffusionBoundaryKind::interior_patch_no_pole;
  policy.phi = dec3d::transport::DiffusionBoundaryKind::periodic;
  return policy;
}

dec3d::state::CanonicalState MakeTwoCellState() {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{2u, 1u, 1u, 2u});
  const double erg_per_keV = dec3d::physics::PhysicsConstantsCGS::erg_per_kev;
  const double gamma_minus_one = 2.0 / 3.0;
  for (std::size_t r = 0; r < 2u; ++r) {
    state.rho(r, 0u, 0u) = 1.0 + static_cast<double>(r);
    state.mom_r(r, 0u, 0u) = 0.0;
    state.mom_theta(r, 0u, 0u) = 0.0;
    state.mom_phi(r, 0u, 0u) = 0.0;
    state.radiation_groups[0](r, 0u, 0u) = 1.0;
    state.radiation_groups[1](r, 0u, 0u) = 2.0;
    const double te = (0.5 + 0.25 * static_cast<double>(r)) * erg_per_keV;
    const double ti = 0.4 * erg_per_keV;
    const double ni = state.rho(r, 0u, 0u) /
                      dec3d::physics::DefaultMeanDTIonMassG();
    const double ne = ni;
    state.e_electron(r, 0u, 0u) = ne * te / gamma_minus_one;
    const double e_ion = ni * ti / gamma_minus_one;
    state.e_fluid_total(r, 0u, 0u) = state.e_electron(r, 0u, 0u) + e_ion;
  }
  return state;
}

dec3d::state::CanonicalState MakeRepeatedPhiState() {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{1u, 1u, 2u, 2u});
  const double erg_per_keV = dec3d::physics::PhysicsConstantsCGS::erg_per_kev;
  const double gamma_minus_one = 2.0 / 3.0;
  for (std::size_t p = 0; p < 2u; ++p) {
    state.rho(0u, 0u, p) = 1.0;
    state.mom_r(0u, 0u, p) = 0.0;
    state.mom_theta(0u, 0u, p) = 0.0;
    state.mom_phi(0u, 0u, p) = 0.0;
    state.radiation_groups[0](0u, 0u, p) = 1.0;
    state.radiation_groups[1](0u, 0u, p) = 2.0;
    const double te = 0.5 * erg_per_keV;
    const double ti = 0.4 * erg_per_keV;
    const double ni = state.rho(0u, 0u, p) /
                      dec3d::physics::DefaultMeanDTIonMassG();
    const double ne = ni;
    state.e_electron(0u, 0u, p) = ne * te / gamma_minus_one;
    const double e_ion = ni * ti / gamma_minus_one;
    state.e_fluid_total(0u, 0u, p) = state.e_electron(0u, 0u, p) + e_ion;
  }
  return state;
}

dec3d::radiation::RadiationGroupLayout MakeTwoGroupLayout() {
  const double kev_to_hz =
      dec3d::physics::PhysicsConstantsCGS::erg_per_kev / PlanckConstantErgS();
  return dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(
      {0.002 * kev_to_hz, 0.02 * kev_to_hz, 0.2 * kev_to_hz});
}

dec3d::radiation::TopsOpacityProviderOptions MakeProviderOptions() {
  dec3d::radiation::TopsOpacityProviderOptions options;
  options.table_root = kTableRoot;
  options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  options.density_clip_policy =
      dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return options;
}

void CheckSuccess(const dec3d::radiation::RadiationCoefficientProviderResult& result) {
  if (!result.success) {
    dec3d::test::Fail("TOPS opacity provider success",
                      __FILE__,
                      __LINE__,
                      result.failure_reason + " | " + result.failure_diagnostics);
  }
}

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    dec3d::test::Fail(label, __FILE__, __LINE__);
  }
}

}  // namespace

int RunTests() {
  {
    const auto table = dec3d::radiation::LoadTopsOpacityTable(MakeProviderOptions());
    DEC3D_CHECK(table.success);
    DEC3D_CHECK_EQ(table.metadata.temperature_count, 66u);
    DEC3D_CHECK_EQ(table.metadata.density_count, 100u);
    DEC3D_CHECK_EQ(table.metadata.photon_energy_count, 32u);
    DEC3D_CHECK_EQ(table.metadata.mean_row_count, 6600u);
    DEC3D_CHECK_EQ(table.metadata.multigroup_row_count, 211200u);
    DEC3D_CHECK_EQ(table.metadata.density_clipping_warning_count, 30u);
    DEC3D_CHECK(table.metadata.material_label == "DT");
    DEC3D_CHECK(table.report_line.find("opacity_data_source=TOPS") != std::string::npos);
    DEC3D_CHECK(table.report_line.find("table_load_count=1") != std::string::npos);
    DEC3D_CHECK(table.report_line.find("lookup_table_resident=true") != std::string::npos);
    DEC3D_CHECK(table.report_line.find("lookup_index_mode=in_memory_coordinate_index") !=
                std::string::npos);
  }

  {
    const auto table = dec3d::radiation::LoadTopsOpacityTable(MakeProviderOptions());
    DEC3D_CHECK(table.success);

    const auto exact = dec3d::radiation::LookupTopsOpacity(
        table.table, 0.001, 0.0001, 0.001, MakeProviderOptions());
    DEC3D_CHECK(exact.success);
    CheckNear(exact.kappaR_mass_cm2_g, 1206.0, 1.0e-10, "exact grid kappaR");
    CheckNear(exact.kappaP_mass_cm2_g, 1238.3, 1.0e-10, "exact grid kappaP");

    const auto interior = dec3d::radiation::LookupTopsOpacity(
        table.table, 0.75, 1.0, 0.05, MakeProviderOptions());
    DEC3D_CHECK(interior.success);
    DEC3D_CHECK(interior.kappaR_mass_cm2_g > 0.0);
    DEC3D_CHECK(interior.kappaP_mass_cm2_g > 0.0);
    DEC3D_CHECK(interior.report_line.find("opacity_interpolation_mode=loglog_trilinear") !=
                std::string::npos);
  }

  {
    const auto table = dec3d::radiation::LoadTopsOpacityTable(MakeProviderOptions());
    DEC3D_CHECK(table.success);
    const auto clipped = dec3d::radiation::LookupTopsOpacity(
        table.table, 0.001, 1000.0, 0.01, MakeProviderOptions());
    DEC3D_CHECK(!clipped.success);
    DEC3D_CHECK(clipped.failure_diagnostics.find("density_clip_policy=hard_fail") !=
                std::string::npos);
  }

  {
    auto state = MakeTwoCellState();
    const auto result = dec3d::radiation::BuildRadiationCoefficientArrays(
        state, MakeTwoCellGeometry(), MakeTwoGroupLayout(), MakeProviderOptions());
    CheckSuccess(result);
    DEC3D_CHECK(dec3d::radiation::ValidateRadiationCoefficientProviderDiagnostics(result));
    DEC3D_CHECK_EQ(result.coefficients.Dbar_cm2_per_s.size(), 2u);
    DEC3D_CHECK_EQ(result.coefficients.kappaP_cm_inv.size(), 2u);
    DEC3D_CHECK_EQ(result.coefficients.B_erg_per_cm3.size(), 2u);
    DEC3D_CHECK(result.min_Dbar_cm2_s > 0.0);
    DEC3D_CHECK(result.max_kappaP_mass_cm2_g > 0.0);
    DEC3D_CHECK(result.max_B_g_erg_cm3 > 0.0);
    DEC3D_CHECK(result.report_line.find("group_edges_source=RadiationGroupLayout") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("opacity_energy_coordinate_source=TOPS photon grid") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("per_cell_csv_io=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("per_lookup_full_table_scan=false") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("double_kB_guard_passed=true") != std::string::npos);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                  static_cast<dec3d::core::AuthoritativeFieldMask>(0));
  }

  {
    auto state = MakeRepeatedPhiState();
    const auto result = dec3d::radiation::BuildRadiationCoefficientArrays(
        state, MakeTwoCellGeometry(), MakeTwoGroupLayout(), MakeProviderOptions());
    CheckSuccess(result);
    DEC3D_CHECK(result.exact_state_cache_hit_count > 0u);
    DEC3D_CHECK(result.exact_state_cache_miss_count > 0u);
    DEC3D_CHECK(result.report_line.find(
                    "exact_state_cache_mode=bitwise_thermodynamic_state") !=
                std::string::npos);
    DEC3D_CHECK(result.report_line.find("exact_state_cache_hit_count=") !=
                std::string::npos);
  }

  {
    auto state = MakeTwoCellState();
    const auto provider = dec3d::radiation::BuildRadiationCoefficientArrays(
        state, MakeTwoCellGeometry(), MakeTwoGroupLayout(), MakeProviderOptions());
    CheckSuccess(provider);

    dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions radiation_options;
    radiation_options.dt_s = 0.0;
    radiation_options.group_layout = MakeTwoGroupLayout();
    radiation_options.coefficients = provider.coefficients;
    radiation_options.boundary_policy = MakePatchBoundary();
    const auto radiation = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state, MakeTwoCellGeometry(), radiation_options);
    if (!radiation.success) {
      dec3d::test::Fail("provider coefficients feed multigroup radiation operator",
                        __FILE__,
                        __LINE__,
                        radiation.failure_reason + " | " +
                            radiation.failure_diagnostics);
    }
    DEC3D_CHECK(radiation.report_line.find(
                    "per_group_coefficients_source=tops_dt_tabulated") !=
                std::string::npos);
    DEC3D_CHECK(radiation.report_line.find("opacity_provider=tops_dt_tabulated") !=
                std::string::npos);
    DEC3D_CHECK(radiation.report_line.find("Bg_source=blackbody_group_integral") !=
                std::string::npos);
    DEC3D_CHECK(radiation.report_line.find("updated_fields=none") !=
                std::string::npos);
  }

  {
    auto state = MakeTwoCellState();
    auto result = dec3d::radiation::BuildRadiationCoefficientArrays(
        state, MakeTwoCellGeometry(), MakeTwoGroupLayout(), MakeProviderOptions());
    CheckSuccess(result);
    const auto token_pos = result.report_line.find("double_kB_guard_passed=true");
    DEC3D_CHECK(token_pos != std::string::npos);
    result.report_line.erase(token_pos, std::string("double_kB_guard_passed=true").size());
    DEC3D_CHECK(!dec3d::radiation::ValidateRadiationCoefficientProviderDiagnostics(result));
  }

  return 0;
}

int main() {
  try {
    return RunTests();
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
