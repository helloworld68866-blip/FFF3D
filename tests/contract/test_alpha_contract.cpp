#include "alpha/alpha_coefficients.hpp"
#include "alpha/alpha_contract.hpp"
#include "core/diagnostics/stage_contracts.hpp"
#include "physics/units/physical_constants.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

namespace {

void CheckNear(double actual, double expected, double tolerance, const std::string& label) {
  if (std::abs(actual - expected) > tolerance) {
    std::ostringstream detail;
    detail << label << ": expected " << expected << ", got " << actual
           << ", tolerance " << tolerance;
    dec3d::test::Fail("near", __FILE__, __LINE__, detail.str());
  }
}

void CheckRelativeNear(
    double actual,
    double expected,
    double relative_tolerance,
    const std::string& label) {
  const double scale = std::max(1.0e-300, std::abs(expected));
  if (std::abs(actual - expected) > relative_tolerance * scale) {
    std::ostringstream detail;
    detail << label << ": expected " << expected << ", got " << actual
           << ", relative tolerance " << relative_tolerance;
    dec3d::test::Fail("relative near", __FILE__, __LINE__, detail.str());
  }
}

dec3d::state::CanonicalState BuildAlphaState() {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{2u, 1u, 1u, 0u});

  const double rho = 10.0;
  const double Te = dec3d::physics::ErgFromKeV(4.0);
  const double Ti = dec3d::physics::ErgFromKeV(3.0);
  const double ni = rho / dec3d::physics::DefaultMeanDTIonMassG();
  const double ne = ni;
  const double gamma_minus_one = 2.0 / 3.0;

  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    state.rho(r, 0, 0) = rho;
    state.mom_r(r, 0, 0) = 0.0;
    state.mom_theta(r, 0, 0) = 0.0;
    state.mom_phi(r, 0, 0) = 0.0;
    state.e_electron(r, 0, 0) = ne * Te / gamma_minus_one;
    state.e_fluid_total(r, 0, 0) = (ne * Te + ni * Ti) / gamma_minus_one;
    state.alpha_state.storage(r, 0, 0) = 1.0e10 * (1.0 + static_cast<double>(r));
  }
  return state;
}

}  // namespace

int main() {
  {
    auto state = BuildAlphaState();
    DEC3D_CHECK(state.alpha_state.has_storage());
    DEC3D_CHECK(state.HasAuthoritativeStorage(dec3d::core::AuthoritativeField::alpha_state));

    const auto result = dec3d::alpha::ValidateAlphaStateContract(state);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.validated_fields ==
                dec3d::core::ToMask(dec3d::core::AuthoritativeField::alpha_state));
    DEC3D_CHECK_EQ(result.updated_fields, static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(result.report_line.find("diagnostic_id=p4.alpha.contract") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("stage_id=A") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("alpha_state_authoritative=alpha_energy_density") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("validated_fields=alpha_state") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("updated_fields=none") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("alpha_energy_unit=erg_per_cm3") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("alpha_transport_model=atzeni_one_group") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("canonical_state_mutated=false") != std::string::npos);
    DEC3D_CHECK(dec3d::alpha::ValidateAlphaContractDiagnostics(result));
  }

  {
    const double eps = 9.0e12;
    const double chi = dec3d::alpha::AlphaPressureScalarFromEnergyDensity(eps);
    const double restored = dec3d::alpha::AlphaEnergyDensityFromPressureScalar(chi);
    CheckNear(restored, eps, 1.0e-10 * eps, "alpha pressure scalar roundtrip");
  }

  {
    auto state = BuildAlphaState();
    dec3d::alpha::AlphaCoefficientProviderOptions options;
    options.composition_model = dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
    options.tau_model = dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
    options.reactivity_model = dec3d::alpha::AlphaReactivityModel::constant_user_supplied;
    options.constant_dt_reactivity_cm3_s = 1.0e-16;
    options.recovery_options.mean_ion_mass_g = dec3d::physics::DefaultMeanDTIonMassG();
    options.recovery_options.zbar = dec3d::physics::DefaultZbar();

    const auto before_alpha_0 = state.alpha_state.storage(0, 0, 0);
    const auto before_alpha_1 = state.alpha_state.storage(1, 0, 0);
    const auto before_electron_0 = state.e_electron(0, 0, 0);
    const auto before_electron_1 = state.e_electron(1, 0, 0);
    const auto before_total_0 = state.e_fluid_total(0, 0, 0);
    const auto before_total_1 = state.e_fluid_total(1, 0, 0);
    const auto before_write_mask = state.last_authoritative_write_mask;
    const auto result = dec3d::alpha::BuildAlphaCoefficientArrays(state, options);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK_EQ(result.cell_count, 2u);
    DEC3D_CHECK(result.coefficients.tau_alphae_s(0, 0, 0) > 0.0);
    DEC3D_CHECK(result.coefficients.lambda_drag_cm(0, 0, 0) > 0.0);
    DEC3D_CHECK(result.coefficients.D_alpha_cm2_s(0, 0, 0) > 0.0);
    DEC3D_CHECK(result.coefficients.birth_source_erg_cm3_s(0, 0, 0) > 0.0);
    DEC3D_CHECK(result.coefficients.lnLambda_alphae_spitzer(0, 0, 0) > 0.0);
    CheckNear(result.coefficients.nD_cm3(0, 0, 0),
              0.5 * result.recovered.cells(0, 0, 0).n_i_cm3,
              1.0e-12 * result.recovered.cells(0, 0, 0).n_i_cm3,
              "equimolar nD");
    CheckNear(result.coefficients.nT_cm3(0, 0, 0),
              result.coefficients.nD_cm3(0, 0, 0),
              1.0e-12 * result.coefficients.nD_cm3(0, 0, 0),
              "equimolar nT");
    DEC3D_CHECK(state.alpha_state.storage(0, 0, 0) == before_alpha_0);
    DEC3D_CHECK(state.alpha_state.storage(1, 0, 0) == before_alpha_1);
    DEC3D_CHECK(state.e_electron(0, 0, 0) == before_electron_0);
    DEC3D_CHECK(state.e_electron(1, 0, 0) == before_electron_1);
    DEC3D_CHECK(state.e_fluid_total(0, 0, 0) == before_total_0);
    DEC3D_CHECK(state.e_fluid_total(1, 0, 0) == before_total_1);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask, before_write_mask);
    DEC3D_CHECK_EQ(state.last_authoritative_write_mask,
                   static_cast<dec3d::core::AuthoritativeFieldMask>(0));
    DEC3D_CHECK(
        result.report_line.find("composition_model=equimolar_dt_from_p2_recovery") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("tau_alphae_model=thesis_spitzer_eq_5_262") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("reactivity_model=constant_user_supplied") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("thesis_reactivity_claim_allowed=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("alpha_transport_model=atzeni_one_group") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("double_kB_guard_passed=true") != std::string::npos);
    DEC3D_CHECK(dec3d::alpha::ValidateAlphaCoefficientDiagnostics(result));
  }

  {
    struct Reference {
      double Ti_keV;
      double reactivity_cm3_s;
    };
    const Reference references[] = {
        {1.0, 6.85688430133071390e-21},
        {3.0, 1.86696927818874975e-18},
        {10.0, 1.13618134482734562e-16},
        {30.0, 6.68102563545722182e-16},
        {100.0, 8.44551910167966373e-16},
    };
    for (const auto& reference : references) {
      const double actual = dec3d::alpha::DtReactivityBoschHaleCm3PerS(reference.Ti_keV);
      CheckRelativeNear(
          actual, reference.reactivity_cm3_s, 1.0e-12, "Bosch-Hale reference value");
    }
  }

  {
    auto state = BuildAlphaState();
    dec3d::alpha::AlphaCoefficientProviderOptions options;
    options.composition_model = dec3d::alpha::AlphaCompositionModel::equimolar_dt_from_p2_recovery;
    options.tau_model = dec3d::alpha::AlphaTauModel::thesis_spitzer_eq_5_262;
    options.reactivity_model = dec3d::alpha::AlphaReactivityModel::bosch_hale_dt;
    options.recovery_options.mean_ion_mass_g = dec3d::physics::DefaultMeanDTIonMassG();
    options.recovery_options.zbar = dec3d::physics::DefaultZbar();

    const auto result = dec3d::alpha::BuildAlphaCoefficientArrays(state, options);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK(result.thesis_reactivity_claim_allowed);
    DEC3D_CHECK(result.bosch_hale_reference_locked_locally);
    DEC3D_CHECK(
        result.report_line.find("reactivity_model_requested=bosch_hale_dt") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("reactivity_model_executed=bosch_hale_dt") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("bosch_hale_temperature_source=Ti_old") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("dt_reactivity_internal_unit=cm3_s") != std::string::npos);
    DEC3D_CHECK(
        result.report_line.find("thesis_reactivity_claim_allowed=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("fallback_used=false") != std::string::npos);

    const double expected = dec3d::alpha::DtReactivityBoschHaleCm3PerS(3.0);
    CheckRelativeNear(result.coefficients.dt_reactivity_cm3_s(0, 0, 0),
                      expected,
                      1.0e-12,
                      "provider uses Ti=3 keV Bosch-Hale");
  }

  {
    const double cm3_s = dec3d::alpha::DtReactivityBoschHaleCm3PerS(10.0);
    const double mistaken_m3_s = cm3_s * 1.0e-6;
    DEC3D_CHECK(cm3_s > mistaken_m3_s);
    CheckRelativeNear(cm3_s / mistaken_m3_s, 1.0e6, 1.0e-15, "Bosch-Hale cgs/SI guard");
  }

  {
    dec3d::alpha::AlphaStateContractResult missing;
    missing.success = true;
    DEC3D_CHECK(!dec3d::alpha::ValidateAlphaContractDiagnostics(missing));
  }

  return 0;
}
