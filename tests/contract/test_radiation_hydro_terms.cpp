#include "radiation/hydro_terms/radiation_hydro_terms.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <exception>
#include <iostream>
#include <string>

namespace {

void CheckNear(double actual, double expected, double tolerance, const char* label) {
  if (std::abs(actual - expected) > tolerance) {
    std::cerr << label << ": actual=" << actual << " expected=" << expected << '\n';
    DEC3D_CHECK(false);
  }
}

}  // namespace

int main() {
  try {
    using dec3d::radiation::BuildRadiationHydroTermBundle;
    using dec3d::radiation::RadiationHydroTermsOptions;
    using dec3d::radiation::RadiationPressureScalarFromUg;
    using dec3d::radiation::UgFromRadiationPressureScalar;
    using dec3d::radiation::ValidateRadiationHydroTermsDiagnostics;

    {
      const double ug = 81.0;
      const double chi = RadiationPressureScalarFromUg(ug);
      CheckNear(chi, std::pow(ug / 3.0, 3.0 / 4.0), 1.0e-14, "chi roundtrip");
      CheckNear(UgFromRadiationPressureScalar(chi), ug, 1.0e-12, "Ug roundtrip");
    }

    {
      const auto layout = dec3d::state::CanonicalStateLayout{2u, 1u, 1u, 2u};
      auto state = dec3d::state::CanonicalState::Create(layout);
      state.radiation_groups[0](0u, 0u, 0u) = 3.0;
      state.radiation_groups[0](1u, 0u, 0u) = 6.0;
      state.radiation_groups[1](0u, 0u, 0u) = 12.0;
      state.radiation_groups[1](1u, 0u, 0u) = 24.0;

      RadiationHydroTermsOptions options;
      options.radiation_energy_floor_erg_per_cm3 = 0.0;
      const auto bundle = BuildRadiationHydroTermBundle(state, options);
      DEC3D_CHECK(bundle.success);
      DEC3D_CHECK_EQ(bundle.chi_rad_by_group.size(), std::size_t{2});
      DEC3D_CHECK(ValidateRadiationHydroTermsDiagnostics(bundle.report_line));
      DEC3D_CHECK(bundle.report_line.find("advected_radiation_scalar=P_g_power_3_over_4") !=
                  std::string::npos);
      DEC3D_CHECK(bundle.report_line.find("passive_Ug_advection=false") !=
                  std::string::npos);
    }

    {
      const auto layout = dec3d::state::CanonicalStateLayout{1u, 1u, 1u, 1u};
      auto state = dec3d::state::CanonicalState::Create(layout);
      state.radiation_groups[0](0u, 0u, 0u) = -1.0;
      const auto bundle = BuildRadiationHydroTermBundle(state, RadiationHydroTermsOptions{});
      DEC3D_CHECK(!bundle.success);
      DEC3D_CHECK(bundle.invalid_radiation_cell_count == 1u);
      DEC3D_CHECK(bundle.failure_diagnostics.find(
                      "diagnostic_id=p3.radiation.hydro_terms.failure") !=
                  std::string::npos);
    }

    {
      const double rho_ratio = 2.0;
      const double expected_ug_ratio = std::pow(rho_ratio, 4.0 / 3.0);
      CheckNear(
          dec3d::radiation::RadiationUgCompressionRatioFromDensityRatio(rho_ratio),
          expected_ug_ratio,
          1.0e-14,
          "uniform compression Ug scaling");
    }

    {
      const double ug_old = 9.0;
      const double rho_old = 1.0;
      const double rho_new = 8.0;
      const double ug_new =
          dec3d::radiation::RadiationUgAfterPressureScalarCompression(
              ug_old, rho_old, rho_new);
      CheckNear(
          ug_new / ug_old,
          std::pow(8.0, 4.0 / 3.0),
          1.0e-12,
          "Ug compression oracle");
    }

    {
      const std::string incomplete =
          "diagnostic_id=p3.radiation.hydro_terms; stage_id=H";
      DEC3D_CHECK(!ValidateRadiationHydroTermsDiagnostics(incomplete));
    }

    {
      const std::string missing =
          "diagnostic_id=p3.radiation.hydro_terms"
          "; stage_id=H"
          "; radiation_hydro_terms_enabled=true"
          "; radiation_advection=enabled"
          "; radiation_pressure_work=enabled"
          "; passive_Ug_advection=false"
          "; radiation_group_count=1";
      DEC3D_CHECK(!ValidateRadiationHydroTermsDiagnostics(missing));
    }

    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
