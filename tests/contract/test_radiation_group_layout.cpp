#include "radiation/groups/radiation_group_layout.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "test_assert.hpp"

#include <limits>
#include <string>

int main() {
  using dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout;
  using dec3d::radiation::MakeGrayFullSpectrumRadiationGroupLayout;
  using dec3d::radiation::RadiationGroupLayout;
  using dec3d::radiation::RadiationGroupMode;
  using dec3d::radiation::ValidateRadiationGroupLayout;
  using dec3d::radiation::ValidateRadiationGroupLayoutDiagnostics;
  using dec3d::radiation::ValidateRadiationGroupStateStorage;

  {
    const auto layout = MakeGrayFullSpectrumRadiationGroupLayout();
    const auto result = ValidateRadiationGroupLayout(layout);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.is_complete());
    DEC3D_CHECK_EQ(result.group_count, 1u);
    DEC3D_CHECK_EQ(result.mode, RadiationGroupMode::gray_full_spectrum);
    DEC3D_CHECK(result.gray_full_spectrum);
    DEC3D_CHECK(!result.frequency_edges_used);
    DEC3D_CHECK(result.report_line.find("diagnostic_id=p3.radiation.group_layout") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("stage_id=R") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("radiation_group_mode=gray_full_spectrum") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("radiation_group_count=1") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("group_count_explicit=true") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("frequency_edges_used=false") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("authoritative_state=radiation_groups") != std::string::npos);
    DEC3D_CHECK(ValidateRadiationGroupLayoutDiagnostics(result));
  }

  {
    RadiationGroupLayout missing;
    const auto result = ValidateRadiationGroupLayout(missing);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("radiation group mode is missing") != std::string::npos);
    DEC3D_CHECK(!ValidateRadiationGroupLayoutDiagnostics(result));
  }

  {
    RadiationGroupLayout bad_gray;
    bad_gray.mode = RadiationGroupMode::gray_full_spectrum;
    bad_gray.group_count = 2u;
    const auto result = ValidateRadiationGroupLayout(bad_gray);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("gray full spectrum requires exactly one group") != std::string::npos);
  }

  {
    const auto finite = MakeExplicitFrequencyRadiationGroupLayout({1.0e15, 2.0e15, 4.0e15});
    const auto result = ValidateRadiationGroupLayout(finite);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK_EQ(result.group_count, 2u);
    DEC3D_CHECK(result.frequency_edges_used);
    DEC3D_CHECK(result.frequency_edges_strictly_increasing);
    DEC3D_CHECK(result.report_line.find("radiation_group_mode=explicit_frequency_groups") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("radiation_group_count=2") != std::string::npos);
    DEC3D_CHECK(result.report_line.find("frequency_edges_used=true") != std::string::npos);
  }

  {
    const auto unsorted = MakeExplicitFrequencyRadiationGroupLayout({1.0e15, 1.0e15, 2.0e15});
    const auto result = ValidateRadiationGroupLayout(unsorted);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("frequency edges must be strictly increasing") != std::string::npos);
  }

  {
    const auto nonpositive = MakeExplicitFrequencyRadiationGroupLayout({0.0, 1.0e15});
    const auto result = ValidateRadiationGroupLayout(nonpositive);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("frequency edges must be finite and positive") != std::string::npos);
  }

  {
    const auto infinite =
        MakeExplicitFrequencyRadiationGroupLayout({1.0e15, std::numeric_limits<double>::infinity()});
    const auto result = ValidateRadiationGroupLayout(infinite);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("frequency edges must be finite and positive") != std::string::npos);
  }

  {
    auto state = dec3d::state::CanonicalState::Create(dec3d::state::CanonicalStateLayout{2u, 2u, 2u, 1u});
    const auto layout = MakeGrayFullSpectrumRadiationGroupLayout();
    const auto result = ValidateRadiationGroupStateStorage(layout, state);
    DEC3D_CHECK(result.success);
    DEC3D_CHECK(result.canonical_group_arrays_present);
    DEC3D_CHECK(result.report_line.find("canonical_group_arrays_present=true") != std::string::npos);
  }

  {
    auto state = dec3d::state::CanonicalState::Create(dec3d::state::CanonicalStateLayout{2u, 2u, 2u, 0u});
    const auto layout = MakeGrayFullSpectrumRadiationGroupLayout();
    const auto result = ValidateRadiationGroupStateStorage(layout, state);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("canonical radiation group count mismatch") != std::string::npos);
  }

  {
    auto state = dec3d::state::CanonicalState::Create(dec3d::state::CanonicalStateLayout{2u, 2u, 2u, 1u});
    state.radiation_groups[0] = dec3d::core::Array3D<double>(1u, 2u, 2u, 0.0);
    const auto layout = MakeGrayFullSpectrumRadiationGroupLayout();
    const auto result = ValidateRadiationGroupStateStorage(layout, state);
    DEC3D_CHECK(!result.success);
    DEC3D_CHECK(result.failure_reason.find("canonical radiation group shape mismatch") != std::string::npos);
  }

  return 0;
}
