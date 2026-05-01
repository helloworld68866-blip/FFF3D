#include "io/radial_profile.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

void WriteText(const std::filesystem::path& path, const std::string& text) {
  std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  out << text;
}

std::filesystem::path TempRoot() {
  auto root = std::filesystem::temp_directory_path() / "dec3d_p5_io_radial_profile_contract";
  std::filesystem::remove_all(root);
  std::filesystem::create_directories(root);
  return root;
}

}  // namespace

int main() {
  const auto root = TempRoot();
  const auto good = root / "good.pro";
  WriteText(good, R"pro(
# good.pro: valid radial profile
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale # header
0.0 10.0 5.0 4.0 -1.0e7 1.0e5 2.0e5 0.0 1.0 # inner
50.0 20.0 3.0 2.0 -5.0e6 0.0 0.0 1.0e9 0.5 # middle
100.0 30.0 1.0 1.0 0.0 0.0 0.0 2.0e9 0.25 # outer
)pro");

  const auto profile = dec3d::io::LoadRadialProfile(good, "um");
  assert(profile.success);
  assert(dec3d::io::ValidateRadialProfileDiagnostics(profile));
  assert(profile.profile.rows.size() == 3u);
  assert(profile.report_line.find("profile_row_count=3") != std::string::npos);
  assert(profile.report_line.find("profile_validation_success=true") != std::string::npos);

  const auto sample = dec3d::io::SampleRadialProfile(profile.profile, 25.0e-4);
  assert(sample.success);
  assert(std::abs(sample.value.rho_g_cm3 - 15.0) < 1.0e-12);
  assert(std::abs(sample.value.Te_keV - 4.0) < 1.0e-12);
  assert(std::abs(sample.value.Ti_keV - 3.0) < 1.0e-12);
  assert(std::abs(sample.value.vr_cm_s + 7.5e6) < 1.0e-6);
  assert(std::abs(sample.value.radiation_scale - 0.75) < 1.0e-12);

  const auto csv = root / "csv.pro";
  WriteText(csv, R"pro(
# csv.pro: comma profile
r_um,rho_g_cm3,Te_keV,Ti_keV,vr_cm_s # header
0,1,2,3,4 # first
1,2,3,4,5 # second
)pro");
  const auto csv_profile = dec3d::io::LoadRadialProfile(csv, "um");
  assert(csv_profile.success);
  assert(csv_profile.profile.rows.size() == 2u);
  assert(csv_profile.profile.rows[0].vt_cm_s == 0.0);
  assert(csv_profile.profile.rows[0].epsilon_alpha_erg_cm3 == 0.0);
  assert(csv_profile.profile.rows[0].radiation_scale == 1.0);

  const auto nonmonotone = root / "nonmonotone.pro";
  WriteText(nonmonotone,
            "r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s\n"
            "0 1 1 1 0\n"
            "0 1 1 1 0\n");
  const auto bad_radius = dec3d::io::LoadRadialProfile(nonmonotone, "um");
  assert(!bad_radius.success);
  assert(bad_radius.failure_reason.find("strictly increasing") != std::string::npos);

  const auto missing = root / "missing.pro";
  WriteText(missing,
            "r_um rho_g_cm3 Te_keV vr_cm_s\n"
            "0 1 1 0\n"
            "1 1 1 0\n");
  const auto missing_columns = dec3d::io::LoadRadialProfile(missing, "um");
  assert(!missing_columns.success);
  assert(missing_columns.failure_reason.find("Ti_keV") != std::string::npos);

  const auto negative = root / "negative.pro";
  WriteText(negative,
            "r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s\n"
            "0 -1 1 1 0\n"
            "1 1 1 1 0\n");
  const auto bad_density = dec3d::io::LoadRadialProfile(negative, "um");
  assert(!bad_density.success);
  assert(bad_density.failure_reason.find("rho") != std::string::npos);

  std::filesystem::remove_all(root);
  return 0;
}
