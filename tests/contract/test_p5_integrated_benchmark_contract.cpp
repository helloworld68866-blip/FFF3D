#include "benchmarks/p5_integrated_benchmark.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <string>

namespace {

bool HasToken(const std::string& text, const char* token) {
  return text.find(token) != std::string::npos;
}

}  // namespace

int main() {
  using namespace dec3d::benchmarks;

  const auto clean = MakeP5IntegratedBenchmarkDescriptor(
      P5IntegratedBenchmarkCase::clean_spherical,
      P5IntegratedBenchmarkVariant::primary_mpi24_ale);
  DEC3D_CHECK(clean.case_id == "p5_clean_spherical");
  DEC3D_CHECK(clean.mode_l == 0);
  DEC3D_CHECK(clean.mode_m == 0);
  DEC3D_CHECK(clean.ale_enabled);
  DEC3D_CHECK(clean.mpi_ranks == 24);
  DEC3D_CHECK(clean.ppm_enabled);
  DEC3D_CHECK(clean.macro_enabled);

  const auto l2 = MakeP5IntegratedBenchmarkDescriptor(
      P5IntegratedBenchmarkCase::legendre_l2_m0,
      P5IntegratedBenchmarkVariant::primary_mpi24_ale);
  DEC3D_CHECK(l2.mode_l == 2);
  DEC3D_CHECK(l2.mode_m == 0);
  DEC3D_CHECK(std::abs(LegendreP(2, 1.0) - 1.0) < 1.0e-14);
  DEC3D_CHECK(std::abs(LegendreP(2, 0.0) + 0.5) < 1.0e-14);

  const auto l4 = MakeP5IntegratedBenchmarkDescriptor(
      P5IntegratedBenchmarkCase::legendre_l4_m0,
      P5IntegratedBenchmarkVariant::mpi12_ale);
  DEC3D_CHECK(l4.mode_l == 4);
  DEC3D_CHECK(l4.mode_m == 0);
  DEC3D_CHECK(l4.mpi_ranks == 12);
  DEC3D_CHECK(l4.ale_enabled);
  DEC3D_CHECK(std::abs(LegendreP(4, 1.0) - 1.0) < 1.0e-14);
  DEC3D_CHECK(std::abs(LegendreP(4, 0.0) - 0.375) < 1.0e-14);

  const double r0 = 24.070450097847356e-4;
  DEC3D_CHECK(std::abs(WooVelocityPerturbationShape(2, r0, r0) - 1.0) <
              1.0e-13);
  DEC3D_CHECK(std::abs(WooVelocityPerturbationShape(4, r0, r0) - 1.0) <
              1.0e-13);

  const auto report = BuildP5IntegratedBenchmarkDiagnosticsLine(l2);
  DEC3D_CHECK(HasToken(report, "diagnostic_id=p5.integrated_benchmark"));
  DEC3D_CHECK(HasToken(report, "case_id=p5_legendre_l2_m0"));
  DEC3D_CHECK(HasToken(report, "variant_id=primary"));
  DEC3D_CHECK(HasToken(report, "perturbation_family=legendre_velocity"));
  DEC3D_CHECK(HasToken(report, "perturbation_applied_to=v_r"));
  DEC3D_CHECK(HasToken(report, "density_angular_perturbation=false"));
  DEC3D_CHECK(HasToken(report, "temperature_angular_perturbation=false"));
  DEC3D_CHECK(HasToken(report, "mode_l=2"));
  DEC3D_CHECK(HasToken(report, "mode_m=0"));
  DEC3D_CHECK(HasToken(report, "stage_order=H,T,E,R,A"));
  DEC3D_CHECK(ValidateP5IntegratedBenchmarkDiagnostics(report));

  auto broken = report;
  const auto pos = broken.find("mode_m=0");
  DEC3D_CHECK(pos != std::string::npos);
  broken.erase(pos, std::string("mode_m=0").size());
  DEC3D_CHECK(!ValidateP5IntegratedBenchmarkDiagnostics(broken));

  const auto manifest = BuildP5IntegratedBenchmarkArtifactManifest(
      l2, "F:/dec3d/analysis/output/p5_integrated_benchmark_contract");
  DEC3D_CHECK(manifest.ContainsRequiredFile("p5_case_descriptor.json"));
  DEC3D_CHECK(manifest.ContainsRequiredFile("p5_stage_reports.json"));
  DEC3D_CHECK(
      manifest.ContainsVisualArtifact("p5_p5_legendre_l2_m0_primary_rz_evolution.gif"));
  DEC3D_CHECK(!manifest.ContainsRequiredFile(
      "p5_p5_legendre_l2_m0_primary_rz_evolution.gif"));

  return 0;
}
