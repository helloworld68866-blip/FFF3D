#include "benchmarks/p5_integrated_benchmark.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace dec3d::benchmarks {
namespace {

[[nodiscard]] bool Contains(const std::string& text, const char* token) noexcept {
  return text.find(token) != std::string::npos;
}

[[nodiscard]] const char* CaseId(P5IntegratedBenchmarkCase case_kind) noexcept {
  switch (case_kind) {
    case P5IntegratedBenchmarkCase::clean_spherical:
      return "p5_clean_spherical";
    case P5IntegratedBenchmarkCase::legendre_l2_m0:
      return "p5_legendre_l2_m0";
    case P5IntegratedBenchmarkCase::legendre_l4_m0:
      return "p5_legendre_l4_m0";
  }
  return "unknown";
}

[[nodiscard]] int ModeL(P5IntegratedBenchmarkCase case_kind) noexcept {
  switch (case_kind) {
    case P5IntegratedBenchmarkCase::clean_spherical:
      return 0;
    case P5IntegratedBenchmarkCase::legendre_l2_m0:
      return 2;
    case P5IntegratedBenchmarkCase::legendre_l4_m0:
      return 4;
  }
  return 0;
}

[[nodiscard]] const char* VariantId(
    P5IntegratedBenchmarkVariant variant) noexcept {
  switch (variant) {
    case P5IntegratedBenchmarkVariant::primary_mpi24_ale:
      return "primary";
    case P5IntegratedBenchmarkVariant::no_ale_mpi24:
      return "no_ale";
    case P5IntegratedBenchmarkVariant::mpi12_ale:
      return "mpi12";
  }
  return "unknown";
}

[[nodiscard]] int VariantRanks(P5IntegratedBenchmarkVariant variant) noexcept {
  switch (variant) {
    case P5IntegratedBenchmarkVariant::primary_mpi24_ale:
    case P5IntegratedBenchmarkVariant::no_ale_mpi24:
      return 24;
    case P5IntegratedBenchmarkVariant::mpi12_ale:
      return 12;
  }
  return 0;
}

[[nodiscard]] bool VariantAleEnabled(
    P5IntegratedBenchmarkVariant variant) noexcept {
  return variant != P5IntegratedBenchmarkVariant::no_ale_mpi24;
}

}  // namespace

bool P5IntegratedBenchmarkArtifactManifest::ContainsRequiredFile(
    const std::string& file_name) const noexcept {
  return std::find(required_files.begin(), required_files.end(), file_name) !=
         required_files.end();
}

bool P5IntegratedBenchmarkArtifactManifest::ContainsVisualArtifact(
    const std::string& file_name) const noexcept {
  return std::find(visual_artifacts.begin(), visual_artifacts.end(), file_name) !=
         visual_artifacts.end();
}

const char* ToString(P5IntegratedBenchmarkCase case_kind) noexcept {
  return CaseId(case_kind);
}

const char* ToString(P5IntegratedBenchmarkVariant variant) noexcept {
  return VariantId(variant);
}

P5IntegratedBenchmarkDescriptor MakeP5IntegratedBenchmarkDescriptor(
    P5IntegratedBenchmarkCase case_kind,
    P5IntegratedBenchmarkVariant variant) noexcept {
  P5IntegratedBenchmarkDescriptor descriptor;
  descriptor.case_kind = case_kind;
  descriptor.variant = variant;
  descriptor.case_id = CaseId(case_kind);
  descriptor.variant_id = VariantId(variant);
  descriptor.mode_l = ModeL(case_kind);
  descriptor.mode_m = 0;
  descriptor.ale_enabled = VariantAleEnabled(variant);
  descriptor.mpi_ranks = VariantRanks(variant);
  return descriptor;
}

double LegendreP(int ell, double mu) noexcept {
  switch (ell) {
    case 0:
      return 1.0;
    case 1:
      return mu;
    case 2:
      return 0.5 * (3.0 * mu * mu - 1.0);
    case 3:
      return 0.5 * (5.0 * mu * mu * mu - 3.0 * mu);
    case 4: {
      const double mu2 = mu * mu;
      return (35.0 * mu2 * mu2 - 30.0 * mu2 + 3.0) / 8.0;
    }
    default:
      return 0.0;
  }
}

double WooVelocityPerturbationShape(
    int ell,
    double r_cm,
    double r0_cm) noexcept {
  if (ell <= 0 || r_cm <= 0.0 || r0_cm <= 0.0) {
    return 0.0;
  }
  const double width = 0.015 * r0_cm;
  const double transition = std::tanh((r_cm - r0_cm) / width);
  const double inner = std::pow(r_cm / r0_cm, ell);
  const double outer = std::pow(r0_cm / r_cm, ell);
  return 0.5 * inner * (1.0 - transition) +
         0.5 * outer * (1.0 + transition);
}

double WooVelocityPerturbationDeltaVr(
    int ell,
    double r_cm,
    double theta_rad,
    double r0_cm,
    double amplitude,
    double reference_velocity_cm_s) noexcept {
  const double mu = std::cos(theta_rad);
  return amplitude * std::abs(reference_velocity_cm_s) *
         WooVelocityPerturbationShape(ell, r_cm, r0_cm) * LegendreP(ell, mu);
}

std::string BuildP5IntegratedBenchmarkDiagnosticsLine(
    const P5IntegratedBenchmarkDescriptor& descriptor) {
  std::ostringstream out;
  out << "diagnostic_id=p5.integrated_benchmark"
      << "; benchmark_id=" << descriptor.benchmark_id
      << "; case_id=" << descriptor.case_id
      << "; variant_id=" << descriptor.variant_id
      << "; stage_order=H,T,E,R,A"
      << "; ppm_enabled=" << (descriptor.ppm_enabled ? "true" : "false")
      << "; macro_enabled=" << (descriptor.macro_enabled ? "true" : "false")
      << "; ale_enabled=" << (descriptor.ale_enabled ? "true" : "false")
      << "; mpi_ranks=" << descriptor.mpi_ranks
      << "; perturbation_family=legendre_velocity"
      << "; perturbation_applied_to=v_r"
      << "; density_angular_perturbation=false"
      << "; temperature_angular_perturbation=false"
      << "; mode_l=" << descriptor.mode_l
      << "; mode_m=" << descriptor.mode_m
      << "; initial_profile_source=image_profile_100ps_hydro_yaml"
      << "; initial_profile_yaml_path=" << descriptor.initial_profile_yaml_path
      << "; initial_profile_csv_path=" << descriptor.initial_profile_csv_path
      << "; perturbation_r0_source=" << descriptor.perturbation_r0_source
      << "; perturbation_r0_um=" << descriptor.perturbation_r0_cm * 1.0e4
      << "; woo_reference_r0_um=" << descriptor.woo_reference_r0_um
      << "; r0_matches_woo_77068="
      << (descriptor.r0_matches_woo_77068 ? "true" : "false")
      << "; parity_claim_allowed=false"
      << "; lilac_reference_available=false"
      << "; one_dimensional_reference_available=false";
  return out.str();
}

bool ValidateP5IntegratedBenchmarkDiagnostics(
    const std::string& report_line) noexcept {
  return Contains(report_line, "diagnostic_id=p5.integrated_benchmark") &&
         Contains(report_line, "case_id=") &&
         Contains(report_line, "variant_id=") &&
         Contains(report_line, "stage_order=H,T,E,R,A") &&
         Contains(report_line, "ppm_enabled=true") &&
         Contains(report_line, "macro_enabled=true") &&
         Contains(report_line, "perturbation_family=legendre_velocity") &&
         Contains(report_line, "perturbation_applied_to=v_r") &&
         Contains(report_line, "density_angular_perturbation=false") &&
         Contains(report_line, "temperature_angular_perturbation=false") &&
         Contains(report_line, "mode_l=") &&
         Contains(report_line, "mode_m=0") &&
         Contains(report_line,
                  "initial_profile_source=image_profile_100ps_hydro_yaml") &&
         Contains(report_line, "perturbation_r0_source=p5_case_descriptor") &&
         Contains(report_line, "r0_matches_woo_77068=false") &&
         Contains(report_line, "parity_claim_allowed=false");
}

P5IntegratedBenchmarkArtifactManifest BuildP5IntegratedBenchmarkArtifactManifest(
    const P5IntegratedBenchmarkDescriptor& descriptor,
    const std::string& output_root) {
  P5IntegratedBenchmarkArtifactManifest manifest;
  manifest.output_directory = output_root;
  manifest.required_files = {
      "p5_case_descriptor.json",
      "p5_stage_reports.json",
      "p5_budget_timeseries.csv",
      "p5_profile_timeseries.csv"};
  manifest.visual_artifacts = {
      "p5_" + descriptor.case_id + "_" + descriptor.variant_id +
          "_rz_evolution.gif",
      "p5_" + descriptor.case_id + "_" + descriptor.variant_id +
          "_profiles.png",
      "p5_" + descriptor.case_id + "_" + descriptor.variant_id +
          "_budgets.png"};
  manifest.comparison_only_files = {
      "p5_" + descriptor.case_id + "_variant_comparison.json"};
  manifest.missing_by_design = {
      "external_lilac_reference.csv",
      "one_dimensional_reference_solution.csv"};
  return manifest;
}

}  // namespace dec3d::benchmarks
