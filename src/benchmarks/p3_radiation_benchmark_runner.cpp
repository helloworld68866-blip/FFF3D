#include "benchmarks/p3_radiation_benchmark_runner.hpp"

#include "mesh/spherical/spherical_mesh.hpp"
#include "physics/units/physical_constants.hpp"
#include "radiation/providers/tops_opacity_provider.hpp"
#include "radiation/transport/multigroup_gray_radiation.hpp"
#include "radiation/transport/one_group_gray_diffusion.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/hydro_state/hydro_view.hpp"
#include "transport/diffusion/generic_diffusion.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace dec3d::benchmarks {
namespace {

[[nodiscard]] std::string FailureLine(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const char* reason) {
  std::ostringstream out;
  out << "diagnostic_id=p3.radiation_benchmark.failure"
      << "; benchmark_id=" << descriptor.case_id
      << "; failure_reason=" << reason
      << "; canonical_state_mutated=false";
  return out.str();
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata OneCellGeometry() {
  dec3d::mesh::SphericalGeometryMetadata geometry;
  geometry.valid = true;
  geometry.radial_faces = {0.1, 0.3};
  geometry.theta_faces = {0.5, 1.2};
  geometry.phi_faces = {0.0, 2.0 * dec3d::physics::PhysicsConstantsCGS::pi};
  const double radial =
      (std::pow(geometry.radial_faces[1], 3) -
       std::pow(geometry.radial_faces[0], 3)) /
      3.0;
  const double polar =
      std::cos(geometry.theta_faces[0]) -
      std::cos(geometry.theta_faces[1]);
  const double azimuth = geometry.phi_faces[1] - geometry.phi_faces[0];
  geometry.cell_volumes = {radial * polar * azimuth};
  geometry.global_volume = geometry.cell_volumes[0];
  return geometry;
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy PatchBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::interior_patch_no_origin,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::interior_patch_no_pole,
      DiffusionBoundaryKind::periodic};
}

[[nodiscard]] dec3d::transport::GenericDiffusionBoundaryPolicy FullSphereBoundary() {
  using dec3d::transport::DiffusionBoundaryKind;
  return dec3d::transport::GenericDiffusionBoundaryPolicy{
      DiffusionBoundaryKind::scalar_origin_remap_required,
      DiffusionBoundaryKind::neumann_zero_flux,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::scalar_pole_remap_required,
      DiffusionBoundaryKind::periodic};
}

[[nodiscard]] dec3d::core::Array3D<double> OneCellArray(double value) {
  return dec3d::core::Array3D<double>(1u, 1u, 1u, value);
}

[[nodiscard]] dec3d::state::CanonicalState MakeB0State(std::size_t groups) {
  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{1u, 1u, 1u, groups});
  state.rho(0, 0, 0) = 1.0;
  state.mom_r(0, 0, 0) = 0.0;
  state.mom_theta(0, 0, 0) = 0.0;
  state.mom_phi(0, 0, 0) = 0.0;
  state.e_electron(0, 0, 0) = 40.0;
  state.e_fluid_total(0, 0, 0) = 90.0;
  for (std::size_t g = 0; g < groups; ++g) {
    state.radiation_groups[g](0, 0, 0) = 2.0 + static_cast<double>(g);
  }
  return state;
}

[[nodiscard]] dec3d::radiation::RadiationGroupLayout MakeB0GroupLayout(
    std::size_t groups) {
  std::vector<double> edges;
  edges.reserve(groups + 1u);
  for (std::size_t i = 0; i <= groups; ++i) {
    edges.push_back((1.0 + static_cast<double>(i)) * 1.0e15);
  }
  return dec3d::radiation::MakeExplicitFrequencyRadiationGroupLayout(edges);
}

[[nodiscard]] P3RadiationBenchmarkRunResult RunB0(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy,
    const P3RadiationBenchmarkRunOptions& run_options) {
  P3RadiationBenchmarkRunResult result;
  const auto validation = ValidateP3RadiationBenchmarkDescriptor(descriptor);
  if (!validation.success) {
    result.failure_reason = validation.failure_reason;
    result.failure_diagnostics = validation.failure_diagnostics;
    return result;
  }

  const std::size_t groups = descriptor.group_count;
  auto state = MakeB0State(groups);
  dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions options;
  options.dt_s = run_options.dt_s;
  options.group_layout = MakeB0GroupLayout(groups);
  options.boundary_policy = PatchBoundary();
  options.radiation_boundary_model =
      descriptor.marshak_enabled
          ? dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum
          : dec3d::radiation::OneGroupGrayRadiationBoundaryModel::
                contract_zero_flux_or_scalar_remap;
  options.coefficients.source = "fixed_user_supplied";
  for (std::size_t g = 0; g < groups; ++g) {
    const double kappa = 1.0e-4 * (1.0 + static_cast<double>(g));
    options.coefficients.Dbar_cm2_per_s.push_back(
        OneCellArray(descriptor.marshak_enabled ? (2.0 + static_cast<double>(g)) : 0.0));
    options.coefficients.kappaP_cm_inv.push_back(OneCellArray(kappa));
    options.coefficients.B_erg_per_cm3.push_back(
        OneCellArray(5.0 + static_cast<double>(g)));
  }

  const auto radiation = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
      state,
      OneCellGeometry(),
      options);
  if (!radiation.success) {
    result.failure_reason = radiation.failure_reason;
    result.failure_diagnostics = radiation.failure_diagnostics;
    return result;
  }

  result.success = true;
  result.operator_oracle_executed = true;
  result.delta_radiation_total = radiation.delta_radiation_total_all_groups;
  result.delta_electron_total = radiation.delta_electron_total;
  result.boundary_leak_total = radiation.boundary_leak_total_all_groups;
  result.global_budget_residual =
      radiation.delta_radiation_total_all_groups +
      radiation.delta_electron_total +
      radiation.boundary_leak_total_all_groups;
  result.budget_closed = std::abs(result.global_budget_residual) < 1.0e-8;
  result.group_count_one_regression_passed =
      descriptor.case_kind ==
          P3RadiationBenchmarkCase::b0_group_count_one_regression
          ? groups == 1u && radiation.group_count == 1u
          : false;
  result.source_oracle_matched =
      descriptor.case_kind !=
      P3RadiationBenchmarkCase::b0_group_count_one_regression;

  const auto base = BuildP3RadiationBenchmarkDiagnosticsLine(
      descriptor,
      reference_policy);
  std::ostringstream out;
  out << std::setprecision(17)
      << base
      << "; operator_oracle_executed=true"
      << "; slab_benchmark_executed=false"
      << "; source_oracle_matched="
      << (result.source_oracle_matched ? "true" : "false")
      << "; group_count_one_regression_passed="
      << (result.group_count_one_regression_passed ? "true" : "false")
      << "; budget_closed=" << (result.budget_closed ? "true" : "false")
      << "; delta_radiation_total=" << result.delta_radiation_total
      << "; delta_electron_total=" << result.delta_electron_total
      << "; boundary_leak_total=" << result.boundary_leak_total
      << "; global_budget_residual=" << result.global_budget_residual
      << "; artifact_gather_used=false"
      << "; profile_artifacts_written=false"
      << "; energy_budget_written=false";
  result.report_line = out.str();
  return result;
}

[[nodiscard]] std::size_t EvenPhiCellCount(std::size_t requested) {
  return requested > 1u && requested % 2u == 0u ? requested : 2u;
}

[[nodiscard]] std::string JoinPath(const std::string& root, const std::string& leaf) {
  return (std::filesystem::path(root) / leaf).generic_string();
}

[[nodiscard]] dec3d::mesh::SphericalGeometryMetadata BuildB1Geometry(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationBenchmarkRunOptions& options) {
  return dec3d::mesh::BuildSphericalGeometry(dec3d::mesh::SphericalMeshDescriptor{
      std::max<std::size_t>(1u, options.radial_cells),
      std::max<std::size_t>(1u, options.theta_cells),
      EvenPhiCellCount(options.phi_cells),
      0.0,
      descriptor.r_max_cm});
}

void FillB1ThermodynamicState(
    dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const P3RadiationBenchmarkDescriptor& descriptor) {
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();

  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    const double radius =
        0.5 * (geometry.radial_faces[r] + geometry.radial_faces[r + 1u]);
    const bool inner = radius < descriptor.r0_cm;
    const double rho = inner ? descriptor.rho_inner_g_cm3
                             : descriptor.rho_outer_g_cm3;
    const double te = dec3d::physics::ErgFromKeV(
        inner ? descriptor.te_inner_keV : descriptor.te_outer_keV);
    const double ti = te;
    const double ni = rho / mean_ion_mass;
    const double ne = dec3d::physics::DefaultZbar() * ni;
    const double e_e = ne * te / gamma_minus_one;
    const double e_i = ni * ti / gamma_minus_one;

    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        state.rho(r, theta, phi) = rho;
        state.mom_r(r, theta, phi) = 0.0;
        state.mom_theta(r, theta, phi) = 0.0;
        state.mom_phi(r, theta, phi) = 0.0;
        state.e_electron(r, theta, phi) = e_e;
        state.e_fluid_total(r, theta, phi) = e_e + e_i;
      }
    }
  }
}

void InitializeRadiationFromBlackbody(
    dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationCoefficientProviderResult& provider) {
  for (std::size_t group = 0u; group < state.radiation_groups.size(); ++group) {
    for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
      for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
        for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
          state.radiation_groups[group](r, theta, phi) =
              provider.coefficients.B_erg_per_cm3[group](r, theta, phi);
        }
      }
    }
  }
}

[[nodiscard]] dec3d::radiation::TopsOpacityProviderOptions MakeB1ProviderOptions(
    const P3RadiationBenchmarkRunOptions& options) {
  dec3d::radiation::TopsOpacityProviderOptions provider_options;
  provider_options.table_root = options.tops_table_root;
  provider_options.opacity_interpolation_mode =
      dec3d::radiation::OpacityInterpolationMode::loglog_trilinear;
  provider_options.lookup_energy_mapping_mode =
      dec3d::radiation::OpacityEnergyMappingMode::geometric_group_energy;
  provider_options.density_clip_policy =
      dec3d::radiation::TopsDensityClipPolicy::hard_fail;
  return provider_options;
}

struct B1ProfileRow {
  double radius_cm{0.0};
  double rho_g_cm3{0.0};
  double te_keV{0.0};
  double sum_Ug_erg_cm3{0.0};
};

[[nodiscard]] std::vector<B1ProfileRow> BuildB1RadialProfile(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  std::vector<B1ProfileRow> rows;
  rows.reserve(state.layout.radial_cells);
  const double gamma_minus_one = dec3d::state::HydroIdealGasGamma() - 1.0;
  const double mean_ion_mass = dec3d::physics::DefaultMeanDTIonMassG();
  for (std::size_t r = 0u; r < state.layout.radial_cells; ++r) {
    B1ProfileRow row;
    row.radius_cm = 0.5 * (geometry.radial_faces[r] + geometry.radial_faces[r + 1u]);
    double rho_sum = 0.0;
    double te_sum = 0.0;
    double ug_sum = 0.0;
    std::size_t count = 0u;
    for (std::size_t theta = 0u; theta < state.layout.theta_cells; ++theta) {
      for (std::size_t phi = 0u; phi < state.layout.phi_cells; ++phi) {
        const double rho = state.rho(r, theta, phi);
        const double ne = dec3d::physics::DefaultZbar() * rho / mean_ion_mass;
        const double te_erg = gamma_minus_one * state.e_electron(r, theta, phi) / ne;
        double sum_ug = 0.0;
        for (const auto& group : state.radiation_groups) {
          sum_ug += group(r, theta, phi);
        }
        rho_sum += rho;
        te_sum += dec3d::physics::KeVFromErg(te_erg);
        ug_sum += sum_ug;
        ++count;
      }
    }
    const double inv_count = count > 0u ? 1.0 / static_cast<double>(count) : 0.0;
    row.rho_g_cm3 = rho_sum * inv_count;
    row.te_keV = te_sum * inv_count;
    row.sum_Ug_erg_cm3 = ug_sum * inv_count;
    rows.push_back(row);
  }
  return rows;
}

[[nodiscard]] bool WriteTextFile(const std::string& path, const std::string& text) {
  std::ofstream output(path, std::ios::binary);
  if (!output) {
    return false;
  }
  output << text;
  return static_cast<bool>(output);
}

void AppendU32(std::vector<unsigned char>& bytes, std::uint32_t value) {
  bytes.push_back(static_cast<unsigned char>((value >> 24u) & 0xffu));
  bytes.push_back(static_cast<unsigned char>((value >> 16u) & 0xffu));
  bytes.push_back(static_cast<unsigned char>((value >> 8u) & 0xffu));
  bytes.push_back(static_cast<unsigned char>(value & 0xffu));
}

[[nodiscard]] std::uint32_t Crc32(
    const unsigned char* data,
    std::size_t size,
    std::uint32_t crc = 0xffffffffu) {
  for (std::size_t i = 0u; i < size; ++i) {
    crc ^= data[i];
    for (int bit = 0; bit < 8; ++bit) {
      crc = (crc >> 1u) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
  }
  return crc;
}

[[nodiscard]] std::uint32_t Adler32(const std::vector<unsigned char>& data) {
  constexpr std::uint32_t mod = 65521u;
  std::uint32_t a = 1u;
  std::uint32_t b = 0u;
  for (const unsigned char byte : data) {
    a = (a + byte) % mod;
    b = (b + a) % mod;
  }
  return (b << 16u) | a;
}

void AppendChunk(
    std::vector<unsigned char>& png,
    const char type[4],
    const std::vector<unsigned char>& data) {
  AppendU32(png, static_cast<std::uint32_t>(data.size()));
  const std::size_t type_offset = png.size();
  png.insert(png.end(), type, type + 4);
  png.insert(png.end(), data.begin(), data.end());
  std::uint32_t crc = Crc32(png.data() + type_offset, png.size() - type_offset);
  crc ^= 0xffffffffu;
  AppendU32(png, crc);
}

[[nodiscard]] std::vector<unsigned char> BuildStoredZlib(
    const std::vector<unsigned char>& raw) {
  std::vector<unsigned char> zlib;
  zlib.push_back(0x78u);
  zlib.push_back(0x01u);
  std::size_t offset = 0u;
  while (offset < raw.size()) {
    const std::size_t chunk = std::min<std::size_t>(65535u, raw.size() - offset);
    const bool final = offset + chunk == raw.size();
    zlib.push_back(final ? 0x01u : 0x00u);
    const auto len = static_cast<std::uint16_t>(chunk);
    const auto nlen = static_cast<std::uint16_t>(~len);
    zlib.push_back(static_cast<unsigned char>(len & 0xffu));
    zlib.push_back(static_cast<unsigned char>((len >> 8u) & 0xffu));
    zlib.push_back(static_cast<unsigned char>(nlen & 0xffu));
    zlib.push_back(static_cast<unsigned char>((nlen >> 8u) & 0xffu));
    zlib.insert(zlib.end(), raw.begin() + static_cast<std::ptrdiff_t>(offset),
                raw.begin() + static_cast<std::ptrdiff_t>(offset + chunk));
    offset += chunk;
  }
  AppendU32(zlib, Adler32(raw));
  return zlib;
}

void SetPixel(
    std::vector<unsigned char>& image,
    int width,
    int height,
    int x,
    int y,
    unsigned char r,
    unsigned char g,
    unsigned char b) {
  if (x < 0 || y < 0 || x >= width || y >= height) {
    return;
  }
  const std::size_t index =
      4u * (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
            static_cast<std::size_t>(x));
  image[index] = r;
  image[index + 1u] = g;
  image[index + 2u] = b;
  image[index + 3u] = 255u;
}

void DrawLine(
    std::vector<unsigned char>& image,
    int width,
    int height,
    int x0,
    int y0,
    int x1,
    int y1,
    unsigned char r,
    unsigned char g,
    unsigned char b) {
  const int dx = std::abs(x1 - x0);
  const int sx = x0 < x1 ? 1 : -1;
  const int dy = -std::abs(y1 - y0);
  const int sy = y0 < y1 ? 1 : -1;
  int err = dx + dy;
  while (true) {
    SetPixel(image, width, height, x0, y0, r, g, b);
    if (x0 == x1 && y0 == y1) {
      break;
    }
    const int e2 = 2 * err;
    if (e2 >= dy) {
      err += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      err += dx;
      y0 += sy;
    }
  }
}

[[nodiscard]] bool WriteLinePlotPng(
    const std::string& path,
    const std::vector<B1ProfileRow>& profile,
    bool plot_temperature) {
  constexpr int width = 720;
  constexpr int height = 440;
  constexpr int left = 56;
  constexpr int right = 24;
  constexpr int top = 24;
  constexpr int bottom = 48;
  std::vector<unsigned char> image(
      static_cast<std::size_t>(width * height * 4),
      255u);
  for (int y = top; y < height - bottom; ++y) {
    SetPixel(image, width, height, left, y, 40, 40, 40);
  }
  for (int x = left; x < width - right; ++x) {
    SetPixel(image, width, height, x, height - bottom, 40, 40, 40);
  }

  double xmin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();
  for (const auto& row : profile) {
    const double y = plot_temperature ? row.te_keV : row.sum_Ug_erg_cm3;
    xmin = std::min(xmin, row.radius_cm);
    xmax = std::max(xmax, row.radius_cm);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
  if (!std::isfinite(xmin) || !std::isfinite(ymin) || xmax <= xmin) {
    return false;
  }
  if (ymax <= ymin) {
    const double pad = std::max(1.0, std::abs(ymin)) * 0.05;
    ymin -= pad;
    ymax += pad;
  }

  auto map_x = [&](double x) {
    const double t = (x - xmin) / (xmax - xmin);
    return left + static_cast<int>(std::llround(
                      t * static_cast<double>(width - left - right - 1)));
  };
  auto map_y = [&](double y) {
    const double t = (y - ymin) / (ymax - ymin);
    return height - bottom - static_cast<int>(std::llround(
                               t * static_cast<double>(height - top - bottom - 1)));
  };

  for (std::size_t i = 1u; i < profile.size(); ++i) {
    const double y0 = plot_temperature ? profile[i - 1u].te_keV
                                       : profile[i - 1u].sum_Ug_erg_cm3;
    const double y1 = plot_temperature ? profile[i].te_keV
                                       : profile[i].sum_Ug_erg_cm3;
    DrawLine(image, width, height,
             map_x(profile[i - 1u].radius_cm), map_y(y0),
             map_x(profile[i].radius_cm), map_y(y1),
             plot_temperature ? 210u : 20u,
             plot_temperature ? 70u : 100u,
             plot_temperature ? 40u : 210u);
  }

  std::vector<unsigned char> raw;
  raw.reserve(static_cast<std::size_t>(height * (1 + width * 4)));
  for (int y = 0; y < height; ++y) {
    raw.push_back(0u);
    const std::size_t offset =
        static_cast<std::size_t>(y) * static_cast<std::size_t>(width) * 4u;
    raw.insert(raw.end(), image.begin() + static_cast<std::ptrdiff_t>(offset),
               image.begin() + static_cast<std::ptrdiff_t>(offset + width * 4));
  }

  std::vector<unsigned char> png = {
      0x89u, 'P', 'N', 'G', '\r', '\n', 0x1au, '\n'};
  std::vector<unsigned char> ihdr;
  AppendU32(ihdr, width);
  AppendU32(ihdr, height);
  ihdr.push_back(8u);
  ihdr.push_back(6u);
  ihdr.push_back(0u);
  ihdr.push_back(0u);
  ihdr.push_back(0u);
  AppendChunk(png, "IHDR", ihdr);
  AppendChunk(png, "IDAT", BuildStoredZlib(raw));
  AppendChunk(png, "IEND", {});

  std::ofstream output(path, std::ios::binary);
  if (!output) {
    return false;
  }
  output.write(reinterpret_cast<const char*>(png.data()),
               static_cast<std::streamsize>(png.size()));
  return static_cast<bool>(output);
}

[[nodiscard]] bool WriteB1Artifacts(
    const std::string& output_directory,
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy,
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::radiation::RadiationCoefficientProviderResult& provider,
    const dec3d::radiation::MultigroupGrayRadiationMatterCouplingResult& radiation,
    const P3RadiationBenchmarkRunResult& result_without_diagnostics) {
  std::filesystem::create_directories(output_directory);
  const auto profile = BuildB1RadialProfile(state, geometry);

  std::ostringstream metadata;
  metadata << std::setprecision(17)
           << "{\n"
           << "  \"diagnostic_id\": \"p3.radiation_benchmark.metadata\",\n"
           << "  \"case_id\": \"" << descriptor.case_id << "\",\n"
           << "  \"benchmark_layer\": \"P3-B1\",\n"
           << "  \"claim_level\": \"" << descriptor.claim_level << "\",\n"
           << "  \"parity_claim_allowed\": false,\n"
           << "  \"lilac_reference_available\": "
           << (reference_policy.lilac_reference_available ? "true" : "false") << ",\n"
           << "  \"radiation_boundary_model\": \"" << descriptor.radiation_boundary_model
           << "\",\n"
           << "  \"group_count\": " << descriptor.group_count << ",\n"
           << "  \"group_edges_source\": \"RadiationGroupLayout\",\n"
           << "  \"opacity_energy_coordinate_source\": \"TOPS photon grid\",\n"
           << "  \"opacity_source_thesis_exact_match\": false,\n"
           << "  \"initial_ti_policy\": \"" << descriptor.initial_ti_policy << "\",\n"
           << "  \"initial_radiation_policy\": \""
           << descriptor.initial_radiation_policy << "\",\n"
           << "  \"radiation_flux_limiter\": \"" << descriptor.radiation_flux_limiter
           << "\"\n"
           << "}\n";
  if (!WriteTextFile(JoinPath(output_directory, "metadata.json"), metadata.str())) {
    return false;
  }

  std::ostringstream budget;
  budget << std::setprecision(17)
         << "diagnostic_id=p3.radiation_benchmark.budget_summary\n"
         << "delta_radiation_total=" << result_without_diagnostics.delta_radiation_total << "\n"
         << "delta_electron_total=" << result_without_diagnostics.delta_electron_total << "\n"
         << "boundary_leak_total=" << result_without_diagnostics.boundary_leak_total << "\n"
         << "global_budget_residual=" << result_without_diagnostics.global_budget_residual << "\n"
         << "budget_closed=" << (result_without_diagnostics.budget_closed ? "true" : "false") << "\n";
  if (!WriteTextFile(JoinPath(output_directory, "budget_summary.txt"), budget.str())) {
    return false;
  }

  std::ostringstream profile_csv;
  profile_csv << std::setprecision(17)
              << "radius_cm,rho_g_cm3,Te_keV,sum_Ug_erg_cm3\n";
  for (const auto& row : profile) {
    profile_csv << row.radius_cm << ',' << row.rho_g_cm3 << ','
                << row.te_keV << ',' << row.sum_Ug_erg_cm3 << '\n';
  }
  if (!WriteTextFile(JoinPath(output_directory, "radial_profiles_step000001.csv"),
                     profile_csv.str())) {
    return false;
  }

  std::ostringstream group_budget;
  group_budget << std::setprecision(17)
               << "group,delta_radiation_total,source_gain_radiation_total,boundary_leak_total\n";
  for (std::size_t group = 0u; group < radiation.group_count; ++group) {
    group_budget << group << ','
                 << radiation.delta_radiation_total_by_group[group] << ','
                 << radiation.source_gain_radiation_total_by_group[group] << ','
                 << radiation.boundary_leak_total_by_group[group] << '\n';
  }
  if (!WriteTextFile(JoinPath(output_directory, "group_energy_budget.csv"),
                     group_budget.str())) {
    return false;
  }

  std::ostringstream provider_diag;
  provider_diag << provider.report_line << '\n';
  if (!WriteTextFile(JoinPath(output_directory, "provider_diagnostics.txt"),
                     provider_diag.str())) {
    return false;
  }

  return WriteLinePlotPng(JoinPath(output_directory, "Te_profile.png"), profile, true) &&
         WriteLinePlotPng(JoinPath(output_directory, "sum_Ug_profile.png"), profile, false);
}

[[nodiscard]] P3RadiationBenchmarkRunResult RunB1(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy,
    const P3RadiationBenchmarkRunOptions& run_options) {
  P3RadiationBenchmarkRunResult result;
  const auto validation = ValidateP3RadiationBenchmarkDescriptor(descriptor);
  if (!validation.success) {
    result.failure_reason = validation.failure_reason;
    result.failure_diagnostics = validation.failure_diagnostics;
    return result;
  }

  const auto group_layout = MakeP3B1TwelveGroupLayout();
  auto geometry = BuildB1Geometry(descriptor, run_options);
  if (!geometry.is_valid()) {
    result.failure_reason = "failed to build P3-B1 spherical slab geometry";
    result.failure_diagnostics = FailureLine(descriptor, result.failure_reason.c_str());
    return result;
  }

  auto state = dec3d::state::CanonicalState::Create(
      dec3d::state::CanonicalStateLayout{
          geometry.radial_faces.size() - 1u,
          geometry.theta_faces.size() - 1u,
          geometry.phi_faces.size() - 1u,
          group_layout.group_count});
  FillB1ThermodynamicState(state, geometry, descriptor);

  const auto provider_options = MakeB1ProviderOptions(run_options);
  auto provider = dec3d::radiation::BuildRadiationCoefficientArrays(
      state,
      geometry,
      group_layout,
      provider_options);
  if (!provider.success ||
      !dec3d::radiation::ValidateRadiationCoefficientProviderDiagnostics(provider)) {
    result.failure_reason = provider.failure_reason.empty()
                                ? "P3-B1 opacity provider failed"
                                : provider.failure_reason;
    result.failure_diagnostics = provider.failure_diagnostics.empty()
                                     ? provider.report_line
                                     : provider.failure_diagnostics;
    return result;
  }

  InitializeRadiationFromBlackbody(state, provider);

  dec3d::radiation::MultigroupGrayRadiationMatterCouplingResult radiation;
  const std::size_t step_count = std::max<std::size_t>(1u, run_options.step_count);
  result.step_wall_seconds.clear();
  result.step_wall_seconds.reserve(step_count);

  const auto total_begin = std::chrono::steady_clock::now();
  for (std::size_t step = 0u; step < step_count; ++step) {
    if (step > 0u) {
      provider = dec3d::radiation::BuildRadiationCoefficientArrays(
          state,
          geometry,
          group_layout,
          provider_options);
      if (!provider.success ||
          !dec3d::radiation::ValidateRadiationCoefficientProviderDiagnostics(provider)) {
        result.failure_reason = provider.failure_reason.empty()
                                    ? "P3-B1 opacity provider failed"
                                    : provider.failure_reason;
        result.failure_diagnostics = provider.failure_diagnostics.empty()
                                         ? provider.report_line
                                         : provider.failure_diagnostics;
        return result;
      }
    }

    dec3d::radiation::MultigroupGrayRadiationMatterCouplingOptions radiation_options;
    radiation_options.dt_s = run_options.dt_s;
    radiation_options.group_layout = group_layout;
    radiation_options.coefficients = provider.coefficients;
    radiation_options.boundary_policy = FullSphereBoundary();
    radiation_options.radiation_boundary_model =
        dec3d::radiation::OneGroupGrayRadiationBoundaryModel::thesis_marshak_vacuum;
    radiation_options.radiation_energy_floor_erg_per_cm3 = 0.0;
    radiation_options.electron_energy_floor_erg_per_cm3 = 0.0;
    radiation_options.ion_energy_floor_erg_per_cm3 = 0.0;
    radiation_options.serial_reference_options.row_limit =
        std::max<std::size_t>(64u,
                              state.layout.radial_cells *
                                  state.layout.theta_cells *
                                  state.layout.phi_cells);
    radiation_options.serial_reference_options.residual_tolerance = 1.0e30;

    const auto step_begin = std::chrono::steady_clock::now();
    radiation = dec3d::radiation::ApplyMultigroupGrayRadiationMatterCoupling(
        state,
        geometry,
        radiation_options);
    const auto step_end = std::chrono::steady_clock::now();
    result.step_wall_seconds.push_back(
        std::chrono::duration<double>(step_end - step_begin).count());

    if (!radiation.success ||
        !dec3d::radiation::ValidateMultigroupGrayRadiationMatterCouplingDiagnostics(
            radiation)) {
      result.failure_reason = radiation.failure_reason.empty()
                                  ? "P3-B1 radiation solve failed"
                                  : radiation.failure_reason;
      result.failure_diagnostics = radiation.failure_diagnostics.empty()
                                       ? radiation.report_line
                                       : radiation.failure_diagnostics;
      return result;
    }

    result.delta_radiation_total += radiation.delta_radiation_total_all_groups;
    result.delta_electron_total += radiation.delta_electron_total;
    result.boundary_leak_total += radiation.boundary_leak_total_all_groups;
    result.global_budget_residual +=
        radiation.global_radiation_plus_electron_residual;
  }
  result.total_wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - total_begin)
          .count();

  result.success = true;
  result.slab_benchmark_executed = true;
  const double budget_scale = std::max(
      1.0,
      std::abs(result.delta_radiation_total) +
          std::abs(result.delta_electron_total) +
          std::abs(result.boundary_leak_total));
  result.budget_closed = result.global_budget_residual <= 1.0e-10 * budget_scale;
  result.output_directory = JoinPath(run_options.output_root, descriptor.case_id);
  if (run_options.write_artifacts) {
    if (!WriteB1Artifacts(result.output_directory,
                          descriptor,
                          reference_policy,
                          state,
                          geometry,
                          provider,
                          radiation,
                          result)) {
      result.success = false;
      result.failure_reason = "failed to write P3-B1 artifact bundle";
      result.failure_diagnostics =
          FailureLine(descriptor, result.failure_reason.c_str());
      return result;
    }
    result.profile_artifacts_written = true;
    result.budget_artifacts_written = true;
    result.metadata_written = true;
    result.diagnostics_written = true;
  }

  const auto base = BuildP3RadiationBenchmarkDiagnosticsLine(
      descriptor,
      reference_policy);
  std::ostringstream out;
  out << std::setprecision(17)
      << base
      << "; operator_oracle_executed=false"
      << "; slab_benchmark_executed=true"
      << "; source_oracle_matched=false"
      << "; group_count_one_regression_passed=false"
      << "; step_count=" << step_count
      << "; effective_radial_cells=" << state.layout.radial_cells
      << "; effective_theta_cells=" << state.layout.theta_cells
      << "; effective_phi_cells=" << state.layout.phi_cells
      << "; provider_report_present=true"
      << "; radiation_report_present=true"
      << "; provider_coefficients_source=" << provider.coefficients.source
      << "; radiation_updated_fields=radiation_groups,e_electron,e_fluid_total"
      << "; budget_closed=" << (result.budget_closed ? "true" : "false")
      << "; delta_radiation_total=" << result.delta_radiation_total
      << "; delta_electron_total=" << result.delta_electron_total
      << "; boundary_leak_total=" << result.boundary_leak_total
      << "; global_budget_residual=" << result.global_budget_residual
      << "; total_wall_seconds=" << result.total_wall_seconds
      << "; step_wall_seconds=";
  for (std::size_t i = 0u; i < result.step_wall_seconds.size(); ++i) {
    if (i > 0u) {
      out << ',';
    }
    out << result.step_wall_seconds[i];
  }
  out
      << "; artifact_gather_used=false"
      << "; output_directory=" << result.output_directory
      << "; profile_artifacts_written="
      << (result.profile_artifacts_written ? "true" : "false")
      << "; energy_budget_written="
      << (result.budget_artifacts_written ? "true" : "false")
      << "; metadata_written=" << (result.metadata_written ? "true" : "false")
      << "; diagnostics_written=" << (result.diagnostics_written ? "true" : "false")
      << "; nested_provider_report={" << provider.report_line << "}"
      << "; nested_radiation_report={" << radiation.report_line << "}";
  result.report_line = out.str();
  if (run_options.write_artifacts &&
      !WriteTextFile(JoinPath(result.output_directory, "diagnostics.txt"),
                     result.report_line + "\n")) {
    result.success = false;
    result.failure_reason = "failed to write P3-B1 diagnostics artifact";
    result.failure_diagnostics = FailureLine(descriptor, result.failure_reason.c_str());
    return result;
  }
  return result;
}

}  // namespace

P3RadiationBenchmarkRunResult RunP3RadiationBenchmarkCase(
    const P3RadiationBenchmarkDescriptor& descriptor,
    const P3RadiationReferencePolicy& reference_policy,
    const P3RadiationBenchmarkRunOptions& options) noexcept {
  try {
    if (descriptor.layer == P3RadiationBenchmarkLayer::b0_operator_budget_oracle) {
      return RunB0(descriptor, reference_policy, options);
    }
    return RunB1(descriptor, reference_policy, options);
  } catch (...) {
    P3RadiationBenchmarkRunResult result;
    result.failure_reason = "unexpected P3 radiation benchmark exception";
    result.failure_diagnostics = FailureLine(descriptor, result.failure_reason.c_str());
    return result;
  }
}

}  // namespace dec3d::benchmarks
