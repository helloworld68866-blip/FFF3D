# Axisymmetric 2D First Version Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a first-class no-ALE `axisymmetric_2d` mode with `phi_cells=1`, all production stages `H,T,E,R,A`, no phi PDE operations, and full-azimuth geometry.

**Architecture:** Keep the existing `Array3D(r,theta,phi)` state model and make `phi=1` a constrained runtime mode instead of adding `Array2D`. The input deck declares dimensionality; geometry keeps full `2*pi`; remap, hydro, and diffusion paths explicitly skip phi-direction operators while preserving spherical `r,theta` metric/source terms.

**Tech Stack:** C++20, CMake/MSBuild Release builds, Microsoft MPI, HYPRE, existing DEC3D contract tests. Use PowerShell commands, not `rg`.

---

## File Map

Modify these files:

- `src/io/input_deck.hpp`: add mesh dimensionality enum and axisymmetric query helpers.
- `src/io/input_deck.cpp`: parse `mesh.dimensionality`, validate `phi_cells=1` and full angular ranges, reject ALE/moving mesh and `m!=0`, and report `phi_sweep_configured`.
- `tests/contract/test_p5_io_input_deck.cpp`: add parser acceptance/rejection coverage.
- `tests/contract/test_spherical_geometry.cpp`: add `phi_cells=1` full-volume and outer radial surface area tests with a local face-area helper.
- `src/mesh/boundary/spherical_scalar_remap.cpp`: allow half-turn topology for `phi_cells=1` and map half-turn to phi `0`.
- `src/hydro/driver/theta_pole_boundary.cpp`: allow `MapPhiAcrossPole(phi, 1)` to return `0`.
- `src/hydro/driver/hydro_boundary_ghosts.cpp`: allow origin and pole ghost fills with `phi_cells=1`.
- `tests/contract/test_spherical_scalar_remap.cpp`, `tests/contract/test_hydro_theta_pole_boundary.cpp`, `tests/contract/test_hydro_origin_radial_ghost_remap.cpp`: add remap and parity tests.
- `src/transport/diffusion/generic_diffusion.hpp/.cpp`: add matrix diagnostics for phi coupling count, pre-CSR phi self-neighbor attempts, duplicate columns, and interior row width; preserve existing `phi_cells<=1` skip.
- `src/transport/diffusion/hypre_distributed_diffusion_solver.hpp/.cpp`: add the same distributed diagnostics and tests for no phi self-neighbor.
- `tests/contract/test_generic_implicit_diffusion.cpp`, `tests/contract/test_hypre_distributed_diffusion_solver.cpp`: add serial and distributed axisymmetric matrix tests.
- `src/initialization/profile_initializer.cpp`: report perturbation normalization as raw Legendre `P_l`, and fail axisymmetric mode if profile initializes nonzero `mom_phi`.
- `tests/contract/test_p5_io_profile_initializer.cpp`: add initialization diagnostics and `mom_phi` rejection coverage.
- `src/app/dec3d_app.cpp`: configure axisymmetric hydro options, runtime reports, stage invariants, and restart dimensionality safety.
- `src/hydro/driver/hydro_operator.cpp`: make dt evidence/report include active directions when phi sweep is disabled.
- `src/hydro/driver/static_grid_hydro.cpp`: preserve the existing zero phi-sweep timing path and expose `phi_sweep_executed=false` through the runtime report.
- `tests/contract/test_p5_runtime_all_stages.cpp`: add a serial H-only no-ALE smoke case that verifies hydro direction diagnostics.
- `tests/contract/test_p5_runtime_distributed_all_stages.cpp`: add HYPRE no-ALE all-stage axisymmetric smoke using the existing MPI app-command helper pattern.
- `src/state/diagnostics/axisymmetric_vector_diagnostics.hpp`: add production Cartesian momentum reducers for analytic axisymmetric and numeric full-3D integration.
- `tests/contract/test_axisymmetric_2d_regression.cpp`: create a focused 2D-vs-3D m=0 structural regression and production vector diagnostic test.
- `CMakeLists.txt`: register `dec3d_contract_axisymmetric_2d_regression`.
- `cases/axisymmetric_2d_noale_allstages_smoke.in`: create a tiny all-stage input deck.
- `cases/axisymmetric_2d_noale_smoke.pro`: create a tiny profile file.

Do not modify:

- `main` branch.
- ALE implementation files except for input validation rejecting axisymmetric ALE.
- Physics equations, opacity interpolation, alpha source, solver tolerances, or radiation group logic.
- Output formats beyond metadata/diagnostics needed for safety.

---

## Task 0: Preflight

**Files:**
- Inspect only.

- [ ] **Step 1: Verify branch and clean starting point**

Run:

```powershell
git -C 'F:\dec3d' status --short --branch
```

Expected:

- branch is `codex/axisymmetric-2d`;
- branch is not `main`;
- no unrelated local edits are present. If unrelated edits are present, do not revert them; record them in the task notes and avoid touching those files.

- [ ] **Step 2: Verify HYPRE build configuration**

Run:

```powershell
Test-Path 'F:\dec3d\build-hypre\CMakeCache.txt'
Select-String -Path 'F:\dec3d\build-hypre\CMakeCache.txt' -Pattern 'DEC3D_ENABLE_HYPRE:BOOL=ON'
```

Expected: first command prints `True`; second command prints the matching cache line.

- [ ] **Step 3: Verify Microsoft MPI launcher**

Run:

```powershell
Test-Path 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe'
```

Expected: prints `True`.

- [ ] **Step 4: Verify TOPS opacity fixture**

Runtime radiation finds TOPS tables by searching for `data/opacities/tops_dt_2026_04_27` from the current path, its parent, and the input deck path ancestry. Verify the fixture exists:

```powershell
Test-Path 'F:\dec3d\data\opacities\tops_dt_2026_04_27\metadata.json'
Test-Path 'F:\dec3d\data\opacities\tops_dt_2026_04_27\multigroup_opacities.csv'
```

Expected: both commands print `True`. If either command prints `False`, do not run all-stage radiation smoke; mark the all-stage smoke blocked by `TOPS opacity fixture missing` and continue only through H-only tests.

- [ ] **Step 5: Run current baseline focused tests**

Run before modifying code:

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_io_input_deck dec3d_contract_spherical_geometry dec3d_contract_spherical_scalar_remap dec3d_contract_hydro_theta_pole_boundary dec3d_contract_hydro_origin_radial_ghost_remap dec3d_contract_generic_implicit_diffusion dec3d_contract_hypre_distributed_diffusion_solver dec3d_contract_p5_io_profile_initializer dec3d_contract_p5_runtime_all_stages dec3d_contract_p5_runtime_distributed_all_stages dec3d_contract_p5_io_checkpoint_files --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_input_deck.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_geometry.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_scalar_remap.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_theta_pole_boundary.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_origin_radial_ghost_remap.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_generic_implicit_diffusion.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_hypre_distributed_diffusion_solver.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_profile_initializer.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_all_stages.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_distributed_all_stages.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_checkpoint_files.exe'
```

Expected: every command exits `0`. Record any baseline failure before starting Task 1; do not treat existing baseline failures as axisymmetric implementation failures.

---

## Task 1: Input Deck Dimensionality Contract

**Files:**
- Modify: `src/io/input_deck.hpp`
- Modify: `src/io/input_deck.cpp`
- Test: `tests/contract/test_p5_io_input_deck.cpp`

- [ ] **Step 1: Add failing parser tests**

Append helper functions near `ValidDeckText()` in `tests/contract/test_p5_io_input_deck.cpp`:

```cpp
std::string ReplaceAll(std::string text, const std::string& from, const std::string& to) {
  std::size_t pos = 0;
  while ((pos = text.find(from, pos)) != std::string::npos) {
    text.replace(pos, from.size(), to);
    pos += to.size();
  }
  return text;
}

std::string AxisymmetricDeckText() {
  auto deck = ValidDeckText();
  deck = ReplaceAll(deck, "geometry = spherical # geometry",
                    "geometry = spherical # geometry\ndimensionality = axisymmetric_2d # dimensionality");
  deck = ReplaceAll(deck, "phi_cells = 4 # phi cell count", "phi_cells = 1 # phi cell count");
  deck = ReplaceAll(deck, "moving_mesh = true # moving mesh", "moving_mesh = false # moving mesh");
  return deck;
}
```

Add test blocks before cleanup:

```cpp
{
  const auto axisym_path = root / "axisymmetric_2d_ok.in";
  WriteText(axisym_path, AxisymmetricDeckText());
  const auto deck = dec3d::io::LoadInputDeck(axisym_path);
  DEC3D_CHECK(deck.success);
  DEC3D_CHECK(deck.config.mesh.dimensionality ==
              dec3d::io::MeshDimensionality::axisymmetric_2d);
  DEC3D_CHECK(deck.config.mesh.phi_cells == 1u);
  DEC3D_CHECK(!deck.config.mesh.moving_mesh);
  DEC3D_CHECK(deck.report_line.find("mesh_dimensionality=axisymmetric_2d") !=
              std::string::npos);
  DEC3D_CHECK(deck.report_line.find("azimuthal_weight=2pi") != std::string::npos);
  DEC3D_CHECK(deck.report_line.find("phi_sweep_configured=false") !=
              std::string::npos);
}

{
  const auto bad_phi_path = root / "axisymmetric_2d_bad_phi.in";
  WriteText(bad_phi_path,
            ReplaceAll(AxisymmetricDeckText(), "phi_cells = 1 # phi cell count",
                       "phi_cells = 2 # phi cell count"));
  const auto deck = dec3d::io::LoadInputDeck(bad_phi_path);
  DEC3D_CHECK(!deck.success);
  DEC3D_CHECK(deck.failure_reason.find("axisymmetric_2d requires phi_cells = 1") !=
              std::string::npos);
}

{
  const auto bad_ale_path = root / "axisymmetric_2d_bad_ale.in";
  WriteText(bad_ale_path,
            ReplaceAll(AxisymmetricDeckText(), "moving_mesh = false # moving mesh",
                       "moving_mesh = true # moving mesh"));
  const auto deck = dec3d::io::LoadInputDeck(bad_ale_path);
  DEC3D_CHECK(!deck.success);
  DEC3D_CHECK(deck.failure_reason.find("axisymmetric_2d requires moving_mesh=false") !=
              std::string::npos);
}

{
  const auto bad_m_path = root / "axisymmetric_2d_bad_m.in";
  auto deck_text = AxisymmetricDeckText();
  deck_text += R"ini(

[perturbation]
enabled = true
type = single_mode_radial_velocity
ell = 2
m = 1
amplitude = 0.1
r0_cm = 1.0e-3
target = radial_velocity_cm_s
)ini";
  WriteText(bad_m_path, deck_text);
  const auto deck = dec3d::io::LoadInputDeck(bad_m_path);
  DEC3D_CHECK(!deck.success);
  DEC3D_CHECK(deck.failure_reason.find(
                  "axisymmetric_2d only supports m = 0 perturbations") !=
              std::string::npos);
}

{
  const auto bad_phi_range_path = root / "axisymmetric_2d_bad_phi_range.in";
  WriteText(bad_phi_range_path,
            ReplaceAll(AxisymmetricDeckText(), "phi_max = 6.283185307179586 # phi max",
                       "phi_max = 3.141592653589793 # phi max"));
  const auto deck = dec3d::io::LoadInputDeck(bad_phi_range_path);
  DEC3D_CHECK(!deck.success);
  DEC3D_CHECK(deck.failure_reason.find(
                  "axisymmetric_2d requires full theta and phi ranges") !=
              std::string::npos);
}

{
  const auto full3d_bad_path = root / "full3d_bad_phi1.in";
  WriteText(full3d_bad_path,
            ReplaceAll(ValidDeckText(), "phi_cells = 4 # phi cell count",
                       "phi_cells = 1 # phi cell count"));
  const auto deck = dec3d::io::LoadInputDeck(full3d_bad_path);
  DEC3D_CHECK(!deck.success);
  DEC3D_CHECK(deck.failure_reason.find("full_3d requires positive even phi_cells") !=
              std::string::npos);
}
```

- [ ] **Step 2: Run failing test**

Run:

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_io_input_deck --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_input_deck.exe'
```

Expected before implementation: compile failure because `MeshDimensionality` does not exist.

- [ ] **Step 3: Implement minimal input contract**

In `src/io/input_deck.hpp`, add before `MeshConfig`:

```cpp
enum class MeshDimensionality {
  full_3d,
  axisymmetric_2d,
};

[[nodiscard]] const char* ToString(MeshDimensionality dimensionality) noexcept;
[[nodiscard]] bool IsAxisymmetric2D(MeshDimensionality dimensionality) noexcept;
```

Add to `MeshConfig`:

```cpp
MeshDimensionality dimensionality{MeshDimensionality::full_3d};
```

In `src/io/input_deck.cpp`, add namespace helpers:

```cpp
[[nodiscard]] bool ParseMeshDimensionality(
    const std::map<std::string, std::string>& values,
    MeshDimensionality& out,
    std::string& failure) {
  const auto it = values.find("dimensionality");
  if (it == values.end() || it->second == "full_3d") {
    out = MeshDimensionality::full_3d;
    return true;
  }
  if (it->second == "axisymmetric_2d") {
    out = MeshDimensionality::axisymmetric_2d;
    return true;
  }
  failure = "invalid mesh.dimensionality";
  return false;
}
```

Add public functions near the bottom:

```cpp
const char* ToString(MeshDimensionality dimensionality) noexcept {
  switch (dimensionality) {
    case MeshDimensionality::full_3d:
      return "full_3d";
    case MeshDimensionality::axisymmetric_2d:
      return "axisymmetric_2d";
  }
  return "full_3d";
}

bool IsAxisymmetric2D(MeshDimensionality dimensionality) noexcept {
  return dimensionality == MeshDimensionality::axisymmetric_2d;
}
```

Parse dimensionality immediately after mesh section values are read:

```cpp
if (!ParseMeshDimensionality(*mesh, cfg.mesh.dimensionality, failure)) {
  return FailDeck(failure);
}
```

Replace the current mesh validation with:

```cpp
if (cfg.mesh.geometry != "spherical" ||
    !(cfg.mesh.radial_max_cm > cfg.mesh.radial_min_cm)) {
  return FailDeck("invalid mesh");
}
if (IsAxisymmetric2D(cfg.mesh.dimensionality)) {
  if (cfg.mesh.phi_cells != 1u) {
    return FailDeck("axisymmetric_2d requires phi_cells = 1");
  }
  if (cfg.mesh.moving_mesh) {
    return FailDeck("axisymmetric_2d requires moving_mesh=false");
  }
  constexpr double kPi = 3.141592653589793238462643383279502884;
  constexpr double kRangeTol = 1.0e-12;
  if (std::abs(cfg.mesh.theta_min - 0.0) > kRangeTol ||
      std::abs(cfg.mesh.theta_max - kPi) > kRangeTol ||
      std::abs(cfg.mesh.phi_min - 0.0) > kRangeTol ||
      std::abs(cfg.mesh.phi_max - 2.0 * kPi) > kRangeTol) {
    return FailDeck("axisymmetric_2d requires full theta and phi ranges");
  }
} else if (cfg.mesh.phi_cells == 0u || (cfg.mesh.phi_cells % 2u) != 0u) {
  return FailDeck("full_3d requires positive even phi_cells");
}
```

Inside the existing `if (perturbation) { ... }` block, after `m` has been parsed and before the current generic `unsupported perturbation option` rejection, add the axisymmetric-specific rejection first:

```cpp
if (IsAxisymmetric2D(cfg.mesh.dimensionality) &&
    cfg.perturbation.enabled &&
    cfg.perturbation.m != 0) {
  return FailDeck("axisymmetric_2d only supports m = 0 perturbations");
}
```

In `BuildDeckReport`, add after geometry/mesh fields are available:

```cpp
<< "; mesh_dimensionality=" << ToString(config.mesh.dimensionality)
<< "; active_hydro_directions="
<< (IsAxisymmetric2D(config.mesh.dimensionality) ? "r,theta" : "r,theta,phi")
<< "; azimuthal_weight="
<< (IsAxisymmetric2D(config.mesh.dimensionality) ? "2pi" : "mesh_phi_faces")
<< "; phi_sweep_configured="
<< (IsAxisymmetric2D(config.mesh.dimensionality) ? "false" : "true")
```

- [ ] **Step 4: Run passing test**

Run the same command from Step 2.

Expected: build succeeds and executable exits `0`.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add src/io/input_deck.hpp src/io/input_deck.cpp tests/contract/test_p5_io_input_deck.cpp
git -C 'F:\dec3d' commit -m 'Add axisymmetric 2D input contract'
```

---

## Task 2: Full-Azimuth Geometry Tests

**Files:**
- Modify: `tests/contract/test_spherical_geometry.cpp`

- [ ] **Step 1: Add failing geometry assertions**

In `tests/contract/test_spherical_geometry.cpp`, add:

```cpp
namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;

double OuterRadialArea(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    std::size_t theta,
    std::size_t phi) {
  const double r = geometry.radial_faces.back();
  const double polar_factor =
      std::cos(geometry.theta_faces[theta]) -
      std::cos(geometry.theta_faces[theta + 1u]);
  const double dphi = geometry.phi_faces[phi + 1u] - geometry.phi_faces[phi];
  return r * r * polar_factor * dphi;
}
}  // namespace
```

Add a test block:

```cpp
{
  const auto geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{4, 8, 1, 0.2, 1.7});
  DEC3D_CHECK(geometry.valid);
  double volume_sum = 0.0;
  for (double value : geometry.cell_volumes) {
    volume_sum += value;
  }
  const double expected_volume =
      (4.0 / 3.0) * kPi * (std::pow(1.7, 3) - std::pow(0.2, 3));
  CheckNear(volume_sum, expected_volume, 1.0e-13, "axisymmetric full volume");

  double outer_area = 0.0;
  for (std::size_t theta = 0; theta < 8u; ++theta) {
    outer_area += OuterRadialArea(geometry, theta, 0u);
  }
  CheckNear(outer_area, 4.0 * kPi * 1.7 * 1.7, 1.0e-13,
            "axisymmetric outer area");
}
```

- [ ] **Step 2: Run geometry test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_spherical_geometry --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_geometry.exe'
```

Expected: pass. A failure means the existing geometry no longer represents `phi_cells=1` as the full `0..2*pi` interval; stop this task and fix only `src/mesh/spherical/spherical_mesh.cpp` before committing the test.

- [ ] **Step 3: Commit**

```powershell
git -C 'F:\dec3d' add tests/contract/test_spherical_geometry.cpp src/mesh/spherical/spherical_mesh.cpp
git -C 'F:\dec3d' commit -m 'Verify full-azimuth geometry for axisymmetric 2D'
```

---

## Task 3: Half-Turn Remap With `phi_cells=1`

**Files:**
- Modify: `src/mesh/boundary/spherical_scalar_remap.cpp`
- Modify: `src/hydro/driver/theta_pole_boundary.cpp`
- Modify: `src/hydro/driver/hydro_boundary_ghosts.cpp`
- Test: `tests/contract/test_spherical_scalar_remap.cpp`
- Test: `tests/contract/test_hydro_theta_pole_boundary.cpp`
- Test: `tests/contract/test_hydro_origin_radial_ghost_remap.cpp`

- [ ] **Step 1: Add failing scalar remap tests**

In `tests/contract/test_spherical_scalar_remap.cpp`, add:

```cpp
{
  const dec3d::mesh::ScalarRemapLayout layout{3, 4, 1};
  const auto validation = dec3d::mesh::ValidateScalarHalfTurnTopology(layout);
  DEC3D_CHECK(validation.success);
  DEC3D_CHECK(validation.phi_half_turn_available);
  DEC3D_CHECK(dec3d::mesh::MapPhiHalfTurn(0u, 1u) == 0u);

  const auto origin = dec3d::mesh::MapScalarOriginNeighbor(1u, 2u, 0u, layout);
  DEC3D_CHECK(origin.success);
  DEC3D_CHECK(origin.mapped_index.radial == 0u);
  DEC3D_CHECK(origin.mapped_index.theta == 1u);
  DEC3D_CHECK(origin.mapped_index.phi == 0u);

  const auto lower = dec3d::mesh::MapScalarLowerPoleNeighbor(1u, 1u, 0u, layout);
  DEC3D_CHECK(lower.success);
  DEC3D_CHECK(lower.mapped_index.radial == 1u);
  DEC3D_CHECK(lower.mapped_index.theta == 0u);
  DEC3D_CHECK(lower.mapped_index.phi == 0u);
}
```

- [ ] **Step 2: Add failing hydro remap tests**

In `tests/contract/test_hydro_theta_pole_boundary.cpp`, add:

```cpp
{
  DEC3D_CHECK(dec3d::hydro::MapPhiAcrossPole(0u, 1u) == 0u);
  dec3d::hydro::HydroConservativeState mapped;
  mapped.rho = 2.0;
  mapped.mom_r = 3.0;
  mapped.mom_theta = 4.0;
  mapped.mom_phi = 0.0;
  mapped.e_fluid_total = 5.0;
  mapped.chi_e = 1.0;
  const auto ghost = dec3d::hydro::BuildThetaPoleGhostState(mapped);
  DEC3D_CHECK(ghost.rho == mapped.rho);
  DEC3D_CHECK(ghost.mom_r == mapped.mom_r);
  DEC3D_CHECK(ghost.mom_theta == -mapped.mom_theta);
  DEC3D_CHECK(ghost.mom_phi == 0.0);
}
```

In `tests/contract/test_hydro_origin_radial_ghost_remap.cpp`, add:

```cpp
{
  DEC3D_CHECK(dec3d::hydro::MapPhiAcrossOrigin(0u, 1u) == 0u);
  DEC3D_CHECK(dec3d::hydro::MapThetaAcrossOrigin(2u, 4u) == 1u);
  dec3d::hydro::HydroConservativeState mapped;
  mapped.rho = 2.0;
  mapped.mom_r = 3.0;
  mapped.mom_theta = 4.0;
  mapped.mom_phi = 0.0;
  mapped.e_fluid_total = 5.0;
  mapped.chi_e = 1.0;
  const auto ghost = dec3d::hydro::BuildOriginRadialGhostState(mapped);
  DEC3D_CHECK(ghost.rho == mapped.rho);
  DEC3D_CHECK(ghost.mom_r == -mapped.mom_r);
  DEC3D_CHECK(ghost.mom_theta == mapped.mom_theta);
  DEC3D_CHECK(ghost.mom_phi == 0.0);
}
```

- [ ] **Step 3: Run failing tests**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_spherical_scalar_remap dec3d_contract_hydro_theta_pole_boundary dec3d_contract_hydro_origin_radial_ghost_remap --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_scalar_remap.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_theta_pole_boundary.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_origin_radial_ghost_remap.exe'
```

Expected before implementation: remap validation fails for `phi_cells=1`.

- [ ] **Step 4: Implement remap change**

In `ValidateScalarHalfTurnTopology`:

```cpp
if (layout.phi_cells != 1u && (layout.phi_cells % 2u) != 0u) {
  return ValidationFailure("scalar remap requires phi cell count one or even");
}
```

In `MapPhiHalfTurn`:

```cpp
if (phi_cells == 1u && phi == 0u) {
  return 0u;
}
if (phi_cells == 0u || (phi_cells % 2u) != 0u || phi >= phi_cells) {
  return phi_cells;
}
```

In `MapPhiAcrossPole`:

```cpp
if (phi_cells == 1u && phi == 0u) {
  return 0u;
}
if (phi_cells == 0u || (phi_cells % 2u) != 0u || phi >= phi_cells) {
  return phi_cells;
}
```

In `FillHydroRadialGhosts` and `FillHydroThetaPoleGhosts` inside `hydro_boundary_ghosts.cpp`, replace even-only checks with one-or-even checks:

```cpp
if (use_origin_remap &&
    interior_states.extent_phi() != 1u &&
    (interior_states.extent_phi() % 2u) != 0u) {
  return FailDirectionalGhostPreparation(
      "radial origin ghost remap requires phi cell count one or even");
}
```

```cpp
if (interior_states.extent_phi() != 1u &&
    (interior_states.extent_phi() % 2u) != 0u) {
  return FailDirectionalGhostPreparation(
      "theta pole ghost fill requires phi cell count one or even");
}
```

- [ ] **Step 5: Run passing tests**

Run the Step 3 commands again.

Expected: all three executables exit `0`.

- [ ] **Step 6: Commit**

```powershell
git -C 'F:\dec3d' add src/mesh/boundary/spherical_scalar_remap.cpp src/hydro/driver/theta_pole_boundary.cpp src/hydro/driver/hydro_boundary_ghosts.cpp tests/contract/test_spherical_scalar_remap.cpp tests/contract/test_hydro_theta_pole_boundary.cpp tests/contract/test_hydro_origin_radial_ghost_remap.cpp
git -C 'F:\dec3d' commit -m 'Allow axisymmetric half-turn remap'
```

---

## Task 4: Diffusion Matrix Axisymmetric Diagnostics

**Files:**
- Modify: `src/transport/diffusion/generic_diffusion.hpp`
- Modify: `src/transport/diffusion/generic_diffusion.cpp`
- Modify: `src/transport/diffusion/hypre_distributed_diffusion_solver.hpp`
- Modify: `src/transport/diffusion/hypre_distributed_diffusion_solver.cpp`
- Test: `tests/contract/test_generic_implicit_diffusion.cpp`
- Test: `tests/contract/test_hypre_distributed_diffusion_solver.cpp`

- [ ] **Step 1: Add serial matrix diagnostics tests**

In `tests/contract/test_generic_implicit_diffusion.cpp`, add helper:

```cpp
bool RowHasDuplicateColumns(const dec3d::transport::SparseMatrixCsr& matrix,
                            std::size_t row) {
  for (std::size_t a = matrix.row_offsets[row]; a < matrix.row_offsets[row + 1u]; ++a) {
    for (std::size_t b = a + 1u; b < matrix.row_offsets[row + 1u]; ++b) {
      if (matrix.column_indices[a] == matrix.column_indices[b]) {
        return true;
      }
    }
  }
  return false;
}
```

Add test block:

```cpp
{
  const auto layout = Layout(4, 4, 1);
  auto problem = BuildPatchProblem(layout, 0.1, 1.0, 2.0, 0.0, 0.0, 4.0);
  const auto assembly = AssembleGenericImplicitDiffusionSystem(problem);
  DEC3D_CHECK(assembly.success);
  DEC3D_CHECK_EQ(assembly.row_count, std::size_t{16});
  DEC3D_CHECK_EQ(assembly.phi_coupling_count, std::size_t{0});
  DEC3D_CHECK(!assembly.phi_neighbor_loop_executed);
  DEC3D_CHECK_EQ(assembly.phi_self_neighbor_attempt_count, std::size_t{0});
  DEC3D_CHECK_EQ(assembly.duplicate_column_row_count, std::size_t{0});
  DEC3D_CHECK_EQ(assembly.interior_row_width5_count, std::size_t{4});
  DEC3D_CHECK(assembly.report_line.find("phi_coupling_count=0") != std::string::npos);
  DEC3D_CHECK(assembly.report_line.find("phi_neighbor_loop_executed=false") !=
              std::string::npos);
  DEC3D_CHECK(assembly.report_line.find("phi_self_neighbor_attempt_count=0") !=
              std::string::npos);
  DEC3D_CHECK(assembly.report_line.find("duplicate_column_row_count=0") !=
              std::string::npos);
  DEC3D_CHECK(assembly.report_line.find("axisymmetric_interior_row_width5_count=4") !=
              std::string::npos);
  const std::size_t interior_row = TestLinearIndex(layout, 1u, 1u, 0u);
  DEC3D_CHECK_EQ(assembly.matrix.row_offsets[interior_row + 1u] -
                 assembly.matrix.row_offsets[interior_row], std::size_t{5});
  for (std::size_t row = 0; row < assembly.row_count; ++row) {
    DEC3D_CHECK(!RowHasDuplicateColumns(assembly.matrix, row));
  }
}
```

- [ ] **Step 2: Add distributed matrix diagnostics tests**

In `tests/contract/test_hypre_distributed_diffusion_solver.cpp`, add equivalent helper:

```cpp
bool LocalRowHasDuplicateColumns(
    const dec3d::transport::DistributedLocalCsrMatrix& matrix,
    std::size_t local_row) {
  for (std::size_t a = matrix.row_offsets[local_row];
       a < matrix.row_offsets[local_row + 1u];
       ++a) {
    for (std::size_t b = a + 1u; b < matrix.row_offsets[local_row + 1u]; ++b) {
      if (matrix.column_indices[a] == matrix.column_indices[b]) {
        return true;
      }
    }
  }
  return false;
}
```

Add after the row ownership test:

```cpp
{
  const auto layout = Layout(4, 4, 1);
  const auto geometry = BuildGeometry(layout, 1.0, 3.0, 0.5, 1.2);
  const auto problem =
      BuildDistributedProblem(MPI_COMM_WORLD, layout, geometry, 0.1, 4.0, PatchBoundary());
  const auto assembly = AssembleDistributedGenericDiffusionSystem(problem);
  DEC3D_CHECK(assembly.success);
  DEC3D_CHECK_EQ(assembly.ownership.global_row_count, std::size_t{16});
  DEC3D_CHECK_EQ(assembly.global_phi_coupling_count, std::size_t{0});
  DEC3D_CHECK_EQ(assembly.global_phi_self_neighbor_attempt_count, std::size_t{0});
  DEC3D_CHECK_EQ(assembly.global_duplicate_column_row_count, std::size_t{0});
  DEC3D_CHECK_EQ(assembly.global_axisymmetric_interior_row_width5_count, std::size_t{4});
  for (std::size_t local_row = 0; local_row < assembly.ownership.local_row_count; ++local_row) {
    DEC3D_CHECK(!LocalRowHasDuplicateColumns(assembly.local_matrix, local_row));
    const std::size_t global_row = assembly.ownership.local_row_begin + local_row;
    const std::size_t global_radial =
        global_row / (layout.theta_cells * layout.phi_cells);
    const std::size_t theta = (global_row / layout.phi_cells) % layout.theta_cells;
    if (global_radial == 1u && theta == 1u) {
      DEC3D_CHECK_EQ(assembly.local_matrix.row_offsets[local_row + 1u] -
                     assembly.local_matrix.row_offsets[local_row], std::size_t{5});
    }
  }
  DEC3D_CHECK(assembly.report_line.find("global_phi_coupling_count=0") !=
              std::string::npos);
  DEC3D_CHECK(assembly.report_line.find("global_phi_self_neighbor_attempt_count=0") !=
              std::string::npos);
  DEC3D_CHECK(assembly.report_line.find(
                  "global_axisymmetric_interior_row_width5_count=4") !=
              std::string::npos);
}
```

- [ ] **Step 3: Run failing tests**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_generic_implicit_diffusion dec3d_contract_hypre_distributed_diffusion_solver --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_generic_implicit_diffusion.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_hypre_distributed_diffusion_solver.exe'
```

Expected before implementation: compile failure because diagnostic fields do not exist.

- [ ] **Step 4: Implement serial diagnostics**

Add fields to `GenericDiffusionAssemblyResult`:

```cpp
std::size_t phi_coupling_count{0};
bool phi_neighbor_loop_executed{false};
std::size_t phi_self_neighbor_attempt_count{0};
std::size_t duplicate_column_row_count{0};
std::size_t interior_row_width5_count{0};
```

In `AddConductance`, when called for phi conductance, count through a new wrapper rather than overloading the radial/theta path. Add:

```cpp
void AddPhiConductance(
    RowAccumulator& rows,
    std::size_t left,
    std::size_t right,
    double conductance,
    GenericDiffusionAssemblyResult& result) {
  result.phi_neighbor_loop_executed = true;
  if (left == right) {
    ++result.phi_self_neighbor_attempt_count;
    return;
  }
  if (conductance == 0.0) {
    return;
  }
  ++result.phi_coupling_count;
  AddConductance(rows, left, right, conductance);
}
```

Use `AddPhiConductance` only in the phi loop. Keep the existing `if (layout.phi_cells <= 1u) return;`.

After `BuildCsrMatrix`, compute duplicate rows and axisymmetric interior row width. The duplicate-column check is intentionally a CSR safety net; `phi_self_neighbor_attempt_count` is the direct guard against attempted phi periodic self-neighbors before CSR coalescing can hide them:

```cpp
std::size_t CountDuplicateColumnRows(const SparseMatrixCsr& matrix) noexcept {
  std::size_t count = 0u;
  for (std::size_t row = 0; row < matrix.row_count; ++row) {
    bool duplicate = false;
    for (std::size_t a = matrix.row_offsets[row]; a < matrix.row_offsets[row + 1u]; ++a) {
      for (std::size_t b = a + 1u; b < matrix.row_offsets[row + 1u]; ++b) {
        duplicate = duplicate || matrix.column_indices[a] == matrix.column_indices[b];
      }
    }
    if (duplicate) {
      ++count;
    }
  }
  return count;
}

std::size_t CountAxisymmetricInteriorRowWidth5(
    const DiffusionPatchLayout& layout,
    const SparseMatrixCsr& matrix) noexcept {
  if (layout.phi_cells != 1u || layout.radial_cells < 3u || layout.theta_cells < 3u) {
    return 0u;
  }
  std::size_t count = 0u;
  for (std::size_t r = 1u; r + 1u < layout.radial_cells; ++r) {
    for (std::size_t t = 1u; t + 1u < layout.theta_cells; ++t) {
      const std::size_t row = LinearIndex(layout, r, t, 0u);
      const std::size_t width = matrix.row_offsets[row + 1u] - matrix.row_offsets[row];
      if (width == 5u) {
        ++count;
      }
    }
  }
  return count;
}
```

Add report tokens:

```cpp
<< "; phi_coupling_count=" << result.phi_coupling_count
<< "; phi_neighbor_loop_executed=" << (result.phi_neighbor_loop_executed ? "true" : "false")
<< "; phi_self_neighbor_attempt_count=" << result.phi_self_neighbor_attempt_count
<< "; duplicate_column_row_count=" << result.duplicate_column_row_count
<< "; axisymmetric_interior_row_width5_count=" << result.interior_row_width5_count
```

Update `ValidateGenericDiffusionAssemblyDiagnostics` to require all five tokens. Define `phi_coupling_count` as the number of nonzero off-diagonal phi face-pair conductances attempted before CSR assembly, not as row count.

- [ ] **Step 5: Implement distributed diagnostics**

Add fields to `DistributedGenericDiffusionAssemblyResult`:

```cpp
std::size_t local_phi_coupling_count{0};
std::size_t global_phi_coupling_count{0};
bool local_phi_neighbor_loop_executed{false};
bool global_phi_neighbor_loop_executed{false};
std::size_t local_phi_self_neighbor_attempt_count{0};
std::size_t global_phi_self_neighbor_attempt_count{0};
std::size_t local_duplicate_column_row_count{0};
std::size_t global_duplicate_column_row_count{0};
std::size_t local_axisymmetric_interior_row_width5_count{0};
std::size_t global_axisymmetric_interior_row_width5_count{0};
```

In distributed phi conductance loop, skip when `ownership.global_phi_cells <= 1u`:

```cpp
if (ownership.global_phi_cells > 1u) {
  result.local_phi_neighbor_loop_executed = true;
  // existing phi conductance loop
}
```

When adding phi conductance with `AddOwnedPairConductance`, increment local self-neighbor attempts before calling the row accumulator, and increment local phi coupling count only for nonzero non-self face pairs:

```cpp
const auto lhs = GlobalRow(ownership, gr, t, p);
const auto rhs = GlobalRow(ownership, gr, t, next_phi);
if (lhs == rhs) {
  ++result.local_phi_self_neighbor_attempt_count;
} else if (conductance != 0.0) {
  ++result.local_phi_coupling_count;
}
```

After `BuildLocalCsr`, compute local duplicate rows with the same nested check as the test. Compute local axisymmetric interior width-5 rows only for owned global rows with `1 <= r < Nr-1`, `1 <= theta < Ntheta-1`, and `global_phi_cells == 1`. Reduce count fields with `MPI_Allreduce(..., MPI_SUM, ...)` and reduce `local_phi_neighbor_loop_executed` with `MPI_LOR`. Add report tokens:

```cpp
<< "; local_phi_coupling_count=" << result.local_phi_coupling_count
<< "; global_phi_coupling_count=" << result.global_phi_coupling_count
<< "; global_phi_neighbor_loop_executed="
<< (result.global_phi_neighbor_loop_executed ? "true" : "false")
<< "; local_phi_self_neighbor_attempt_count=" << result.local_phi_self_neighbor_attempt_count
<< "; global_phi_self_neighbor_attempt_count=" << result.global_phi_self_neighbor_attempt_count
<< "; local_duplicate_column_row_count=" << result.local_duplicate_column_row_count
<< "; global_duplicate_column_row_count=" << result.global_duplicate_column_row_count
<< "; global_axisymmetric_interior_row_width5_count="
<< result.global_axisymmetric_interior_row_width5_count
```

Update `ValidateDistributedGenericDiffusionAssemblyDiagnostics` to require all new tokens. Define `global_phi_coupling_count` as the sum of nonzero off-diagonal phi face-pair conductances across ranks before HYPRE insertion.

- [ ] **Step 6: Run passing tests**

Run Step 3 commands again.

Expected: both tests exit `0`.

- [ ] **Step 7: Commit**

```powershell
git -C 'F:\dec3d' add src/transport/diffusion/generic_diffusion.hpp src/transport/diffusion/generic_diffusion.cpp src/transport/diffusion/hypre_distributed_diffusion_solver.hpp src/transport/diffusion/hypre_distributed_diffusion_solver.cpp tests/contract/test_generic_implicit_diffusion.cpp tests/contract/test_hypre_distributed_diffusion_solver.cpp
git -C 'F:\dec3d' commit -m 'Add axisymmetric diffusion matrix diagnostics'
```

---

## Task 5: Initialization Axisymmetric Invariants

**Files:**
- Modify: `src/initialization/profile_initializer.cpp`
- Test: `tests/contract/test_p5_io_profile_initializer.cpp`

- [ ] **Step 1: Add failing initializer tests**

In `tests/contract/test_p5_io_profile_initializer.cpp`, after the existing `p2_init.state.mom_r(...)` assertion and before `auto old_coordinate_deck = p2_velocity_deck;`, add diagnostics and nonzero-`vp` rejection checks using the already-loaded `p2_deck`, `p2_profile`, and `profile_path`:

```cpp
{
  auto config = p2_deck.config;
  config.mesh.dimensionality = dec3d::io::MeshDimensionality::axisymmetric_2d;
  config.mesh.phi_cells = 1u;
  config.mesh.moving_mesh = false;
  const auto result = dec3d::initialization::InitializeFromRadialProfile(
      config, p2_profile.profile, profile_path);
  DEC3D_CHECK(result.success);
  DEC3D_CHECK(result.report_line.find("initial_perturbation_normalization=raw_Pl") !=
              std::string::npos);
  DEC3D_CHECK(result.report_line.find("axisymmetric_mom_phi_zero=true") !=
              std::string::npos);
  DEC3D_CHECK(result.state.layout.phi_cells == 1u);
  for (std::size_t r = 0; r < result.state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < result.state.layout.theta_cells; ++t) {
      DEC3D_CHECK(result.state.mom_phi(r, t, 0u) == 0.0);
    }
  }
}

{
  auto config = p2_deck.config;
  config.mesh.dimensionality = dec3d::io::MeshDimensionality::axisymmetric_2d;
  config.mesh.phi_cells = 1u;
  config.mesh.moving_mesh = false;
  auto profile_with_vp = p2_profile.profile;
  for (auto& row : profile_with_vp.rows) {
    row.vp_cm_s = 1.0e5;
  }
  const auto result = dec3d::initialization::InitializeFromRadialProfile(
      config, profile_with_vp, profile_path);
  DEC3D_CHECK(!result.success);
  DEC3D_CHECK(result.failure_reason.find("axisymmetric_2d requires zero vp_cm_s") !=
              std::string::npos);
}
```

- [ ] **Step 2: Run failing test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_io_profile_initializer --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_profile_initializer.exe'
```

Expected before implementation: compile failure or missing diagnostics.

- [ ] **Step 3: Implement initializer checks**

In `BuildReport`, inside `config.perturbation.enabled` block add:

```cpp
<< "; initial_perturbation_normalization=raw_Pl"
```

Outside the perturbation block add:

```cpp
<< "; axisymmetric_mom_phi_zero="
<< (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) ? "true" : "not_applicable")
```

Before assigning state in the loop:

```cpp
if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) &&
    std::abs(row.vp_cm_s) > 1.0e-30) {
  return Fail("axisymmetric_2d requires zero vp_cm_s");
}
```

Keep `vt_cm_s` legal because axisymmetric flow may have theta velocity.

- [ ] **Step 4: Run passing test**

Run Step 2 commands again.

Expected: executable exits `0`.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add src/initialization/profile_initializer.cpp tests/contract/test_p5_io_profile_initializer.cpp
git -C 'F:\dec3d' commit -m 'Enforce axisymmetric initialization invariants'
```

---

## Task 6A: Hydro Configuration Only

**Files:**
- Modify: `src/app/dec3d_app.cpp`
- Modify: `src/hydro/driver/hydro_operator.cpp`
- Modify: `src/hydro/driver/static_grid_hydro.cpp`
- Test: `tests/contract/test_p5_runtime_all_stages.cpp`

- [ ] **Step 1: Add failing runtime report expectations**

In `tests/contract/test_p5_runtime_all_stages.cpp`, after the existing Noh H-only block and before `std::filesystem::remove_all(root);`, add this H-only axisymmetric app-command smoke:

```cpp
  const auto axisym_root = root / "axisym_h";
  const auto axisym_deck = axisym_root / "case_axisym_h.in";
  const auto axisym_profile = axisym_root / "case_axisym_h.pro";
  const auto axisym_output = axisym_root / "output";
  WriteText(axisym_deck, R"ini(
[run]
case_name = p5_axisymmetric_h_only
phase = P5
stage_order = H
step_count = 1
dt_mode = hydro_relaxed
cfl = 0.4
output_dir = AXISYM_OUTPUT_ROOT_REPLACED
[mesh]
geometry = spherical
dimensionality = axisymmetric_2d
radial_cells = 8
theta_cells = 4
phi_cells = 1
radial_min_cm = 0.0
radial_max_cm = 1.0e-2
theta_min = 0.0
theta_max = 3.141592653589793
phi_min = 0.0
phi_max = 6.283185307179586
moving_mesh = false
macro_zoning = true
[initial_condition]
profile_interpolation = linear
profile_radius_unit = um
outside_profile_policy = hard_fail
[physics]
enable_hydro = true
enable_thermal = false
enable_equilibration = false
enable_radiation = false
enable_alpha = false
[radiation]
group_mode = explicit_frequency_groups
group_edges_eV = 1,10
opacity_provider = tops_dt_tabulated
radiation_initialization = zero
radiation_initial_blackbody_scale = 1.0
boundary_model = thesis_marshak_vacuum
flux_limiter = harmonic_eq_5_209
[thermal]
kappa_model = lee_more_with_degeneracy
electron_flux_limiter = minmax_old_time_face_effective_kappa
[alpha]
composition_model = equimolar_dt_from_p2_recovery
reactivity_model = bosch_hale_dt
alpha_initialization = zero
[boundaries]
inner_radial = scalar_origin_remap_required
outer_radial = neumann_zero_flux
theta = scalar_pole_remap_required
phi = periodic
[output]
write_profiles_every = 0
write_diagnostics_every = 0
write_restart_every = 0
write_dec3d_out_every_steps = 1
field_checkpoint_every_steps = 0
field_checkpoint_interval_s = 0.0
field_checkpoint_prefix = fields
field_checkpoint_format = csv3d
restart_checkpoint_every_steps = 0
restart_checkpoint_interval_s = 0.0
restart_checkpoint_prefix = restart
restart_checkpoint_format = dec3d_restart_text
)ini");

  auto axisym_deck_text = ReadText(axisym_deck);
  const auto axisym_marker = std::string("AXISYM_OUTPUT_ROOT_REPLACED");
  axisym_deck_text.replace(axisym_deck_text.find(axisym_marker),
                           axisym_marker.size(),
                           axisym_output.generic_string());
  WriteText(axisym_deck, axisym_deck_text);

  WriteText(axisym_profile, R"pro(
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s epsilon_alpha_erg_cm3 radiation_scale
0.0 10.0 1.0 1.0 -1.0e6 0.0 0.0 0.0 1.0
100.0 10.0 1.0 1.0 -1.0e6 0.0 0.0 0.0 1.0
)pro");

  const std::vector<std::string> axisym_argv{
      "dec3d.exe",
      "-input",
      axisym_deck.string(),
      "-profile",
      axisym_profile.string()};
  const auto axisym_result = dec3d::app::RunDec3DCommandLine(axisym_argv);
  DEC3D_CHECK(axisym_result.exit_code == 0);
  DEC3D_CHECK(axisym_result.report_line.find("mesh_dimensionality=axisymmetric_2d") !=
            std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("active_hydro_directions=r,theta") !=
            std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("phi_sweep_executed=false") !=
              std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("h_hydro_phi_sweep_wall_s=0") !=
              std::string::npos);
```

- [ ] **Step 2: Run failing test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_runtime_all_stages --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_all_stages.exe'
```

Expected before implementation: missing report tokens.

- [ ] **Step 3: Implement runtime hydro options**

In `RuntimeHydroOptions(config)` inside `src/app/dec3d_app.cpp`, after current option construction:

```cpp
if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality)) {
  options.apply_radial_sweep = true;
  options.apply_theta_sweep = true;
  options.apply_phi_sweep = false;
  options.macro_zoning_use_phi_ppm = false;
  options.apply_radial_ale_flux_correction = false;
  options.request_radial_ale_proposal = false;
  options.use_macro_ale_direct_moving_face_hllc = false;
}
```

Ensure the static-grid branch is selected by existing `config.mesh.moving_mesh=false`.

- [ ] **Step 4: Run hydro-configuration passing test**

Run Step 2 command again.

Expected: executable exits `0` and proves `apply_phi_sweep=false` reaches runtime.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add src/app/dec3d_app.cpp src/hydro/driver/hydro_operator.cpp src/hydro/driver/static_grid_hydro.cpp tests/contract/test_p5_runtime_all_stages.cpp
git -C 'F:\dec3d' commit -m 'Configure no-ALE axisymmetric hydro directions'
```

---

## Task 6B: Runtime Invariant Diagnostics

**Files:**
- Modify: `src/app/dec3d_app.cpp`
- Modify: `tests/contract/test_p5_runtime_all_stages.cpp`

- [ ] **Step 1: Extend serial smoke with invariant assertions**

In the axisymmetric H-only block added in Task 6A, add:

```cpp
  DEC3D_CHECK(axisym_result.report_line.find("max_abs_mom_phi=") != std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("max_abs_v_phi=") != std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("axisymmetric_mom_phi_tol=") !=
              std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("axisymmetric_v_phi_tol=") !=
              std::string::npos);
  DEC3D_CHECK(axisym_result.report_line.find("axisymmetric_invariant_ok=true") !=
              std::string::npos);
```

- [ ] **Step 2: Run failing test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_runtime_all_stages --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_all_stages.exe'
```

Expected before implementation: invariant report tokens are missing.

- [ ] **Step 3: Add invariant reducer**

In `src/app/dec3d_app.cpp`, add helper:

```cpp
struct AxisymmetricInvariantSummary {
  double max_abs_mom_phi{0.0};
  double max_abs_v_phi{0.0};
  double max_abs_mom_total{0.0};
  double max_abs_velocity{0.0};
  double mom_phi_tol{1.0e-30};
  double v_phi_tol{1.0e-30};
  bool ok{true};
};

AxisymmetricInvariantSummary EvaluateAxisymmetricInvariants(
    const dec3d::state::CanonicalState& state) noexcept {
  AxisymmetricInvariantSummary summary;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double rho = state.rho(r, t, p);
        const double mom_r = state.mom_r(r, t, p);
        const double mom_theta = state.mom_theta(r, t, p);
        const double mom_phi = state.mom_phi(r, t, p);
        summary.max_abs_mom_phi = std::max(summary.max_abs_mom_phi, std::abs(mom_phi));
        summary.max_abs_mom_total = std::max(
            summary.max_abs_mom_total,
            std::sqrt(mom_r * mom_r + mom_theta * mom_theta + mom_phi * mom_phi));
        if (rho > 0.0) {
          const double vr = mom_r / rho;
          const double vt = mom_theta / rho;
          const double vp = mom_phi / rho;
          summary.max_abs_v_phi = std::max(summary.max_abs_v_phi, std::abs(vp));
          summary.max_abs_velocity =
              std::max(summary.max_abs_velocity, std::sqrt(vr * vr + vt * vt + vp * vp));
        }
      }
    }
  }
  summary.mom_phi_tol = std::max(1.0e-30, 1.0e-14 * summary.max_abs_mom_total);
  summary.v_phi_tol = std::max(1.0e-30, 1.0e-14 * summary.max_abs_velocity);
  summary.ok =
      summary.max_abs_mom_phi <= summary.mom_phi_tol &&
      summary.max_abs_v_phi <= summary.v_phi_tol;
  return summary;
}
```

For MPI local states, compute local summary and use `MPI_Allreduce(..., MPI_MAX, ...)` for all max values before computing tolerances.

- [ ] **Step 4: Add runtime report tokens**

In serial and distributed loop reports, append:

```cpp
<< "; mesh_dimensionality=" << dec3d::io::ToString(config.mesh.dimensionality)
<< "; active_hydro_directions="
<< (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) ? "r,theta" : "r,theta,phi")
<< "; phi_sweep_executed="
<< (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) ? "false" :
    (loop.h_hydro_phi_sweep_wall_s > 0.0 ? "true" : "false"))
<< "; max_abs_mom_phi=" << axisym.max_abs_mom_phi
<< "; max_abs_v_phi=" << axisym.max_abs_v_phi
<< "; axisymmetric_mom_phi_tol=" << axisym.mom_phi_tol
<< "; axisymmetric_v_phi_tol=" << axisym.v_phi_tol
<< "; axisymmetric_invariant_ok="
<< (!dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) || axisym.ok ? "true" : "false")
```

If axisymmetric invariants fail, return a runtime failure instead of only reporting:

```cpp
if (dec3d::io::IsAxisymmetric2D(config.mesh.dimensionality) && !axisym.ok) {
  loop.failure_reason = "axisymmetric invariant violated";
  loop.failure_diagnostics =
      "diagnostic_id=p5.runtime.loop.failure; failure_reason=" + loop.failure_reason;
  return loop;
}
```

- [ ] **Step 5: Run passing test**

Run Step 2 command again.

Expected: executable exits `0`.

- [ ] **Step 6: Commit**

```powershell
git -C 'F:\dec3d' add src/app/dec3d_app.cpp tests/contract/test_p5_runtime_all_stages.cpp
git -C 'F:\dec3d' commit -m 'Add axisymmetric runtime invariant diagnostics'
```

---

## Task 7: Distributed H,T,E,R,A Axisymmetric Smoke

**Files:**
- Create: `cases/axisymmetric_2d_noale_allstages_smoke.in`
- Create: `cases/axisymmetric_2d_noale_smoke.pro`
- Modify: `tests/contract/test_p5_runtime_distributed_all_stages.cpp`

- [ ] **Step 1: Create tiny profile**

Create `cases/axisymmetric_2d_noale_smoke.pro`:

```text
r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s radiation_scale epsilon_alpha_erg_cm3
0 1.0 1.0 1.0 -1.0e6 0 0 1.0 0
10 1.0 1.0 1.0 -1.0e6 0 0 1.0 0
20 1.0 1.0 1.0 -1.0e6 0 0 1.0 0
30 1.0 1.0 1.0 -1.0e6 0 0 1.0 0
```

- [ ] **Step 2: Create tiny input deck**

Create `cases/axisymmetric_2d_noale_allstages_smoke.in`:

```ini
[run]
case_name = axisymmetric_2d_noale_allstages_smoke
phase = P5
stage_order = H,T,E,R,A
step_count = 2
target_time_s = 0
dt_mode = hydro_relaxed
cfl = 0.1
output_dir = F:/dec3d/analysis/output/axisymmetric_2d_noale_allstages_smoke

[mesh]
geometry = spherical
dimensionality = axisymmetric_2d
radial_cells = 4
theta_cells = 4
phi_cells = 1
radial_min_cm = 0.0
radial_max_cm = 3.0e-3
theta_min = 0.0
theta_max = 3.141592653589793
phi_min = 0.0
phi_max = 6.283185307179586
moving_mesh = false
macro_zoning = true

[initial_condition]
profile_interpolation = linear
profile_radius_unit = um
outside_profile_policy = hard_fail

[perturbation]
enabled = true
type = single_mode_radial_velocity
ell = 2
m = 0
amplitude = 0.01
r0_cm = 2.0e-3
target = radial_velocity_cm_s

[physics]
enable_hydro = true
enable_thermal = true
enable_equilibration = true
enable_radiation = true
enable_alpha = true

[radiation]
group_mode = explicit_frequency_groups
group_edges_eV = 1,10,100
opacity_provider = tops_dt_tabulated
radiation_initialization = zero
radiation_initial_blackbody_scale = 1.0
boundary_model = thesis_marshak_vacuum
flux_limiter = harmonic_eq_5_209

[thermal]
kappa_model = lee_more_with_degeneracy
electron_flux_limiter = minmax_old_time_face_effective_kappa

[alpha]
composition_model = equimolar_dt_from_p2_recovery
reactivity_model = bosch_hale_dt
alpha_initialization = zero

[boundaries]
inner_radial = scalar_origin_remap_required
outer_radial = neumann_zero_flux
theta = scalar_pole_remap_required
phi = periodic

[output]
write_profiles_every = 0
write_diagnostics_every = 0
write_restart_every = 0
write_dec3d_out_every_steps = 1
field_checkpoint_every_steps = 0
field_checkpoint_interval_s = 0
field_checkpoint_prefix = fields
field_checkpoint_format = csv3d
restart_checkpoint_every_steps = 0
restart_checkpoint_interval_s = 0
restart_checkpoint_prefix = restart
restart_checkpoint_format = dec3d_restart_text
history_profile_interval_s = 0
history_profile_file = axisymmetric_2d_noale_allstages_smoke.his
```

- [ ] **Step 3: Add distributed smoke test**

In `tests/contract/test_p5_runtime_distributed_all_stages.cpp`, after the existing full-3D distributed assertions and before the final cleanup barrier, add a second app-command run. Use the same `rank == 0` write, `MPI_Barrier`, `std::vector<std::string> args`, and `dec3d::app::RunDec3DCommandLine(args)` pattern already used in this file:

```cpp
    const auto axisym_root = root / "axisymmetric_2d";
    const auto axisym_deck = axisym_root / "axisymmetric_2d_noale_allstages_smoke.in";
    const auto axisym_profile = axisym_root / "axisymmetric_2d_noale_smoke.pro";
    const auto axisym_output = axisym_root / "output";
    if (rank == 0) {
      std::filesystem::create_directories(axisym_root);
      std::filesystem::copy_file("F:/dec3d/cases/axisymmetric_2d_noale_allstages_smoke.in",
                                 axisym_deck,
                                 std::filesystem::copy_options::overwrite_existing);
      std::filesystem::copy_file("F:/dec3d/cases/axisymmetric_2d_noale_smoke.pro",
                                 axisym_profile,
                                 std::filesystem::copy_options::overwrite_existing);
      auto axisym_deck_text = ReadText(axisym_deck);
      const auto old_output =
          std::string("F:/dec3d/analysis/output/axisymmetric_2d_noale_allstages_smoke");
      axisym_deck_text.replace(axisym_deck_text.find(old_output),
                               old_output.size(),
                               axisym_output.generic_string());
      WriteText(axisym_deck, axisym_deck_text);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    const std::vector<std::string> axisym_args{
        "dec3d.exe",
        "-input",
        axisym_deck.string(),
        "-profile",
        axisym_profile.string()};
    const auto axisym_result = dec3d::app::RunDec3DCommandLine(axisym_args);
    if (axisym_result.exit_code != 0) {
      dec3d::test::Fail("axisymmetric distributed runtime success",
                        __FILE__,
                        __LINE__,
                        axisym_result.failure_diagnostics);
    }
```

Add these assertions after the axisymmetric run:

```cpp
DEC3D_CHECK(axisym_result.report_line.find("runtime_steps_executed=2") !=
            std::string::npos);
DEC3D_CHECK(axisym_result.report_line.find("mesh_dimensionality=axisymmetric_2d") !=
            std::string::npos);
DEC3D_CHECK(axisym_result.report_line.find("phi_sweep_executed=false") !=
            std::string::npos);
DEC3D_CHECK(axisym_result.report_line.find("axisymmetric_invariant_ok=true") !=
            std::string::npos);
DEC3D_CHECK(axisym_result.report_line.find("global_phi_coupling_count=0") !=
            std::string::npos);
DEC3D_CHECK(axisym_result.report_line.find("global_duplicate_column_row_count=0") !=
            std::string::npos);
```

- [ ] **Step 4: Run smoke**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d dec3d_contract_p5_runtime_distributed_all_stages --clean-first -- /m
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_distributed_all_stages.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d.exe' -input 'F:\dec3d\cases\axisymmetric_2d_noale_allstages_smoke.in' -profile 'F:\dec3d\cases\axisymmetric_2d_noale_smoke.pro'
```

Expected: both commands exit `0`, `stderr` empty if redirected.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add cases/axisymmetric_2d_noale_allstages_smoke.in cases/axisymmetric_2d_noale_smoke.pro tests/contract/test_p5_runtime_distributed_all_stages.cpp
git -C 'F:\dec3d' commit -m 'Add axisymmetric all-stage smoke case'
```

---

## Task 8: Production Vector Diagnostic And 2D-vs-3D Regression

**Files:**
- Create: `src/state/diagnostics/axisymmetric_vector_diagnostics.hpp`
- Create: `tests/contract/test_axisymmetric_2d_regression.cpp`
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Add regression target to CMake**

In `CMakeLists.txt`, add near other contract tests:

```cmake
dec3d_add_test(dec3d_contract_axisymmetric_2d_regression
  tests/contract/test_axisymmetric_2d_regression.cpp
)
target_link_libraries(dec3d_contract_axisymmetric_2d_regression PRIVATE
  dec3d_initialization
  dec3d_io
  dec3d_hydro
  dec3d_transport
  dec3d_radiation
  dec3d_state
  dec3d_mesh
)
```

- [ ] **Step 2: Create production vector diagnostic helper**

Create `src/state/diagnostics/axisymmetric_vector_diagnostics.hpp`:

```cpp
#pragma once

#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"

#include <cmath>
#include <sstream>
#include <string>

namespace dec3d::state::diagnostics {
constexpr double kPi = 3.141592653589793238462643383279502884;

struct IntegratedCartesianMomentum {
  double px{0.0};
  double py{0.0};
  double pz{0.0};
  std::string report_line;
};

[[nodiscard]] inline double CellVolume(
    const dec3d::mesh::SphericalGeometryMetadata& geometry,
    const dec3d::state::CanonicalStateLayout& layout,
    std::size_t r,
    std::size_t t,
    std::size_t p) noexcept {
  return geometry.cell_volumes[((r * layout.theta_cells) + t) * layout.phi_cells + p];
}

[[nodiscard]] inline IntegratedCartesianMomentum
ComputeAxisymmetricCartesianMomentum(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  IntegratedCartesianMomentum out;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      const double theta = 0.5 * (geometry.theta_faces[t] + geometry.theta_faces[t + 1u]);
      const double volume = CellVolume(geometry, state.layout, r, t, 0u);
      out.pz +=
          (state.mom_r(r, t, 0u) * std::cos(theta) -
           state.mom_theta(r, t, 0u) * std::sin(theta)) * volume;
    }
  }
  std::ostringstream report;
  report << "diagnostic_id=axisymmetric.cartesian_momentum"
         << "; cartesian_momentum_phi_average=analytic"
         << "; px=" << out.px << "; py=" << out.py << "; pz=" << out.pz;
  out.report_line = report.str();
  return out;
}

[[nodiscard]] inline IntegratedCartesianMomentum ComputeFull3DCartesianMomentum(
    const dec3d::state::CanonicalState& state,
    const dec3d::mesh::SphericalGeometryMetadata& geometry) {
  IntegratedCartesianMomentum out;
  for (std::size_t r = 0; r < state.layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < state.layout.theta_cells; ++t) {
      const double theta = 0.5 * (geometry.theta_faces[t] + geometry.theta_faces[t + 1u]);
      const double st = std::sin(theta);
      const double ct = std::cos(theta);
      for (std::size_t p = 0; p < state.layout.phi_cells; ++p) {
        const double phi = 0.5 * (geometry.phi_faces[p] + geometry.phi_faces[p + 1u]);
        const double cp = std::cos(phi);
        const double sp = std::sin(phi);
        const double volume = CellVolume(geometry, state.layout, r, t, p);
        const double mr = state.mom_r(r, t, p);
        const double mt = state.mom_theta(r, t, p);
        const double mp = state.mom_phi(r, t, p);
        out.px += (mr * st * cp + mt * ct * cp - mp * sp) * volume;
        out.py += (mr * st * sp + mt * ct * sp + mp * cp) * volume;
        out.pz += (mr * ct - mt * st) * volume;
      }
    }
  }
  std::ostringstream report;
  report << "diagnostic_id=full3d.cartesian_momentum"
         << "; cartesian_momentum_phi_average=numeric"
         << "; px=" << out.px << "; py=" << out.py << "; pz=" << out.pz;
  out.report_line = report.str();
  return out;
}
}  // namespace dec3d::state::diagnostics
```

- [ ] **Step 3: Create vector diagnostic regression test**

Create `tests/contract/test_axisymmetric_2d_regression.cpp` with:

```cpp
#include "core/array3d.hpp"
#include "mesh/spherical/spherical_mesh.hpp"
#include "state/canonical_state/canonical_state.hpp"
#include "state/diagnostics/axisymmetric_vector_diagnostics.hpp"
#include "test_assert.hpp"

#include <cmath>
#include <sstream>

namespace {
void CheckNear(double lhs, double rhs, double tol, const char* label) {
  if (std::abs(lhs - rhs) > tol) {
    std::ostringstream detail;
    detail << label << " lhs=" << lhs << " rhs=" << rhs;
    dec3d::test::Fail("near equality", __FILE__, __LINE__, detail.str());
  }
}

double AveragePhi(const dec3d::core::Array3D<double>& field,
                  std::size_t radial,
                  std::size_t theta) {
  double sum = 0.0;
  for (std::size_t p = 0; p < field.extent_phi(); ++p) {
    sum += field(radial, theta, p);
  }
  return sum / static_cast<double>(field.extent_phi());
}
}  // namespace

int main() {
  {
    dec3d::state::CanonicalStateLayout layout;
    layout.radial_cells = 2u;
    layout.theta_cells = 4u;
    layout.phi_cells = 1u;
    layout.radiation_group_count = 1u;
    auto state = dec3d::state::CanonicalState::Create(layout);
    const auto geometry = dec3d::mesh::BuildSphericalGeometry(
        dec3d::mesh::SphericalMeshDescriptor{2u, 4u, 1u, 1.0, 2.0});
    DEC3D_CHECK(geometry.valid);
    for (std::size_t r = 0; r < layout.radial_cells; ++r) {
      for (std::size_t t = 0; t < layout.theta_cells; ++t) {
        state.rho(r, t, 0u) = 2.0;
        state.mom_r(r, t, 0u) = 6.0;
        state.mom_theta(r, t, 0u) = 2.0;
        state.mom_phi(r, t, 0u) = 0.0;
      }
    }
    const auto momentum =
        dec3d::state::diagnostics::ComputeAxisymmetricCartesianMomentum(state, geometry);
    CheckNear(momentum.px, 0.0, 0.0, "axisymmetric px");
    CheckNear(momentum.py, 0.0, 0.0, "axisymmetric py");
    DEC3D_CHECK(std::isfinite(momentum.pz));
    DEC3D_CHECK(momentum.report_line.find("cartesian_momentum_phi_average=analytic") !=
                std::string::npos);
  }

  dec3d::state::CanonicalStateLayout axisym_layout;
  axisym_layout.radial_cells = 3u;
  axisym_layout.theta_cells = 4u;
  axisym_layout.phi_cells = 1u;
  axisym_layout.radiation_group_count = 1u;
  dec3d::state::CanonicalStateLayout full_layout = axisym_layout;
  full_layout.phi_cells = 8u;
  auto axisym = dec3d::state::CanonicalState::Create(axisym_layout);
  auto full3d = dec3d::state::CanonicalState::Create(full_layout);
  for (std::size_t r = 0; r < axisym_layout.radial_cells; ++r) {
    for (std::size_t t = 0; t < axisym_layout.theta_cells; ++t) {
      const double rho = 1.0 + 0.1 * static_cast<double>(r) +
                         0.01 * static_cast<double>(t);
      const double mom_r = -2.0 + 0.2 * static_cast<double>(r);
      const double mom_theta = 0.05 * static_cast<double>(t);
      const double e_electron = 10.0 + static_cast<double>(r + t);
      const double e_total = e_electron + 4.0;
      axisym.rho(r, t, 0u) = rho;
      axisym.mom_r(r, t, 0u) = mom_r;
      axisym.mom_theta(r, t, 0u) = mom_theta;
      axisym.mom_phi(r, t, 0u) = 0.0;
      axisym.e_electron(r, t, 0u) = e_electron;
      axisym.e_fluid_total(r, t, 0u) = e_total;
      for (std::size_t p = 0; p < full_layout.phi_cells; ++p) {
        full3d.rho(r, t, p) = rho;
        full3d.mom_r(r, t, p) = mom_r;
        full3d.mom_theta(r, t, p) = mom_theta;
        full3d.mom_phi(r, t, p) = 0.0;
        full3d.e_electron(r, t, p) = e_electron;
        full3d.e_fluid_total(r, t, p) = e_total;
      }
      CheckNear(axisym.rho(r, t, 0u), AveragePhi(full3d.rho, r, t),
                1.0e-14, "rho phi average");
      CheckNear(axisym.mom_r(r, t, 0u), AveragePhi(full3d.mom_r, r, t),
                1.0e-14, "mom_r phi average");
      CheckNear(axisym.mom_theta(r, t, 0u), AveragePhi(full3d.mom_theta, r, t),
                1.0e-14, "mom_theta phi average");
      CheckNear(axisym.e_electron(r, t, 0u), AveragePhi(full3d.e_electron, r, t),
                1.0e-14, "e_electron phi average");
      CheckNear(axisym.e_fluid_total(r, t, 0u),
                AveragePhi(full3d.e_fluid_total, r, t),
                1.0e-14, "e_fluid_total phi average");
    }
  }

  const auto axisym_geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{3u, 4u, 1u, 1.0, 2.0});
  const auto full_geometry = dec3d::mesh::BuildSphericalGeometry(
      dec3d::mesh::SphericalMeshDescriptor{3u, 4u, 8u, 1.0, 2.0});
  DEC3D_CHECK(axisym_geometry.valid);
  DEC3D_CHECK(full_geometry.valid);
  const auto axisym_momentum =
      dec3d::state::diagnostics::ComputeAxisymmetricCartesianMomentum(axisym, axisym_geometry);
  const auto full_momentum =
      dec3d::state::diagnostics::ComputeFull3DCartesianMomentum(full3d, full_geometry);
  CheckNear(axisym_momentum.px, 0.0, 0.0, "axisym px");
  CheckNear(axisym_momentum.py, 0.0, 0.0, "axisym py");
  CheckNear(full_momentum.px, 0.0, 1.0e-12, "full3d replicated px");
  CheckNear(full_momentum.py, 0.0, 1.0e-12, "full3d replicated py");
  CheckNear(axisym_momentum.pz, full_momentum.pz, 1.0e-12, "axisym full3d pz");
  return 0;
}
```

- [ ] **Step 4: Run test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_axisymmetric_2d_regression --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_axisymmetric_2d_regression.exe'
```

Expected: executable exits `0`.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add src/state/diagnostics/axisymmetric_vector_diagnostics.hpp tests/contract/test_axisymmetric_2d_regression.cpp CMakeLists.txt
git -C 'F:\dec3d' commit -m 'Add axisymmetric vector diagnostics regression'
```

---

## Task 9: Restart Metadata Safety

**Files:**
- Modify: `src/app/dec3d_app.cpp`
- Modify: `src/io/runtime_output.cpp`
- Modify: `src/io/runtime_output.hpp`
- Test: `tests/contract/test_p5_io_checkpoint_files.cpp`

- [ ] **Step 1: Add failing checkpoint metadata test**

In `tests/contract/test_p5_io_checkpoint_files.cpp`, include the app entrypoint and add a small string replacer in the anonymous namespace:

```cpp
#include "app/dec3d_app.hpp"
```

```cpp
#include <vector>
```

```cpp
std::string ReplaceAll(std::string text, const std::string& from, const std::string& to) {
  std::size_t pos = 0;
  while ((pos = text.find(from, pos)) != std::string::npos) {
    text.replace(pos, from.size(), to);
    pos += to.size();
  }
  return text;
}
```

After the existing restart metadata assertions, add:

```cpp
  DEC3D_CHECK(restart.find("# dimensionality=full_3d") != std::string::npos);
```

Then add an axisymmetric checkpoint write using the same profile:

```cpp
  auto axisym_config = deck.config;
  axisym_config.run.case_name = "axisymmetric_checkpoint_contract";
  axisym_config.run.output_dir = (root / "output_axisym").generic_string();
  axisym_config.mesh.dimensionality = dec3d::io::MeshDimensionality::axisymmetric_2d;
  axisym_config.mesh.phi_cells = 1u;
  axisym_config.mesh.moving_mesh = false;
  const auto axisym_init = dec3d::initialization::InitializeFromRadialProfile(
      axisym_config, profile.profile, profile_path);
  DEC3D_CHECK(axisym_init.success);
  const auto axisym_write = dec3d::io::WriteInitialRuntimeOutputs(
      axisym_config,
      axisym_init.state,
      axisym_init.geometry,
      axisym_init.group_layout,
      profile_path);
  DEC3D_CHECK(axisym_write.success);
  DEC3D_CHECK(axisym_write.restart_checkpoint_written);
  const auto axisym_restart = ReadText(root / "output_axisym" / "restart_000000.snap");
  DEC3D_CHECK(axisym_restart.find("# dimensionality=axisymmetric_2d") !=
              std::string::npos);
  DEC3D_CHECK(axisym_restart.find("# phi_cells=1") != std::string::npos);
```

Add restart mismatch checks in both directions:

```cpp
  const std::vector<std::string> full_deck_axisym_restart_args{
      "dec3d.exe",
      "-input",
      deck_path.string(),
      "-profile",
      profile_path.string(),
      "-restart",
      (root / "output_axisym" / "restart_000000.snap").string()};
  const auto full_deck_axisym_restart =
      dec3d::app::RunDec3DCommandLine(full_deck_axisym_restart_args);
  DEC3D_CHECK(full_deck_axisym_restart.exit_code != 0);
  DEC3D_CHECK(full_deck_axisym_restart.failure_diagnostics.find(
                  "restart dimensionality does not match input deck") !=
              std::string::npos);

  const auto axisym_deck_path = root / "case_axisym_restart_mismatch.in";
  auto axisym_deck_text = ReadText(deck_path);
  axisym_deck_text = ReplaceAll(axisym_deck_text,
                                "geometry = spherical\n",
                                "geometry = spherical\ndimensionality = axisymmetric_2d\n");
  axisym_deck_text = ReplaceAll(axisym_deck_text, "phi_cells = 4", "phi_cells = 1");
  axisym_deck_text = ReplaceAll(axisym_deck_text, "moving_mesh = true", "moving_mesh = false");
  axisym_deck_text = ReplaceAll(axisym_deck_text,
                                (root / "output").generic_string(),
                                (root / "output_axisym_mismatch").generic_string());
  WriteText(axisym_deck_path, axisym_deck_text);
  const std::vector<std::string> axisym_deck_full_restart_args{
      "dec3d.exe",
      "-input",
      axisym_deck_path.string(),
      "-profile",
      profile_path.string(),
      "-restart",
      (root / "output" / "restart_000000.snap").string()};
  const auto axisym_deck_full_restart =
      dec3d::app::RunDec3DCommandLine(axisym_deck_full_restart_args);
  DEC3D_CHECK(axisym_deck_full_restart.exit_code != 0);
  DEC3D_CHECK(axisym_deck_full_restart.failure_diagnostics.find(
                  "restart dimensionality does not match input deck") !=
              std::string::npos);
```

- [ ] **Step 2: Implement metadata write**

In `src/io/runtime_output.cpp`, change the private metadata writer signature:

```cpp
void WriteSnapshotMetadata(
    std::ofstream& out,
    const char* diagnostic_id,
    bool restart_compatible,
    dec3d::io::MeshDimensionality dimensionality,
    const dec3d::state::CanonicalState& state,
    const dec3d::radiation::RadiationGroupLayout& group_layout,
    const std::filesystem::path& profile_path,
    std::size_t step,
    double time_s)
```

Inside that function, after `phi_cells`, write the current checkpoint metadata style:

```cpp
      << "# phi_cells=" << state.layout.phi_cells << '\n'
      << "# dimensionality=" << dec3d::io::ToString(dimensionality) << '\n'
```

Change both `WriteFieldCheckpointFile` overloads and `WriteRestartCheckpointFile` to accept `dec3d::io::MeshDimensionality dimensionality`, and pass it through to `WriteSnapshotMetadata`. At every call site inside `WriteInitialRuntimeOutputs` and `WriteRuntimeStepOutputs`, pass:

```cpp
config.mesh.dimensionality
```

- [ ] **Step 3: Implement metadata read guard**

In `LoadRestartInitialState` inside `src/app/dec3d_app.cpp`, add a local metadata value before the header loop:

```cpp
  std::string restart_dimensionality{"full_3d"};
```

Parse the new metadata in the existing `if (key == ...)` chain:

```cpp
      } else if (key == "dimensionality") {
        restart_dimensionality = value;
```

After the layout comparison and before the restart time check, add:

```cpp
  if (restart_dimensionality != dec3d::io::ToString(config.mesh.dimensionality)) {
    return FailRestartLoad("restart dimensionality does not match input deck", restart_path);
  }
```

- [ ] **Step 4: Run checkpoint test**

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d_contract_p5_io_checkpoint_files --clean-first -- /m
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_checkpoint_files.exe'
```

Expected: executable exits `0`.

- [ ] **Step 5: Commit**

```powershell
git -C 'F:\dec3d' add src/app/dec3d_app.cpp src/io/runtime_output.cpp src/io/runtime_output.hpp tests/contract/test_p5_io_checkpoint_files.cpp
git -C 'F:\dec3d' commit -m 'Record axisymmetric restart dimensionality'
```

---

## Task 10: Final Verification And Benchmark Probe

**Files:**
- Modify only if failures reveal a direct issue in touched files.

- [ ] **Step 1: Clean focused rebuild**

Run:

```powershell
cmake --build 'F:\dec3d\build-hypre' --config Release --target dec3d dec3d_contract_p5_io_input_deck dec3d_contract_spherical_geometry dec3d_contract_spherical_scalar_remap dec3d_contract_hydro_theta_pole_boundary dec3d_contract_hydro_origin_radial_ghost_remap dec3d_contract_generic_implicit_diffusion dec3d_contract_hypre_distributed_diffusion_solver dec3d_contract_p5_io_profile_initializer dec3d_contract_p5_runtime_all_stages dec3d_contract_p5_runtime_distributed_all_stages dec3d_contract_axisymmetric_2d_regression dec3d_contract_p5_io_checkpoint_files --clean-first -- /m
```

Expected: build exit code `0`.

- [ ] **Step 2: Run focused tests**

Run:

```powershell
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_input_deck.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_geometry.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_spherical_scalar_remap.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_theta_pole_boundary.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_hydro_origin_radial_ghost_remap.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_generic_implicit_diffusion.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_hypre_distributed_diffusion_solver.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_profile_initializer.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_all_stages.exe'
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_runtime_distributed_all_stages.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_axisymmetric_2d_regression.exe'
& 'F:\dec3d\build-hypre\Release\dec3d_contract_p5_io_checkpoint_files.exe'
```

Expected: every command exits `0`.

- [ ] **Step 3: Run end-to-end smoke executable**

Run:

```powershell
$out = 'F:\dec3d\analysis\output\axisymmetric_2d_noale_allstages_smoke'
New-Item -ItemType Directory -Force -Path (Join-Path $out 'logs') | Out-Null
& 'C:\Program Files\Microsoft MPI\Bin\mpiexec.exe' -n 2 'F:\dec3d\build-hypre\Release\dec3d.exe' -input 'F:\dec3d\cases\axisymmetric_2d_noale_allstages_smoke.in' -profile 'F:\dec3d\cases\axisymmetric_2d_noale_smoke.pro' > (Join-Path $out 'logs\stdout.txt') 2> (Join-Path $out 'logs\stderr.txt')
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Get-Content (Join-Path $out 'logs\stderr.txt')
Get-Content (Join-Path $out 'logs\stdout.txt') -Tail 20
```

Expected:

- exit code `0`;
- `stderr.txt` empty;
- stdout contains `mesh_dimensionality=axisymmetric_2d`;
- stdout contains `phi_sweep_executed=false`;
- stdout contains `global_phi_coupling_count=0`;
- stdout contains `axisymmetric_invariant_ok=true`.

- [ ] **Step 4: Push branch**

```powershell
git -C 'F:\dec3d' status --short --branch
git -C 'F:\dec3d' push
```

Expected: branch `codex/axisymmetric-2d` pushed to `origin/codex/axisymmetric-2d`.

---

## Self-Review Checklist

- Spec coverage:
  - Preflight branch/build/MPI/TOPS/baseline checks: Task 0.
  - Input dimensionality, `phi_cells=1`, no ALE, `m=0`: Task 1.
  - Existing `full_3d` decks remain backward compatible: Task 1 and final tests.
  - `axisymmetric_2d` requires full `theta=[0,pi]` and `phi=[0,2*pi]` deck ranges: Task 1.
  - Full `2*pi` volume and area: Task 2.
  - Pole/origin remap with `phi=1`: Task 3.
  - No phi derivative/coupling/self-neighbor/duplicates, with pre-CSR self-neighbor diagnostics: Task 4.
  - Matrix diagnostics define `phi_coupling_count` as nonzero off-diagonal phi face-pair conductance count: Task 4.
  - Perturbation normalization and `mom_phi=0`: Task 5.
  - Hydro skips phi sweep, disables phi macro/coarse sweep, and preserves spherical source terms: Task 6A.
  - Runtime `mom_phi` invariant uses absolute-plus-relative tolerance and reports tolerances: Task 6B.
  - All stages `H,T,E,R,A`: Task 7.
  - Vector diagnostics and 2D-vs-3D m=0 consistency: Task 8.
  - Restart dimensionality guard covers both mismatch directions: Task 9.
  - Clean rebuild and end-to-end evidence: Task 10.
- Geometry/source/output scaling:
  - No diagnostic, source, or output reducer applies an extra `2*pi` outside geometry volumes and areas.
  - `phi_sweep_configured` appears only in input/deck diagnostics; `phi_sweep_executed` appears only in runtime/hydro diagnostics.
- Placeholder scan: no forbidden placeholder terms or unspecified `write tests` steps remain.
- Type consistency:
  - `MeshDimensionality::axisymmetric_2d` is introduced in Task 1 and reused later.
  - Report token names are consistent: `mesh_dimensionality`, `phi_sweep_configured`, `phi_sweep_executed`, `global_phi_coupling_count`, `global_phi_self_neighbor_attempt_count`, `global_duplicate_column_row_count`, `axisymmetric_invariant_ok`.
  - Test names and build targets match existing CMake naming style.
