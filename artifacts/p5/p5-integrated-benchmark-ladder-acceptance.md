# P5 Integrated Benchmark Ladder Acceptance

Date: 2026-04-29

## Scope

This slice adds the P5 integrated benchmark ladder harness for the existing
`H -> T -> E -> R -> A` production chain. It covers:

- clean spherical descriptor support;
- `l=2,m=0` and `l=4,m=0` Legendre radial-velocity perturbation descriptors;
- PPM + macro + ALE + MPI variant metadata;
- MPI12 and no-ALE comparison metadata;
- initial radial profile preparation from
  `F:/python_project/cases/image_profile_100ps_hydro.yaml`;
- CSV/JSON benchmark artifacts and Python PNG/GIF visualization artifacts.

## Initial Profile Source

The P5 initial profile source is:

```text
F:/python_project/cases/image_profile_100ps_hydro.yaml
```

The YAML resolves to:

```text
F:/python_project/intialState/extracted_initial_profiles.npz
```

The generated C++-friendly CSV is:

```text
F:/dec3d/analysis/output/p5_initial_profile_from_image_profile_100ps_hydro.csv
```

The supplied profile spans `0..50 um`. Therefore the original Woo shot-77068
`r0=68 um` is not silently reused. The first P5 descriptor uses
`perturbation_r0_um=24.070450097847356`, the density-gradient feature extracted
from the supplied profile, and diagnostics report:

```text
r0_matches_woo_77068=false
```

## Verification

Focused checks executed:

```text
cmake --build F:/dec3d/build --config Debug --target dec3d_contract_p5_integrated_benchmark_contract --clean-first -- /m:4
ctest --test-dir F:/dec3d/build -C Debug -R dec3d_contract_p5_integrated_benchmark_contract --output-on-failure

cmake --build F:/dec3d/build-hypre --config Debug --target dec3d_contract_p5_integrated_benchmark_hypre --clean-first -- /m:4
ctest --test-dir F:/dec3d/build-hypre -C Debug -R dec3d_contract_p5_integrated_benchmark_hypre --output-on-failure
```

Both focused contract tests passed during implementation.

## Generated P5 Artifacts

The HYPRE/MPI contract writes artifacts to:

```text
F:/dec3d/analysis/output/p5_integrated_benchmark_hypre_contract
```

The current artifact set includes:

- `p5_case_descriptor.json`
- `p5_stage_reports.json`
- `p5_radial_profiles.csv`
- `p5_mode_amplitudes.csv`
- `p5_budget_history.csv`
- `p5_stage_timing.csv`
- `p5_variant_comparison_summary.json`
- `p5_manifest.json`
- `p5_radial_profiles.png`
- `p5_mode_amplitudes.png`
- `p5_budget_residuals.png`
- `p5_p5_legendre_l2_m0_primary_rz_evolution.gif`

## Claim Boundary

This slice is an integrated benchmark harness and artifact pipeline. It does
not add new P1-P4 physics.

Allowed claim:

```text
P5 integrated benchmark ladder harness exists and can execute the P4
production orchestration path with P5 descriptors, initial-profile provenance,
variant metadata, diagnostics, and visual artifacts.
```

Forbidden claims:

```text
P5 full MPI24 long-run benchmark acceptance is complete.
P5 provides external LILAC or 1-D parity.
P5 adds fuel depletion, radiation force feedback, or new alpha/radiation physics.
```

## Remaining Work

- Run the full intended MPI24 primary variants and MPI12/no-ALE comparisons as
  long benchmark jobs outside the focused contract smoke.
- Replace the current compact smoke artifact generation with full-resolution
  physical profile sampling once long-run runtime cost is acceptable.
- Promote the P5 visual artifacts from smoke evidence to final benchmark
  evidence only after the MPI24/MPI12/no-ALE benchmark matrix is run.
