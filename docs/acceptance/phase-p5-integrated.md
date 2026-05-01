# Phase P5 Integrated Acceptance

Document version: v0.2  
Last updated: 2026-04-29

## Goal

`P5` integrates the already-real `H/T/E/R/A` stages into the thesis-style orchestration flow. It does not redefine earlier phase completion.

## Required Stages

- `H`
- `T`
- `E`
- `R`
- `A`

## Authoritative Write-Set

`P5` does not create new physics write authority. It validates that the combined stage sequence preserves the write-set contracts already established in `P1` through `P4`.

## Acceptance Artifacts

- clean rebuild report
- integrated failing-first regression evidence
- full `H/T/E/R/A` orchestration stage report
- combined write-set conformance summary
- integrated numerical regression summary
- diagnostics completeness summary
- assumption-ledger delta

## Current P5 Benchmark Ladder Slice

The current P5 implementation adds an integrated benchmark ladder harness for:

- `p5_clean_spherical`
- `p5_legendre_l2_m0`
- `p5_legendre_l4_m0`

The benchmark initial radial profile is sourced from:

```text
F:/python_project/cases/image_profile_100ps_hydro.yaml
```

That profile resolves to:

```text
F:/python_project/intialState/extracted_initial_profiles.npz
```

The supplied profile spans `0..50 um`; therefore the original Woo shot-77068
`r0=68 um` is not silently reused. P5 diagnostics must report:

```text
perturbation_r0_source=p5_case_descriptor
r0_matches_woo_77068=false
```

Focused acceptance evidence is recorded in:

```text
F:/dec3d/artifacts/p5/p5-integrated-benchmark-ladder-acceptance.md
```

The focused HYPRE/MPI contract generates the visual and tabular artifact set in:

```text
F:/dec3d/analysis/output/p5_integrated_benchmark_hypre_contract
```

This focused contract is a benchmark harness smoke, not the final MPI24 long-run
benchmark matrix.

## Failure Conditions

- any earlier phase write-set contract is weakened or rewritten
- any stage is bypassed by an integration-only solver
- orchestration cannot prove the required stage order
- integrated regression is used to excuse missing lower-phase evidence
- document thresholds or done definitions changed without versioned reason
