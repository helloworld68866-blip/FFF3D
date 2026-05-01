# Phase P2 Thermal Acceptance

Document version: v0.2  
Last updated: 2026-04-27

## Goal

`P2` adds two-temperature thermal conduction and electron-ion equilibration on top of the completed hydro path.

## Required Stages

- `H`
- `T`
- `E`

Stages `R` and `A` must remain unregistered in `P2`.

## Authoritative Write-Set

`T` must write:

- `E_electron`
- `E_fluid_total`

`E` must write:

- `E_electron`

`E` must preserve:

- `E_fluid_total`, up to accepted floating-point tolerance

## Acceptance Artifacts

- clean rebuild report
- failing-first thermal and equilibration test evidence
- `H/T/E` stage reports with execution evidence
- matrix build and solve evidence for thermal paths
- thermal numerical check summary
- write-set conformance summary
- assumption-ledger delta
- P2-8 benchmark harness acceptance summary:
  `F:\dec3d\artifacts\p2\p2-8-acceptance-summary.md`

## P2-8 Benchmark Harness

P2-8 is a benchmark and harness closeout layer. It does not add thermal
physics and does not claim LILAC parity without authoritative external LILAC
reference data.

Accepted harness families:

- Woo 5.5.4 slab-style thermal benchmark descriptors:
  `spitzer_pure_thermal`, `lee_more_pure_thermal`,
  `lee_more_thermal_plus_ei`
- constant-kappa exact operator checks:
  `exact_l1m1_constant_kappa`, `exact_l0_radial_constant_kappa`
- production-path smoke/stress checks:
  `flux_limiter_stress`, `ei_0d_thesis_spitzer`,
  `distributed_hte_smoke`

Required reference policy when no external LILAC data is present:

- `lilac_reference_available=false`
- `parity_claim_allowed=false`

## ALE H-to-T Geometry Contract

If hydro stage `H` commits an ALE mesh update, thermal stage `T` must consume
the committed post-hydro geometry. The production orchestrator must fail before
`E` if `T` reports a stale geometry epoch.

Required evidence:

- stale `T` geometry epoch fails with `stage_id=T`
- matching `T` geometry epoch succeeds
- successful report contains `ale_geometry_committed=true`
- successful report contains `thermal_geometry_epoch=<committed_epoch>`

## Failure Conditions

- `T` or `E` not registered
- thermal solve path lacks matrix or solve evidence
- equilibration silently changes `E_fluid_total`
- stale recovered thermodynamic caches used across stages
- stale geometry used by `T` after an ALE `H` commit
- diagnostics missing for any required stage
