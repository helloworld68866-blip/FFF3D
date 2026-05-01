# P4-B Alpha Benchmark Ladder Acceptance

Date: 2026-04-29

## Scope

P4-B is a standalone frozen-fluid alpha benchmark ladder. It verifies the
P4-1/P4-2 Atzeni one-group alpha operator and distributed alpha path; it does
not add new alpha physics.

Implemented layers:

- `P4-B0`: operator and budget oracles.
- `P4-B1`: frozen-fluid alpha hotspot, serial-vs-distributed cross-validation.
- `P4-B3`: distributed timing artifact smoke.

Intentionally out of scope:

- Fuel depletion or separate authoritative D/T species.
- Multigroup alpha.
- External alpha parity claims.
- Production `H/T/E/R/A` registration changes.

## Boundary And Model Contract

The benchmark uses the thesis-consistent P4 alpha model already implemented in
P4-1/P4-2:

- `alpha_transport_model=atzeni_one_group`
- `reactivity_model_executed=bosch_hale_dt`
- `composition_model=equimolar_dt_from_p2_recovery`
- `fuel_depletion_enabled=false`
- `separate_dt_species_authoritative=false`

The 1-D spherical benchmark starts at the origin and uses the scalar origin
remap boundary contract:

- `boundary_model=benchmark_1d_spherical_scalar_origin_remap_outer_zero_flux`
- inner radial boundary: `scalar_origin_remap_required`
- outer radial boundary: `neumann_zero_flux`
- even `phi` count for scalar origin remap topology

## Acceptance Evidence

Focused default contract:

```text
cmake --build F:\dec3d\build --config Debug --clean-first --target dec3d_contract_p4_alpha_benchmark_contract -- /m:4
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_p4_alpha_benchmark_contract" --output-on-failure
```

Focused HYPRE/MPI contract:

```text
cmake --build F:\dec3d\build-hypre --config Debug --clean-first --target dec3d_contract_p4_alpha_benchmark_hypre -- /m:4
ctest --test-dir F:\dec3d\build-hypre -C Debug -R "dec3d_contract_p4_alpha_benchmark_hypre" --output-on-failure
```

Generated artifact directory:

```text
F:\dec3d\analysis\output\p4_alpha_benchmark_hypre_contract
```

Required artifacts:

- `p4_b0_operator_oracles.json`
- `p4_b1_alpha_hotspot_case_descriptor.json`
- `p4_b1_serial_profiles.csv`
- `p4_b1_distributed_profiles.csv`
- `p4_b1_budget_summary.csv`
- `p4_b1_serial_vs_distributed_summary.json`
- `p4_b3_distributed_performance.csv`
- `p4_benchmark_manifest.json`

Visual artifacts:

- `p4_b1_profiles.png`
- `p4_b1_budget_closure.png`
- `p4_b1_serial_distributed_error.png`
- `p4_b3_timing_breakdown.png`

## Numerical Evidence

The B0 oracle report closes the alpha/electron budget:

```text
global_alpha_plus_electron_budget_residual=0
```

The B1 serial-vs-distributed comparison is code-internal cross-validation, not
external parity:

```text
serial_distributed_cross_validation_present=true
rank0_gather_solve_used=false
epsilon_alpha_relative_linf < 1e-8
Te_linf < 1e-10
```

The B3 timing CSV is a timing artifact smoke, not a strong/weak scaling claim.

## Claim Boundary

Allowed claim:

```text
P4-B standalone alpha benchmark ladder is implemented for B0/B1/B3: operator
oracles, frozen-fluid serial-vs-distributed alpha hotspot cross-validation, and
distributed timing artifact smoke.
```

Forbidden claims:

```text
Alpha benchmark parity passed.
Fuel depletion is implemented.
Multigroup alpha is implemented.
P4-B proves production H/T/E/R/A behavior.
P4-B3 is a completed scaling benchmark.
```
