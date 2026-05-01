# P3-3 Multigroup Gray Radiation-Matter Acceptance

Date: 2026-04-27

## Scope

P3-3 adds a callable frozen-hydro multigroup gray radiation-matter coupling
operator. It solves every explicit frequency group independently with fixed
user-supplied coefficients, then applies a single conservative electron and
total-fluid writeback accumulated over all group implicit source terms.

## Required Diagnostics

- `diagnostic_id=p3.radiation.multigroup_gray_matter_coupling`
- `stage_id=R`
- `radiation_group_mode=explicit_frequency_groups`
- `frequency_edges_unit=Hz`
- `group_layout_report_present=true`
- `per_group_solve_count`
- `all_groups_updated`
- `per_group_matrix_report_present=true`
- `per_group_solve_report_present=true`
- `first_failing_group_index`
- `delta_radiation_total_all_groups`
- `delta_electron_total`
- `source_gain_radiation_total_all_groups`
- `boundary_leak_total_all_groups`
- `radiation_electron_exchange_residual`
- `global_radiation_plus_electron_residual`

All global budgets are volume-integrated first per group, then summed over
groups. Unweighted cell sums are not accepted as conservation evidence.

## Contract Tests

- two-group implicit source oracle;
- `group_count=1` explicit-frequency regression against P3-2;
- all groups update, not only group 0;
- zero-flux all-group `U + e_electron` conservation;
- Marshak all-group boundary leak accounting;
- later group failure with zero publish;
- `dt_s=0` no-change diagnostics;
- missing required diagnostics fail.

## Non-Claims

This slice does not implement:

- opacity or `Bg(Te)` providers;
- radiation flux limiter;
- distributed/HYPRE multigroup radiation;
- group-level parallelism;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` registration;
- Woo 5.4 benchmark or LILAC parity.

## Verification

Focused default-build P3-3 contract:

```text
cmake --build F:\dec3d\build --config Debug --target dec3d_contract_radiation_multigroup_gray_coupling --clean-first -- /m:4
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_radiation_multigroup_gray_coupling" --output-on-failure
```

Result: `1/1` passed.

Focused P3 radiation contracts:

```text
cmake --build F:\dec3d\build --config Debug --clean-first --target dec3d_contract_radiation_group_layout dec3d_contract_radiation_one_group_gray_diffusion dec3d_contract_radiation_multigroup_gray_coupling -- /m:4
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_radiation_group_layout|dec3d_contract_radiation_one_group_gray_diffusion|dec3d_contract_radiation_multigroup_gray_coupling" --output-on-failure
```

Result: `3/3` passed.

Focused P2/P3 regression:

```text
cmake --build F:\dec3d\build --config Debug --clean-first --target dec3d_contract_thermodynamic_recovery dec3d_contract_electron_ion_equilibration dec3d_contract_thermal_conduction_constant_kappa dec3d_contract_radiation_one_group_gray_diffusion dec3d_contract_radiation_multigroup_gray_coupling -- /m:4
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_thermodynamic_recovery|dec3d_contract_electron_ion_equilibration|dec3d_contract_thermal_conduction_constant_kappa|dec3d_contract_radiation_one_group_gray_diffusion|dec3d_contract_radiation_multigroup_gray_coupling" --output-on-failure
```

Result: `5/5` passed.

Default clean-first full suite:

```text
cmake --build F:\dec3d\build --config Debug --clean-first -- /m:4
ctest --test-dir F:\dec3d\build -C Debug --output-on-failure
```

Result: `99/99` passed.

HYPRE clean-first full suite:

```text
cmake --build F:\dec3d\build-hypre --config Debug --clean-first -- /m:4
ctest --test-dir F:\dec3d\build-hypre -C Debug --output-on-failure
```

Result: `105/105` passed.

The P3-3 contract still reports `backend_executed=serial_dense_reference`;
the HYPRE build regression only proves the new slice does not break the
HYPRE-enabled build and test suite.
