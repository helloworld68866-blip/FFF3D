# P3-2 Radiation-Matter Coupling Acceptance

Date: 2026-04-27

## Scope

P3-2 adds frozen-hydro one-group gray radiation-electron coupling. The operator
updates `Ug`, `e_electron`, and `e_fluid_total` in one staged transaction.

The coupling uses the same implicit source term used by the radiation solve:

```text
source_gain_radiation = dt * c * kappaP * (Bgray - U_new)
delta_electron = -source_gain_radiation
```

## Required Diagnostics

- `diagnostic_id=p3.radiation.one_group_gray_matter_coupling`
- `stage_id=R`
- `radiation_electron_coupling_enabled=true`
- `backend_requested`
- `backend_executed`
- `radiation_report_present=true`
- `thermodynamic_recovery_report_present=true`
- `delta_radiation_total`
- `delta_electron_total`
- `source_gain_radiation_total`
- `boundary_leak_total`
- `radiation_electron_exchange_residual`
- `global_radiation_plus_electron_residual`
- `updated_fields=radiation_groups,e_electron,e_fluid_total`

All global exchange and leak budgets are volume-integrated quantities.
Unweighted cell sums are not accepted as conservation evidence.

## Contract Tests

- one-cell implicit absorption/emission source oracle;
- zero-flux volume-integrated `U + e_electron` conservation;
- Marshak boundary leak accounting separate from matter exchange;
- exact write-set check;
- `dt_s=0` no-change diagnostics;
- electron floor failure with zero publish;
- ion floor failure with zero publish;
- missing required diagnostics fail.

## Non-Claims

This slice does not implement:

- multigroup transport;
- opacity, Planck/Rosseland, or `Bgray` providers;
- radiation flux limiting;
- radiation advection;
- radiation pressure work;
- production `H/T/E/R` registration;
- Woo 5.4 benchmark or LILAC parity.

## Verification

Focused implementation verification:

```text
cmake --build F:\dec3d\build --config Debug --target dec3d_contract_radiation_one_group_gray_diffusion --clean-first
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_radiation_one_group_gray_diffusion" --output-on-failure
```

Result:

```text
1/1 passed
```

Focused P3 contract verification:

```text
cmake --build F:\dec3d\build --config Debug --clean-first --target dec3d_contract_radiation_group_layout dec3d_contract_radiation_one_group_gray_diffusion
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_radiation_group_layout|dec3d_contract_radiation_one_group_gray_diffusion" --output-on-failure
```

Result:

```text
2/2 passed
```

Focused P2/P3 regression verification:

```text
cmake --build F:\dec3d\build --config Debug --clean-first --target dec3d_contract_thermodynamic_recovery dec3d_contract_electron_ion_equilibration dec3d_contract_thermal_conduction_constant_kappa dec3d_contract_radiation_one_group_gray_diffusion
ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_thermodynamic_recovery|dec3d_contract_electron_ion_equilibration|dec3d_contract_thermal_conduction_constant_kappa|dec3d_contract_radiation_one_group_gray_diffusion" --output-on-failure
```

Result:

```text
4/4 passed
```

Full default verification:

```text
cmake --build F:\dec3d\build --config Debug --clean-first -- /m:1
ctest --test-dir F:\dec3d\build -C Debug --output-on-failure
```

Result:

```text
98/98 passed
```

Full HYPRE-gated verification:

```text
cmake --build F:\dec3d\build-hypre --config Debug --clean-first -- /m:1
ctest --test-dir F:\dec3d\build-hypre -C Debug --output-on-failure
```

Result:

```text
104/104 passed
```

Observed non-blocking warnings:

- existing manual Noh benchmark `C4127` constant-condition warnings;
- existing HYPRE/MSVC runtime-library `LNK4098` warnings.
