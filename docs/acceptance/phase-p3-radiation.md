# Phase P3 Radiation Acceptance

Document version: v0.4  
Last updated: 2026-04-27

## Goal

`P3` adds multigroup radiation transport as a real implicit stage.

## Required Stages

- `H`
- `T`
- `E`
- `R`

Stage `A` must remain unregistered in `P3`.

## Authoritative Write-Set

`R` must write:

- at least one `U_rad[g]`
- `E_electron`
- `E_fluid_total` where coupling requires it

## Acceptance Artifacts

- clean rebuild report
- failing-first radiation test evidence
- `H/T/E/R` stage reports with execution evidence
- radiation matrix build and solve evidence
- radiation-group field update summary
- radiation numerical check summary
- assumption-ledger delta

## Failure Conditions

- `R` not registered
- no radiation group shows real update evidence
- radiation path executes without complete solve diagnostics
- radiation stage writes fields outside its contract
- placeholder radiation registration used instead of real transport

## P3-0 / P3-1a Status

Implemented scope:

- explicit radiation group/state contract;
- one-group gray full-spectrum radiation layout;
- frozen-hydro one-group gray radiation diffusion operator;
- radiation update writes only `radiation_groups`;
- P3-1a diagnostics explicitly report `marshak_enabled=false` and
  `parity_claim_allowed=false`.

Not yet claimed:

- Woo 5.4 radiation benchmark parity;
- multigroup radiation;
- radiation-matter electron energy writeback;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` runtime registration.

## P3-1b Status

Implemented scope:

- thesis Marshak/P1 vacuum boundary as an explicit one-group radiation option;
- generic diffusion `outer_radial=radiation_marshak_vacuum` Robin boundary support;
- Marshak outer boundary contributes a real diagonal leakage term to the CSR matrix;
- one-cell Marshak loss oracle;
- one-cell Marshak plus absorption/emission oracle;
- `dt_s=0 + Marshak` no-change regression;
- radiation diagnostics report `radiation_boundary_model=thesis_marshak_vacuum`,
  `marshak_enabled=true`, and `parity_claim_allowed=false`.

Not yet claimed:

- Woo 5.4 radiation benchmark parity;
- multigroup radiation;
- radiation-matter electron energy writeback;
- radiation flux limiter;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` runtime registration.

## P3-2 Status

Implemented scope:

- frozen-hydro one-group gray radiation-matter coupling;
- electron writeback from the same implicit `c kappaP (Bgray - U_new)` source
  used by the radiation solve;
- authoritative write set:
  `radiation_groups,e_electron,e_fluid_total`;
- volume-integrated radiation/electron exchange diagnostics;
- zero-flux/scalar-remap `U + e_electron` conservation check;
- Marshak boundary leakage separated from radiation-electron exchange;
- staged writeback with zero canonical publish on thermodynamic failure;
- diagnostics include `stage_id=R`, backend request/execution evidence, nested
  radiation report presence, and thermodynamic recovery report presence.

Not yet claimed:

- Woo 5.4 radiation benchmark parity;
- multigroup radiation;
- opacity, Planck/Rosseland, or `Bgray` provider;
- radiation flux limiter;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` runtime registration.

## P3-3 Status

Implemented scope:

- explicit-frequency multigroup gray radiation-matter coupling;
- one serial reference diffusion solve per group;
- `group_count=1` explicit-frequency regression to P3-2 semantics;
- group-indexed fixed user-supplied `Dbar/kappaP/Bg` coefficient arrays;
- one staged electron and total-fluid writeback after all groups solve;
- volume-integrated all-group radiation/electron exchange diagnostics;
- Marshak boundary leakage summed over groups and separated from matter exchange;
- zero canonical publish on any group failure.

Not yet claimed:

- Woo 5.4 radiation benchmark parity;
- opacity, Planck/Rosseland, or `Bg(Te)` provider;
- radiation flux limiter;
- distributed/HYPRE multigroup radiation solve;
- group-parallel execution;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` runtime registration.
