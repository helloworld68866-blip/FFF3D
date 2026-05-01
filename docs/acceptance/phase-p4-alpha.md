# Phase P4 Alpha Acceptance

Document version: v0.2  
Last updated: 2026-04-29

## Goal

`P4` adds alpha-particle transport and deposition as a distinct physical stage.

## Required Stages

- `H`
- `T`
- `E`
- `R`
- `A`

## Authoritative Write-Set

`A` must write:

- `AuthoritativeAlphaState`
- `E_electron`
- `E_fluid_total` where deposition coupling requires it

## Acceptance Artifacts

- clean rebuild report
- failing-first alpha test evidence
- `H/T/E/R/A` stage reports with execution evidence
- alpha state update summary
- alpha deposition evidence
- alpha numerical check summary
- assumption-ledger delta
- `artifacts/p4/p4-0-alpha-state-provider-contract-acceptance.md`
- `artifacts/p4/p4-1-bosch-hale-dt-reference.md`
- `artifacts/p4/p4-1-one-group-alpha-operator-acceptance.md`
- `artifacts/p4/p4-2-distributed-alpha-solve-acceptance.md`
- `artifacts/p4/p4-3-alpha-hydro-side-advection-acceptance.md`
- `artifacts/p4/p4-4-production-hte-r-a-integration-acceptance.md`
- `artifacts/p4/p4-b-alpha-benchmark-ladder-acceptance.md`

## Current P4 Slice Status

### P4-0

Completed as an alpha state / units / composition / provider contract slice.

Implemented:

- authoritative `alpha_state` contract
- `P_alpha=(2/3)epsilon_alpha` and `chi_alpha=P_alpha^(3/5)` roundtrip
- equimolar DT composition closure
- thesis-Spitzer `tau_alphae` and Atzeni one-group drag coefficient provider

### P4-1

Completed as a single-rank callable Atzeni one-group alpha diffusion/source operator.

Implemented:

- Bosch-Hale DT reactivity path with local reference values and `Ti_old` input
- generic diffusion mapping `A=1`, `D=D_alpha`, `C=-1/tau_alphae`, `B=birth_source`
- implicit alpha-electron drag deposition using `epsilon_alpha_new`
- transactional writeback to `alpha_state`, `e_electron`, and `e_fluid_total`
- volume-integrated single-rank budget diagnostics

Not yet implemented:

- distributed alpha HYPRE
- production `H/T/E/R/A`
- alpha hydro-side scalar advection
- fuel depletion
- alpha benchmark parity
- multigroup alpha

### P4-2

Completed as a callable distributed provider-fed Atzeni one-group alpha solve.

Implemented:

- provider-fed old-time alpha coefficients with Bosch-Hale DT reactivity
- distributed generic diffusion assembly for `epsilon_alpha`
- real distributed HYPRE ParCSR GMRES + BoomerAMG solve
- owned-slab coefficient build and owned-slab writeback
- seam `D_alpha` halo / off-rank conductance evidence
- implicit alpha-electron drag deposition using `epsilon_alpha_new`
- two-stage MPI transaction and zero publish on global failure
- MPI allreduced volume-integrated budget diagnostics

Still not implemented:

- production `H/T/E/R/A`
- alpha hydro-side scalar advection
- fuel depletion or separate D/T authoritative species
- alpha benchmark parity
- multigroup alpha

### P4-3

Completed as alpha hydro-side scalar advection in `H`.

Implemented:

- `chi_alpha=P_alpha^(3/5)` derived from authoritative `alpha_state`
- H-stage scalar transport of `chi_alpha`
- H writeback recovery of authoritative `alpha_state`
- controlled compression/expansion oracle with `5/3` exponent
- static H-path evidence that alpha travels through the hydro scalar bundle

Still not implemented:

- production `H/T/E/R/A`
- fuel depletion or separate D/T authoritative species
- alpha benchmark parity
- multigroup alpha
- alpha momentum-force feedback beyond scalar advection

### P4-4

Completed as production `H/T/E/R/A` integration, non-benchmark.

Implemented:

- backend-neutral P4 production orchestrator hook layer
- `phase_id=p4` runtime registration for `H/T/E/R/A`
- `H -> T -> E -> R -> A` execution order with shared `dt_s`
- P3-8 radiation hydro-side evidence required from stage diagnostics
- P4-3 alpha hydro-side evidence required from H-stage diagnostics
- A-stage epoch checks:
  `post-H` ALE geometry, `post-H` alpha state, and `post-R` thermal state
- A-stage write-set validation for `alpha_state/e_electron/e_fluid_total`
- nested alpha failure diagnostics preservation
- stage-local atomic semantics with no whole-step rollback claim

Still not implemented:

- fuel depletion or separate D/T authoritative species
- alpha benchmark parity
- multigroup alpha
- radiation or alpha momentum-force feedback beyond scalar advection

### P4-B

Completed as a standalone frozen-fluid alpha benchmark ladder for B0/B1/B3.

Implemented:

- B0 operator and budget oracles
- B1 frozen-fluid alpha hotspot serial-vs-distributed cross-validation
- B3 distributed timing artifact smoke
- scalar-origin-remap benchmark boundary from `r=0`
- Python-generated visual artifacts for profiles, budget residuals,
  serial-vs-distributed pointwise error, and timing breakdown

Still not implemented:

- external alpha benchmark parity
- fuel depletion or separate D/T authoritative species
- multigroup alpha
- full strong/weak scaling benchmark campaign

## Failure Conditions

- `A` not registered
- alpha behavior hidden inside thermal or radiation source logic
- alpha stage lacks independent execution evidence
- alpha authoritative state not owned by canonical state
- diagnostics missing for alpha transport or deposition
- P4 production chain accepts missing P3-8 radiation hydro-side diagnostics
- P4 production chain accepts missing P4-3 alpha hydro-side diagnostics
- `A` writes fields outside `alpha_state/e_electron/e_fluid_total`
