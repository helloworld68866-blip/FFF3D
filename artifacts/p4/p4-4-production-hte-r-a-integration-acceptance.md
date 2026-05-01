# P4-4 Production H/T/E/R/A Integration Acceptance

Date: 2026-04-29

## Scope

P4-4 integrates the alpha `A` stage into the production split chain:

```text
H -> T -> E -> R -> A
```

This is a non-benchmark runtime-integration slice. It does not add new alpha
physics beyond the already-built P4-2 distributed Atzeni one-group alpha solve
and P4-3 alpha hydro-side scalar advection.

## Implemented

- Added a backend-neutral P4 production orchestrator hook layer.
- Registered and executed `H/T/E/R/A` in order for `phase_id=p4`.
- Passed the same hydro-relaxed `dt_s` to all five stages.
- Required P3-8 radiation hydro-side evidence from the `H/R` stage reports:
  `radiation_hydro_coupling_mode=hydro_hllc_species_scalar_radiation_pressure`.
- Required P4-3 alpha hydro-side evidence from the `H` stage report:
  `alpha_hydro_coupling_mode=hydro_hllc_species_scalar_alpha_pressure`.
- Required `A` to consume:
  `alpha_geometry_epoch=post_H_committed_ALE_geometry`,
  `alpha_state_epoch=post_H_committed`, and
  `alpha_thermal_state_epoch=post_R_committed`.
- Required `A` old-time coefficients to mean lagged to `A` stage entry,
  after `H/T/E/R` have committed.
- Enforced the `A` write set:
  `alpha_state,e_electron,e_fluid_total`.
- Preserved nested alpha failure diagnostics when the `A` hook fails.
- Preserved stage-local atomic semantics:
  `step_atomic_across_H_T_E_R_A=false`.

## Failure Semantics

P4-4 does not implement whole-step rollback. If `A` fails:

```text
success=false
advance_to_next_timestep=false
committed_prior_stages=H,T,E,R
a_stage_published=false
```

The failed `A` stage must publish nothing. Already committed prior stages remain
committed and the caller must not advance to the next timestep.

## Not Implemented

- Fuel depletion or separate D/T authoritative species.
- Multigroup alpha.
- Alpha benchmark parity.
- Radiation or alpha momentum-force feedback beyond scalar hydro-side terms.
- Direct HYPRE dependency from the backend-neutral runtime orchestrator.

## Required Verification

```text
dec3d_contract_p4_production_orchestrator
dec3d_contract_p3_production_orchestrator
```

The P4 contract must demonstrate:

```text
H/T/E/R/A success order
missing A hook preflight failure
H failure skips T/E/R/A
T failure skips E/R/A
E failure skips R/A
R failure skips A
A failure preserves committed H/T/E/R and nested alpha diagnostics
old frozen radiation hydro mode rejected
missing alpha hydro terms rejected
stale A geometry/state/thermal epochs rejected
wrong A write set rejected
dt_s=0 alpha no-op success
runtime phase_id=p4 registration and ingestion
phase_id=p3 rejected by P4 orchestrator
P3 orchestrator still rejects alpha registration
```

Full clean-first default and HYPRE regressions must pass before P4-4 closure.
