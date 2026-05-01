# P4-2 Distributed Alpha Solve Acceptance

Date: 2026-04-28

## Scope

P4-2 promotes the P4-1 single-rank alpha operator to a callable distributed
Atzeni one-group alpha transport solve.

## Implemented

- Provider-fed old-time alpha coefficients from P4-0/P4-1:
  `tau_alphae`, `D_alpha`, `birth_source`, and Bosch-Hale DT reactivity.
- Distributed generic diffusion assembly for `epsilon_alpha`.
- Real HYPRE ParCSR GMRES + BoomerAMG solve.
- Owned-slab coefficient build and owned-slab writeback.
- Seam `D_alpha` halo evidence for off-rank radial face conductance.
- Implicit alpha-electron drag deposition using `epsilon_alpha_new`.
- Two-stage MPI transaction:
  `local_stage_ok -> global_stage_ok -> local_publishable -> global_publish_ok`.
- Transactional writeback to `alpha_state`, `e_electron`, and `e_fluid_total`.
- MPI allreduced volume-integrated alpha/electron budget diagnostics.

## Not Implemented

- Production `H/T/E/R/A` runtime registration.
- Alpha hydro-side scalar advection in `H`.
- Fuel depletion or separate D/T authoritative species.
- Alpha benchmark parity.
- Multigroup alpha transport.

## Required Verification

```text
dec3d_contract_alpha_distributed_one_group_operator
```

The focused HYPRE-gated contract must demonstrate:

```text
dt_s=0 no-change
two-rank HYPRE execution
distributed-vs-serial parity
provider failure -> global no publish
publish preflight failure -> global no publish
missing diagnostics fail
rank0_gather_solve_used=false
serial_dense_fallback_used=false
```

Full clean-first default and HYPRE build regressions must pass before P4-2
closure.
