# P4-1 One-Group Alpha Operator Acceptance

Date: 2026-04-28

## Scope

P4-1 implements a single-rank callable Atzeni one-group alpha
diffusion/source operator.

## Implemented

- `bosch_hale_dt` DT reactivity provider with local reference values.
- Generic diffusion mapping: `A=1`, `D=D_alpha`, `C=-1/tau_alphae`, `B=birth_source`.
- Implicit drag deposition to `e_electron` using `epsilon_alpha_new`.
- Transactional writeback to `alpha_state`, `e_electron`, and `e_fluid_total`.
- Volume-integrated single-rank budget diagnostics.

## Not Implemented

- Distributed alpha HYPRE.
- Production `H/T/E/R/A`.
- Alpha hydro-side scalar advection.
- Fuel depletion.
- Alpha benchmark parity.
- Multigroup alpha.

## Required Verification

```text
dec3d_contract_alpha_contract
dec3d_contract_alpha_one_group_operator
```

Full clean-first default and HYPRE build regressions must pass before closure.
