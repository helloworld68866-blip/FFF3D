# P4-3 Alpha Hydro-Side Scalar Advection Acceptance

Date: 2026-04-28

## Scope

P4-3 implements alpha hydro-side scalar advection in the `H` stage.

## Implemented

- Derived `chi_alpha=P_alpha^(3/5)` from authoritative `alpha_state`.
- Transported `chi_alpha` through the hydro scalar bundle.
- Recovered `alpha_state=epsilon_alpha=(3/2)chi_alpha^(5/3)` during H writeback.
- Added H-stage write evidence for `alpha_state`.
- Added diagnostics proving `passive_epsilon_alpha_advection=false`.
- Verified the controlled compression/expansion oracle:
  `epsilon_alpha_new/epsilon_alpha_old=(rho_new/rho_old)^(5/3)`.
- Verified the static H path carries alpha through the same scalar bundle plumbing.

## Not Implemented

- Production `H/T/E/R/A` registration.
- Fuel depletion.
- Alpha benchmark parity.
- Multigroup alpha.
- Alpha momentum-force feedback beyond scalar advection.

## Required Verification

```text
dec3d_contract_alpha_hydro_terms
dec3d_contract_hydro_view
dec3d_contract_ppm_reconstruction
dec3d_contract_hydro_macro_ale_single_rank
dec3d_contract_alpha_contract
dec3d_contract_alpha_one_group_operator
dec3d_contract_alpha_distributed_one_group_operator
```

Full clean-first default and HYPRE regressions must pass before P4-3 closure.
