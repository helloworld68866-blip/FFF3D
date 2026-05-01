# P3-1b Marshak Vacuum Boundary Acceptance

Date: 2026-04-27

## Scope

P3-1b implements thesis Marshak/P1 vacuum radiation boundary support as an
explicit outer-radial Robin boundary in the generic implicit diffusion
substrate. One-group gray radiation can opt into this boundary model.

## Runtime Evidence Required

Required diagnostics include:

- `outer_boundary_policy=radiation_marshak_vacuum`
- `marshak_boundary_used=true`
- `marshak_nonzero_diagonal_loss_count>0`
- `radiation_boundary_model=thesis_marshak_vacuum`
- `marshak_enabled=true`
- `parity_claim_allowed=false`

## Numerical Contract

The Marshak boundary adds a diagonal leakage term:

```text
F_rad dot n = c U / 2
G_face = (h * (D_face / d)) / (h + (D_face / d))
h = c / 2
d = r_outer_face - r_outer_cell_center
diag += A_face * G_face
```

It is not implemented as a ghost-state shortcut, post-solve subtraction, or a
cell-centered source.

## Acceptance Tests

Focused tests must cover:

- generic CSR diagonal loss for Marshak;
- one-cell Marshak loss oracle;
- one-cell Marshak plus absorption/emission oracle;
- `dt_s=0 + Marshak` no-change behavior;
- Marshak only allowed on true outer radial boundary;
- default zero-flux thermal/generic behavior unchanged;
- one-group gray radiation reports Marshak enabled when requested;
- missing Marshak diagnostics fail validation.

## Non-Claims

This slice does not claim:

- Woo 5.4 radiation benchmark parity;
- LILAC parity;
- multigroup radiation;
- electron radiation-matter writeback;
- radiation flux limiter;
- radiation advection or `Pg div(v)` work;
- production `H/T/E/R` runtime registration.

