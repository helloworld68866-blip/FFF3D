# DEC3D Hydro Design

Document version: v0.1  
Last updated: 2026-04-22  
Primary phase: `P1`

## Scope

`P1` delivers the complete explicit hydro core described at thesis level. It is not limited to static-grid Eulerian hydro.

`P1` includes:

- conservative fluid state update
- spherical-coordinate geometric terms
- reconstruction
- Riemann solve
- flux divergence
- timestep estimation
- radial ALE
- macro-zoning with real coarse/fine coupling
- radial MPI hydro acceptance

If implementation pacing requires smaller internal milestones, they remain subparts of the same phase:

- `P1a`: static-grid hydro core
- `P1b`: radial ALE
- `P1c`: macro-zoning and radial MPI hydro acceptance

## Canonical Inputs And Outputs

At stage entry, hydro reads authoritative fine-grid state from `canonical_state`.

Hydro constructs the local conservative package:

- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `chi_e`

`chi_e` is recovered from the canonical electron thermodynamic state through the thesis-style hydro view of electron pressure. It is never stored as global authoritative truth.

At stage exit, hydro writes back:

- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `E_electron`, through `chi_e -> Pe -> E_electron`

Hydro must not write:

- `U_rad[g]`
- `AuthoritativeAlphaState`

## Hydro-Only Infrastructure

The following mesh features are hydro-only infrastructure even if they live under `src/mesh/`:

- radial ALE
- macro-zoning

Implicit modules may not reuse them as generic mesh services.

## ALE Boundary

Radial ALE is part of the hydro core, not a later integration feature.

Hydro may:

- compute an ALE mesh-motion proposal
- include ALE-dependent geometry in hydro updates
- return a `mesh_update_proposal` in its `StageResult`

Hydro may not:

- directly mutate canonical mesh ownership
- bypass runtime when committing mesh updates

## Macro-Zoning Boundary

Macro-zoning exists only to support explicit hydro and relax the hydro CFL burden. It must not become a second authoritative state.

Required rules:

- coarse representations are short-lived hydro work views
- fine-grid canonical truth remains authoritative
- coarse/fine transfers are part of hydro stage logic
- macro-zoning completion requires real coupling, not interface-only placeholders

## Acceptance Intent

`P1` is complete only when the thesis-style hydro path is demonstrably real.

That means:

- hydro runtime evidence proves the new path executed
- ALE evidence proves mesh motion was proposed and committed through runtime
- macro-zoning evidence proves real coarse/fine coupling occurred
- numerical checks prove the hydro path behaves acceptably on minimal validation problems

## Assumption Policy

Any thesis-external hydro completion, including limiter constants or hydro-specific numerical details not spelled out in the thesis, must be recorded in `docs/architecture/assumption-ledger.md`.
