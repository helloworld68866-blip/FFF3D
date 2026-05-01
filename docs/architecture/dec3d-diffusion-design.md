# DEC3D Diffusion And Implicit Solver Design

Document version: v0.1  
Last updated: 2026-04-22  
Primary phases: `P2` through `P4`

## Scope

This document defines the shared numerical substrate used by the implicit operators:

- thermal conduction
- electron-ion equilibration support where matrix infrastructure is needed
- multigroup radiation transport
- alpha transport and deposition support where matrix infrastructure is needed

The design goal is shared numerical infrastructure without shared hidden physics semantics.

## Layer Boundary

`src/diffusion/` is a numerical foundation layer parallel to `mesh/state`, not a replacement for the physics modules.

It may provide:

- sparse matrix helpers
- vector assembly helpers
- boundary-condition support utilities
- Hypre solver integration
- reporting hooks for matrix and solve evidence

It may not provide:

- a generic physics backend that decides operator meaning
- a hidden legacy shell
- a substitute for explicit operator-owned assembly logic

## Ownership Rules

Each implicit physics module owns:

- its governing discrete operator
- its coefficients
- its source terms
- its authoritative write-set

`diffusion infrastructure` owns:

- reusable sparse assembly mechanics
- reusable Hypre bridge code
- reusable solve configuration surfaces

This keeps the physics meaning in the operator and the solver mechanics in the infrastructure.

## Hypre Contract

The Hypre integration layer must expose enough structure for runtime evidence and gate checks.

Every implicit solve path must be able to report:

- matrix build occurred
- right-hand side build occurred
- solver invocation occurred
- solve status
- iteration count or equivalent convergence evidence

Missing solve diagnostics are gate failures.

## Module-Specific Expectations

### Thermal

Thermal writes:

- `E_electron`
- `E_fluid_total`

Thermal recovers:

- `Te`
- `Ti`
- conductivity

Thermal must not mutate radiation-group or alpha authoritative state.

### Equilibration

Equilibration transfers energy between electron and ion thermal reservoirs.

Equilibration writes:

- `E_electron`

Equilibration must preserve `E_fluid_total` up to acceptable floating-point tolerance.

### Radiation

Radiation owns multigroup transport semantics.

Radiation writes:

- `U_rad[g]`
- `E_electron`
- `E_fluid_total` where coupling requires it

Radiation completion requires at least one group to show real participation in state update.

### Alpha

Alpha owns alpha transport and deposition semantics.

Alpha writes:

- `AuthoritativeAlphaState`
- `E_electron`
- `E_fluid_total` where deposition coupling requires it

Alpha may share numerical infrastructure with other implicit modules, but it may not be folded into thermal or radiation semantics.

## Runtime Interaction

Implicit modules must use the same runtime-facing contract as hydro:

- `bind(...)`
- `estimate_dt(...) -> DtAdvice`
- `advance(...) -> StageResult`

Implicit modules may return soft timestep advice, but they do not own the final global timestep.

## Assumption Policy

Any thesis-external decision about solver tolerances, preconditioning choices, group closures, or coefficient regularization must be recorded in `docs/architecture/assumption-ledger.md`.
