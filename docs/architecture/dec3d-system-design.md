# DEC3D System Design

Document version: v0.1  
Last updated: 2026-04-22  
Applies to: `P0` through `P5`

## Goal

This project reconstructs the DEC3D program described by Woo's thesis as a new `C++/MPI/Hypre` codebase. The goal is thesis-level structural fidelity, not source-code archaeology.

That means:

- preserve the thesis-level architecture
- preserve the thesis-level physics split and runtime organization
- allow only minimal engineering completion where the thesis is silent
- record every thesis-external completion in `docs/architecture/assumption-ledger.md`

## Phase Model

The project is organized as six phases, each adding exactly one class of real capability.

- `P0`: foundation contracts and runtime substrate
- `P1`: complete explicit hydro core
- `P2`: two-temperature thermal conduction and electron-ion equilibration
- `P3`: multigroup radiation transport
- `P4`: alpha-particle transport and deposition
- `P5`: integrated `H/T/E/R/A` orchestration and regression

Each phase must satisfy the same completion contract:

- code exists for the new path
- a failing test exists first
- the new path later passes
- runtime evidence proves the new path executed
- missing diagnostics are gate failures
- validation uses a clean rebuild, not incremental artifacts

## Architecture Layers

The repository is organized into four major layers plus top-level orchestration.

### `core`

`core` contains only generic infrastructure:

- `constants`
- `units`
- `array`
- `logging`
- `diagnostics`

It may not depend on physics modules.

### `mesh/state`

`mesh` and `state` provide the physical substrate shared by operators.

`mesh` contains two kinds of capability:

- project-wide mesh infrastructure: spherical geometry, ownership, ghost topology
- hydro-only mesh infrastructure: radial ALE and macro-zoning

Hydro-only mesh infrastructure may live under `src/mesh/`, but it is not a generic service for implicit modules.

`state` owns the unique fine-grid authoritative truth. Physics modules may construct local views, but they may not own long-lived second truths.

### `diffusion infrastructure`

`src/diffusion/common`, `src/diffusion/matrix`, and `src/diffusion/hypre` form a shared numerical substrate for implicit operators.

This layer provides:

- common sparse operator utilities
- matrix assembly support
- Hypre integration

This layer does not own physics meaning and must not become a generic hidden backend that replaces explicit operator logic.

### `physics operators`

Physics is implemented as separate operators:

- `hydro`
- `thermal`
- `equilibration`
- `radiation`
- `alpha`

Operators own their own discrete updates. They do not own global time orchestration and they do not call one another directly.

### `runtime/init/io`

`runtime` is the only layer that knows the global `H/T/E/R/A` stage order and the final timestep decision. `init` constructs initial conditions. `io` handles snapshots and reports only.

## Dependency Rules

The intended dependency picture is:

- `core`
- `mesh/state` depends on `core`
- `diffusion infrastructure` depends on `core`
- `hydro` depends on `core + mesh/state`
- `thermal` depends on `core + mesh/state + diffusion infrastructure`
- `equilibration` depends on `core + state`
- `radiation` depends on `core + mesh/state + diffusion infrastructure`
- `alpha` depends on `core + mesh/state + diffusion infrastructure`
- `runtime/init/io` may depend on all lower layers

No physics operator may depend on `runtime`.

## Authoritative State Model

`canonical_state` is the unique fine-grid authoritative physical state.

Its authoritative fields are:

- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `E_electron`
- `U_rad[g]` for each radiation group
- `AuthoritativeAlphaState`

`E_fluid_total` excludes radiation-group energy and alpha-state stored energy.

Recovered or cached quantities are not authoritative:

- `Te`
- `Ti`
- `Pe`
- `Pi`
- `chi_e`
- sound speed
- conductivity
- opacity

Hydro must construct its local conservative package as:

- `[rho, mom_r, mom_theta, mom_phi, E_fluid_total, chi_e]`

`chi_e` is operator-local. It is built from the canonical electron thermodynamic state at hydro stage entry and written back through `Pe -> E_electron` at hydro stage exit.

After any authoritative write:

- all dependent caches are invalid
- the next operator must recover fresh thermodynamic views

Checkpoint and restart rules follow the same boundary:

- authoritative fields may be checkpointed as restart truth
- recovered or cached fields are not restart truth
- if recovered or cached fields are emitted for debug or inspection, restart must mark them invalid and rebuild them before the next operator executes

Canonical truth must remain thermodynamically feasible after every stage commit:

- `E_electron >= E_e_min`
- `E_ion = E_fluid_total - E_kinetic - E_electron >= E_i_min`

## Operator Contract

All physics operators must implement the same runtime-facing contract:

- `bind(context, mesh, state, diagnostics_sink)`
- `estimate_dt(stage_view) -> DtAdvice`
- `advance(stage_span, state) -> StageResult`

`DtAdvice` must distinguish:

- `hard_cap_dt`
- `soft_advice_dt`
- `reason`
- `evidence`

`StageResult` must include:

- success or failure
- authoritative write-set
- execution evidence
- diagnostics payload
- optional `mesh_update_proposal` for hydro

Diagnostics are not optional side effects. They are part of stage completion.

## Runtime Orchestration

`runtime` is the sole owner of:

- stage registration
- final timestep selection
- canonical mesh commit
- authoritative write-set tracking
- cache invalidation propagation
- stage ordering

The canonical stage order is:

- `H`
- `T`
- `E`
- `R`
- `A`

The runtime loop is:

1. build or load fine-grid canonical state
2. collect timestep advice
3. execute registered stages in order
4. validate feasibility, diagnostics, and execution evidence after each stage
5. commit canonical mesh only through runtime
6. emit reports through `io`

Hydro may propose a mesh update for radial ALE, but it may not directly commit canonical mesh changes.

Stages that are not required by the active phase must remain unregistered. Placeholder success stages are forbidden.

## Validation Philosophy

All verification is failure-first and phase-scoped.

Every phase must provide:

- recorded failing-first evidence
- clean rebuild verification
- execution evidence for every required new path
- write-set conformance
- complete diagnostics
- numerical acceptance checks
- an assumption-ledger delta

Integrated regression may validate orchestration, but it may not weaken earlier phase done definitions.

## Document Map

This system design is the project-wide authority for architecture and phase boundaries. Related documents are:

- `docs/architecture/dec3d-hydro-design.md`
- `docs/architecture/dec3d-diffusion-design.md`
- `docs/architecture/assumption-ledger.md`
- `docs/acceptance/dec3d-gates.md`
- `docs/acceptance/phase-p0.md`
- `docs/acceptance/phase-p1-hydro.md`
- `docs/acceptance/phase-p2-thermal.md`
- `docs/acceptance/phase-p3-radiation.md`
- `docs/acceptance/phase-p4-alpha.md`
- `docs/acceptance/phase-p5-integrated.md`
