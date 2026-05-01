# P2-6a Scalar Origin/Pole Contract Ledger

Date: 2026-04-27

## Current Code Reality

The current `generic_diffusion` full-sphere scalar origin/pole path reports:

- `origin_remap_coupling_mode=metric_regular_origin`
- `pole_remap_coupling_mode=metric_regular_axis`

The current contract test also asserts that the inspected origin/pole opposite-side
CSR entries are zero, not finite nonzero mapped-neighbor off-diagonal entries.

Evidence:

- `src/transport/diffusion/generic_diffusion.cpp`
  - `AssembleScalarRemapConductances(...)` is currently a no-op.
  - diagnostics report `metric_regular_origin` / `metric_regular_axis`.
- `tests/contract/test_generic_implicit_diffusion.cpp`
  - origin/pole inspected CSR values are expected to be `0.0`.

## Drift From Earlier P2-6a Spec

The original P2-6a design/implementation plan described a different acceptance
contract:

- full-sphere scalar origin/pole remap should create real mapped-neighbor CSR
  coupling;
- mapped origin/pole pairs should have finite nonzero off-diagonal evidence.

That wording no longer matches the current implementation.

## Current Acceptance Wording

Until this is explicitly redesigned, do not describe current P2-6a as:

- `mapped-neighbor nonzero CSR coupling`
- `finite nonzero origin/pole mapped off-diagonal closeout`

Use the current implementation wording instead:

- metric-regular scalar origin/pole closure in generic diffusion;
- no artificial half-turn / opposite-phi CSR edge is added at origin or pole;
- full-sphere scalar boundary behavior is accepted by the current tests under
  the metric-regular closure semantics.

## P3-1b Relationship

This ledger does not block P3-1b.

P3-1b concerns the true outer radial physical boundary and adds a Marshak/P1
vacuum Robin diagonal loss there. It does not rely on origin/pole mapped
off-diagonal coupling.

## Follow-Up Decision

Before claiming thesis-final full-sphere scalar diffusion stencil parity, choose
one explicit contract and update spec/tests accordingly:

1. Keep the current metric-regular closure and revise old P2-6a/P2-6b wording.
2. Restore a mapped-neighbor nonzero CSR coupling design and rework tests/code.

No benchmark or production parity claim should depend on the old nonzero
mapped-neighbor wording until this decision is made.
