# Phase P0 Acceptance

Document version: v0.1  
Last updated: 2026-04-22

## Goal

`P0` establishes the foundation contracts and runtime substrate. It does not claim completion of any real hydro, thermal, radiation, or alpha physics solve.

## Required Stages

- no physics stage is required to be registered in `P0`

## Authoritative Write-Set

`P0` may validate state creation and runtime plumbing, but it must not claim a physics-stage write-set.

Allowed writes are limited to:

- canonical state initialization
- diagnostics scaffolding outputs
- runtime metadata

Checkpoint and restart in `P0` must treat authoritative state as the only restart truth. Recovered or cached fields may be omitted from checkpoints or written only as invalid-on-restart debug payloads.

## Acceptance Artifacts

- clean rebuild report
- failing-first evidence for foundation contract tests
- runtime substrate report
- checkpoint and restart contract report
- diagnostics completeness report
- assumption-ledger delta

## Failure Conditions

- any physics stage registered as a placeholder success
- missing diagnostics for runtime substrate checks
- hidden delegation to legacy hydro or legacy backends
- inability to prove that new runtime plumbing executed
- recovered or cached fields treated as authoritative restart truth
