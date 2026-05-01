# DEC3D Project Gate Rules

Document version: v0.1  
Last updated: 2026-04-22  
Applies to: all phases

## Global Gate Philosophy

DEC3D validation is failure-first.

Green status is invalid unless all of the following are true:

- the required stage is registered for the active phase
- the new path actually executed
- diagnostics are complete
- authoritative writes conform to the phase contract
- the reported numerical checks passed

## Mandatory Gate Checks

Every phase gate must check:

- clean rebuild verification using `-CleanFirst`
- recorded failing-first evidence
- runtime execution evidence
- authoritative write-set conformance
- numerical acceptance checks
- complete diagnostics
- thermodynamic feasibility after commit
- assumption-ledger delta presence

## Registration Rules

For any stage required by the active phase contract, missing registration is a gate failure.

For any stage outside the active phase contract:

- the stage must be absent
- placeholder registration is forbidden
- placeholder success is forbidden

## Diagnostics Rules

The following are immediate gate failures:

- missing diagnostics
- default-zero diagnostics that hide missing data
- NaNs
- stale-cache reads after authoritative writes
- unproven runtime-path execution
- authoritative write-set violations
- infeasible recovered thermodynamic states

## Failing-First Evidence

Each new phase capability must preserve failing-first evidence as an auditable artifact:

- failing test case identifier
- failing log or `ctest` summary
- fix patch or equivalent implementation identifier
- passing rerun result

Claims without retained evidence do not satisfy the phase gate.

## Acceptance Artifacts

Every phase must retain or generate:

- build report
- test report
- stage reports
- field update summary
- numerical check summary
- assumption-ledger delta

## Change Control

Any modification to this document or to any `phase-p*.md` file that changes:

- done definition
- thresholds
- required stages
- authoritative write-set contract
- assumption-ledger status rules

must record a document version change and a reason.
