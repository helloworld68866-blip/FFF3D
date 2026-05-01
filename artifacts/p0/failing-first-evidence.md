# P0 Failing-First Evidence

Date: 2026-04-22
Phase: `P0`
Status: recorded

## Representative Red-To-Green Evidence

- `dec3d_contract_stage_result`
  Initial failure: `result.is_semantically_complete()` was false because diagnostics and execution evidence were incomplete.
  Green condition: `StageResult` semantic-completeness contract implemented and test passed.

- `dec3d_contract_cache_invalidation`
  Initial failure: cache invalidation mask did not invalidate dependent temperature fields.
  Green condition: authoritative write-set invalidation mapping implemented and test passed.

- `dec3d_contract_spherical_geometry`
  Initial failure: spherical face monotonicity / geometry contract was incomplete.
  Green condition: geometry metadata and positive-volume contract implemented and test passed.

- `dec3d_contract_radial_ownership`
  Initial failure: ownership slice initialization / continuity contract failed.
  Green condition: contiguous radial ownership decomposition implemented and test passed.

- `dec3d_contract_runtime_stage_result_ingestion`
  Initial failure: unregistered, forbidden, duplicate, or semantically incomplete stage results were not fully rejected.
  Green condition: phase-aware `TryIngestStageResult(...)` contract implemented and test passed.

- `dec3d_contract_checkpoint_contract`
  Initial failure 1: clean rebuild failed because `src/runtime/checkpoint_contract.cpp` did not exist.
  Initial failure 2: broken runtime snapshot still allowed a checkpoint payload to appear complete.
  Green condition: authoritative-only checkpoint payload, restart validation, runtime-snapshot fail path, and restart invalidation digest implemented and test passed.

- `dec3d_acceptance_phase_p0`
  Initial failure 1: clean rebuild failed because `src/runtime/p0_acceptance.cpp` did not exist.
  Initial failure 2: acceptance summary did not expose `runtime_ingested_stage_mask`, `runtime_missing_stage_report_mask`, or `runtime_authoritative_write_set_digest`.
  Green condition: `P0AcceptanceSummary` closeout contract implemented and test passed.

## Passing Rerun

```powershell
ctest --test-dir build -C Debug --output-on-failure
```

Result: `20/20` tests passed.
