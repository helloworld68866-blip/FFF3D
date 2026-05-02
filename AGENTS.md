# Global Instructions

On this machine, `rg` and `rg --files` are not usable and must not be used.

When searching for files or file contents, use PowerShell native commands instead.

Preferred alternatives:
- File search: `Get-ChildItem -Recurse -File`
- Content search: `Get-ChildItem -Recurse -File | Select-String -Pattern '...'`

Python on this machine:
- Use `C:\Users\Administrator\anaconda3\envs\spyder55-pip\python.exe` for Python scripts, plotting, and post-processing.
- In PowerShell, invoke it as `& 'C:\Users\Administrator\anaconda3\envs\spyder55-pip\python.exe' ...`.
- Do not use `C:\Users\Administrator\AppData\Local\Microsoft\WindowsApps\python.exe`; it is the Microsoft Store alias/stub and may fail with "Python was not found".

For build and test verification in this repository, always prefer a clean rebuild first:
- Use `-CleanFirst` before verification runs.
- Do not rely on incremental build results when validating changes, because stale build artifacts have already caused false failures and misleading behavior in this project.

# AGENTS.md

## Project
This repo is a **from-scratch DEC3D implementation**.

## Source of truth
Follow docs in this order:

"F:\dec3d\docs\Woo博士论文.md"
"F:\dec3d\docs\DEC3D物理方程与常数说明.md"
"F:\dec3d\docs\notes.md"

If code and docs disagree, **docs win**.

## Hard rules
1. **No legacy delegation**
   - Do not import, wrap, or call previous hydro/radiation implementations.
   - A backend shell that forwards into old code is forbidden.
   
2. **No fake completion**
   - A task is not done unless the runtime actually executes the new implementation.
   
3. **Missing diagnostics must fail**
   - Required metrics may not default to zero or silently disappear.
   
4. **One authoritative state**
   - Views/workspaces/ghost buffers are allowed.
   - A second authoritative simulation state is not allowed.
   
5. **No phase skipping**
   - Only implement the active phase.
   - Do not add future-phase physics early.
   
6. **No fake bug fixing**
   - 以通用的方法修复bug.
   
   - 严禁做case specific or case-aware的bug修复
7. **同一个 build 目录里多次单 target --clean-first 会反复清掉别的测试 exe，所以建议改成一次 clean-first 构建多个 focused targets，避免互相清理。**

## Development style
Work in **small vertical slices**:
1. write/update failing test
2. implement minimum code
3. run relevant tests
4. prove runtime uses the new path
5. then mark done

For simple, low-risk tasks such as documentation edits, plotting-only changes, file organization, command output inspection, or small configuration tweaks, do not force the red-test/TDD workflow; handle the task directly and use only appropriately lightweight verification.

## Efficiency optimization rules
1. Use the `128x32x32` case as the primary efficiency optimization benchmark unless the user explicitly changes the benchmark.
2. Keep the benchmark protocol fixed when comparing performance: same MPI rank count, Release build, input profile, radiation groups, output settings, step count, machine state, and command line.
3. Separate initialization cost from per-step cost. Track `init_wall_s`, first-step wall time, and steady per-step wall time independently.
4. Optimize the most expensive stage first. Within a stage, start from the most expensive measured sub-step.
5. Optimize one bottleneck at a time. Do not mix unrelated hydro, thermal, radiation, alpha, I/O, or initialization changes in the same optimization slice unless the coupling is proven necessary.
6. Before changing a target sub-step, write down its theoretical efficiency limit. Analyze the limit in three layers:
   - algorithmic lower bound: required cell visits, stencil operations, reductions, halo exchanges, solver iterations, and I/O;
   - hardware lower bound: memory bandwidth, FLOP throughput, MPI latency/bandwidth, and cache behavior;
   - engineering lower bound: best likely result without large unrelated rewrites.
7. For a numeric operation such as `c = a + b`, estimate the ideal operation count and memory traffic first, then explain why the current implementation misses that limit before optimizing it.
8. Classify every bottleneck as compute-bound, memory-bound, communication-bound, synchronization-bound, solver-bound, allocation-bound, or I/O-bound before choosing an optimization.
9. Every optimization must pass a numerical consistency gate. Compare key diagnostics such as density, pressure, velocity, electron/ion temperature, total energy where available, history profiles, and relevant norms. Efficiency changes must not change physics except for justified floating-point roundoff.
10. Use repeated measurements for timing-sensitive conclusions. Prefer at least two to three short probes and report mean and spread when practical.
11. Record Amdahl's-law impact before optimizing a sub-step: sub-step fraction of total time, best possible local speedup, and maximum possible end-to-end speedup.
12. Use clear stop conditions. Stop optimizing a sub-step when it is no longer a dominant cost, additional speedup is below measurement noise or about 5 percent end-to-end, the implementation is within roughly 2x of the justified practical limit, or further work requires a larger approved redesign.
13. Document every optimization in `optimize.md`, including target step, baseline time, theoretical limit, bottleneck diagnosis, implementation technique, changed files, verification command, new timing, speedup, numerical differences, and next bottleneck.
14. Do not make case-aware optimizations. The `128x32x32` case is a benchmark, not a license to special-case a grid, profile, mode, output directory, or input deck.

## Forbidden shortcuts
- delegating execution to old code
- placeholder kernels presented as final
- default-zero diagnostics for required metrics
- silently relaxing thresholds

## Final principle
This repo succeeds only by delivering a **real executable numerical core** with **real diagnostics** and **real acceptance evidence**.
