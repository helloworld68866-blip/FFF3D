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
Do not patch, wrap, or delegate to any older codebase.

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

## Completion rule
A task is complete only if:
- it matches current docs
- tests fail first and pass after
- runtime executes the new path
- required diagnostics are emitted
- no forbidden shortcut was used

## Forbidden shortcuts
- delegating execution to old code
- placeholder kernels presented as final
- fake HLLC / fake PPM
- copying velocity into momentum buffers and calling it conservative
- default-zero diagnostics for required metrics
- silently relaxing thresholds

## Phase policy
Phases:
- `P0` foundation
- `P1` hydro
- `P2` thermal + e-i
- `P3` radiation
- `P4` alpha
- `P5` integrated DEC3D

Do not move to the next phase until the current phase acceptance gates are green.

## Final principle
This repo succeeds only by delivering a **real executable numerical core** with **real diagnostics** and **real acceptance evidence**.
