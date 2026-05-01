# DEC3D

This repository is the clean-room reconstruction of the DEC3D program described in Woo's thesis, implemented in `C++`, `MPI`, and `Hypre`.

The current implementation target is `P0`: foundation contracts, runtime substrate, authoritative-state scaffolding, and failure-first test harnesses.

## Local Build

On this machine, the most reliable build path is Visual Studio 2022 via CMake:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Debug --clean-first
ctest --test-dir build -C Debug --output-on-failure
```

## Phase References

Project design and acceptance rules live in:

- `docs/architecture/dec3d-system-design.md`
- `docs/architecture/phase-p0-implementation-plan.md`
- `docs/acceptance/dec3d-gates.md`
- `docs/acceptance/phase-p0.md`
