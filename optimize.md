# DEC3D Optimization Results

Record each completed optimization here.

## Entries

### 2026-05-02 - Baseline, 128x32x32, mpi24, all stages

Case: `artificial_profile_accel2e17_tauE100_p2m0p20_noale_baseline50_128x32x32_mpi24_allstages`

Executable: `F:\dec3d\build-hypre\Release\dec3d.exe`

Input: `F:\dec3d\cases\artificial_profile_accel2e17_tauE100_p2m0p20_noale_baseline50_128x32x32_mpi24_allstages.in`

Profile: `F:\dec3d\cases\artificial_profile_accel2e17_tauE100_1024.pro`

Output: `F:\dec3d\analysis\output\artificial_profile_accel2e17_tauE100_p2m0p20_noale_baseline50_128x32x32_mpi24_allstages`

Result:
- Exit code: 0.
- Steps: 50.
- End physical time: 0.078723 ps; final dt: 0.001574 ps.
- Time-loop wall: 45.308 s; efficiency: 0.906 s/step.
- End-to-end wall: 53.403 s; efficiency: 1.068 s/step.
- Physical throughput: 6.26 ps/hour by time-loop wall, 5.31 ps/hour end-to-end.
- 100 ps estimate: 16.0 h by time-loop wall, 18.8 h end-to-end.

Stage timing:

| Stage | Wall s | Share |
|---|---:|---:|
| R radiation | 36.980 | 81.6% |
| H hydro | 4.124 | 9.1% |
| T thermal | 2.115 | 4.7% |
| A alpha | 1.170 | 2.6% |
| E equilibration | 0.701 | 1.5% |

Current bottleneck:
- First target: `R` stage, especially `r_hypre_setup_wall_s = 16.495 s` (36.4% total, 44.6% of R).
- Other major R costs: `r_coefficient_provider_wall_s = 6.255 s`, `r_hypre_solve_wall_s = 3.777 s`, `r_assembly_wall_s = 2.759 s`, `r_flux_limiter_wall_s = 0.563 s`.

### 2026-05-02 - Optimization 1, BoomerAMG 3D strong threshold

Target: `r_hypre_setup_wall_s`.

Theoretical limit:
- The 50-step benchmark does 500 radiation HYPRE setups: 50 steps times 10 groups.
- With the current GMRES+BoomerAMG backend and group-by-group matrices, the setup cannot be removed without reusing/changing the preconditioner. The local lower bound is one AMG hierarchy construction per group matrix, over about 907264 global matrix nonzeros.
- The optimization lever is therefore AMG hierarchy complexity, not the radiation discretization.

Change:
- Set `HYPRE_BoomerAMGSetStrongThreshold(..., 0.5)` in the generic HYPRE solver path.
- Kept GMRES tolerance, max iterations, matrix assembly, and radiation equations unchanged.
- Trial `max_coarse_size=256` was rejected because it worsened setup and solve time.
- Trial `strong_threshold=0.6` was rejected because the R stage failed the distributed solve gate.

Result:
- Case: `artificial_profile_accel2e17_tauE100_p2m0p20_noale_opt1_strong0p5_128x32x32_mpi24_allstages`.
- Exit code: 0.
- Time-loop wall: 45.308 s -> 44.643 s, improvement 0.665 s (1.5%).
- End-to-end wall: 53.403 s -> 51.174 s, improvement 2.229 s (4.2%).
- `R` stage: 36.980 s -> 36.425 s, improvement 0.556 s (1.5%).
- `r_hypre_setup_wall_s`: 16.495 s -> 15.322 s, improvement 1.173 s (7.1%).
- `r_hypre_solve_wall_s`: 3.777 s -> 4.374 s, regression 0.598 s.
- `r_solver_iterations`: 2839 -> 2888.
- Time-loop throughput: 6.26 ps/hour -> 6.35 ps/hour.
- 100 ps estimate by time-loop wall: 16.0 h -> 15.8 h.

Verification:
- Clean Release rebuild of `dec3d`.
- `dec3d_contract_hypre_diffusion_solver`: passed.
- `dec3d_contract_hypre_distributed_diffusion_solver` with `mpiexec -n 2`: passed.
- `dec3d_contract_radiation_distributed_provider_fed_multigroup` with `mpiexec -n 2`: passed.
- Final benchmark summaries stayed close to baseline; largest final scalar-summary difference was `rho_max` by about `8.4e-7`.

### 2026-05-02 - Optimization 2, batch HYPRE IJ matrix insertion

Target: low-risk HYPRE input assembly overhead around the radiation solve.

Theoretical limit:
- The matrix entries are already assembled in local CSR form. The previous path repacked and inserted each row separately, causing one HYPRE insertion call and two temporary vectors per local row.
- The low-risk lower bound is one linear CSR-to-HYPRE conversion plus one batched `HYPRE_IJMatrixSetValues` call per local matrix.

Change:
- Replaced row-by-row `HYPRE_IJMatrixSetValues` calls with one batched local CSR insertion.
- Matrix entries, row order, solver settings, and convergence checks are unchanged.

Result:
- Case: `artificial_profile_accel2e17_tauE100_p2m0p20_noale_opt2_batchsetvalues_128x32x32_mpi24_allstages`.
- Exit code: 0.
- Time-loop wall: 44.643 s -> 44.363 s, improvement 0.280 s (0.6%) from Optimization 1.
- `R` stage: 36.425 s -> 36.309 s, improvement 0.116 s.
- Time-loop throughput: 6.35 ps/hour -> 6.39 ps/hour.
- 100 ps estimate by time-loop wall: 15.75 h -> 15.65 h.

Verification:
- Clean Release rebuild of `dec3d`.
- `dec3d_contract_hypre_diffusion_solver`: passed.
- `dec3d_contract_hypre_distributed_diffusion_solver` with `mpiexec -n 2`: passed.
- `dec3d_contract_radiation_distributed_provider_fed_multigroup` with `mpiexec -n 2`: passed.
- Final benchmark scalar summaries matched Optimization 1 at printed precision.

### 2026-05-02 - Optimization 3, low-risk hot-path lookup and storage access cleanup

Target: `R` stage coefficient provider and repeated regular-grid storage access, following `docs/优化点.md`.

Theoretical limit:
- The TOPS provider must still visit each local cell for each radiation group and perform the same opacity interpolation and blackbody evaluation.
- The low-risk lower bound removes avoidable overhead around that required work: no bounds-checked `Array3D` access in Release hot loops, linear storage sweeps for regular grids, and average O(1) TOPS table row lookup instead of tree lookup.

Change:
- `Array3D::operator()` now uses unchecked storage access in Release, with explicit `checked()` preserved for tests/debug access.
- Radiation provider and distributed radiation writeback use linear `Array3D::storage()` sweeps where the layout is already validated.
- TOPS opacity table row index changed from `std::map` to reserved `std::unordered_map`; table values and interpolation logic are unchanged.

Result:
- Case: `artificial_profile_accel2e17_tauE100_p2m0p20_noale_opt5_tops_hash_lookup_128x32x32_mpi24_allstages`.
- Exit code: 0 for two 50-step probes.
- Time-loop wall: 44.363 s -> 42.60 s and 42.39 s, mean 42.50 s; improvement 1.87 s (4.2%) from Optimization 2.
- End-to-end wall: 50.921 s -> 49.341 s and 48.714 s, mean 49.03 s; improvement 1.89 s (3.7%) from Optimization 2.
- Mean `R` stage across ranks: 36.257 s -> 34.744 s and 34.559 s, mean 34.652 s.
- Mean `r_coefficient_provider_wall_s` across ranks: 6.328 s -> 5.629 s and 5.574 s, mean 5.601 s.
- Final benchmark scalar summaries matched Optimization 2 at printed precision.

Verification:
- Clean Release rebuild of `dec3d`.
- `dec3d_unit_array3d`: passed.
- `dec3d_contract_radiation_tops_opacity_provider`: passed.
- `dec3d_contract_radiation_distributed_provider_fed_multigroup` with `mpiexec -n 2`: passed.
- `dec3d_contract_hypre_distributed_diffusion_solver` with `mpiexec -n 2`: passed.

### 2026-05-02 - Optimization 4, lagged BoomerAMG setup reuse

Target: `r_hypre_setup_wall_s`.

Theoretical limit:
- The 50-step benchmark does 500 radiation group solves. With full rebuild, each solve performs one BoomerAMG setup.
- If the AMG hierarchy could be reused for every 4 solves without extra iterations, the local setup lower bound would be roughly 25% of the previous setup cost.
- Practical lower bound is higher because the matrix-change gate and first solve per group still require rebuilds, and a reused preconditioner can increase solve time.

Change:
- Added a per-radiation-group lagged BoomerAMG cache in the distributed runtime.
- Reuse is only attempted when the matrix shape matches, the max relative matrix change against the last AMG-setup matrix is at most `0.1`, and the cache has been reused fewer than `rebuild_every - 1` times.
- Current runtime setting: `rebuild_every=4`, `max_matrix_relative_change=0.1`, `max_iteration_growth=1.5`.
- Reuse still runs GMRES setup/solve with the cached BoomerAMG hierarchy, then accepts only if the recomputed residual gate passes and iterations do not exceed the growth gate.
- Rejected reuse falls back to a fresh BoomerAMG setup. No radiation equation, matrix assembly, RHS, or convergence tolerance was changed.

Result:
- Case: `artificial_profile_accel2e17_tauE100_p2m0p20_noale_opt6_lagged_amg_128x32x32_mpi24_allstages`.
- Exit code: 0 for two 50-step probes; stderr empty.
- Time-loop wall: Optimization 3 mean `42.50 s` -> `39.29 s` and `40.31 s`, mean `39.80 s`; improvement `2.70 s` (6.4%).
- End-to-end wall: Optimization 3 mean `49.03 s` -> `46.01 s` and `47.12 s`, mean `46.57 s`; improvement `2.46 s` (5.0%).
- Mean `R` stage: Optimization 3 mean `34.652 s` -> `31.319 s` and `32.093 s`, mean `31.706 s`; improvement `2.946 s` (8.5%).
- `r_hypre_setup_wall_s`: saved Optimization 3 probe `15.241 s` -> `11.532 s` and `11.788 s`, mean `11.660 s`; improvement about 23.5%.
- `r_hypre_solve_wall_s`: saved Optimization 3 probe `4.205 s` -> `4.431 s` and `4.570 s`, mean `4.501 s`; expected regression from reusing older preconditioners.
- `r_solver_iterations`: stayed at `2888`.
- Lagged AMG counts per 50-step run: `169` accepted reuses, `331` rebuilds, `0` fallback rebuilds.

Numerical consistency:
- Compared `dec3d.out` against Optimization 3 over all 50 history rows.
- Max relative differences: `rho_max=3.7e-11`, `Te_min_keV=1.9e-8`, `Te_max_keV=5.0e-11`, `Ti_min_keV=1.7e-9`, `Ti_max_keV=1.4e-11`.
- Final scalar summaries are unchanged to the intended benchmark precision.

Verification:
- Clean Release rebuild of `dec3d`, `dec3d_contract_hypre_diffusion_solver`, `dec3d_contract_hypre_distributed_diffusion_solver`, and `dec3d_contract_radiation_distributed_provider_fed_multigroup`: passed.
- `dec3d_contract_hypre_diffusion_solver`: passed.
- `dec3d_contract_hypre_distributed_diffusion_solver` with `mpiexec -n 2`: passed.
- `dec3d_contract_radiation_distributed_provider_fed_multigroup` with `mpiexec -n 2`: passed.
