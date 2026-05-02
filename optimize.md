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
