# Phase P1 Hydro Acceptance

Document version: v0.7  
Last updated: 2026-04-23  
Revision reason: aligned pre-PPM source acceptance with the implemented old-time explicit geometric source step, retained the theta/pole closeout gate, clarified that full directional static-grid hydro alone does not yet satisfy thesis-level spherical hydro, and added the closed single-rank/full-MPI PPM artifact requirements now present in the repository.

## Goal

`P1` completes the explicit hydro core at thesis level, including radial ALE, real macro-zoning coupling, and auditable visualization of a true 3D hydro-only case.

## Required Stages

- `H`

Stages `T`, `E`, `R`, and `A` must remain unregistered in `P1`.

## Authoritative Write-Set

`H` must write:

- `rho`
- `mom_r`
- `mom_theta`
- `mom_phi`
- `E_fluid_total`
- `E_electron`

`H` must not write:

- `U_rad[g]`
- `AuthoritativeAlphaState`

## Acceptance Artifacts

- clean rebuild report
- failing-first hydro test evidence
- hydro stage report with execution evidence
- runtime evidence showing executed `r/theta/phi` directional hydro updates
- runtime evidence showing executed spherical old-time geometric source step
- theta/pole boundary contract evidence before PPM closeout
- ALE mesh proposal and runtime commit evidence
- macro-zoning coarse/fine coupling evidence
- radial MPI parity and seam evidence
- real multi-rank radial-wave parity/seam analysis outputs under `analysis/output/case2_radial_wave_mpi_real/`
- real multi-rank radial-wave budget outputs under `analysis/output/case2_radial_wave_mpi_real/` including `single_rank_budget_residual.txt`, `multi_rank_budget_residual.txt`, and `budget_residual_comparison.txt`
- real multi-rank radial-wave full-PPM outputs under `analysis/output/case2_radial_wave_mpi_real_ppm/` showing `ng=3` ghost-consuming radial high-order parity/seam sanity and budget comparison against the single-rank full-PPM baseline
- real multi-rank true-3D low-mode full-PPM outputs under `analysis/output/case3_true3d_low_mode_mpi_real_ppm/` showing `ng=3` ghost-consuming full-direction high-order parity sanity and budget comparison against the single-rank full-PPM baseline
- those full-PPM analysis directories must record `ppm_executed=true`, `reconstruction_ghost_layers=3`, and an auditable reconstruction mode in `summary.txt`, `mpi_parity_diagnostics.txt`, and `case_manifest.txt`
- hydro numerical check summary
- budget / conservation residual summary for `rho` and `E_fluid_total`
- minimal true-3D static-grid hydro check before ALE and macro-zoning closeout
- true 3D hydro case manifest
- `case3_true3d_low_mode/budget_residual.txt` and `case4_pole_seam_stress/budget_residual.txt`
- real multi-rank `case3_true3d_low_mode_mpi_real/` and `case4_pole_seam_stress_mpi_real/` outputs including `mpi_parity_diagnostics.txt` and `budget_residual_comparison.txt`
- runnable case decks for hydro-only `P1` cases
- `txt` field outputs plus Python-generated plot artifacts for `rho`, recovered hydro-only `Te`, and velocity field or velocity magnitude under `analysis/output/<case-folder>/`
- authoritative field update summary
- `Workstream 6a` proposal/commit contract evidence showing that hydro may propose radial mesh motion, runtime owns canonical mesh commit, and missing commit is surfaced through auditable runtime diagnostics
- `Workstream 6b` radial ALE evidence showing that the active radial hydro path consumes `F_r^{ALE} = F_r^{HLLC} - w_{face} U_*`, selects the moving-interface HLLC branch against `w_face`, updates the ALE path in conservative extensive form `(VQ)`, and emits auditable ALE flux diagnostics
- `Workstream 7a` macro-zoning substrate evidence showing thesis-style `r Δtheta` / `r sin(theta) Δphi` detection, conservative transfer over the hydro package, and explicit rejection of coarse buffers as authoritative truth
- `Workstream 7b` single-rank macro-zoning evidence showing the active hydro path consumes coarse/fine coupling on angular sweeps and emits auditable macro-zoning diagnostics
- `Workstream 7c` real-MPI macro-zoning evidence showing `case2` radial-wave parity / seam sanity, budget comparison, and transaction safety with macro-zoning enabled
- `Workstream 7c` real-MPI macro-zoning evidence showing `case3` true-3D low-mode parity / budget artifacts with macro-zoning enabled
- `Workstream 7c` real-MPI macro-zoning evidence showing `case4` pole / seam stress parity, budget artifacts, and explicit angular-momentum seam diagnostics with macro-zoning enabled
- `artifacts/p1/*.md` generated through the `P1` acceptance harness, with `p1-closeout-summary.md` explicitly reporting `workstream7_complete=true` and `full_p1_closeout_ready=false`
- assumption-ledger delta

## Failure Conditions

- `H` not registered
- `H` registered but unable to prove real execution
- static-grid-only completion presented as full `P1`
- radial-only static-grid completion presented as full directional hydro
- full directional static-grid hydro presented as thesis-level spherical hydro before geometric source step closeout
- direct theta boundary physical-flux substitution presented as a valid pole-ready angular boundary contract
- single-rank-only completion presented as full `P1` when radial-MPI acceptance is still missing
- real MPI seam evidence bypassing later integrated runtime/MPI closeout presented as full `P1` completion
- ALE implemented outside runtime mesh-commit rules
- proposal/commit-only ALE presented as numerically complete `P1b` without proof that the active radial flux path consumes the thesis-style ALE correction with `w_face`-based branch selection and extensive conservation
- macro-zoning substrate-only completion presented as full coarse/fine hydro coupling
- single-rank first-order macro-zoning coupling presented as full `P1c` closeout before MPI seam/parity and later high-order macro-zoning closure
- `case2` real-MPI macro-zoning parity presented as full `P1c` closeout before broader true-3D macro-on MPI checks and later higher-order macro-zoning closure
- `case3` real-MPI true-3D macro-zoning parity presented as full `P1c` closeout before broader macro-on MPI coverage and later higher-order macro-zoning closure
- `case4` real-MPI pole/seam macro-zoning parity presented as full `P1c` closeout before broader macro-on/off and later higher-order macro-zoning closure
- macro-zoning left as interface-only scaffolding
- no auditable true-3D hydro visualization outputs
- unauthorized writes to radiation or alpha authoritative fields
