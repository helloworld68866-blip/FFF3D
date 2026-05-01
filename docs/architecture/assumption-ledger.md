# DEC3D Assumption Ledger

Document version: v0.4  
Last updated: 2026-04-28  
Scope: repository-wide thesis-external engineering assumptions for the DEC3D reconstruction effort

## Purpose

This ledger records every engineering completion that is not explicitly specified by Woo's thesis or the supporting notes in this repository. No phase may close without an assumption-ledger delta, even if the delta is empty.

The governing rule is:

- thesis-external decisions may fill implementation gaps
- those decisions must not change the thesis-level architecture
- every such decision must remain auditable

## Entry Format

Every assumption entry must contain:

- `assumption_id`
- `phase`
- `module`
- `thesis_gap`
- `decision`
- `status`
- `impact_on_gates`
- `impact_on_result_interpretation`
- `owner`
- `date`

Allowed `status` values:

- `temporary`
- `promoted_to_contract`

## Review Rules

- A new assumption entry is required whenever implementation completes a thesis gap.
- Promoting an assumption into a formal contract requires an explicit document version bump in the affected design or acceptance document.
- If a phase introduces no new assumptions, that phase must still record an empty delta in its acceptance artifacts.
- Assumptions that alter result interpretation must be called out in phase reports and integrated regression summaries.

## Initial Ledger

### A-P1-001

- `assumption_id`: `A-P1-001`
- `phase`: `P1`
- `module`: `hydro/riemann`
- `thesis_gap`: the thesis locks `chi_e = P_e^{3/5}` to HLLC-consistent scalar advection but does not spell out a unique star-state construction formula for `chi_e` inside the HLLC solver
- `decision`: the current implementation scales `chi_e` in the left/right star states with the same density-ratio factor used for the HLLC star density and advects `chi_e` through the HLLC-selected interface flux
- `status`: `temporary`
- `impact_on_gates`: `P1` numerical checks may validate the current HLLC path, but final hydro closeout must keep this assumption visible rather than presenting it as thesis-directly specified
- `impact_on_result_interpretation`: current `chi_e` transport is a reasonable engineering completion for the scaffold and static-grid hydro stages, but its exact star-state form remains an auditable assumption until benchmarked or tightened
- `owner`: `Codex`
- `date`: `2026-04-22`

### A-P1-002

- `assumption_id`: `A-P1-002`
- `phase`: `P1`
- `module`: `hydro/driver`
- `thesis_gap`: the thesis text fixes the outer radial boundary to zero inflow, but the currently available repository material does not provide a line-by-line ghost-state formula for all hydro variables at the outer shell boundary in the present first-order static-grid implementation
- `decision`: the current implementation constructs the upper radial ghost state by copying density, pressure, tangential velocity, and `chi_e` from the outer interior cell and clamping inward radial velocity to zero before sending the ghost/interior pair through HLLC; the lower radial boundary remains an explicit shell-style ghost copy until origin-side 3D mapping is implemented as part of the later ghost substrate work
- `status`: `temporary`
- `impact_on_gates`: `P1a` contract and numerical checks may validate the current radial ghost-state boundary contract, but shock-tube and later ALE/MPI acceptance must continue to treat this as an auditable engineering completion rather than thesis-directly specified final boundary logic
- `impact_on_result_interpretation`: outer radial inflow suppression now follows an explicit, test-covered contract instead of a placeholder physical-flux substitution, but the exact ghost-state construction remains a temporary assumption pending the fuller hydro ghost substrate
- `owner`: `Codex`
- `date`: `2026-04-22`

### A-P1-003

- `assumption_id`: `A-P1-003`
- `phase`: `P1`
- `module`: `hydro/driver`
- `thesis_gap`: the thesis and repository design lock pole and periodic angular relationships for 3D spherical hydro, but the current pre-PPM codebase still lacks a fully specified `ng=3` source-level hydro ghost-fill API and an origin-side mapping contract for all directions
- `decision`: the current implementation introduces a temporary directional ghost substrate that fills `theta` pole ghosts through `phi + Nphi/2` remap plus `mom_theta/mom_phi` sign flips, fills `phi` ghosts by periodic copy, and fills radial ghost layers by repeating the current one-layer radial boundary contract so the first-order hydro path can execute through a unified directional ghost entry before true PPM reconstruction is added
- `status`: `temporary`
- `impact_on_gates`: `P1a` may use the new directional ghost substrate as the pre-PPM hydro boundary contract, but `PPM` closeout still requires the future `ng=3` hydro ghost fill to be treated as a reconstruction-ready substrate rather than assuming this temporary first-order shell is already the final high-order boundary system
- `impact_on_result_interpretation`: current `theta/phi` ghost preparation is sufficient to make the static-grid first-order hydro path auditable and non-placeholder, but it remains an engineering completion that must stay visible until the full PPM-ready hydro ghost substrate is implemented and benchmarked
- `owner`: `Codex`
- `date`: `2026-04-22`

### A-P1-004

- `assumption_id`: `A-P1-004`
- `phase`: `P1`
- `module`: `hydro/driver`
- `thesis_gap`: the current repository docs require hydro-only `Te` visualization outputs in `P1`, but the thesis-aligned implementation at this stage does not yet carry a full EOS/units-backed electron temperature recovery contract for plotting
- `decision`: pre-PPM numerical checks and analysis outputs recover hydro-only electron temperature as `Te ~ Pe / rho`, with `Pe` taken from the current gamma-law electron channel and unit gas-constant normalization, and surface that quantity only as an auditable visualization proxy
- `status`: `temporary`
- `impact_on_gates`: `P1a` and pre-PPM numerical checks may emit `hydro-only Te` text outputs and plots using this proxy, but final `P1` closeout must continue to treat the plotted temperature as a hydro-only visualization quantity rather than a fully closed EOS-backed thermodynamic observable
- `impact_on_result_interpretation`: current `Te` shell maps and polar slices are useful for comparing hydro paths and spotting angular artifacts, but their magnitude should be interpreted as a plotting proxy until the later thermodynamic phases tighten the temperature recovery contract
- `owner`: `Codex`
- `date`: `2026-04-22`

### A-P1-005

- `assumption_id`: `A-P1-005`
- `phase`: `P1`
- `module`: `hydro/driver`
- `thesis_gap`: the current repository still lacks a real runtime-owned MPI hydro execution path, but the pre-PPM acceptance work now needs radial seam and parity evidence before full `PPM` closeout
- `decision`: the current `case2_radial_wave_mpi_parity` evidence is produced by simulating split-rank radial ownership on one process, reusing the same first-order static-grid hydro kernel with ownership-driven seam ghost overrides and comparing the assembled result against the single-rank baseline
- `status`: `temporary`
- `impact_on_gates`: `P1a` may use this split-rank parity/seam sanity case as stronger numerical evidence that the current ghost substrate and radial seam contract are not single-rank-only scaffolding, but `P1` closeout must continue to treat it as pre-MPI evidence rather than as full radial-MPI acceptance
- `impact_on_result_interpretation`: current parity/seam outputs are valuable for checking that the radial ownership decomposition does not inject seam artifacts into the one-step hydro update, but they do not yet prove that the runtime, halo exchange, and diagnostics all work under a real multi-process MPI launch
- `owner`: `Codex`
- `date`: `2026-04-22`

### A-P1-006

- `assumption_id`: `A-P1-006`
- `phase`: `P1`
- `module`: `mesh/ghost` + `hydro/driver`
- `thesis_gap`: the repository now includes a real multi-process MPI parity/seam sanity case, but the current pre-PPM hydro core still only consumes one exchanged radial ghost layer and runs the directional hydro update through the static-grid helper rather than through the later full runtime/MPI orchestration
- `decision`: the current `case2_radial_wave_mpi_real` launches under real `mpiexec`, performs ownership-local radial halo exchange for `rho/momentum/E` and reuses the first-order static-grid hydro kernel on each rank before gathering the result back to root for seam/parity comparison against the single-rank baseline
- `status`: `temporary`
- `impact_on_gates`: this real MPI case is now valid pre-PPM evidence that radial halo exchange and seam handling do not break the current one-step hydro update, but it must still be treated as a targeted parity/seam sanity case rather than as the final integrated MPI acceptance for `P1`
- `impact_on_result_interpretation`: the current real MPI outputs demonstrate that the active first-order radial-wave case remains parity-consistent across a true multi-process decomposition, but they do not yet prove higher-order ghost consumption, ALE, macro-zoning, or full runtime-owned MPI execution
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-007

- `assumption_id`: `A-P1-007`
- `phase`: `P1`
- `module`: `mesh/ghost` + `hydro/driver`
- `thesis_gap`: the repository now exercises true-3D low-mode and pole/seam stress cases under real multi-process MPI, but these checks still reuse the current first-order static-grid hydro kernel, a single exchanged radial ghost layer, and root-side parity comparison against a single-rank baseline rather than the later fully integrated runtime-owned MPI orchestration
- `decision`: the current `case3_true3d_low_mode_mpi_real` and `case4_pole_seam_stress_mpi_real` launch under real `mpiexec`, perform ownership-local radial halo exchange, gather the post-step states back to root, compare them against the single-rank full-directional baseline with `mpi_parity_diagnostics.txt`, and emit `budget_residual_comparison.txt` as pre-PPM stronger numerical evidence
- `status`: `temporary`
- `impact_on_gates`: these real MPI case3/case4 outputs are now valid pre-PPM evidence that the current ghost substrate and directional hydro path remain parity-consistent on true-3D and pole-adjacent stress cases, but they must still be treated as targeted stronger checks rather than as the final integrated MPI closeout for `P1`
- `impact_on_result_interpretation`: the current outputs demonstrate that the active first-order hydro path preserves true-3D low-mode and pole/seam stress behavior across a real multi-process radial decomposition, but they do not yet prove `ng=3` reconstruction consumption, ALE, macro-zoning, or the later fully integrated runtime-owned MPI execution path
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-008

- `assumption_id`: `A-P1-008`
- `phase`: `P1`
- `module`: `hydro/reconstruction`
- `thesis_gap`: the repository now has a single-rank direction-local characteristic/traced-interface PPM path, a real-MPI radial-only high-order sanity harness, and a real-MPI full-direction true-3D high-order sanity harness, but broader MPI high-order closeout and later tuning items such as limiter refinement remain outside the currently validated scope
- `decision`: the current implementation upgrades `Workstream 4a` from a narrow face-value shell to a direction-local characteristic PPM path that consumes `ng=3` directional ghost layers, builds parabolic cell profiles, traces left/right interface states with `dt/dx`, downgrades troubled interfaces to first-order HLLC when basis construction, tracing, monotonicity, or physicality checks fail, proves real-MPI radial high-order consumption through `case2_radial_wave_mpi_real_ppm`, and proves real-MPI full-direction high-order ghost consumption through `case3_true3d_low_mode_mpi_real_ppm`
- `status`: `temporary`
- `impact_on_gates`: `Workstream 4a` can now claim that `HydroOperator` truly executes a characteristic / traced-interface PPM reconstruction consuming `ng=3` directional ghosts on single-rank runs, that the radial-only high-order path remains parity-consistent under real MPI, and that a true-3D full-direction high-order path also remains parity-consistent under real MPI, but final `P1` PPM closeout still requires broader MPI/high-order closeout and later tuning before the phase can treat the reconstruction path as fully closed
- `impact_on_result_interpretation`: current PPM diagnostics demonstrate that a higher-order characteristic/traced-interface reconstruction path exists, is runtime-executed on top of the unified ghost substrate, survives a real-MPI radial parity/seam sanity check, and survives a real-MPI full-direction true-3D low-mode sanity check, but its downgrade heuristics and the lack of validated broader MPI/high-order closeout mean the path should still be interpreted as a staged engineering closeout rather than as the final thesis-level MPI-ready PPM endpoint
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-009

- `assumption_id`: `A-P1-009`
- `phase`: `P1`
- `module`: `hydro/driver`
- `thesis_gap`: the accepted `P1` documents require hydro to provide an explicit `DtAdvice` hard cap from the active spherical hydro geometry, but the current repository documents do not lock a single thesis-direct CFL prefactor or a line-by-line formula for combining the directional signal speeds and effective widths in the present pre-ALE implementation
- `decision`: the current implementation computes `DtAdvice` with a conservative explicit CFL prefactor of `0.5`, local signal speeds `max(|v_r| + c_s, |v_theta| + c_s, |v_phi| + c_s)` over the active hydro directions, and direction-local effective widths `dr`, `r Δtheta`, and `r sin(theta) Δphi`, then returns the minimum positive hard cap as both `hard_cap_dt` and `soft_advice_dt` while leaving final timestep ownership in runtime
- `status`: `temporary`
- `impact_on_gates`: `Workstream 5` may now claim a real non-placeholder hydro `DtAdvice` path and failing-first CFL sensitivity tests, but later ALE, macro-zoning, and integrated `P1` closeout must continue to treat the exact CFL prefactor and combination rule as an auditable engineering completion until they are benchmarked or promoted into a tighter contract
- `impact_on_result_interpretation`: current explicit timestep evidence is now derived from the real full-direction spherical hydro geometry rather than scaffold defaults, but absolute timestep magnitudes should still be interpreted as a conservative engineering choice rather than as a thesis-directly fixed final limiter constant
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-010

- `assumption_id`: `A-P1-010`
- `phase`: `P1`
- `module`: `mesh/ale` + `runtime`
- `thesis_gap`: the thesis fixes the radial moving-mesh form to the upwind moving-interface flux `F_r^{ALE} = F_r^{HLLC} - w_{face} U_*`, conservative `(VQ)` update, and an outer-boundary-driven radial face motion, but it does not give a repository-ready line-by-line prescription for the exact outer-face scaling parameter `beta`, broader ALE high-order MPI closure, or later macro-zoning coupling
- `decision`: the current implementation now introduces a runtime-auditable `mesh_update_proposal` carrying homologously predicted proposed radial faces and face speeds, computes the outer face speed as `beta * v_CM^shell(t^n)` with the current default `beta = 1.0`, keeps `HydroOperator` on a const canonical mesh snapshot, allows only runtime-owned commit of the proposed radial faces, selects the ALE HLLC branch against `w_face` rather than `x/t = 0`, and updates the active ALE path in conservative extensive form `(VQ)` before recovering new cell averages against the proposal geometry
- `status`: `temporary`
- `impact_on_gates`: `Workstream 6` may now claim a thesis-consistent radial ALE proposal/commit contract plus active radial moving-mesh flux entry with failing-first tests for moving-interface branch selection and extensive conservation, but broader ALE closeout must still cover later high-order MPI consumption, macro-zoning coupling, and tuning before `P1b/P1c` can be treated as fully complete
- `impact_on_result_interpretation`: current ALE evidence proves ownership and auditability of mesh motion proposals, thesis-style `w_face` branch selection, and conservative extensive `(VQ)` ALE update on the active radial path, but results should still be interpreted as a radial ALE closure rather than as final moving-mesh hydro or later integrated `P1` completion
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-011

- `assumption_id`: `A-P1-011`
- `phase`: `P1`
- `module`: `mesh/macro_zoning`
- `thesis_gap`: the thesis fixes macro-zoning as a three-step detect / coarse-mesh generation / prolong-restrict algorithm within each radial layer, but it does not provide a repository-ready object model for representing coarse theta-bands, variable `phi` coarse factors by band, or an explicit guardrail API for enforcing that coarse hydro state remains a short-lived work view rather than canonical truth
- `decision`: `Workstream 7a` now introduces a substrate-only `MacroZoneMap` plus conservative fine/coarse transfer API under `src/mesh/macro_zoning/`; detection is driven by `ΔS_theta = r Δtheta` and `ΔS_phi = r sin(theta) Δphi` against `Δr/2`, `phi` coarse factors may vary by theta-band, conservative transfer covers `rho`, `mom_r`, `mom_theta`, `mom_phi`, `E_fluid_total`, and `chi_e`, and the coarse package is explicitly marked as a short-lived non-authoritative work view
- `status`: `temporary`
- `impact_on_gates`: `Workstream 7a` may now claim a real macro-zoning substrate with failing-first tests for thesis-style detection, conservative transfer, and non-authoritative coarse buffers, but `P1c` remains incomplete until active hydro coarse/fine coupling, MPI seam/parity, and macro-zoning acceptance artifacts are closed
- `impact_on_result_interpretation`: current macro-zoning evidence proves that the repository has a thesis-aligned detect/map/transfer substrate rather than interface-only placeholders, but it does not yet prove that active hydro consumes coarse states or that radial MPI plus macro-zoning are numerically closed together
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-012

- `assumption_id`: `A-P1-012`
- `phase`: `P1`
- `module`: `mesh/macro_zoning` + `hydro/driver`
- `thesis_gap`: the thesis fixes binary recombination and coarse/fine transfer within each radial layer, but it does not provide a repository-ready statement that the current implementation must assume uniform angular spacing, require power-of-two fine `theta/phi` counts, or use zeroth-order constant-fill prolongation in the first active coupling slice
- `decision`: the current `Workstream 7` implementation assumes uniform angular spacing through the active spherical mesh metadata, requires power-of-two fine `theta/phi` resolutions to support binary recombination, and uses zeroth-order constant-fill prolongation after coarse angular hydro updates; the first active coupling slice applies macro-zoning only to single-rank first-order `theta/phi` sweeps while radial sweep remains on the fine mesh
- `status`: `temporary`
- `impact_on_gates`: `Workstream 7b` may now claim a real single-rank macro-zoned angular hydro path with failing-first tests for active consumption and transactional failure behavior, but later work must still close high-order macro-zoning consumption, MPI parity/seam with macro-zoning enabled, and improved prolongation quality before `P1c` can be treated as complete
- `impact_on_result_interpretation`: current macro-zoning results should be interpreted as a first usable single-rank coarse/fine coupling slice rather than as the final macro-zoning numerical quality; future non-uniform angular meshes, broader MPI coupling, and higher-order prolongation may require revisiting these assumptions
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-013

- `assumption_id`: `A-P1-013`
- `phase`: `P1`
- `module`: `mesh/macro_zoning` + `mesh/ghost` + `hydro/driver`
- `thesis_gap`: the thesis requires macro-zoning plus radial MPI hydro acceptance, but it does not give a repository-ready staged closeout plan for how much real-MPI evidence is enough before broader true-3D macro-on cases, macro-on/off comparison harnesses, and higher-order macro-zoning are implemented
- `decision`: the current `Workstream 7c` slice is intentionally narrow: real MPI now proves macro-zoning-enabled `case2` radial-wave parity / seam sanity, `case3` true-3D low-mode parity / budget artifacts, and `case4` pole/seam stress parity / budget artifacts against macro-zoning-enabled single-rank baselines, while higher-order macro-zoning and richer macro-on/off artifacts remain for later slices
- `status`: `temporary`
- `impact_on_gates`: `Workstream 7c` may now claim real-MPI macro-zoning footholds for `case2`, `case3`, and `case4`, but `P1c` remains incomplete until later slices close broader macro-on/off coverage and acceptance artifacts
- `impact_on_result_interpretation`: current MPI macro-zoning outputs should be interpreted as proof that the active first-order macro-zoned hydro path can survive real radial decomposition on the radial-wave, true-3D low-mode, and pole/seam stress cases, not as proof that all macro-zoning + MPI combinations are numerically closed
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-014

- `assumption_id`: `A-P1-014`
- `phase`: `P1`
- `module`: `hydro benchmark / Sedov`
- `thesis_gap`: the current benchmark design draft proposes a Sedov baseline with `gamma = 1.4`, but the active hydro implementation still fixes the plasma ideal-gas ratio to `gamma = 5/3`; the repository also does not yet include a full analytic Sedov profile evaluator, only a first benchmark slice centered on shock-radius scaling and shell symmetry
- `decision`: `Sedov 1a` is implemented against the current code-fixed `gamma = 5/3`, uses a spherical shock-radius reference coefficient `C_3 ≈ 1.15` for the first benchmark slice, and accepts a baseline shock-radius relative error tolerance of `15%` while the benchmark remains restricted to `single-rank / ALE off / macro off / PPM on`
- `status`: `temporary`
- `impact_on_gates`: Sedov baseline evidence may now be claimed as the first hydro benchmark validation slice, but later benchmark closeout must still revisit configurable `gamma`, richer analytic/profile comparison, and broader ALE/macro/MPI benchmark variants
- `impact_on_result_interpretation`: current Sedov outputs prove that the active spherical hydro path can propagate a strong spherical blast with reasonable shock-radius scaling, shell symmetry, and quiet tangential momentum under the code's present `gamma = 5/3` model, but they are not yet the final benchmark-grade validation of all Sedov observables
- `owner`: `Codex`
- `date`: `2026-04-23`

### A-P1-015

- `assumption_id`: `A-P1-015`
- `phase`: `P1`
- `module`: `mesh/ale`
- `thesis_gap`: A1 now makes radial ALE proposals global face-indexed objects shared by single-rank and radial-MPI paths, but Woo's thesis does not provide a repository-ready formula for the exact shell velocity averaging rule that feeds the current outer-face radial motion predictor
- `decision`: the current A1/A2 shell velocity truth inherits the existing single-rank ALE definition: for each radial shell, compute per-cell `mom_r / rho` and use the arithmetic mean over angular cells; if later A2/A3 evidence requires shell-center-of-mass or mass-weighted shell velocity, the replacement must be made before the shared global proposal builder so single-rank and multi-rank paths consume the same `global_shell_velocity`
- `status`: `temporary`
- `impact_on_gates`: A1 may close because the acceptance target is proposal parity, seam uniqueness, and commit atomicity rather than final physical optimality of shell velocity averaging; A2 may reuse this predictor while focusing on single-rank macro/ALE compatibility
- `impact_on_result_interpretation`: current ALE benchmark and compatibility-mode results should be interpreted as using the inherited arithmetic shell-velocity predictor, not as proof that this predictor is the final thesis-optimal moving-mesh velocity definition
- `owner`: `Codex`
- `date`: `2026-04-25`

### A-P3-001

- `assumption_id`: `A-P3-001`
- `phase`: `P3`
- `module`: `radiation/hydro coupling`
- `thesis_gap`: the visible Woo thesis material gives the multigroup radiation energy equation, derives the radiation moment equation used to obtain the P1/Fick diffusion closure, and specifies the hydro-step transport of radiation pressure through `P_photon^{3/4}`, but it does not expose a repository-ready fluid momentum source term that feeds radiation force or radiation momentum back into `rho v`
- `decision`: the current P3 implementation treats radiation force / momentum feedback as out of scope; P3 claims only thesis-consistent radiation hydro-side energy transport via `P_g^{3/4}`, radiation diffusion/source updates, and radiation-electron energy exchange, not a direct radiation acceleration of the fluid momentum
- `status`: `temporary`
- `impact_on_gates`: P3 acceptance gates must not require a radiation-force source in the hydro momentum equations, and P3 closeout language must not claim radiation force / momentum feedback unless a later thesis-backed or explicitly documented extension adds and tests that source term
- `impact_on_result_interpretation`: current radiation-enabled results may include radiative heating/cooling, Marshak leakage, flux-limited diffusion, and radiation energy response to compression/expansion, but they should not be interpreted as including direct radiation pressure or radiation flux forces that push the fluid momentum
- `owner`: `Codex`
- `date`: `2026-04-28`

### A-P4-001

- `assumption_id`: `A-P4-001`
- `phase`: `P4`
- `module`: `alpha/composition`
- `thesis_gap`: the current repository state does not carry separate authoritative deuterium and tritium number-density fields even though the alpha birth source needs `n_D n_T`
- `decision`: P4-0 uses the existing single-fluid DT closure and P2 thermodynamic recovery with `composition_model=equimolar_dt_from_p2_recovery`, setting `n_D=n_T=0.5*n_i`
- `status`: `temporary`
- `impact_on_gates`: P4-0/P4-1 contract tests may validate alpha source coefficients under equimolar DT, but later species-resolved composition support must replace this assumption before claiming species-resolved alpha yields
- `impact_on_result_interpretation`: alpha birth source results are DT-equimolar model results, not arbitrary-mixture alpha yields
- `owner`: `Codex`
- `date`: `2026-04-28`

### A-P4-002

- `assumption_id`: `A-P4-002`
- `phase`: `P4`
- `module`: `alpha/reactivity`
- `thesis_gap`: the thesis states DT reactivity comes from the Bosch-Hale model, and P4-1 now uses a project-owner supplied coefficient reference to lock the local implementation form
- `decision`: P4-1 promotes `bosch_hale_dt` from a hard-fail target to a locally locked provider contract using the `project_owner_2026_04_28` coefficients and reference values recorded in `artifacts/p4/p4-1-bosch-hale-dt-reference.md`; `constant_user_supplied` remains only an oracle/test mode
- `status`: `promoted_to_contract`
- `impact_on_gates`: P4-1 acceptance requires at least one successful Bosch-Hale physical-path operator run and must not close on constant reactivity alone
- `impact_on_result_interpretation`: Bosch-Hale alpha birth source results use the locally locked DT reactivity formula with `Ti_old` and cgs `cm^3/s` units, while composition remains equimolar DT under `A-P4-001`
- `owner`: `Codex`
- `date`: `2026-04-28`

### A-P4-003

- `assumption_id`: `A-P4-003`
- `phase`: `P4`
- `module`: `alpha/coefficients`
- `thesis_gap`: the thesis states alpha birth energy and initial speed but does not give a repository-ready explicit formula for `v_alpha0`
- `decision`: P4-0 computes `v_alpha0=sqrt(2E_alpha0/m_alpha)` from nonrelativistic birth energy and reports `v_alpha0_source=nonrelativistic_birth_energy`
- `status`: `temporary`
- `impact_on_gates`: P4-0 coefficient tests can validate `lambda_drag` and `D_alpha` units and positivity, but later benchmark closeout must keep this assumption visible
- `impact_on_result_interpretation`: `D_alpha` is based on a nonrelativistic birth-speed completion rather than a separately thesis-tabulated velocity
- `owner`: `Codex`
- `date`: `2026-04-28`
