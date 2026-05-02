# Axisymmetric 2D First Version Design

## Goal

Add a first-class `axisymmetric_2d` runtime mode for pure axisymmetric perturbation studies.

The first version is deliberately scoped to:

- `phi_cells = 1`.
- no ALE or moving mesh.
- all production stages enabled: `H,T,E,R,A`.
- existing `Array3D(r, theta, phi)` state storage, with `phi` extent equal to one.

The first version does not introduce `Array2D`, a new hydro solver family, a dedicated 2D checkpoint format, ALE support, or special case-specific physics.

## Physical Semantics

The single phi cell represents the full azimuthal interval `0..2*pi`, not a wedge. Cell volumes, face areas, total mass, total energy, radiation energy, and alpha quantities therefore continue to represent the complete rotated spherical body.

The PDE degrees of freedom are two-dimensional:

- fields vary in `r` and `theta`;
- fields are independent of `phi`;
- all phi-direction PDE operations are disabled;
- `mom_phi` and `v_phi` must remain zero to roundoff.

This mode corresponds to the `m = 0` two-dimensional axisymmetric family discussed in Woo-style single-mode studies, not to a replicated 3D `m = 0` calculation.

## Input Contract

Add a mesh dimensionality option:

```ini
[mesh]
dimensionality = axisymmetric_2d
radial_cells = 256
theta_cells = 128
phi_cells = 1
ale_enabled = false
```

Validation rules:

- `axisymmetric_2d` requires spherical geometry.
- `axisymmetric_2d` requires `phi_cells = 1`.
- `axisymmetric_2d` requires no ALE in the first version.
- `axisymmetric_2d` allows only `m = 0` perturbations.
- `full_3d` keeps the existing even positive `phi_cells` requirement.

Runtime diagnostics must report:

- `mesh_dimensionality=axisymmetric_2d`;
- `active_hydro_directions=r,theta`;
- `phi_cells=1`;
- `azimuthal_weight=2pi`;
- `phi_sweep_executed=false`.

## Geometry

`BuildSphericalGeometry` should keep the existing phi faces over `0..2*pi`. With `phi_cells=1`, the single phi interval naturally has width `2*pi`.

Add a geometry contract test:

```text
sum(cell_volume over all r,theta,phi=0) ~= 4/3*pi*(r_outer^3-r_inner^3)
```

This test proves that the single phi cell represents the complete azimuthal rotation.

## Boundary Remap

Existing pole and origin remap code uses a half-turn phi mapping. In `axisymmetric_2d`, the half-turn maps back to the only phi cell:

```cpp
if (phi_cells == 1u) {
  return 0u;
}
```

Do not disable pole or origin remap. Scalar and vector parity rules still apply:

- theta pole remap flips `mom_theta` and `mom_phi`;
- origin remap flips `mom_r` and `mom_phi`;
- with `mom_phi = 0`, these transforms preserve axisymmetry.

Affected areas include hydro boundary ghosts and scalar remap used by thermal, radiation, and alpha diffusion assembly.

## Initialization

The profile initializer keeps the existing axisymmetric Legendre perturbation form:

```text
delta v_r = amplitude * f(r) * P_l(cos(theta))
```

For `axisymmetric_2d`:

- reject `m != 0`;
- require initialized `mom_phi` to be zero within a strict tolerance;
- report the perturbation formula and `m=0`.

Large nonzero `mom_phi` should fail rather than being silently cleared.

## Hydro

The runtime should configure hydro as:

```text
apply_radial_sweep = true
apply_theta_sweep = true
apply_phi_sweep = false
```

Macro zoning may remain enabled in `r` and `theta`, but the phi macro extent must be one and no phi coarse sweep may run.

The CFL estimate should use only active directions in this mode. Because the first version rejects ALE, only the static-grid hydro path is required.

## Implicit Diffusion Modules

Thermal conduction, electron-ion equilibration, radiation, and alpha should not get separate 2D implementations.

The shared distributed diffusion assembly should reduce naturally to `r,theta` by omitting phi neighbor coupling when `global_phi_cells == 1`. The expected matrix has:

- center;
- radial minus and plus;
- theta minus and plus;
- no phi minus or plus entries.

Required diagnostics:

- global unknown count equals `Nr * Ntheta`;
- phi coupling count is zero in axisymmetric mode.

## Radiation And Alpha

Radiation group loops and alpha provider loops continue to visit the local cells. Since `phi_cells=1`, the total unknown count and provider work drop accordingly.

No group physics, opacity interpolation, blackbody source, alpha deposition, residual tolerance, or solver convergence gate changes are part of this design.

## Output And Restart

The first version keeps the existing output and checkpoint style with `phi_cells=1`. A dedicated 2D CSV format is out of scope for this first slice.

Restart safety:

- axisymmetric restart metadata must record `dimensionality=axisymmetric_2d`;
- reading an axisymmetric checkpoint with a `full_3d` deck should fail unless a future explicit expansion tool is added.

## Invariants

After major stages in axisymmetric mode, diagnostics should include:

- `max_abs_mom_phi`;
- `max_abs_v_phi`;
- `phi_sweep_executed`;
- `hypre_phi_coupling_count`;
- `global_scalar_unknowns`;
- total mass and total energy where already available.

The first implementation may warn or fail on nonzero `mom_phi`, but it must not silently erase a large violation.

## Test Plan

Add focused tests in this order:

1. Input parser accepts `axisymmetric_2d` with `phi_cells=1`.
2. Input parser rejects `axisymmetric_2d` with `phi_cells != 1`.
3. Input parser rejects `axisymmetric_2d` with `m != 0`.
4. Geometry volume sum with `phi_cells=1` equals the full spherical shell volume.
5. Pole and origin remap handle `phi_cells=1`.
6. Distributed diffusion assembly has `Nr*Ntheta` rows and zero phi couplings.
7. H-only `axisymmetric_2d` smoke case runs several steps with `phi_sweep_executed=false`.
8. H,T,E,R,A `axisymmetric_2d` smoke case runs with no ALE.
9. Small regression compares `axisymmetric_2d phi=1` against `full_3d phi=8, m=0` phi averages for density, radial velocity, pressure or temperature histories.

## Acceptance Criteria

The first version is accepted when:

- no-ALE all-stage axisymmetric case runs from input deck without case-specific hacks;
- `phi_cells=1` covers full `2*pi` geometry;
- hydro reports no phi sweep;
- implicit matrices report no phi coupling;
- `mom_phi` remains zero to tolerance;
- small 2D vs replicated-3D `m=0` regression agrees within a documented numerical tolerance;
- clean Release rebuild and focused tests pass.
