# P4-1 Bosch-Hale DT Reactivity Local Reference

Date: 2026-04-28

## Scope

This artifact locally locks the project-owner supplied Bosch-Hale DT
reactivity coefficient reference for P4-1. It is used only for DT
`<sigma v>` in the alpha birth source:

```text
birth_source = n_D n_T <sigma v>_DT E_alpha0
```

## Formula

Input:

```text
T = Ti_keV
b_G = 34.3827
mrc2 = 1124656.0
```

Formula:

```text
num = 1.51361e-2 T + 4.60643e-3 T^2 - 1.06750e-4 T^3
den = 1 + 7.51886e-2 T + 1.35000e-2 T^2 + 1.3600e-5 T^3
theta = T / (1 - num / den)
zeta = (b_G^2 / (4 theta))^(1/3)

<sigma v>_DT_cm3_s =
  1.17302e-9 theta sqrt(zeta / (mrc2 T^3)) exp(-3 zeta)
```

The project-owner Python helper returns `m^3/s` by multiplying the cgs formula
by `1e-6`. DEC3D stores and consumes `cm^3/s` because densities are `cm^-3`.

## Reference Values

| Ti_keV | `<sigma v>` cm^3/s |
|---:|---:|
| 1.0 | 6.85688430133071390e-21 |
| 2.0 | 2.97742663781859648e-19 |
| 3.0 | 1.86696927818874975e-18 |
| 5.0 | 1.36578648247469899e-17 |
| 10.0 | 1.13618134482734562e-16 |
| 20.0 | 4.33035875744403697e-16 |
| 30.0 | 6.68102563545722182e-16 |
| 50.0 | 8.64897375739629370e-16 |
| 100.0 | 8.44551910167966373e-16 |

## Claim Boundary

`thesis_reactivity_claim_allowed=true` is permitted for `bosch_hale_dt` only
when these values are locked by contract tests. This does not model fuel
depletion or species-resolved DT composition.
