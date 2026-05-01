# TOPS DT Opacity Table

Source: F:\dec3d\docs\TOPS Opacities - Summary.html

Material label: DT

TOPS date: Apr 27, 2026

Units: opacities in cm^2/g, temperature in keV, density in g/cm^3, photon energy in keV.

Generated files:

- `metadata.json`: grids, units, composition, counts, and caveats.
- `mean_opacities.csv`: Rosseland/Planck mean opacity and free-electron columns for each T/rho point.
- `multigroup_opacities.csv`: photon-energy-indexed Rosseland/Planck multigroup opacity rows for each T/rho point.
- `density_clipping_warnings.csv`: TOPS off-table density warnings.

Important caveats:

- 30 requested density points were clipped by TOPS; use `density_was_clipped` and `density_used_g_cm3` when consuming tables.
- The source labels the material as DT, but the normalized composition row is Z=1 / H / Mat ID=4525. Preserve the material label and material id together.
- The photon grid is exported as `photon_energy_keV`; this converter does not claim those values are P3 group edges.
