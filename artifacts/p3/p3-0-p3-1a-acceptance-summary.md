# P3-0 / P3-1a Acceptance Summary

Date: 2026-04-27

Implemented:

- P3-0 radiation group/state contract.
- P3-1a frozen-hydro one-group gray radiation diffusion.

Runtime evidence:

- `dec3d_contract_radiation_group_layout`
- `dec3d_contract_radiation_one_group_gray_diffusion`

Required diagnostics:

- `diagnostic_id=p3.radiation.group_layout`
- `diagnostic_id=p3.radiation.one_group_gray_diffusion`
- `radiation_group_mode=gray_full_spectrum`
- `frequency_edges_used=false`
- `radiation_energy_floor_erg_per_cm3=...`
- `Bgray_source=fixed_user_supplied`
- `marshak_enabled=false`
- `parity_claim_allowed=false`
- `electron_writeback=false`
- `updated_fields=radiation_groups` on nontrivial successful radiation updates.
- `updated_fields=none`, `metadata_written=false`, and
  `canonical_state_mutated=false` on `dt_s=0` or no-op success paths.

Forbidden claims:

- Woo 5.4 benchmark parity.
- Marshak vacuum boundary.
- Multigroup radiation.
- Radiation-matter electron writeback.
- Production `H/T/E/R` registration.
