# P4-0 Alpha State / Provider Contract Acceptance

Status: implemented and verified

Scope:
- Authoritative alpha energy density state contract.
- Alpha pressure scalar roundtrip for future H-stage `P_alpha^(3/5)` transport.
- Equimolar DT composition provider using P2 thermodynamic recovery.
- Thesis-Spitzer alpha-electron relaxation time provider for Eq. 5.262 form.
- Atzeni one-group drag coefficient arrays.
- Constant user-supplied DT reactivity contract-test path.

Non-claims:
- No alpha diffusion solve.
- No distributed alpha HYPRE backend.
- No production `A` stage registration.
- No Bosch-Hale thesis reactivity claim.
- No alpha benchmark parity.

Assumption ledger delta:
- `A-P4-001`: equimolar DT composition from P2 recovery.
- `A-P4-002`: Bosch-Hale coefficients are not yet locally locked; `constant_user_supplied` is contract-test only.
- `A-P4-003`: `v_alpha0` uses nonrelativistic birth-energy completion.

Focused verification:
- `cmake --build F:\dec3d\build --config Debug --target dec3d_contract_alpha_contract --clean-first -- /m:4`
- `ctest --test-dir F:\dec3d\build -C Debug -R "dec3d_contract_alpha_contract" --output-on-failure`
- Result: `1/1` focused test passed.

Full regression:
- `cmake --build F:\dec3d\build --config Debug --clean-first -- /m:4`
- `ctest --test-dir F:\dec3d\build -C Debug --output-on-failure`
- Result: `105/105` default build tests passed.
- Note: existing manual Noh benchmark C4127 warnings are still emitted by unrelated manual benchmark translation units.

HYPRE regression:
- `cmake --build F:\dec3d\build-hypre --config Debug --clean-first -- /m:4`
- `ctest --test-dir F:\dec3d\build-hypre -C Debug --output-on-failure`
- Result: `112/112` HYPRE build tests passed.
- Note: existing HYPRE-link `LNK4098` warnings and existing manual Noh benchmark C4127 warnings remain unrelated to P4-0.
