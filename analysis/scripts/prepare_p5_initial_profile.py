import argparse
from pathlib import Path

import numpy as np
import yaml


REQUIRED_FIELDS = ["r_um", "rho_g_cm3", "te_kev", "ti_kev", "vr_cm_s"]


def resolve_profile_path(yaml_path: Path, profile_path: str) -> Path:
    candidate = Path(profile_path)
    if not candidate.is_absolute():
        candidate = (yaml_path.parent / candidate).resolve()
    return candidate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--yaml", required=True)
    parser.add_argument("--output-csv", required=True)
    parser.add_argument("--summary-json", required=True)
    args = parser.parse_args()

    yaml_path = Path(args.yaml)
    with yaml_path.open("r", encoding="utf-8") as handle:
        deck = yaml.safe_load(handle)

    profile_path = resolve_profile_path(
        yaml_path, deck["initial_profile"]["profile_path"]
    )
    data = np.load(profile_path)
    missing = [field for field in REQUIRED_FIELDS if field not in data.files]
    if missing:
        raise RuntimeError(f"missing required fields: {missing}")

    arrays = {field: np.asarray(data[field], dtype=float) for field in REQUIRED_FIELDS}
    n = len(arrays["r_um"])
    for field, values in arrays.items():
        if len(values) != n:
            raise RuntimeError(f"field {field} has length {len(values)} not {n}")
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"field {field} contains non-finite values")

    output_csv = Path(args.output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    table = np.column_stack(
        [
            arrays["r_um"] * 1.0e-4,
            arrays["r_um"],
            arrays["rho_g_cm3"],
            arrays["te_kev"],
            arrays["ti_kev"],
            arrays["vr_cm_s"],
        ]
    )
    header = "r_cm,r_um,rho_g_cm3,te_kev,ti_kev,vr_cm_s"
    np.savetxt(output_csv, table, delimiter=",", header=header, comments="")

    # The supplied profile spans 0..50 um. The original Woo 77068 text reports
    # r0=68 um, so the P5 descriptor must supply a profile-compatible r0.
    grad = np.abs(np.gradient(arrays["rho_g_cm3"], arrays["r_um"]))
    density_gradient_r0_um = float(arrays["r_um"][int(np.nanargmax(grad))])

    summary_json = Path(args.summary_json)
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(
        "{\n"
        '  "diagnostic_id": "p5.initial_profile.prepare",\n'
        f'  "source_yaml": "{yaml_path.as_posix()}",\n'
        f'  "profile_path": "{profile_path.as_posix()}",\n'
        f'  "output_csv": "{output_csv.as_posix()}",\n'
        '  "profile_fields": "r_um,rho_g_cm3,te_kev,ti_kev,vr_cm_s",\n'
        f'  "point_count": {n},\n'
        f'  "r_min_um": {float(arrays["r_um"].min())},\n'
        f'  "r_max_um": {float(arrays["r_um"].max())},\n'
        '  "woo_reference_r0_um": 68.0,\n'
        f'  "suggested_density_gradient_r0_um": {density_gradient_r0_um},\n'
        '  "r0_matches_woo_77068": false\n'
        "}\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
