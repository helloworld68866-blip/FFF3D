from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


@dataclass(frozen=True)
class FieldSpec:
    key: str
    title: str
    colorbar: str
    scale: float = 1.0
    signed: bool = False


FIELDS = [
    FieldSpec("rho", "Density", "g/cm^3"),
    FieldSpec("Te_keV", "Electron temperature", "keV"),
    FieldSpec("Ti_keV", "Ion temperature", "keV"),
    FieldSpec("vr_cm_s", "Radial velocity", "10^7 cm/s", 1.0e-7, True),
    FieldSpec("vtheta_cm_s", "Polar velocity", "10^7 cm/s", 1.0e-7, True),
    FieldSpec("vphi_cm_s", "Azimuthal velocity", "10^7 cm/s", 1.0e-7, True),
]


def parse_headers(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            text = line[1:].strip()
            if "=" in text:
                key, value = text.split("=", 1)
                headers[key.strip()] = value.strip()
    return headers


def robust_limits(values: np.ndarray, signed: bool) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    if signed:
        bound = float(np.nanpercentile(np.abs(finite), 99.0))
        if not math.isfinite(bound) or bound <= 0.0:
            bound = float(np.nanmax(np.abs(finite)))
        if not math.isfinite(bound) or bound <= 0.0:
            bound = 1.0
        return -bound, bound
    lo = float(np.nanpercentile(finite, 1.0))
    hi = float(np.nanpercentile(finite, 99.0))
    if not math.isfinite(lo):
        lo = float(np.nanmin(finite))
    if not math.isfinite(hi):
        hi = float(np.nanmax(finite))
    if hi <= lo:
        pad = max(abs(hi), 1.0) * 1.0e-6
        lo -= pad
        hi += pad
    return lo, hi


def load_rz_plane(
    path: Path,
    radial_min_cm: float,
    radial_max_cm: float,
    theta_min: float,
    theta_max: float,
    phi_index: int,
) -> tuple[dict[str, str], dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]]:
    headers = parse_headers(path)
    radial_cells = int(headers["radial_cells"])
    theta_cells = int(headers["theta_cells"])
    phi_cells = int(headers["phi_cells"])
    opposite_phi = (phi_index + phi_cells // 2) % phi_cells
    selected_phi = {phi_index % phi_cells: 1.0, opposite_phi: -1.0}

    dr_um = (radial_max_cm - radial_min_cm) * 1.0e4 / radial_cells
    r0_um = radial_min_cm * 1.0e4
    dtheta = (theta_max - theta_min) / theta_cells

    buffers: dict[str, list[tuple[float, float, float]]] = {field.key: [] for field in FIELDS}
    wanted = set(buffers)
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#") or line.startswith("field,"):
                continue
            parts = line.strip().split(",")
            if len(parts) != 6:
                continue
            field = parts[0]
            if field not in wanted or int(parts[1]) != -1:
                continue
            radial_index = int(parts[2])
            theta_index = int(parts[3])
            phi = int(parts[4])
            side = selected_phi.get(phi)
            if side is None:
                continue
            radius_um = r0_um + (radial_index + 0.5) * dr_um
            theta = theta_min + (theta_index + 0.5) * dtheta
            x_um = side * radius_um * math.sin(theta)
            z_um = radius_um * math.cos(theta)
            buffers[field].append((x_um, z_um, float(parts[5])))

    arrays: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for key, rows in buffers.items():
        if not rows:
            arrays[key] = (np.array([]), np.array([]), np.array([]))
            continue
        data = np.asarray(rows, dtype=float)
        arrays[key] = (data[:, 0], data[:, 1], data[:, 2])
    return headers, arrays


def plot_snap(args: argparse.Namespace) -> None:
    snap = Path(args.snap)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    headers, arrays = load_rz_plane(
        snap,
        args.radial_min_cm,
        args.radial_max_cm,
        args.theta_min,
        args.theta_max,
        args.phi_index,
    )
    time_ps = float(headers.get("time_s", "0")) * 1.0e12
    step = headers.get("step", "?")

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.2), constrained_layout=True)
    for ax, spec in zip(axes.ravel(), FIELDS):
        x_um, z_um, values = arrays[spec.key]
        values = values * spec.scale
        vmin, vmax = robust_limits(values, spec.signed)
        levels = np.linspace(vmin, vmax, args.levels)
        triangulation = mtri.Triangulation(x_um, z_um)
        image = ax.tricontourf(
            triangulation,
            values,
            levels=levels,
            cmap="coolwarm" if spec.signed else "viridis",
            extend="both",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("r (um)")
        ax.set_ylabel("z (um)")
        ax.set_title(spec.title)
        ax.grid(True, alpha=0.18, linewidth=0.5)
        fig.colorbar(image, ax=ax, label=spec.colorbar, shrink=0.88)

    fig.suptitle(f"{snap.stem}: step {step}, t = {time_ps:.4g} ps")
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snap", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--radial-min-cm", type=float, default=0.0)
    parser.add_argument("--radial-max-cm", type=float, required=True)
    parser.add_argument("--theta-min", type=float, default=0.0)
    parser.add_argument("--theta-max", type=float, default=math.pi)
    parser.add_argument("--phi-index", type=int, default=0)
    parser.add_argument("--levels", type=int, default=64)
    parser.add_argument("--dpi", type=int, default=170)
    plot_snap(parser.parse_args())


if __name__ == "__main__":
    main()
