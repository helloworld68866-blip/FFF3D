from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


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


def load_velocity_slice(
    path: Path,
    radial_min_cm: float,
    radial_max_cm: float,
    theta_min: float,
    theta_max: float,
    phi_index: int,
) -> tuple[
    dict[str, str],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    headers = parse_headers(path)
    radial_cells = int(headers["radial_cells"])
    theta_cells = int(headers["theta_cells"])
    phi_cells = int(headers["phi_cells"])
    opposite_phi = (phi_index + phi_cells // 2) % phi_cells
    selected_phi = {phi_index % phi_cells: 1.0, opposite_phi: -1.0}

    fields: dict[str, dict[tuple[int, int, int], float]] = {
        "vr_cm_s": {},
        "vtheta_cm_s": {},
    }
    wanted = set(fields)
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
            phi = int(parts[4])
            if phi not in selected_phi:
                continue
            fields[field][(int(parts[2]), int(parts[3]), phi)] = float(parts[5])

    dr_um = (radial_max_cm - radial_min_cm) * 1.0e4 / radial_cells
    r0_um = radial_min_cm * 1.0e4
    dtheta = (theta_max - theta_min) / theta_cells
    x_values: list[float] = []
    z_values: list[float] = []
    radius_values: list[float] = []
    radial_indices: list[int] = []
    theta_indices: list[int] = []
    vx_values: list[float] = []
    vz_values: list[float] = []

    for radial in range(radial_cells):
        radius_um = r0_um + (radial + 0.5) * dr_um
        for theta_index in range(theta_cells):
            theta = theta_min + (theta_index + 0.5) * dtheta
            sin_theta = math.sin(theta)
            cos_theta = math.cos(theta)
            for phi, side in selected_phi.items():
                key = (radial, theta_index, phi)
                vr = fields["vr_cm_s"].get(key)
                vtheta = fields["vtheta_cm_s"].get(key)
                if vr is None or vtheta is None:
                    continue
                x_values.append(side * radius_um * sin_theta)
                z_values.append(radius_um * cos_theta)
                radius_values.append(radius_um)
                radial_indices.append(radial)
                theta_indices.append(theta_index)
                vx_values.append(side * (vr * sin_theta + vtheta * cos_theta) * 1.0e-7)
                vz_values.append((vr * cos_theta - vtheta * sin_theta) * 1.0e-7)

    return (
        headers,
        np.asarray(x_values, dtype=float),
        np.asarray(z_values, dtype=float),
        np.asarray(radius_values, dtype=float),
        np.asarray(radial_indices, dtype=int),
        np.asarray(theta_indices, dtype=int),
        np.asarray(vx_values, dtype=float),
        np.asarray(vz_values, dtype=float),
    )


def plot_velocity_arrows(args: argparse.Namespace) -> None:
    snap = Path(args.snap)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    (
        headers,
        x_values,
        z_values,
        radius_values,
        radial_indices,
        theta_indices,
        vx_values,
        vz_values,
    ) = load_velocity_slice(
        snap,
        args.radial_min_cm,
        args.radial_max_cm,
        args.theta_min,
        args.theta_max,
        args.phi_index,
    )
    radius_um = args.radial_max_cm * 1.0e4
    speed = np.hypot(vx_values, vz_values)
    finite_speed = speed[np.isfinite(speed)]
    if finite_speed.size == 0:
        vmax = 1.0
    else:
        vmax = float(np.nanpercentile(finite_speed, 99.0))
        if not math.isfinite(vmax) or vmax <= 0.0:
            vmax = float(np.nanmax(finite_speed))
        if not math.isfinite(vmax) or vmax <= 0.0:
            vmax = 1.0

    time_ps = float(headers.get("time_s", "0")) * 1.0e12
    step = headers.get("step", "?")
    fig, ax = plt.subplots(figsize=(7.2, 6.8), constrained_layout=True)
    levels = np.linspace(0.0, vmax, args.levels)
    triangulation = mtri.Triangulation(x_values, z_values)
    image = ax.tricontourf(
        triangulation,
        speed,
        levels=levels,
        cmap="magma",
        extend="max",
    )
    display = radius_values >= args.min_radius_um
    if args.radial_stride > 1:
        display &= (radial_indices % args.radial_stride) == (args.radial_stride // 2)
    if args.theta_stride > 1:
        display &= (theta_indices % args.theta_stride) == 0
    if not np.any(display):
        display = np.ones_like(radius_values, dtype=bool)

    x_arrow = x_values[display]
    z_arrow = z_values[display]
    vx_arrow = vx_values[display]
    vz_arrow = vz_values[display]
    speed_arrow = np.hypot(vx_arrow, vz_arrow)
    nonzero = speed_arrow > 0.0
    arrow_u = np.zeros_like(vx_arrow)
    arrow_w = np.zeros_like(vz_arrow)
    arrow_u[nonzero] = vx_arrow[nonzero] / speed_arrow[nonzero] * args.arrow_length_um
    arrow_w[nonzero] = vz_arrow[nonzero] / speed_arrow[nonzero] * args.arrow_length_um
    quiver_kwargs = {
        "angles": "xy",
        "scale_units": "xy",
        "scale": 1.0,
        "pivot": "mid",
        "headwidth": 5.2,
        "headlength": 6.2,
        "headaxislength": 5.4,
        "minshaft": 1.2,
        "minlength": 0.0,
    }
    ax.quiver(
        x_arrow,
        z_arrow,
        arrow_u,
        arrow_w,
        color="black",
        alpha=0.72,
        width=args.arrow_width * 1.9,
        zorder=4,
        **quiver_kwargs,
    )
    ax.quiver(
        x_arrow,
        z_arrow,
        arrow_u,
        arrow_w,
        color="white",
        alpha=0.95,
        width=args.arrow_width,
        zorder=5,
        **quiver_kwargs,
    )
    boundary = plt.Circle((0.0, 0.0), radius_um, fill=False, color="black", linewidth=0.9)
    ax.add_patch(boundary)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-radius_um, radius_um)
    ax.set_ylim(-radius_um, radius_um)
    ax.set_xlabel("r (um)")
    ax.set_ylabel("z (um)")
    ax.set_title(f"{snap.stem}: velocity direction arrows, step {step}, t = {time_ps:.4g} ps")
    ax.grid(True, alpha=0.18, linewidth=0.5)
    fig.colorbar(image, ax=ax, label="|v_rz| (10^7 cm/s)", shrink=0.88)
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
    parser.add_argument("--arrow-length-um", type=float, default=2.0)
    parser.add_argument("--arrow-width", type=float, default=0.003)
    parser.add_argument("--radial-stride", type=int, default=8)
    parser.add_argument("--theta-stride", type=int, default=2)
    parser.add_argument("--min-radius-um", type=float, default=4.0)
    parser.add_argument("--levels", type=int, default=80)
    parser.add_argument("--dpi", type=int, default=180)
    plot_velocity_arrows(parser.parse_args())


if __name__ == "__main__":
    main()
