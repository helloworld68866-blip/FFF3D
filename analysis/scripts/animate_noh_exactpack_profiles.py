from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.noh.noh1 import Noh
from PIL import Image

from plot_hydro_fields import read_headers, read_numeric_rows


@dataclass(frozen=True)
class Profile:
    step_index: int
    time_s: float
    radius: np.ndarray
    rho: np.ndarray
    radial_velocity: np.ndarray


@dataclass(frozen=True)
class FrameData:
    frame_index: int
    time_s: float
    radius: np.ndarray
    numerical_rho: np.ndarray
    numerical_velocity: np.ndarray
    exact_rho: np.ndarray
    exact_velocity: np.ndarray
    exact_shock_radius: float


def parse_profile(path: Path) -> Profile:
    headers = read_headers(path)
    rows = read_numeric_rows(path)
    if not rows:
        raise ValueError(f"profile has no numeric rows: {path}")
    return Profile(
        step_index=int(headers.get("step_index", "0")),
        time_s=float(headers.get("time_s", "0")),
        radius=np.array([float(row[1]) for row in rows], dtype=float),
        rho=np.array([float(row[2]) for row in rows], dtype=float),
        radial_velocity=np.array([float(row[3]) for row in rows], dtype=float),
    )


def load_profiles(case_dir: Path) -> list[Profile]:
    profiles = [parse_profile(path) for path in sorted(case_dir.glob("radial_profile_step*.txt"))]
    if len(profiles) < 2:
        raise ValueError(f"need at least two radial_profile_step*.txt files in {case_dir}")
    profiles.sort(key=lambda profile: profile.time_s)

    radius = profiles[0].radius
    for profile in profiles[1:]:
        if profile.radius.shape != radius.shape or not np.allclose(profile.radius, radius):
            raise ValueError("radial profile grids are not consistent across time")
    return profiles


def interpolate_profile(profiles: list[Profile], target_time_s: float) -> Profile:
    times = np.array([profile.time_s for profile in profiles], dtype=float)
    if target_time_s < times[0] - 1.0e-15 or target_time_s > times[-1] + 1.0e-15:
        raise ValueError(f"target time {target_time_s:.17g} lies outside profile range")

    upper = int(np.searchsorted(times, target_time_s, side="left"))
    if upper < len(profiles) and abs(times[upper] - target_time_s) <= 1.0e-15:
        return profiles[upper]
    if upper == 0:
        return profiles[0]
    if upper >= len(profiles):
        return profiles[-1]

    lower = upper - 1
    denominator = times[upper] - times[lower]
    weight = 0.0 if denominator <= 0.0 else (target_time_s - times[lower]) / denominator
    return Profile(
        step_index=-1,
        time_s=target_time_s,
        radius=profiles[lower].radius,
        rho=(1.0 - weight) * profiles[lower].rho + weight * profiles[upper].rho,
        radial_velocity=(1.0 - weight) * profiles[lower].radial_velocity +
        weight * profiles[upper].radial_velocity,
    )


def finite_range(values: list[np.ndarray], padding_fraction: float = 0.05) -> tuple[float, float]:
    finite_parts = [value[np.isfinite(value)] for value in values if value.size > 0]
    finite_values = np.concatenate(finite_parts) if finite_parts else np.array([], dtype=float)
    if finite_values.size == 0:
        return 0.0, 1.0
    minimum = float(np.min(finite_values))
    maximum = float(np.max(finite_values))
    if minimum == maximum:
        padding = max(abs(minimum) * padding_fraction, 1.0)
        return minimum - padding, maximum + padding
    padding = (maximum - minimum) * padding_fraction
    return minimum - padding, maximum + padding


def build_frame_data(
    profiles: list[Profile],
    solver: Noh,
    target_times: np.ndarray,
) -> list[FrameData]:
    frames: list[FrameData] = []
    for frame_index, target_time_s in enumerate(target_times, start=1):
        profile = interpolate_profile(profiles, float(target_time_s))
        exact = solver(profile.radius, float(target_time_s))
        exact_rho = np.asarray(exact["density"], dtype=float)
        exact_velocity = np.asarray(exact["velocity"], dtype=float)
        exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else target_time_s / 3.0
        frames.append(
            FrameData(
                frame_index=frame_index,
                time_s=float(target_time_s),
                radius=profile.radius,
                numerical_rho=profile.rho,
                numerical_velocity=profile.radial_velocity,
                exact_rho=exact_rho,
                exact_velocity=exact_velocity,
                exact_shock_radius=exact_shock_radius,
            )
        )
    return frames


def render_frame(
    frame: FrameData,
    output_path: Path,
    x_limit: tuple[float, float],
    rho_limit: tuple[float, float],
    velocity_limit: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), sharex=False)

    axes[0].plot(frame.radius, frame.exact_rho, label="ExactPack", linewidth=1.8, zorder=2)
    axes[0].plot(
        frame.radius,
        frame.numerical_rho,
        label="DEC3D numerical",
        linewidth=1.8,
        linestyle="--",
        zorder=3,
    )
    axes[0].axvline(frame.exact_shock_radius, color="black", linestyle="--", linewidth=1.0)
    axes[0].set_xlim(*x_limit)
    axes[0].set_ylim(*rho_limit)
    axes[0].set_xlabel("radius")
    axes[0].set_ylabel("density")
    axes[0].set_title("rho-r")
    axes[0].grid(alpha=0.25)

    axes[1].plot(frame.radius, frame.exact_velocity, label="ExactPack", linewidth=1.8, zorder=2)
    axes[1].plot(
        frame.radius,
        frame.numerical_velocity,
        label="DEC3D numerical",
        linewidth=1.8,
        linestyle="--",
        zorder=3,
    )
    axes[1].axvline(frame.exact_shock_radius, color="black", linestyle="--", linewidth=1.0, label="ExactPack shock")
    axes[1].set_xlim(*x_limit)
    axes[1].set_ylim(*velocity_limit)
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("radial velocity")
    axes[1].set_title("velocity-r")
    axes[1].grid(alpha=0.25)

    fig.suptitle(f"Spherical Noh profile comparison, t={frame.time_s:.3f} s")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94))
    fig.savefig(output_path, dpi=140)
    plt.close(fig)


def write_gif(frame_paths: list[Path], gif_path: Path, duration_ms: int) -> None:
    images = [Image.open(path) for path in frame_paths]
    try:
        palette_images = [image.convert("P", palette=Image.Palette.ADAPTIVE) for image in images]
        palette_images[0].save(
            gif_path,
            save_all=True,
            append_images=palette_images[1:],
            duration=duration_ms,
            loop=0,
            optimize=False,
        )
    finally:
        for image in images:
            image.close()


def compact_decimal(value: float) -> str:
    text = f"{value:g}".replace("-", "m").replace(".", "p")
    return text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", required=True)
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--time-interval", type=float, default=0.01)
    parser.add_argument("--x-max", type=float, default=0.4)
    parser.add_argument("--duration-ms", type=int, default=180)
    args = parser.parse_args()

    case_json = Path(args.case_json)
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else (case_dir / "profile_animation")
    frames_dir = output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    vr0 = float(config["init"]["vr0"])
    final_time_s = float(config["runtime"]["t_end_s"])

    profiles = load_profiles(case_dir)
    interval = float(args.time_interval)
    frame_count = int(round(final_time_s / interval))
    target_times = np.array([interval * index for index in range(0, frame_count + 1)], dtype=float)
    solver = Noh(geometry=3, gamma=gamma, rho0=rho0, u0=vr0)
    frames = build_frame_data(profiles, solver, target_times)

    x_limit = (0.0, float(args.x_max))
    x_masks = [(frame.radius >= x_limit[0]) & (frame.radius <= x_limit[1]) for frame in frames]
    rho_limit = finite_range(
        [frame.numerical_rho[mask] for frame, mask in zip(frames, x_masks)] +
        [frame.exact_rho[mask] for frame, mask in zip(frames, x_masks)]
    )
    velocity_limit = finite_range(
        [frame.numerical_velocity[mask] for frame, mask in zip(frames, x_masks)] +
        [frame.exact_velocity[mask] for frame, mask in zip(frames, x_masks)]
    )
    rho_limit = (min(0.0, rho_limit[0]), rho_limit[1])

    frame_paths: list[Path] = []
    for frame in frames:
        frame_path = frames_dir / f"noh_profile_compare_frame_{frame.frame_index:03d}_t{frame.time_s:.3f}.png"
        render_frame(frame, frame_path, x_limit, rho_limit, velocity_limit)
        frame_paths.append(frame_path)

    gif_path = output_dir / (
        f"noh_rho_velocity_exactpack_compare_dt{compact_decimal(interval)}_"
        f"xmax{compact_decimal(float(args.x_max))}.gif"
    )
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    manifest = [
        "# kind=noh_profile_animation",
        f"case_dir={case_dir}",
        f"case_json={case_json}",
        "time_start_s=0",
        f"time_end_s={final_time_s:.17g}",
        f"time_interval_s={interval:.17g}",
        f"frame_count={len(frame_paths)}",
        f"x_min={x_limit[0]:.17g}",
        f"x_max={x_limit[1]:.17g}",
        f"rho_y_min={rho_limit[0]:.17g}",
        f"rho_y_max={rho_limit[1]:.17g}",
        f"velocity_y_min={velocity_limit[0]:.17g}",
        f"velocity_y_max={velocity_limit[1]:.17g}",
        "numerical_time_sampling=linear interpolation between radial_profile_step*.txt snapshots",
        "exact_solution=ExactPack Noh evaluated at each target frame time",
        f"frames_dir={frames_dir}",
        f"gif={gif_path}",
    ]
    (output_dir / "manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(gif_path)
    print(output_dir / "manifest.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
