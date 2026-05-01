from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.sedov.sedov import Sedov
from PIL import Image

from plot_hydro_fields import read_headers, read_numeric_rows


@dataclass(frozen=True)
class Profile:
    step_index: int
    time_s: float
    radius: np.ndarray
    rho: np.ndarray
    mom_r: np.ndarray


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
        mom_r=np.array([float(row[5]) for row in rows], dtype=float),
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
        mom_r=(1.0 - weight) * profiles[lower].mom_r + weight * profiles[upper].mom_r,
    )


def finite_range(values: list[np.ndarray], padding_fraction: float = 0.05) -> tuple[float, float]:
    finite_values = np.concatenate([value[np.isfinite(value)] for value in values if value.size > 0])
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
    solver: Sedov,
    target_times: np.ndarray,
) -> list[FrameData]:
    frames: list[FrameData] = []
    for frame_index, target_time_s in enumerate(target_times, start=1):
        profile = interpolate_profile(profiles, float(target_time_s))
        numerical_velocity = np.divide(
            profile.mom_r,
            profile.rho,
            out=np.zeros_like(profile.mom_r),
            where=profile.rho > 0.0,
        )
        if target_time_s > 0.0:
            exact = solver(profile.radius, float(target_time_s))
            exact_rho = np.asarray(exact["density"], dtype=float)
            exact_velocity = np.asarray(exact["velocity"], dtype=float)
            exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")
        else:
            exact_rho = np.full_like(profile.rho, np.nan)
            exact_velocity = np.full_like(profile.rho, np.nan)
            exact_shock_radius = float("nan")
        frames.append(
            FrameData(
                frame_index=frame_index,
                time_s=float(target_time_s),
                radius=profile.radius,
                numerical_rho=profile.rho,
                numerical_velocity=numerical_velocity,
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

    if frame.time_s > 0.0:
        axes[0].plot(frame.radius, frame.exact_rho, label="ExactPack", linewidth=1.8, zorder=2)
        axes[0].axvline(frame.exact_shock_radius, color="black", linestyle="--", linewidth=1.0)
    else:
        axes[0].text(
            0.03,
            0.9,
            "ExactPack point-blast is singular at t=0",
            transform=axes[0].transAxes,
            fontsize=9,
        )
    axes[0].plot(
        frame.radius,
        frame.numerical_rho,
        label="DEC3D numerical",
        linewidth=1.8,
        linestyle="--",
        zorder=3,
    )
    axes[0].set_xlim(*x_limit)
    axes[0].set_ylim(*rho_limit)
    axes[0].set_xlabel("radius")
    axes[0].set_ylabel("density")
    axes[0].set_title("rho-r")
    axes[0].grid(alpha=0.25)

    if frame.time_s > 0.0:
        axes[1].plot(frame.radius, frame.exact_velocity, label="ExactPack", linewidth=1.8, zorder=2)
        axes[1].axvline(frame.exact_shock_radius, color="black", linestyle="--", linewidth=1.0, label="ExactPack shock")
    else:
        axes[1].text(
            0.03,
            0.9,
            "ExactPack not plotted at t=0",
            transform=axes[1].transAxes,
            fontsize=9,
        )
    axes[1].plot(
        frame.radius,
        frame.numerical_velocity,
        label="DEC3D numerical",
        linewidth=1.8,
        linestyle="--",
        zorder=3,
    )
    axes[1].set_xlim(*x_limit)
    axes[1].set_ylim(*velocity_limit)
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("radial velocity")
    axes[1].set_title("velocity-r")
    axes[1].grid(alpha=0.25)

    fig.suptitle(f"Sedov profile comparison, t={frame.time_s:.6f} s")
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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", required=True)
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--time-interval", type=float, default=1.0e-4)
    parser.add_argument("--early-end", type=float)
    parser.add_argument("--early-interval", type=float)
    parser.add_argument("--late-interval", type=float)
    parser.add_argument("--include-zero", action="store_true")
    parser.add_argument("--x-max", type=float, default=0.15)
    parser.add_argument("--duration-ms", type=int, default=350)
    args = parser.parse_args()

    case_json = Path(args.case_json)
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else (case_dir / "profile_animation")
    frames_dir = output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    eblast = float(config["init"]["blast_energy"])
    final_time_s = float(config["runtime"]["t_end_s"])

    profiles = load_profiles(case_dir)
    interval = float(args.time_interval)
    if args.early_end is not None or args.early_interval is not None or args.late_interval is not None:
        if args.early_end is None or args.early_interval is None or args.late_interval is None:
            raise ValueError("--early-end, --early-interval, and --late-interval must be provided together")
        early_end = float(args.early_end)
        early_interval = float(args.early_interval)
        late_interval = float(args.late_interval)
        times: list[float] = [0.0] if args.include_zero else []
        early_count = int(round(early_end / early_interval))
        times.extend(early_interval * index for index in range(1, early_count + 1))
        late_count = int(round((final_time_s - early_end) / late_interval))
        times.extend(early_end + late_interval * index for index in range(1, late_count + 1))
        target_times = np.array(times, dtype=float)
    else:
        frame_count = int(round(final_time_s / interval))
        start_index = 0 if args.include_zero else 1
        target_times = np.array([interval * index for index in range(start_index, frame_count + 1)], dtype=float)
    solver = Sedov(geometry=3, gamma=gamma, rho0=rho0, omega=0.0, eblast=eblast)
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
    velocity_limit = (min(0.0, velocity_limit[0]), velocity_limit[1])

    frame_paths: list[Path] = []
    for frame in frames:
        frame_path = frames_dir / f"profile_compare_frame_{frame.frame_index:03d}_t{frame.time_s:.6f}.png"
        render_frame(frame, frame_path, x_limit, rho_limit, velocity_limit)
        frame_paths.append(frame_path)

    if args.early_end is not None:
        gif_name = "sedov_rho_velocity_exactpack_compare_t0_dt000025_dt00005_xmax015.gif"
    else:
        gif_name = "sedov_rho_velocity_exactpack_compare_dt0001_xmax015.gif"
    gif_path = output_dir / gif_name
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    manifest = [
        "# kind=sedov_profile_animation",
        f"case_dir={case_dir}",
        f"case_json={case_json}",
        f"time_start_s={float(target_times[0]):.17g}",
        f"time_end_s={final_time_s:.17g}",
        f"time_interval_s={'piecewise' if args.early_end is not None else f'{interval:.17g}'}",
        f"early_end_s={args.early_end if args.early_end is not None else ''}",
        f"early_interval_s={args.early_interval if args.early_interval is not None else ''}",
        f"late_interval_s={args.late_interval if args.late_interval is not None else ''}",
        f"frame_count={len(frame_paths)}",
        f"x_min={x_limit[0]:.17g}",
        f"x_max={x_limit[1]:.17g}",
        f"rho_y_min={rho_limit[0]:.17g}",
        f"rho_y_max={rho_limit[1]:.17g}",
        f"velocity_y_min={velocity_limit[0]:.17g}",
        f"velocity_y_max={velocity_limit[1]:.17g}",
        "numerical_time_sampling=linear interpolation between radial_profile_step*.txt snapshots",
        "exact_solution=ExactPack Sedov evaluated at each target frame time",
        "t0_note=ExactPack point-blast solution is singular at t=0; if included, only the numerical profile is plotted on that frame",
        f"frames_dir={frames_dir}",
        f"gif={gif_path}",
    ]
    (output_dir / "manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(gif_path)
    print(output_dir / "manifest.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
