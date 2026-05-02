from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.sedov.sedov import Sedov
from PIL import Image


@dataclass(frozen=True)
class HistoryProfile:
    step_index: int
    time_s: float
    radius: np.ndarray
    rho: np.ndarray
    radial_velocity: np.ndarray


@dataclass(frozen=True)
class FrameData:
    frame_index: int
    profile: HistoryProfile
    exact_rho: np.ndarray
    exact_velocity: np.ndarray
    exact_shock_radius: float


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


def compact_decimal(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def read_history_profiles(path: Path) -> list[HistoryProfile]:
    grouped: dict[tuple[int, float], list[dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(rows)
        required = {"step", "time_s", "r_cm", "rho_g_cm3", "vr_cm_s"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"history profile is missing required columns: {path}")

        for row in reader:
            step = int(row["step"])
            time_s = float(row["time_s"])
            grouped.setdefault((step, time_s), []).append(row)

    profiles: list[HistoryProfile] = []
    for (step, time_s), rows in grouped.items():
        rows.sort(key=lambda item: int(item["r_index"]))
        profiles.append(
            HistoryProfile(
                step_index=step,
                time_s=time_s,
                radius=np.array([float(row["r_cm"]) for row in rows], dtype=float),
                rho=np.array([float(row["rho_g_cm3"]) for row in rows], dtype=float),
                radial_velocity=np.array([float(row["vr_cm_s"]) for row in rows], dtype=float),
            )
        )

    if len(profiles) < 2:
        raise ValueError(f"need at least two history profile snapshots in {path}")
    profiles.sort(key=lambda profile: (profile.time_s, profile.step_index))

    radius = profiles[0].radius
    for profile in profiles[1:]:
        if profile.radius.shape != radius.shape or not np.allclose(profile.radius, radius):
            raise ValueError("history profile radial grids are not consistent across time")
    return profiles


def select_raw_profiles(
    profiles: list[HistoryProfile],
    target_time_interval_s: float,
    include_zero: bool,
) -> list[HistoryProfile]:
    if target_time_interval_s <= 0.0:
        return profiles if include_zero else [profile for profile in profiles if profile.time_s > 0.0]

    start_time = 0.0 if include_zero else target_time_interval_s
    final_time = profiles[-1].time_s
    target_times = np.arange(start_time, final_time + 0.5 * target_time_interval_s, target_time_interval_s)
    selected: list[HistoryProfile] = []
    used_steps: set[int] = set()
    profile_times = np.array([profile.time_s for profile in profiles], dtype=float)

    for target_time in target_times:
        nearest_index = int(np.argmin(np.abs(profile_times - target_time)))
        profile = profiles[nearest_index]
        if profile.step_index in used_steps:
            continue
        if not include_zero and profile.time_s <= 0.0:
            continue
        selected.append(profile)
        used_steps.add(profile.step_index)

    if len(selected) < 2:
        raise ValueError("raw profile selection produced fewer than two frames")
    return selected


def build_frame_data(profiles: list[HistoryProfile], solver: Sedov) -> list[FrameData]:
    frames: list[FrameData] = []
    for frame_index, profile in enumerate(profiles, start=1):
        if profile.time_s > 0.0:
            exact = solver(profile.radius, profile.time_s)
            exact_rho = np.asarray(exact["density"], dtype=float)
            exact_velocity = np.asarray(exact["velocity"], dtype=float)
            exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")
        else:
            exact_rho = np.full_like(profile.rho, np.nan)
            exact_velocity = np.full_like(profile.radial_velocity, np.nan)
            exact_shock_radius = float("nan")

        frames.append(
            FrameData(
                frame_index=frame_index,
                profile=profile,
                exact_rho=exact_rho,
                exact_velocity=exact_velocity,
                exact_shock_radius=exact_shock_radius,
            )
        )
    return frames


def plot_exact_or_proxy(ax: plt.Axes, radius: np.ndarray, values: np.ndarray, label: str) -> None:
    if np.any(np.isfinite(values)):
        ax.plot(radius, values, color="black", linewidth=2.4, label=label)
    else:
        ax.plot([], [], color="black", linewidth=2.4, label=label)


def plot_shock_or_proxy(ax: plt.Axes, shock_radius: float) -> None:
    if np.isfinite(shock_radius):
        ax.axvline(shock_radius, color="0.25", linestyle=":", linewidth=1.5, label="Exact shock")
    else:
        ax.plot([], [], color="0.25", linestyle=":", linewidth=1.5, label="Exact shock")


def render_frame(
    frame: FrameData,
    output_path: Path,
    case_label: str,
    x_limit: tuple[float, float],
    rho_limit: tuple[float, float],
    velocity_limit: tuple[float, float],
) -> None:
    profile = frame.profile
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 12.0), sharex=True)

    fig.suptitle(
        f"{case_label}, raw step {profile.step_index}, t={profile.time_s:.6e} s",
        fontsize=14,
        y=0.985,
    )

    plot_exact_or_proxy(axes[0], profile.radius, frame.exact_rho, "ExactPack")
    axes[0].plot(
        profile.radius,
        profile.rho,
        color="#1f77b4",
        linestyle="--",
        linewidth=2.0,
        label="DEC3D no-ALE raw",
    )
    plot_shock_or_proxy(axes[0], frame.exact_shock_radius)
    axes[0].set_xlim(*x_limit)
    axes[0].set_ylim(*rho_limit)
    axes[0].set_ylabel("density")
    axes[0].set_title("Spherical Sedov density")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="upper right")

    plot_exact_or_proxy(axes[1], profile.radius, frame.exact_velocity, "ExactPack")
    axes[1].plot(
        profile.radius,
        profile.radial_velocity,
        color="#1f77b4",
        linestyle="--",
        linewidth=2.0,
        label="DEC3D no-ALE raw",
    )
    plot_shock_or_proxy(axes[1], frame.exact_shock_radius)
    axes[1].set_xlim(*x_limit)
    axes[1].set_ylim(*velocity_limit)
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("radial velocity")
    axes[1].set_title("Spherical Sedov radial velocity")
    axes[1].grid(alpha=0.25)
    axes[1].legend(loc="lower right")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    fig.savefig(output_path, dpi=100)
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
    parser.add_argument("--history-profile", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--case-label", default="Sedov nr1024 PPM+Macro MPI24 no-ALE")
    parser.add_argument("--target-time-interval", type=float, default=2.0e-5)
    parser.add_argument("--x-max", type=float, default=0.12)
    parser.add_argument("--duration-ms", type=int, default=220)
    parser.add_argument("--skip-zero", action="store_true")
    args = parser.parse_args()

    case_json = Path(args.case_json)
    history_profile = Path(args.history_profile)
    case_dir = history_profile.parent
    output_dir = Path(args.output_dir) if args.output_dir else (case_dir / "profile_animation_no_interp_latest")
    frames_dir = output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    eblast = float(config["init"]["blast_energy"])

    available_profiles = read_history_profiles(history_profile)
    selected_profiles = select_raw_profiles(
        available_profiles,
        float(args.target_time_interval),
        include_zero=not bool(args.skip_zero),
    )

    solver = Sedov(geometry=3, gamma=gamma, rho0=rho0, omega=0.0, eblast=eblast)
    frames = build_frame_data(selected_profiles, solver)

    x_limit = (0.0, float(args.x_max))
    x_masks = [(frame.profile.radius >= x_limit[0]) & (frame.profile.radius <= x_limit[1]) for frame in frames]
    rho_limit = finite_range(
        [frame.profile.rho[mask] for frame, mask in zip(frames, x_masks)]
        + [frame.exact_rho[mask] for frame, mask in zip(frames, x_masks)]
    )
    velocity_limit = finite_range(
        [frame.profile.radial_velocity[mask] for frame, mask in zip(frames, x_masks)]
        + [frame.exact_velocity[mask] for frame, mask in zip(frames, x_masks)]
    )
    rho_limit = (0.0, rho_limit[1])
    velocity_limit = (min(0.0, velocity_limit[0]), velocity_limit[1])

    frame_paths: list[Path] = []
    for frame in frames:
        frame_path = (
            frames_dir
            / f"sedov_noale_latest_frame_{frame.frame_index:04d}_step{frame.profile.step_index:06d}.png"
        )
        render_frame(frame, frame_path, args.case_label, x_limit, rho_limit, velocity_limit)
        frame_paths.append(frame_path)

    gif_path = output_dir / (
        f"sedov_noale_latest_raw_profiles_dt{compact_decimal(float(args.target_time_interval))}.gif"
    )
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    manifest = [
        "# kind=sedov_noale_latest_raw_profile_animation",
        f"case_dir={case_dir}",
        f"history_profile={history_profile}",
        f"case={args.case_label}",
        f"exact_solution=ExactPack Sedov geometry=3 gamma={gamma:.17g} rho0={rho0:.17g} eblast={eblast:.17g}",
        "sampling=no interpolation; each frame uses one raw history_profile snapshot",
        f"available_profile_count={len(available_profiles)}",
        f"frame_count={len(frame_paths)}",
        f"time_start_s={frames[0].profile.time_s:.17g}",
        f"time_end_s={frames[-1].profile.time_s:.17g}",
        f"target_time_interval_s={float(args.target_time_interval):.17g}",
        f"x_min={x_limit[0]:.17g}",
        f"x_max={x_limit[1]:.17g}",
        f"rho_y_min={rho_limit[0]:.17g}",
        f"rho_y_max={rho_limit[1]:.17g}",
        f"velocity_y_min={velocity_limit[0]:.17g}",
        f"velocity_y_max={velocity_limit[1]:.17g}",
        f"frames_dir={frames_dir}",
        f"gif={gif_path}",
    ]
    (output_dir / "manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(gif_path)
    print(output_dir / "manifest.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
