from __future__ import annotations

import argparse
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
class Frame:
    index: int
    target_time_s: float
    no_ale_time_s: float
    ale_time_s: float
    no_ale_radius: np.ndarray
    ale_radius: np.ndarray
    exact_radius: np.ndarray
    no_ale_rho: np.ndarray
    ale_rho: np.ndarray
    exact_rho: np.ndarray
    no_ale_velocity: np.ndarray
    ale_velocity: np.ndarray
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
    profile_paths = sorted(case_dir.glob("radial_profile_step*.txt"))
    t1_path = case_dir / "radial_profile_t1.txt"
    if t1_path.exists():
        profile_paths.append(t1_path)
    profiles = [parse_profile(path) for path in profile_paths]
    if len(profiles) < 2:
        raise ValueError(f"need at least two radial_profile_step*.txt files in {case_dir}")
    profiles.sort(key=lambda profile: profile.time_s)
    unique_profiles: list[Profile] = []
    for profile in profiles:
        if unique_profiles and abs(unique_profiles[-1].time_s - profile.time_s) <= 1.0e-14:
            unique_profiles[-1] = profile
        else:
            unique_profiles.append(profile)
    return unique_profiles


def nearest_profile(profiles: list[Profile], target_time_s: float) -> Profile:
    return min(profiles, key=lambda profile: abs(profile.time_s - target_time_s))


def exact_noh_profile(
    solver: Noh,
    radius: np.ndarray,
    time_s: float,
    rho0: float,
    vr0: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    if time_s <= 0.0:
        return (
            np.full_like(radius, rho0, dtype=float),
            np.full_like(radius, vr0, dtype=float),
            0.0,
        )
    eval_radius = radius.copy()
    if eval_radius.size > 0 and eval_radius[0] <= 0.0:
        positive = eval_radius[eval_radius > 0.0]
        eval_radius[0] = 0.5 * positive[0] if positive.size > 0 else 1.0e-12
    exact = solver(eval_radius, time_s)
    exact_rho = np.asarray(exact["density"], dtype=float)
    exact_velocity = np.asarray(exact["velocity"], dtype=float)
    exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else time_s / 3.0
    return exact_rho, exact_velocity, exact_shock_radius


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


def render_frame(
    frame: Frame,
    output_path: Path,
    x_limit: tuple[float, float],
    rho_limit: tuple[float, float],
    velocity_limit: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 8.0), sharex=True)

    axes[0].plot(frame.exact_radius, frame.exact_rho, color="black", label="ExactPack", linewidth=2.0)
    axes[0].plot(frame.no_ale_radius, frame.no_ale_rho, label="PPM+Macro MPI24", linewidth=1.8, linestyle="--")
    axes[0].plot(frame.ale_radius, frame.ale_rho, label="PPM+Macro MPI24 + ALE", linewidth=1.8, linestyle="-.")
    axes[0].axvline(frame.exact_shock_radius, color="0.25", linestyle=":", linewidth=1.2)
    axes[0].set_xlim(*x_limit)
    axes[0].set_ylim(*rho_limit)
    axes[0].set_ylabel("density")
    axes[0].set_title("Spherical Noh density")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="upper right")

    axes[1].plot(frame.exact_radius, frame.exact_velocity, color="black", label="ExactPack", linewidth=2.0)
    axes[1].plot(
        frame.no_ale_radius,
        frame.no_ale_velocity,
        label="PPM+Macro MPI24",
        linewidth=1.8,
        linestyle="--",
    )
    axes[1].plot(
        frame.ale_radius,
        frame.ale_velocity,
        label="PPM+Macro MPI24 + ALE",
        linewidth=1.8,
        linestyle="-.",
    )
    axes[1].axvline(frame.exact_shock_radius, color="0.25", linestyle=":", linewidth=1.2, label="Exact shock")
    axes[1].set_xlim(*x_limit)
    axes[1].set_ylim(*velocity_limit)
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("radial velocity")
    axes[1].set_title("Spherical Noh radial velocity")
    axes[1].grid(alpha=0.25)
    axes[1].legend(loc="lower right")

    fig.suptitle(
        "Noh radial profiles, "
        f"target t={frame.target_time_s:.5f} s "
        f"(no-ALE {frame.no_ale_time_s:.5f}, ALE {frame.ale_time_s:.5f})"
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
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
    return f"{value:g}".replace("-", "m").replace(".", "p")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-ale-case-dir", required=True)
    parser.add_argument("--ale-case-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--rho0", type=float, default=1.0)
    parser.add_argument("--vr0", type=float, default=-1.0)
    parser.add_argument("--time-interval", type=float, default=0.005)
    parser.add_argument("--x-max", type=float, default=0.15)
    parser.add_argument("--radial-samples", type=int, default=800)
    parser.add_argument("--duration-ms", type=int, default=160)
    args = parser.parse_args()

    no_ale_case_dir = Path(args.no_ale_case_dir)
    ale_case_dir = Path(args.ale_case_dir)
    output_dir = Path(args.output_dir)
    frames_dir = output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    no_ale_profiles = load_profiles(no_ale_case_dir)
    ale_profiles = load_profiles(ale_case_dir)
    final_time_s = min(no_ale_profiles[-1].time_s, ale_profiles[-1].time_s)
    interval = float(args.time_interval)
    frame_count = int(round(final_time_s / interval))
    target_times = np.array([interval * index for index in range(frame_count + 1)], dtype=float)
    if target_times[-1] < final_time_s - 1.0e-12:
        target_times = np.append(target_times, final_time_s)
    target_times[-1] = final_time_s

    x_limit = (0.0, float(args.x_max))
    exact_radius = np.linspace(x_limit[0], x_limit[1], int(args.radial_samples), dtype=float)
    solver = Noh(geometry=3, gamma=float(args.gamma), rho0=float(args.rho0), u0=float(args.vr0))

    frames: list[Frame] = []
    for index, time_s in enumerate(target_times, start=1):
        no_ale_profile = nearest_profile(no_ale_profiles, float(time_s))
        ale_profile = nearest_profile(ale_profiles, float(time_s))
        exact_rho, exact_velocity, exact_shock_radius = exact_noh_profile(
            solver,
            exact_radius,
            float(time_s),
            float(args.rho0),
            float(args.vr0),
        )
        frames.append(
            Frame(
                index=index,
                target_time_s=float(time_s),
                no_ale_time_s=no_ale_profile.time_s,
                ale_time_s=ale_profile.time_s,
                no_ale_radius=no_ale_profile.radius,
                ale_radius=ale_profile.radius,
                exact_radius=exact_radius,
                no_ale_rho=no_ale_profile.rho,
                ale_rho=ale_profile.rho,
                exact_rho=exact_rho,
                no_ale_velocity=no_ale_profile.radial_velocity,
                ale_velocity=ale_profile.radial_velocity,
                exact_velocity=exact_velocity,
                exact_shock_radius=exact_shock_radius,
            )
        )

    rho_limit = finite_range(
        [frame.exact_rho for frame in frames] +
        [frame.no_ale_rho[(frame.no_ale_radius >= x_limit[0]) & (frame.no_ale_radius <= x_limit[1])]
         for frame in frames] +
        [frame.ale_rho[(frame.ale_radius >= x_limit[0]) & (frame.ale_radius <= x_limit[1])]
         for frame in frames]
    )
    rho_limit = (0.0, rho_limit[1])
    velocity_limit = finite_range(
        [frame.exact_velocity for frame in frames] +
        [frame.no_ale_velocity[(frame.no_ale_radius >= x_limit[0]) & (frame.no_ale_radius <= x_limit[1])]
         for frame in frames] +
        [frame.ale_velocity[(frame.ale_radius >= x_limit[0]) & (frame.ale_radius <= x_limit[1])]
         for frame in frames]
    )

    frame_paths: list[Path] = []
    for frame in frames:
        frame_path = frames_dir / f"noh_ale_noale_exact_frame_{frame.index:03d}_t{frame.target_time_s:.5f}.png"
        render_frame(frame, frame_path, x_limit, rho_limit, velocity_limit)
        frame_paths.append(frame_path)

    gif_path = output_dir / (
        f"noh_ale_noale_exact_rho_velocity_dt{compact_decimal(interval)}_"
        f"xmax{compact_decimal(float(args.x_max))}.gif"
    )
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    manifest = [
        "# kind=noh_ale_noale_exact_profile_animation",
        f"no_ale_case_dir={no_ale_case_dir}",
        f"ale_case_dir={ale_case_dir}",
        "exact_solution=ExactPack Noh geometry=3",
        "sampling=no interpolation; each numerical curve uses its nearest raw radial_profile snapshot and native radius grid",
        f"time_start_s={target_times[0]:.17g}",
        f"time_end_s={target_times[-1]:.17g}",
        f"time_interval_s={interval:.17g}",
        f"frame_count={len(frame_paths)}",
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
