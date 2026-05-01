from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.noh2.noh2 import Noh2
from PIL import Image

from plot_hydro_fields import read_headers, read_numeric_rows


@dataclass(frozen=True)
class Profile:
    step_index: int
    time_s: float
    radius: np.ndarray
    rho: np.ndarray
    radial_velocity: np.ndarray
    e_fluid_total: np.ndarray


@dataclass(frozen=True)
class Frame:
    index: int
    profile: Profile
    exact_rho: np.ndarray
    exact_radial_velocity: np.ndarray
    exact_e_fluid_total: np.ndarray


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
        e_fluid_total=np.array([float(row[6]) for row in rows], dtype=float),
    )


def load_profiles(case_dir: Path, stride: int) -> list[Profile]:
    profiles = [parse_profile(path) for path in sorted(case_dir.glob("radial_profile_step*.txt"))]
    if not profiles:
        raise ValueError(f"no radial_profile_step*.txt files found in {case_dir}")
    profiles.sort(key=lambda profile: profile.time_s)

    unique_profiles: list[Profile] = []
    for profile in profiles:
        if unique_profiles and abs(unique_profiles[-1].time_s - profile.time_s) <= 1.0e-14:
            unique_profiles[-1] = profile
        else:
            unique_profiles.append(profile)

    selected = unique_profiles[::max(1, stride)]
    if selected[-1].step_index != unique_profiles[-1].step_index:
        selected.append(unique_profiles[-1])
    return selected


def exact_noh2_on_profile(
    solver: Noh2,
    profile: Profile,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    exact = solver(profile.radius, profile.time_s)
    exact_rho = np.asarray(exact["density"], dtype=float)
    exact_velocity = np.asarray(exact["velocity"], dtype=float)
    exact_sie = np.asarray(exact["specific_internal_energy"], dtype=float)
    exact_e_fluid_total = exact_rho * exact_sie + 0.5 * exact_rho * exact_velocity * exact_velocity
    return exact_rho, exact_velocity, exact_e_fluid_total


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
    energy_limit: tuple[float, float],
) -> None:
    profile = frame.profile
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10.0), sharex=True)

    axes[0].plot(profile.radius, frame.exact_rho, color="black", linewidth=2.0, label="ExactPack Noh2")
    axes[0].plot(profile.radius, profile.rho, color="#1f77b4", linewidth=1.8, linestyle="--", label="DEC3D")
    axes[0].set_ylabel("density")
    axes[0].set_ylim(*rho_limit)
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="best")

    axes[1].plot(profile.radius, frame.exact_radial_velocity, color="black", linewidth=2.0)
    axes[1].plot(profile.radius, profile.radial_velocity, color="#d62728", linewidth=1.8, linestyle="--")
    axes[1].set_ylabel("radial velocity")
    axes[1].set_ylim(*velocity_limit)
    axes[1].grid(alpha=0.25)

    axes[2].plot(profile.radius, frame.exact_e_fluid_total, color="black", linewidth=2.0)
    axes[2].plot(profile.radius, profile.e_fluid_total, color="#2ca02c", linewidth=1.8, linestyle="--")
    axes[2].set_ylabel("e_fluid_total")
    axes[2].set_xlabel("radius")
    axes[2].set_ylim(*energy_limit)
    axes[2].set_xlim(*x_limit)
    axes[2].grid(alpha=0.25)

    fig.suptitle(
        "Noh2 spherical profile comparison, "
        f"step={profile.step_index}, t={profile.time_s:.6f} s"
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
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--rho0", type=float, default=1.0)
    parser.add_argument("--e0", type=float, default=1.0)
    parser.add_argument("--x-max", type=float, default=1.0)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--duration-ms", type=int, default=120)
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else case_dir / "noh2_profile_animation"
    frames_dir = output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    profiles = load_profiles(case_dir, int(args.stride))
    solver = Noh2(geometry=3, gamma=float(args.gamma), rho0=float(args.rho0), e0=float(args.e0))

    frames: list[Frame] = []
    for index, profile in enumerate(profiles, start=1):
        exact_rho, exact_velocity, exact_energy = exact_noh2_on_profile(solver, profile)
        frames.append(
            Frame(
                index=index,
                profile=profile,
                exact_rho=exact_rho,
                exact_radial_velocity=exact_velocity,
                exact_e_fluid_total=exact_energy,
            )
        )

    x_limit = (0.0, float(args.x_max))
    masks = [(frame.profile.radius >= x_limit[0]) & (frame.profile.radius <= x_limit[1]) for frame in frames]
    rho_limit = finite_range(
        [frame.profile.rho[mask] for frame, mask in zip(frames, masks)] +
        [frame.exact_rho[mask] for frame, mask in zip(frames, masks)]
    )
    velocity_limit = finite_range(
        [frame.profile.radial_velocity[mask] for frame, mask in zip(frames, masks)] +
        [frame.exact_radial_velocity[mask] for frame, mask in zip(frames, masks)]
    )
    energy_limit = finite_range(
        [frame.profile.e_fluid_total[mask] for frame, mask in zip(frames, masks)] +
        [frame.exact_e_fluid_total[mask] for frame, mask in zip(frames, masks)]
    )

    frame_paths: list[Path] = []
    for frame in frames:
        frame_path = frames_dir / f"noh2_profile_frame_{frame.index:03d}_step{frame.profile.step_index:06d}.png"
        render_frame(frame, frame_path, x_limit, rho_limit, velocity_limit, energy_limit)
        frame_paths.append(frame_path)

    gif_path = output_dir / f"noh2_exactpack_profile_xmax{compact_decimal(float(args.x_max))}.gif"
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    manifest = [
        "# kind=noh2_exactpack_profile_animation",
        f"case_dir={case_dir}",
        "exact_solution=ExactPack Noh2 geometry=3",
        "sampling=no interpolation; each frame uses one raw radial_profile_step*.txt snapshot",
        "exact_sampling=ExactPack evaluated on the same native numerical radius grid for each frame",
        f"time_start_s={frames[0].profile.time_s:.17g}",
        f"time_end_s={frames[-1].profile.time_s:.17g}",
        f"frame_count={len(frame_paths)}",
        f"x_min={x_limit[0]:.17g}",
        f"x_max={x_limit[1]:.17g}",
        f"rho_y_min={rho_limit[0]:.17g}",
        f"rho_y_max={rho_limit[1]:.17g}",
        f"velocity_y_min={velocity_limit[0]:.17g}",
        f"velocity_y_max={velocity_limit[1]:.17g}",
        f"e_fluid_total_y_min={energy_limit[0]:.17g}",
        f"e_fluid_total_y_max={energy_limit[1]:.17g}",
        f"frames_dir={frames_dir}",
        f"gif={gif_path}",
    ]
    (output_dir / "manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(gif_path)
    print(output_dir / "manifest.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
