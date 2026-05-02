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


ERG_PER_KEV = 1.602176634e-9
DEUTERON_MASS_G = 3.3435837724e-24
TRITON_MASS_G = 5.0073567446e-24
MEAN_DT_ION_MASS_G = 0.5 * (DEUTERON_MASS_G + TRITON_MASS_G)
DEFAULT_ZBAR = 1.0


@dataclass(frozen=True)
class HistoryProfile:
    step_index: int
    time_s: float
    radius: np.ndarray
    rho: np.ndarray
    te_kev: np.ndarray
    ti_kev: np.ndarray
    radial_velocity: np.ndarray


@dataclass(frozen=True)
class ReferenceProfile:
    radius: np.ndarray
    rho: np.ndarray
    te_kev: np.ndarray
    ti_kev: np.ndarray
    radial_velocity: np.ndarray
    shock_radius: float


def compact_decimal(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


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


def positive_range(values: list[np.ndarray], padding_factor: float = 1.2) -> tuple[float, float]:
    finite_parts = [value[np.isfinite(value) & (value > 0.0)] for value in values if value.size > 0]
    finite_values = np.concatenate(finite_parts) if finite_parts else np.array([], dtype=float)
    if finite_values.size == 0:
        return 1.0e-30, 1.0
    minimum = float(np.min(finite_values))
    maximum = float(np.max(finite_values))
    if minimum == maximum:
        return minimum / padding_factor, maximum * padding_factor
    return minimum / padding_factor, maximum * padding_factor


def read_history_profiles(path: Path) -> list[HistoryProfile]:
    grouped: dict[tuple[int, float], list[dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(rows)
        required = {"step", "time_s", "r_index", "r_cm", "rho_g_cm3", "Te_keV", "Ti_keV", "vr_cm_s"}
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
                te_kev=np.array([float(row["Te_keV"]) for row in rows], dtype=float),
                ti_kev=np.array([float(row["Ti_keV"]) for row in rows], dtype=float),
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
    profile_times = np.array([profile.time_s for profile in profiles], dtype=float)
    selected: list[HistoryProfile] = []
    used_steps: set[int] = set()

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


def sedov_reference_profile(
    profile: HistoryProfile,
    initial: HistoryProfile,
    solver: Sedov,
    ambient_pressure: float,
    electron_energy_fraction: float,
) -> ReferenceProfile:
    if profile.time_s <= 0.0:
        return ReferenceProfile(
            radius=initial.radius,
            rho=initial.rho,
            te_kev=initial.te_kev,
            ti_kev=initial.ti_kev,
            radial_velocity=initial.radial_velocity,
            shock_radius=0.0,
        )

    exact = solver(profile.radius, profile.time_s)
    rho = np.asarray(exact["density"], dtype=float)
    pressure = np.asarray(exact["pressure"], dtype=float)
    radial_velocity = np.asarray(exact["velocity"], dtype=float)
    shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")

    rho = np.where(np.isfinite(rho) & (rho > 0.0), rho, 1.0)
    pressure = np.where(np.isfinite(pressure) & (pressure > ambient_pressure), pressure, ambient_pressure)
    radial_velocity = np.where(np.isfinite(radial_velocity), radial_velocity, 0.0)

    ion_number_density = rho / MEAN_DT_ION_MASS_G
    electron_number_density = DEFAULT_ZBAR * ion_number_density
    electron_pressure = electron_energy_fraction * pressure
    ion_pressure = (1.0 - electron_energy_fraction) * pressure
    te_kev = electron_pressure / electron_number_density / ERG_PER_KEV
    ti_kev = ion_pressure / ion_number_density / ERG_PER_KEV

    return ReferenceProfile(
        radius=profile.radius,
        rho=rho,
        te_kev=te_kev,
        ti_kev=ti_kev,
        radial_velocity=radial_velocity,
        shock_radius=shock_radius,
    )


def masked_values(profile: HistoryProfile, x_limit: tuple[float, float]) -> dict[str, np.ndarray]:
    mask = (profile.radius >= x_limit[0]) & (profile.radius <= x_limit[1])
    return {
        "rho": profile.rho[mask],
        "te_kev": profile.te_kev[mask],
        "ti_kev": profile.ti_kev[mask],
        "radial_velocity": profile.radial_velocity[mask],
    }


def masked_reference_values(profile: ReferenceProfile, x_limit: tuple[float, float]) -> dict[str, np.ndarray]:
    mask = (profile.radius >= x_limit[0]) & (profile.radius <= x_limit[1])
    return {
        "rho": profile.rho[mask],
        "te_kev": profile.te_kev[mask],
        "ti_kev": profile.ti_kev[mask],
        "radial_velocity": profile.radial_velocity[mask],
    }


def render_panel(
    ax: plt.Axes,
    radius: np.ndarray,
    numerical: np.ndarray,
    reference: np.ndarray,
    ylabel: str,
    x_limit: tuple[float, float],
    y_limit: tuple[float, float],
    *,
    log_y: bool = False,
    show_legend: bool = False,
) -> None:
    ax.plot(radius, numerical, color="#d62728", linewidth=2.0, label="DEC3D numerical")
    ax.plot(radius, reference, color="black", linestyle="--", linewidth=1.8, label="Two-temperature Sedov reference")
    ax.set_xlim(*x_limit)
    if log_y:
        ax.set_yscale("log")
    ax.set_ylim(*y_limit)
    ax.set_xlabel("r (cm)")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    if show_legend:
        ax.legend(loc="upper right")


def render_frame(
    current: HistoryProfile,
    reference: ReferenceProfile,
    output_path: Path,
    title_prefix: str,
    x_limit: tuple[float, float],
    rho_limit: tuple[float, float],
    velocity_limit: tuple[float, float],
    te_limit: tuple[float, float],
    ti_limit: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14.0, 9.0))
    time_text = f"{current.time_s:.3e}" if current.time_s < 1.0e-2 else f"{current.time_s:.3f}"
    fig.suptitle(
        f"{title_prefix}, t={time_text} s",
        fontsize=16,
        y=0.98,
    )

    render_panel(
        axes[0, 0],
        current.radius,
        current.rho,
        reference.rho,
        r"Density $\rho$ (g/cm$^3$)",
        x_limit,
        rho_limit,
        show_legend=True,
    )
    render_panel(
        axes[0, 1],
        current.radius,
        current.radial_velocity,
        reference.radial_velocity,
        r"Radial velocity $v_r$ (cm/s)",
        x_limit,
        velocity_limit,
    )
    render_panel(
        axes[1, 0],
        current.radius,
        current.te_kev,
        reference.te_kev,
        r"Electron temperature $T_e$ (keV)",
        x_limit,
        te_limit,
        log_y=True,
    )
    render_panel(
        axes[1, 1],
        current.radius,
        current.ti_kev,
        reference.ti_kev,
        r"Ion temperature $T_i$ (keV)",
        x_limit,
        ti_limit,
        log_y=True,
    )

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    fig.savefig(output_path, dpi=120)
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
    parser.add_argument("--history-profile", required=True)
    parser.add_argument("--case-json")
    parser.add_argument("--output-dir")
    parser.add_argument("--output-basename")
    parser.add_argument("--title-prefix", default="Spherical Sedov H-only, no-ALE exact reference")
    parser.add_argument("--target-time-interval", type=float, default=2.0e-5)
    parser.add_argument("--x-max", type=float, default=0.12)
    parser.add_argument("--duration-ms", type=int, default=220)
    parser.add_argument("--skip-zero", action="store_true")
    args = parser.parse_args()

    history_profile = Path(args.history_profile)
    case_dir = history_profile.parent
    output_dir = Path(args.output_dir) if args.output_dir else case_dir
    output_basename = args.output_basename if args.output_basename else f"{history_profile.stem}_profiles_with_exact"
    frames_dir = output_dir / "profiles_with_exact_frames"
    frames_dir.mkdir(parents=True, exist_ok=True)

    gamma = 5.0 / 3.0
    rho0 = 1.0
    eblast = 1.0
    ambient_pressure = 1.0e-8
    electron_energy_fraction = 0.5
    if args.case_json:
        config = json.loads(Path(args.case_json).read_text(encoding="utf-8"))
        gamma = float(config["hydro"]["gamma"])
        rho0 = float(config["init"]["rho0"])
        eblast = float(config["init"]["blast_energy"])
        ambient_pressure = float(config["init"]["p0"])
        electron_energy_fraction = float(config["init"].get("electron_energy_fraction", 0.5))

    available_profiles = read_history_profiles(history_profile)
    selected_profiles = select_raw_profiles(
        available_profiles,
        float(args.target_time_interval),
        include_zero=not bool(args.skip_zero),
    )
    initial = available_profiles[0]
    solver = Sedov(geometry=3, gamma=gamma, rho0=rho0, omega=0.0, eblast=eblast)
    references = [
        sedov_reference_profile(
            profile,
            initial,
            solver,
            ambient_pressure,
            electron_energy_fraction,
        )
        for profile in selected_profiles
    ]

    x_limit = (0.0, float(args.x_max))
    selected_masked = [masked_values(profile, x_limit) for profile in selected_profiles]
    reference_masked = [masked_reference_values(profile, x_limit) for profile in references]
    rho_limit = finite_range(
        [values["rho"] for values in selected_masked] + [values["rho"] for values in reference_masked]
    )
    te_limit = positive_range(
        [values["te_kev"] for values in selected_masked] + [values["te_kev"] for values in reference_masked]
    )
    ti_limit = positive_range(
        [values["ti_kev"] for values in selected_masked] + [values["ti_kev"] for values in reference_masked]
    )
    velocity_limit = finite_range(
        [values["radial_velocity"] for values in selected_masked]
        + [values["radial_velocity"] for values in reference_masked]
    )
    rho_limit = (min(0.0, rho_limit[0]), rho_limit[1])
    velocity_limit = (min(0.0, velocity_limit[0]), velocity_limit[1])

    frame_paths: list[Path] = []
    for index, (profile, reference) in enumerate(zip(selected_profiles, references), start=1):
        frame_path = frames_dir / f"frame_{index:04d}_step{profile.step_index:06d}.png"
        render_frame(
            profile,
            reference,
            frame_path,
            args.title_prefix,
            x_limit,
            rho_limit,
            velocity_limit,
            te_limit,
            ti_limit,
        )
        frame_paths.append(frame_path)

    gif_path = output_dir / f"{output_basename}.gif"
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    final_png = output_dir / f"{output_basename}.png"
    render_frame(
        selected_profiles[-1],
        references[-1],
        final_png,
        args.title_prefix,
        x_limit,
        rho_limit,
        velocity_limit,
        te_limit,
        ti_limit,
    )

    summary_path = output_dir / f"{output_basename}_summary.json"
    summary = {
        "case": history_profile.stem,
        "history_file": str(history_profile),
        "frame_count": len(frame_paths),
        "final_time_s": selected_profiles[-1].time_s,
        "shock_radius_cm": references[-1].shock_radius,
        "final_png": str(final_png),
        "gif": str(gif_path),
        "reference": "ExactPack Sedov with two-temperature pressure split",
        "gamma": gamma,
        "rho0": rho0,
        "eblast": eblast,
        "ambient_pressure_floor": ambient_pressure,
        "electron_energy_fraction": electron_energy_fraction,
        "panel_layout": "2x2",
        "panels": ["rho_g_cm3", "vr_cm_s", "Te_keV", "Ti_keV"],
        "legend_labels": ["DEC3D numerical", "Two-temperature Sedov reference"],
        "frames_dir": str(frames_dir),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(gif_path)
    print(final_png)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
