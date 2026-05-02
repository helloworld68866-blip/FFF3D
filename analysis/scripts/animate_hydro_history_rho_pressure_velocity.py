from __future__ import annotations

import argparse
import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.noh.noh1 import Noh
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
    radius: Iterable[float]
    rho: Iterable[float]
    te_kev: Iterable[float]
    ti_kev: Iterable[float]
    radial_velocity: Iterable[float]


@dataclass(frozen=True)
class ReferenceProfile:
    radius: Iterable[float]
    rho: Iterable[float]
    pressure: Iterable[float]
    radial_velocity: Iterable[float]
    shock_radius: float


@dataclass(frozen=True)
class CaseSettings:
    gamma: float
    rho0: float
    p0: float
    vr0: float
    eblast: float


def as_array(values: Iterable[float]) -> np.ndarray:
    return np.asarray(list(values) if not isinstance(values, np.ndarray) else values, dtype=float)


def compute_total_pressure(
    rho: Iterable[float],
    te_kev: Iterable[float],
    ti_kev: Iterable[float],
    *,
    ion_mass_g: float = MEAN_DT_ION_MASS_G,
    zbar: float = DEFAULT_ZBAR,
) -> np.ndarray:
    rho_array = as_array(rho)
    te_array = as_array(te_kev)
    ti_array = as_array(ti_kev)
    ion_number_density = rho_array / ion_mass_g
    return ion_number_density * (zbar * te_array + ti_array) * ERG_PER_KEV


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


def nonnegative_limit(limit: tuple[float, float]) -> tuple[float, float]:
    upper = limit[1] if limit[1] > 0.0 else 1.0
    return 0.0, upper


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

    radius = as_array(profiles[0].radius)
    for profile in profiles[1:]:
        profile_radius = as_array(profile.radius)
        if profile_radius.shape != radius.shape or not np.allclose(profile_radius, radius):
            raise ValueError("history profile radial grids are not consistent across time")
    return profiles


def select_history_profiles(
    profiles: list[HistoryProfile],
    target_time_interval_s: float,
    *,
    include_zero: bool = True,
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

    final_profile = profiles[-1]
    if final_profile.step_index not in used_steps:
        selected.append(final_profile)

    if len(selected) < 2:
        raise ValueError("history profile selection produced fewer than two frames")
    selected.sort(key=lambda profile: (profile.time_s, profile.step_index))
    return selected


def read_case_settings(case_json: Path | None, defaults: CaseSettings) -> CaseSettings:
    if case_json is None:
        return defaults

    config = json.loads(case_json.read_text(encoding="utf-8"))
    hydro = config.get("hydro", {})
    init = config.get("init", {})
    return CaseSettings(
        gamma=float(hydro.get("gamma", defaults.gamma)),
        rho0=float(init.get("rho0", defaults.rho0)),
        p0=float(init.get("p0", defaults.p0)),
        vr0=float(init.get("vr0", init.get("v_r0", defaults.vr0))),
        eblast=float(init.get("blast_energy", defaults.eblast)),
    )


def initial_reference(profile: HistoryProfile, pressure_value: float | None = None) -> ReferenceProfile:
    radius = as_array(profile.radius)
    if pressure_value is None:
        pressure = compute_total_pressure(profile.rho, profile.te_kev, profile.ti_kev)
    else:
        pressure = np.full_like(radius, pressure_value, dtype=float)
    return ReferenceProfile(
        radius=radius,
        rho=as_array(profile.rho),
        pressure=pressure,
        radial_velocity=as_array(profile.radial_velocity),
        shock_radius=0.0,
    )


def sedov_reference_profile(
    profile: HistoryProfile,
    initial: HistoryProfile,
    solver: Sedov,
    ambient_pressure: float,
) -> ReferenceProfile:
    if profile.time_s <= 0.0:
        return initial_reference(initial, ambient_pressure)

    radius = as_array(profile.radius)
    exact = solver(radius, profile.time_s)
    rho = np.asarray(exact["density"], dtype=float)
    pressure = np.asarray(exact["pressure"], dtype=float)
    radial_velocity = np.asarray(exact["velocity"], dtype=float)
    shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")

    rho = np.where(np.isfinite(rho) & (rho > 0.0), rho, as_array(initial.rho))
    pressure = np.where(np.isfinite(pressure), pressure, ambient_pressure)
    radial_velocity = np.where(np.isfinite(radial_velocity), radial_velocity, 0.0)

    return ReferenceProfile(
        radius=radius,
        rho=rho,
        pressure=pressure,
        radial_velocity=radial_velocity,
        shock_radius=shock_radius,
    )


def noh_reference_profile(
    profile: HistoryProfile,
    initial: HistoryProfile,
    solver: Noh,
    ambient_pressure: float,
) -> ReferenceProfile:
    if profile.time_s <= 0.0:
        return initial_reference(initial, ambient_pressure)

    radius = as_array(profile.radius)
    exact = solver(radius, profile.time_s)
    rho = np.asarray(exact["density"], dtype=float)
    pressure = np.asarray(exact["pressure"], dtype=float)
    radial_velocity = np.asarray(exact["velocity"], dtype=float)
    shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")

    rho = np.where(np.isfinite(rho) & (rho > 0.0), rho, as_array(initial.rho))
    pressure = np.where(np.isfinite(pressure), pressure, ambient_pressure)
    radial_velocity = np.where(np.isfinite(radial_velocity), radial_velocity, as_array(initial.radial_velocity))

    return ReferenceProfile(
        radius=radius,
        rho=rho,
        pressure=pressure,
        radial_velocity=radial_velocity,
        shock_radius=shock_radius,
    )


def masked(array: Iterable[float], radius: Iterable[float], x_limit: tuple[float, float]) -> np.ndarray:
    radius_array = as_array(radius)
    values = as_array(array)
    mask = (radius_array >= x_limit[0]) & (radius_array <= x_limit[1])
    return values[mask]


def compute_plot_limits(
    frames: list[HistoryProfile],
    references: list[ReferenceProfile],
    x_limit: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
    rho_values: list[np.ndarray] = []
    pressure_values: list[np.ndarray] = []
    velocity_values: list[np.ndarray] = []
    has_evolved_frame = any(frame.time_s > 0.0 for frame in frames)

    for frame, reference in zip(frames, references):
        radius = as_array(frame.radius)
        mask = (radius >= x_limit[0]) & (radius <= x_limit[1])
        rho_values.extend([as_array(frame.rho)[mask], as_array(reference.rho)[mask]])
        if not has_evolved_frame or frame.time_s > 0.0:
            pressure_values.extend(
                [
                    compute_total_pressure(frame.rho, frame.te_kev, frame.ti_kev)[mask],
                    as_array(reference.pressure)[mask],
                ]
            )
        velocity_values.extend([as_array(frame.radial_velocity)[mask], as_array(reference.radial_velocity)[mask]])

    rho_limit = nonnegative_limit(finite_range(rho_values))
    pressure_limit = nonnegative_limit(finite_range(pressure_values))
    velocity_limit = finite_range(velocity_values)
    return rho_limit, pressure_limit, velocity_limit


def render_panel(
    ax: plt.Axes,
    radius: Iterable[float],
    numerical: Iterable[float],
    reference: Iterable[float],
    ylabel: str,
    x_limit: tuple[float, float],
    y_limit: tuple[float, float],
    reference_label: str,
    *,
    show_legend: bool = False,
) -> None:
    ax.plot(as_array(radius), as_array(numerical), color="#d62728", linewidth=2.0, label="DEC3D numerical")
    ax.plot(
        as_array(radius),
        as_array(reference),
        color="black",
        linestyle="--",
        linewidth=1.8,
        label=reference_label,
    )
    ax.set_xlim(*x_limit)
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
    reference_label: str,
    x_limit: tuple[float, float],
    rho_limit: tuple[float, float],
    pressure_limit: tuple[float, float],
    velocity_limit: tuple[float, float],
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.6))
    time_text = f"{current.time_s:.3e}" if current.time_s < 1.0e-2 else f"{current.time_s:.3f}"
    fig.suptitle(f"{title_prefix}, t={time_text} s", fontsize=15, y=0.98)

    pressure = compute_total_pressure(current.rho, current.te_kev, current.ti_kev)
    render_panel(
        axes[0],
        current.radius,
        current.rho,
        reference.rho,
        r"Density $\rho$ (g/cm$^3$)",
        x_limit,
        rho_limit,
        reference_label,
        show_legend=True,
    )
    render_panel(
        axes[1],
        current.radius,
        pressure,
        reference.pressure,
        r"Pressure $P$ (erg/cm$^3$)",
        x_limit,
        pressure_limit,
        reference_label,
    )
    render_panel(
        axes[2],
        current.radius,
        current.radial_velocity,
        reference.radial_velocity,
        r"Radial velocity $v_r$ (cm/s)",
        x_limit,
        velocity_limit,
        reference_label,
    )

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
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


def build_references(
    case_kind: str,
    frames: list[HistoryProfile],
    settings: CaseSettings,
) -> tuple[list[ReferenceProfile], str]:
    initial = frames[0]
    if case_kind == "sedov":
        solver = Sedov(geometry=3, gamma=settings.gamma, rho0=settings.rho0, omega=0.0, eblast=settings.eblast)
        return (
            [sedov_reference_profile(frame, initial, solver, settings.p0) for frame in frames],
            "ExactPack Sedov",
        )
    if case_kind == "noh":
        solver = Noh(geometry=3, gamma=settings.gamma, rho0=settings.rho0, u0=settings.vr0)
        return (
            [noh_reference_profile(frame, initial, solver, settings.p0) for frame in frames],
            "ExactPack Noh",
        )
    raise ValueError(f"unsupported case kind: {case_kind}")


def clean_frame_dir(frames_dir: Path) -> None:
    frames_dir.mkdir(parents=True, exist_ok=True)
    for path in frames_dir.glob("*.png"):
        path.unlink()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-kind", choices=["sedov", "noh"], required=True)
    parser.add_argument("--history-profile", required=True)
    parser.add_argument("--case-json")
    parser.add_argument("--output-dir")
    parser.add_argument("--output-basename")
    parser.add_argument("--title-prefix")
    parser.add_argument("--target-time-interval", type=float, default=0.0)
    parser.add_argument("--x-min", type=float, default=0.0)
    parser.add_argument("--x-max", type=float)
    parser.add_argument("--duration-ms", type=int, default=220)
    parser.add_argument("--limits-from", choices=["final", "all-evolved", "all"], default="final")
    parser.add_argument("--gamma", type=float, default=5.0 / 3.0)
    parser.add_argument("--rho0", type=float, default=1.0)
    parser.add_argument("--p0", type=float, default=1.0e-8)
    parser.add_argument("--vr0", type=float, default=-1.0)
    parser.add_argument("--eblast", type=float, default=1.0)
    parser.add_argument("--reference-label")
    args = parser.parse_args()

    history_profile = Path(args.history_profile)
    output_dir = Path(args.output_dir) if args.output_dir else history_profile.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    output_basename = args.output_basename or f"{history_profile.stem}_rho_pressure_velocity_with_exact"
    frames_dir = output_dir / "rho_pressure_velocity_frames"
    clean_frame_dir(frames_dir)

    defaults = CaseSettings(
        gamma=float(args.gamma),
        rho0=float(args.rho0),
        p0=float(args.p0),
        vr0=float(args.vr0),
        eblast=float(args.eblast),
    )
    settings = read_case_settings(Path(args.case_json) if args.case_json else None, defaults)

    profiles = read_history_profiles(history_profile)
    frames = select_history_profiles(profiles, float(args.target_time_interval), include_zero=True)
    references, default_reference_label = build_references(args.case_kind, frames, settings)
    reference_label = args.reference_label or default_reference_label

    radius = as_array(frames[0].radius)
    x_limit = (float(args.x_min), float(args.x_max) if args.x_max is not None else float(np.max(radius)))
    if args.limits_from == "final":
        limit_frames = [frames[-1]]
        limit_references = [references[-1]]
    elif args.limits_from == "all-evolved":
        limit_pairs = [(frame, reference) for frame, reference in zip(frames, references) if frame.time_s > 0.0]
        if not limit_pairs:
            limit_pairs = list(zip(frames, references))
        limit_frames = [frame for frame, _ in limit_pairs]
        limit_references = [reference for _, reference in limit_pairs]
    else:
        limit_frames = frames
        limit_references = references
    rho_limit, pressure_limit, velocity_limit = compute_plot_limits(limit_frames, limit_references, x_limit)
    title_prefix = args.title_prefix or (
        "Spherical Sedov H-only, no-ALE exact reference"
        if args.case_kind == "sedov"
        else "Spherical Noh H-only, exact inflow fixed"
    )

    frame_paths: list[Path] = []
    for frame_index, (frame, reference) in enumerate(zip(frames, references), start=1):
        frame_path = frames_dir / (
            f"frame_{frame_index:03d}_step{frame.step_index:06d}_t{compact_decimal(frame.time_s)}.png"
        )
        render_frame(
            current=frame,
            reference=reference,
            output_path=frame_path,
            title_prefix=title_prefix,
            reference_label=reference_label,
            x_limit=x_limit,
            rho_limit=rho_limit,
            pressure_limit=pressure_limit,
            velocity_limit=velocity_limit,
        )
        frame_paths.append(frame_path)

    png_path = output_dir / f"{output_basename}.png"
    gif_path = output_dir / f"{output_basename}.gif"
    summary_path = output_dir / f"{output_basename}_summary.json"
    shutil.copyfile(frame_paths[-1], png_path)
    write_gif(frame_paths, gif_path, int(args.duration_ms))

    summary = {
        "case": history_profile.stem,
        "case_kind": args.case_kind,
        "history_file": str(history_profile),
        "frame_count": len(frame_paths),
        "final_time_s": frames[-1].time_s,
        "shock_radius_cm": references[-1].shock_radius,
        "final_png": str(png_path),
        "gif": str(gif_path),
        "panels": ["rho_g_cm3", "pressure_erg_cm3", "vr_cm_s"],
        "reference_label": reference_label,
        "limits_from": args.limits_from,
        "x_limit_cm": [x_limit[0], x_limit[1]],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
