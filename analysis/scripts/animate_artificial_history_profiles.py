from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


@dataclass(frozen=True)
class HistoryProfile:
    step: int
    time_s: float
    r_um: np.ndarray
    rho: np.ndarray
    te_kev: np.ndarray
    ti_kev: np.ndarray
    vr_1e7: np.ndarray
    vtheta_1e7: np.ndarray
    vphi_1e7: np.ndarray


def finite_limit(values: list[np.ndarray], padding_fraction: float = 0.06) -> tuple[float, float]:
    finite = [value[np.isfinite(value)] for value in values if value.size > 0]
    merged = np.concatenate(finite) if finite else np.array([], dtype=float)
    if merged.size == 0:
        return 0.0, 1.0
    lo = float(np.min(merged))
    hi = float(np.max(merged))
    if lo == hi:
        pad = max(abs(lo), 1.0) * padding_fraction
        return lo - pad, hi + pad
    pad = (hi - lo) * padding_fraction
    return lo - pad, hi + pad


def symmetric_limit(values: list[np.ndarray], padding_fraction: float = 0.08) -> tuple[float, float]:
    finite = [value[np.isfinite(value)] for value in values if value.size > 0]
    merged = np.concatenate(finite) if finite else np.array([], dtype=float)
    if merged.size == 0:
        return -1.0, 1.0
    bound = float(np.max(np.abs(merged)))
    if bound <= 0.0 or not np.isfinite(bound):
        bound = 1.0
    bound *= 1.0 + padding_fraction
    return -bound, bound


def read_profiles(path: Path) -> list[HistoryProfile]:
    grouped: dict[tuple[int, float], list[dict[str, str]]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if line.strip() and not line.startswith("#"))
        reader = csv.DictReader(rows)
        required = {
            "step",
            "time_s",
            "r_index",
            "r_um",
            "rho_g_cm3",
            "Te_keV",
            "Ti_keV",
            "vr_cm_s",
            "vtheta_cm_s",
            "vphi_cm_s",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"history profile missing required columns: {path}")
        for row in reader:
            grouped.setdefault((int(row["step"]), float(row["time_s"])), []).append(row)

    profiles: list[HistoryProfile] = []
    for (step, time_s), rows in grouped.items():
        rows.sort(key=lambda item: int(item["r_index"]))
        profiles.append(
            HistoryProfile(
                step=step,
                time_s=time_s,
                r_um=np.array([float(row["r_um"]) for row in rows], dtype=float),
                rho=np.array([float(row["rho_g_cm3"]) for row in rows], dtype=float),
                te_kev=np.array([float(row["Te_keV"]) for row in rows], dtype=float),
                ti_kev=np.array([float(row["Ti_keV"]) for row in rows], dtype=float),
                vr_1e7=np.array([float(row["vr_cm_s"]) * 1.0e-7 for row in rows], dtype=float),
                vtheta_1e7=np.array([float(row["vtheta_cm_s"]) * 1.0e-7 for row in rows], dtype=float),
                vphi_1e7=np.array([float(row["vphi_cm_s"]) * 1.0e-7 for row in rows], dtype=float),
            )
        )
    profiles.sort(key=lambda profile: (profile.time_s, profile.step))
    if len(profiles) < 2:
        raise ValueError("need at least two history profiles for animation")
    return profiles


def select_profiles(profiles: list[HistoryProfile], interval_ps: float) -> list[HistoryProfile]:
    if interval_ps <= 0.0:
        return profiles
    times_ps = np.array([profile.time_s * 1.0e12 for profile in profiles], dtype=float)
    target_times = np.arange(0.0, times_ps[-1] + 0.5 * interval_ps, interval_ps)
    selected: list[HistoryProfile] = []
    used_steps: set[int] = set()
    for target_time in target_times:
        index = int(np.argmin(np.abs(times_ps - target_time)))
        profile = profiles[index]
        if profile.step in used_steps:
            continue
        selected.append(profile)
        used_steps.add(profile.step)
    if profiles[-1].step not in used_steps:
        selected.append(profiles[-1])
    selected.sort(key=lambda profile: (profile.time_s, profile.step))
    return selected


def render_frame(
    profile: HistoryProfile,
    output_path: Path,
    x_limit: tuple[float, float],
    limits: dict[str, tuple[float, float]],
) -> None:
    panels = [
        ("rho", "Density", "rho (g/cm^3)", profile.rho, "#1f77b4"),
        ("te", "Electron temperature", "Te (keV)", profile.te_kev, "#d62728"),
        ("ti", "Ion temperature", "Ti (keV)", profile.ti_kev, "#ff7f0e"),
        ("vr", "Radial velocity", "v_r (10^7 cm/s)", profile.vr_1e7, "#2ca02c"),
        ("vtheta", "Polar velocity", "v_theta (10^7 cm/s)", profile.vtheta_1e7, "#9467bd"),
        ("vphi", "Azimuthal velocity", "v_phi (10^7 cm/s)", profile.vphi_1e7, "#8c564b"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.4), constrained_layout=True)
    for ax, (key, title, ylabel, values, color) in zip(axes.ravel(), panels):
        ax.plot(profile.r_um, values, color=color, linewidth=2.0)
        ax.set_xlim(*x_limit)
        ax.set_ylim(*limits[key])
        ax.set_xlabel("r (um)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, alpha=0.25)
    fig.suptitle(f"Radial profile evolution, step {profile.step}, t = {profile.time_s * 1.0e12:.3f} ps")
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def make_gif(frame_paths: list[Path], output_path: Path, duration_ms: int) -> None:
    images: list[Image.Image] = []
    for frame_path in frame_paths:
        with Image.open(frame_path) as image:
            images.append(image.convert("P", palette=Image.Palette.ADAPTIVE).copy())
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--history", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--frame-interval-ps", type=float, default=1.0)
    parser.add_argument("--duration-ms", type=int, default=180)
    parser.add_argument("--keep-frames", action="store_true")
    args = parser.parse_args()

    history_path = Path(args.history)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profiles = read_profiles(history_path)
    selected = select_profiles(profiles, args.frame_interval_ps)

    x_limit = finite_limit([profile.r_um for profile in selected], padding_fraction=0.0)
    limits = {
        "rho": finite_limit([profile.rho for profile in selected]),
        "te": finite_limit([profile.te_kev for profile in selected]),
        "ti": finite_limit([profile.ti_kev for profile in selected]),
        "vr": finite_limit([profile.vr_1e7 for profile in selected]),
        "vtheta": symmetric_limit([profile.vtheta_1e7 for profile in selected]),
        "vphi": symmetric_limit([profile.vphi_1e7 for profile in selected]),
    }

    frame_dir = output_path.with_suffix("")
    frame_dir = frame_dir.parent / f"{frame_dir.name}_frames"
    if frame_dir.exists():
        shutil.rmtree(frame_dir)
    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[Path] = []
    for index, profile in enumerate(selected):
        frame_path = frame_dir / f"frame_{index:04d}_step_{profile.step:06d}.png"
        render_frame(profile, frame_path, x_limit, limits)
        frame_paths.append(frame_path)

    make_gif(frame_paths, output_path, args.duration_ms)
    metadata = (
        f"history={history_path}\n"
        f"output={output_path}\n"
        f"raw_profile_count={len(profiles)}\n"
        f"frame_count={len(selected)}\n"
        f"first_time_ps={selected[0].time_s * 1.0e12:.12g}\n"
        f"last_time_ps={selected[-1].time_s * 1.0e12:.12g}\n"
        f"frame_interval_ps={args.frame_interval_ps}\n"
    )
    output_path.with_suffix(".txt").write_text(metadata, encoding="utf-8")
    if not args.keep_frames:
        shutil.rmtree(frame_dir)
    print(metadata, end="")


if __name__ == "__main__":
    main()
