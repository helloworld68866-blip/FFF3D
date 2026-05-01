import argparse
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image


def read_profile(path: Path):
    numeric_rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.replace(",", " ").split()
            try:
                numeric_rows.append([float(part) for part in parts])
            except ValueError:
                continue
    return numeric_rows


def read_exact_profile(path: Path):
    time_s = 0.0
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("# time_s"):
                parts = stripped.split()
                if len(parts) >= 3:
                    time_s = float(parts[2])
                continue
            if stripped.startswith("#"):
                continue
            parts = stripped.replace(",", " ").split()
            if len(parts) < 5:
                continue
            rows.append([float(part) for part in parts[:5]])
    return {
        "time_s": time_s,
        "r_cm": [row[0] for row in rows],
        "te_numeric_keV": [row[1] for row in rows],
        "te_exact_keV": [row[2] for row in rows],
        "ti_numeric_keV": [row[3] for row in rows],
        "ti_exact_keV": [row[4] for row in rows],
    }


def read_exact_slice(path: Path):
    time_s = 0.0
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("# time_s"):
                parts = stripped.split()
                if len(parts) >= 3:
                    time_s = float(parts[2])
                continue
            if stripped.startswith("#"):
                continue
            parts = stripped.replace(",", " ").split()
            if len(parts) < 6:
                continue
            rows.append([float(part) for part in parts[:6]])
    return {
        "time_s": time_s,
        "x_um": [row[0] for row in rows],
        "z_um": [row[1] for row in rows],
        "te_numeric_keV": [row[2] for row in rows],
        "te_exact_keV": [row[3] for row in rows],
        "ti_numeric_keV": [row[4] for row in rows],
        "ti_exact_keV": [row[5] for row in rows],
    }


def bounds(values, pad_fraction=0.06):
    finite = [value for value in values if value == value]
    if not finite:
        return 0.0, 1.0
    lo = min(finite)
    hi = max(finite)
    if lo == hi:
        pad = max(1.0e-12, abs(lo) * pad_fraction)
    else:
        pad = (hi - lo) * pad_fraction
    return lo - pad, hi + pad


def read_key_value_diagnostics(path: Path):
    values = {}
    if not path.exists():
        return values
    text = path.read_text(encoding="utf-8")
    for token in text.replace("\n", ";").split(";"):
        stripped = token.strip()
        if "=" not in stripped:
            continue
        key, value = stripped.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def fixed_temperature_bounds_from_diagnostics(values, bar_key, amp_key, fallback_values):
    try:
        center = float(values[bar_key])
        amplitude = abs(float(values[amp_key]))
    except (KeyError, ValueError):
        return bounds(fallback_values)
    pad = max(0.1 * amplitude, 1.0e-9)
    return center - amplitude - pad, center + amplitude + pad


def radial_axis_bounds_um(r_cm_values):
    finite = sorted({value for value in r_cm_values if value == value})
    if not finite:
        return 0.0, 1.0
    if len(finite) == 1:
        half_width = max(abs(finite[0]) * 0.02, 1.0e-6)
        lo_cm = max(0.0, finite[0] - half_width)
        hi_cm = finite[0] + half_width
    else:
        positive_diffs = [
            finite[i + 1] - finite[i]
            for i in range(len(finite) - 1)
            if finite[i + 1] > finite[i]
        ]
        dr_cm = min(positive_diffs) if positive_diffs else 0.0
        lo_cm = max(0.0, finite[0] - 0.5 * dr_cm)
        hi_cm = finite[-1] + 0.5 * dr_cm
    return lo_cm * 1.0e4, hi_cm * 1.0e4


def write_exact_gif(case_dir: Path, exact_profiles):
    frames = [read_exact_profile(path) for path in exact_profiles]
    frames = [frame for frame in frames if frame["r_cm"]]
    if not frames:
        return False

    all_r = []
    all_te = []
    all_ti = []
    for frame in frames:
        all_r.extend(frame["r_cm"])
        all_te.extend(frame["te_numeric_keV"])
        all_te.extend(frame["te_exact_keV"])
        all_ti.extend(frame["ti_numeric_keV"])
        all_ti.extend(frame["ti_exact_keV"])

    diagnostic_values = read_key_value_diagnostics(case_dir / "benchmark_diagnostics.txt")
    xlim = radial_axis_bounds_um(all_r)
    te_ylim = fixed_temperature_bounds_from_diagnostics(
        diagnostic_values,
        "exact_te_bar_keV",
        "exact_te_amp_keV",
        all_te,
    )
    ti_ylim = fixed_temperature_bounds_from_diagnostics(
        diagnostic_values,
        "exact_ti_bar_keV",
        "exact_ti_amp_keV",
        all_ti,
    )

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.0), sharex=True)
    axes[0].set_ylabel("T_e (keV)")
    axes[1].set_ylabel("T_i (keV)")
    axes[1].set_xlabel("r (um)")
    axes[0].set_xlim(*xlim)
    axes[0].set_ylim(*te_ylim)
    axes[1].set_ylim(*ti_ylim)
    axes[0].set_autoscale_on(False)
    axes[1].set_autoscale_on(False)
    axes[0].grid(alpha=0.25)
    axes[1].grid(alpha=0.25)

    exact_te_line, = axes[0].plot([], [], color="black", linewidth=2.0, label="Exact Te")
    numeric_te_line, = axes[0].plot(
        [], [], color="#d95f02", linestyle="--", linewidth=1.8, label="DEC3D Te"
    )
    exact_ti_line, = axes[1].plot([], [], color="black", linewidth=2.0, label="Exact Ti")
    numeric_ti_line, = axes[1].plot(
        [], [], color="#1b9e77", linestyle="--", linewidth=1.8, label="DEC3D Ti"
    )
    title = fig.suptitle("")
    axes[0].legend(loc="best")
    axes[1].legend(loc="best")

    def update(index):
        frame = frames[index]
        r_um = [value * 1.0e4 for value in frame["r_cm"]]
        exact_te_line.set_data(r_um, frame["te_exact_keV"])
        numeric_te_line.set_data(r_um, frame["te_numeric_keV"])
        exact_ti_line.set_data(r_um, frame["ti_exact_keV"])
        numeric_ti_line.set_data(r_um, frame["ti_numeric_keV"])
        title.set_text(
            "P2 exact l=1,m=1 thermal diffusion, "
            f"t={frame['time_s']:.6e} s"
        )
        return exact_te_line, numeric_te_line, exact_ti_line, numeric_ti_line, title

    update(0)
    fig.tight_layout()
    fig.savefig(case_dir / "exact_te_ti_profiles_final.png", dpi=160)

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=len(frames),
        interval=550,
        blit=False,
        repeat=True,
    )
    anim.save(case_dir / "exact_te_ti_evolution.gif", writer=animation.PillowWriter(fps=2))
    plt.close(fig)
    return True


def write_exact_slice_gif(case_dir: Path, exact_slices):
    frames = [read_exact_slice(path) for path in exact_slices]
    frames = [frame for frame in frames if frame["x_um"]]
    if not frames:
        return False

    diagnostic_values = read_key_value_diagnostics(case_dir / "benchmark_diagnostics.txt")
    all_x = []
    all_z = []
    all_te = []
    all_ti = []
    all_te_error = []
    all_ti_error = []
    for frame in frames:
        all_x.extend(frame["x_um"])
        all_z.extend(frame["z_um"])
        all_te.extend(frame["te_numeric_keV"])
        all_te.extend(frame["te_exact_keV"])
        all_ti.extend(frame["ti_numeric_keV"])
        all_ti.extend(frame["ti_exact_keV"])
        all_te_error.extend(
            numeric - exact
            for numeric, exact in zip(frame["te_numeric_keV"], frame["te_exact_keV"])
        )
        all_ti_error.extend(
            numeric - exact
            for numeric, exact in zip(frame["ti_numeric_keV"], frame["ti_exact_keV"])
        )

    axis_bound = max(max(abs(value) for value in all_x), max(abs(value) for value in all_z))
    te_bounds = fixed_temperature_bounds_from_diagnostics(
        diagnostic_values,
        "exact_te_bar_keV",
        "exact_te_amp_keV",
        all_te,
    )
    ti_bounds = fixed_temperature_bounds_from_diagnostics(
        diagnostic_values,
        "exact_ti_bar_keV",
        "exact_ti_amp_keV",
        all_ti,
    )
    te_error_bound = max(max(abs(value) for value in all_te_error), 1.0e-12)
    ti_error_bound = max(max(abs(value) for value in all_ti_error), 1.0e-12)

    frame_dir = case_dir / "_exact_xz_contour_frames"
    frame_dir.mkdir(exist_ok=True)
    for stale_frame in frame_dir.glob("frame_*.png"):
        stale_frame.unlink()

    te_levels = np.linspace(te_bounds[0], te_bounds[1], 42)
    ti_levels = np.linspace(ti_bounds[0], ti_bounds[1], 42)
    te_error_levels = np.linspace(-te_error_bound, te_error_bound, 42)
    ti_error_levels = np.linspace(-ti_error_bound, ti_error_bound, 42)
    panel_specs = [
        ("DEC3D Te", "te_numeric_keV", te_levels, "inferno", "T_e (keV)"),
        ("Exact Te", "te_exact_keV", te_levels, "inferno", "T_e (keV)"),
        ("Te error", "te_error", te_error_levels, "coolwarm", "num - exact (keV)"),
        ("DEC3D Ti", "ti_numeric_keV", ti_levels, "viridis", "T_i (keV)"),
        ("Exact Ti", "ti_exact_keV", ti_levels, "viridis", "T_i (keV)"),
        ("Ti error", "ti_error", ti_error_levels, "coolwarm", "num - exact (keV)"),
    ]
    png_paths = []
    for index, frame in enumerate(frames):
        x = np.asarray(frame["x_um"])
        z = np.asarray(frame["z_um"])
        triangulation = mtri.Triangulation(x, z)
        values_by_name = {
            "te_numeric_keV": np.asarray(frame["te_numeric_keV"]),
            "te_exact_keV": np.asarray(frame["te_exact_keV"]),
            "te_error": np.asarray(frame["te_numeric_keV"]) - np.asarray(frame["te_exact_keV"]),
            "ti_numeric_keV": np.asarray(frame["ti_numeric_keV"]),
            "ti_exact_keV": np.asarray(frame["ti_exact_keV"]),
            "ti_error": np.asarray(frame["ti_numeric_keV"]) - np.asarray(frame["ti_exact_keV"]),
        }
        fig, axes = plt.subplots(2, 3, figsize=(13.4, 7.8), sharex=True, sharey=True)
        for ax, (title, value_key, levels, cmap, colorbar_label) in zip(axes.flat, panel_specs):
            contour = ax.tricontourf(
                triangulation,
                values_by_name[value_key],
                levels=levels,
                cmap=cmap,
                extend="both",
            )
            ax.set_title(title)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(-axis_bound, axis_bound)
            ax.set_ylim(-axis_bound, axis_bound)
            ax.set_xlabel("x (um)")
            ax.set_ylabel("z (um)")
            ax.set_facecolor("#f7f7f4")
            circle = plt.Circle(
                (0.0, 0.0),
                axis_bound,
                color="white",
                fill=False,
                linewidth=0.8,
                alpha=0.75,
            )
            ax.add_patch(circle)
            colorbar = fig.colorbar(contour, ax=ax, fraction=0.046, pad=0.035)
            colorbar.set_label(colorbar_label)
        fig.suptitle(
            "P2 exact l=1,m=1 x-z slice, "
            f"t={frame['time_s']:.6e} s",
            y=0.995,
        )
        fig.tight_layout()
        png_path = frame_dir / f"frame_{index:06d}.png"
        fig.savefig(png_path, dpi=140)
        if index == len(frames) - 1:
            fig.savefig(case_dir / "exact_xz_te_ti_comparison_final.png", dpi=160)
        plt.close(fig)
        png_paths.append(png_path)

    images = [Image.open(path).convert("P", palette=Image.ADAPTIVE) for path in png_paths]
    images[0].save(
        case_dir / "exact_xz_te_ti_comparison.gif",
        save_all=True,
        append_images=images[1:],
        duration=360,
        loop=0,
    )
    for image in images:
        image.close()
    return True


def write_native_profile_png(case_dir: Path, profile: Path):
    rows = read_profile(profile)
    if not rows:
        return False
    r_cm = [row[0] for row in rows if len(row) >= 4]
    te_keV = [row[2] for row in rows if len(row) >= 4]
    ti_keV = [row[3] for row in rows if len(row) >= 4]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("P2 thermal benchmark profiles")
    ax.set_xlabel("r (cm)")
    ax.set_ylabel("T (keV)")
    if r_cm and len(r_cm) == len(te_keV) == len(ti_keV):
        ax.plot(r_cm, te_keV, label="DEC3D Te")
        ax.plot(r_cm, ti_keV, label="DEC3D Ti")
        ax.legend()
    else:
        ax.text(0.05, 0.5, "native profile artifact", transform=ax.transAxes)
    fig.savefig(case_dir / "te_ti_profiles.png", dpi=160)
    plt.close(fig)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--lilac-reference", default="")
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    summary = case_dir / "postprocess_summary.txt"
    exact_profiles = sorted(case_dir.glob("profile_exact_vs_numeric_step*.txt"))
    exact_slices = sorted(case_dir.glob("slice_xz_exact_vs_numeric_step*.txt"))
    profile = case_dir / "profile_step000001.txt"

    lilac_reference_available = bool(args.lilac_reference)
    parity_claim_allowed = lilac_reference_available
    exact_gif_written = False
    exact_slice_gif_written = False
    native_png_written = False

    if exact_profiles:
        exact_gif_written = write_exact_gif(case_dir, exact_profiles)
        if exact_slices:
            exact_slice_gif_written = write_exact_slice_gif(case_dir, exact_slices)
    elif profile.exists():
        native_png_written = write_native_profile_png(case_dir, profile)
    else:
        summary.write_text(
            "diagnostic_id=p2.thermal_benchmark.postprocess; success=false; "
            "failure_reason=profile_artifact_missing\n",
            encoding="utf-8",
        )
        return 1

    if exact_profiles and not exact_gif_written:
        summary.write_text(
            "diagnostic_id=p2.thermal_benchmark.postprocess; success=false; "
            "failure_reason=exact_profiles_contain_no_numeric_rows\n",
            encoding="utf-8",
        )
        return 1
    if profile.exists() and not exact_profiles and not native_png_written:
        summary.write_text(
            "diagnostic_id=p2.thermal_benchmark.postprocess; success=false; "
            "failure_reason=profile_contains_no_numeric_rows\n",
            encoding="utf-8",
        )
        return 1

    summary.write_text(
        "diagnostic_id=p2.thermal_benchmark.postprocess; success=true; "
        f"lilac_reference_available={str(lilac_reference_available).lower()}; "
        f"parity_claim_allowed={str(parity_claim_allowed).lower()}; "
        f"exact_profiles_detected={str(bool(exact_profiles)).lower()}; "
        f"exact_slices_detected={str(bool(exact_slices)).lower()}; "
        f"exact_evolution_gif_written={str(exact_gif_written).lower()}; "
        f"exact_slice_evolution_gif_written={str(exact_slice_gif_written).lower()}; "
        f"te_ti_profiles_written={str(native_png_written or exact_gif_written).lower()}; "
        f"axis_x_unit={'um' if exact_profiles else 'cm'}; "
        f"fixed_axis_ranges={str(bool(exact_profiles)).lower()}\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
