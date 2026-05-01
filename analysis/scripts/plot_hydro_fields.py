from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_headers(path: Path) -> dict[str, str]:
    headers: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                continue
            text = line[1:].strip()
            if "=" in text:
                key, value = text.split("=", 1)
                headers[key.strip()] = value.strip()
    return headers


def parse_float_list(text: str) -> list[float]:
    if not text:
        return []
    return [float(item) for item in text.split(",") if item]


def read_numeric_rows(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            rows.append(stripped.split())
    return rows


def write_text(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def compute_pchip_derivatives(x: list[float], y: list[float]) -> list[float]:
    count = len(x)
    if count < 2:
        return [0.0 for _ in x]
    if count == 2:
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return [slope, slope]

    h = [x[index + 1] - x[index] for index in range(count - 1)]
    delta = [(y[index + 1] - y[index]) / h[index] for index in range(count - 1)]
    derivatives = [0.0 for _ in x]

    for index in range(1, count - 1):
        if delta[index - 1] == 0.0 or delta[index] == 0.0 or sign(delta[index - 1]) != sign(delta[index]):
            derivatives[index] = 0.0
        else:
            w1 = 2.0 * h[index] + h[index - 1]
            w2 = h[index] + 2.0 * h[index - 1]
            derivatives[index] = (w1 + w2) / ((w1 / delta[index - 1]) + (w2 / delta[index]))

    d0 = ((2.0 * h[0] + h[1]) * delta[0] - h[0] * delta[1]) / (h[0] + h[1])
    if sign(d0) != sign(delta[0]):
        d0 = 0.0
    elif sign(delta[0]) != sign(delta[1]) and abs(d0) > abs(3.0 * delta[0]):
        d0 = 3.0 * delta[0]
    derivatives[0] = d0

    dn = ((2.0 * h[-1] + h[-2]) * delta[-1] - h[-1] * delta[-2]) / (h[-1] + h[-2])
    if sign(dn) != sign(delta[-1]):
        dn = 0.0
    elif sign(delta[-1]) != sign(delta[-2]) and abs(dn) > abs(3.0 * delta[-1]):
        dn = 3.0 * delta[-1]
    derivatives[-1] = dn

    return derivatives


def hermite_derivative(
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    m0: float,
    m1: float,
    x_value: float,
) -> float:
    interval = x1 - x0
    if interval <= 0.0:
        return 0.0

    t = (x_value - x0) / interval
    dh00 = 6.0 * t * t - 6.0 * t
    dh10 = 3.0 * t * t - 4.0 * t + 1.0
    dh01 = -6.0 * t * t + 6.0 * t
    dh11 = 3.0 * t * t - 2.0 * t
    return (
        dh00 * y0 / interval
        + dh10 * m0
        + dh01 * y1 / interval
        + dh11 * m1
    )


@dataclass(frozen=True)
class ShockProfileFeatures:
    shock_radius: float
    inner_rise_radius: float
    density_peak_radius: float


def detect_shock_features_from_profile(radial_centers: list[float], rho_values: list[float]) -> ShockProfileFeatures:
    if len(radial_centers) < 3 or len(radial_centers) != len(rho_values):
        return ShockProfileFeatures(0.0, 0.0, 0.0)

    derivatives = compute_pchip_derivatives(radial_centers, rho_values)
    outer_shock_radius = radial_centers[0]
    inner_rise_radius = radial_centers[0]
    density_peak_radius = radial_centers[0]
    strongest_outer_fall = float("inf")
    strongest_inner_rise = float("-inf")
    peak_density = float("-inf")
    samples_per_interval = 64

    for radius, density in zip(radial_centers, rho_values):
        if density > peak_density:
            peak_density = density
            density_peak_radius = radius

    for index in range(len(radial_centers) - 1):
        x0 = radial_centers[index]
        x1 = radial_centers[index + 1]
        y0 = rho_values[index]
        y1 = rho_values[index + 1]
        m0 = derivatives[index]
        m1 = derivatives[index + 1]
        interval = x1 - x0
        if interval <= 0.0:
            continue

        for sample_index in range(samples_per_interval + 1):
            fraction = sample_index / samples_per_interval
            radius = x0 + fraction * interval
            gradient = hermite_derivative(x0, x1, y0, y1, m0, m1, radius)
            if gradient < strongest_outer_fall:
                strongest_outer_fall = gradient
                outer_shock_radius = radius
            if gradient > strongest_inner_rise:
                strongest_inner_rise = gradient
                inner_rise_radius = radius

    return ShockProfileFeatures(outer_shock_radius, inner_rise_radius, density_peak_radius)


def recompute_shock_radius_from_profile(radial_centers: list[float], rho_values: list[float]) -> float:
    return detect_shock_features_from_profile(radial_centers, rho_values).shock_radius


def plot_shell_map(path: Path, output_path: Path) -> None:
    headers = read_headers(path)
    rows = read_numeric_rows(path)
    theta_cells = int(headers["theta_cells"])
    phi_cells = int(headers["phi_cells"])
    data = [[0.0 for _ in range(phi_cells)] for _ in range(theta_cells)]

    for row in rows:
        theta = int(row[0])
        phi = int(row[1])
        value = float(row[4])
        data[theta][phi] = value

    fig, ax = plt.subplots(figsize=(7, 4.5))
    image = ax.imshow(data, origin="lower", aspect="auto")
    ax.set_title(f"{headers.get('field', 'field')} @ {headers.get('time', 't?')}")
    ax.set_xlabel("phi index")
    ax.set_ylabel("theta index")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_mode_projection(path: Path, output_path: Path) -> None:
    rows = read_numeric_rows(path)
    labels = [row[0] for row in rows]
    amplitudes = [float(row[1]) for row in rows]
    phases = [float(row[2]) for row in rows]

    fig, (ax_amp, ax_phase) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
    ax_amp.plot(labels, amplitudes, marker="o")
    ax_amp.set_ylabel("amplitude")
    ax_amp.set_title("mode projection")
    ax_phase.plot(labels, phases, marker="o")
    ax_phase.set_ylabel("phase")
    ax_phase.set_xlabel("snapshot")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_shock_radius_history(path: Path, output_path: Path) -> None:
    rows = read_numeric_rows(path)
    time_s = [float(row[0]) for row in rows]
    numerical = [float(row[1]) for row in rows]
    reference = [float(row[2]) for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(time_s, numerical, marker="o", label="numerical")
    ax.plot(time_s, reference, marker="s", label="reference")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("outer shock radius")
    ax.set_title("outer shock radius vs time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_radial_profile(path: Path, output_path: Path) -> None:
    rows = read_numeric_rows(path)
    radial_center = [float(row[1]) for row in rows]
    rho = [float(row[2]) for row in rows]
    te = [float(row[3]) for row in rows]
    velocity = [float(row[4]) for row in rows]
    mom_r = [float(row[5]) for row in rows] if rows and len(rows[0]) > 5 else []
    e_fluid_total = [float(row[6]) for row in rows] if rows and len(rows[0]) > 6 else []

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(radial_center, rho, label="rho")
    ax.plot(radial_center, te, label="hydro-only Te")
    ax.plot(radial_center, velocity, label="|v|")
    if mom_r:
        ax.plot(radial_center, mom_r, label="mom_r")
    if e_fluid_total:
        ax.plot(radial_center, e_fluid_total, label="E_fluid_total")
    ax.set_xlabel("radial center")
    ax.set_title(path.stem)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_radial_time_map(path: Path, output_path: Path) -> None:
    headers = read_headers(path)
    rows = read_numeric_rows(path)
    radial_centers = parse_float_list(headers.get("radial_centers", ""))
    if not rows or not radial_centers:
        return

    time_s = [float(row[0]) for row in rows]
    reference_shock = [float(row[1]) for row in rows]
    values = [[float(value) for value in row[3:]] for row in rows]
    if not values or len(values[0]) != len(radial_centers):
        return

    fig, ax = plt.subplots(figsize=(8, 4.8))
    image = ax.imshow(
        values,
        origin="lower",
        aspect="auto",
        extent=[radial_centers[0], radial_centers[-1], time_s[0], time_s[-1]],
    )
    ax.plot(reference_shock, time_s, color="white", linewidth=1.5, label="analytic shock radius")
    ax.set_xlabel("radius")
    ax.set_ylabel("time (s)")
    ax.set_title(f"{headers.get('field', 'field')} r-t map")
    ax.legend(loc="upper right")
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def synthesize_radial_time_maps(case_dir: Path) -> list[Path]:
    profile_paths = sorted(case_dir.glob("radial_profile_step*.txt"))
    if not profile_paths:
        return []

    samples: list[dict[str, object]] = []
    for path in profile_paths:
        headers = read_headers(path)
        rows = read_numeric_rows(path)
        if not rows:
            continue
        radial_centers = [float(row[1]) for row in rows]
        samples.append(
            {
                "time_s": float(headers.get("time_s", "0")),
                "reference_shock_radius": float(headers.get("reference_shock_radius", "0")),
                "radial_centers": radial_centers,
                "rho": [float(row[2]) for row in rows],
                "hydro_only_te": [float(row[3]) for row in rows],
                "velocity_magnitude": [float(row[4]) for row in rows],
                "mom_r": [float(row[5]) for row in rows],
                "e_fluid_total": [float(row[6]) for row in rows],
            }
        )

    if not samples:
        return []

    for sample in samples:
        features = detect_shock_features_from_profile(
            sample["radial_centers"],
            sample["rho"],
        )
        sample["numerical_shock_radius"] = features.shock_radius
        sample["inner_rise_radius"] = features.inner_rise_radius
        sample["density_peak_radius"] = features.density_peak_radius

    radial_centers = samples[0]["radial_centers"]
    field_map = {
        "rho_rt_map.txt": "rho",
        "hydro_only_te_rt_map.txt": "hydro_only_te",
        "velocity_magnitude_rt_map.txt": "velocity_magnitude",
        "mom_r_rt_map.txt": "mom_r",
        "e_fluid_total_rt_map.txt": "e_fluid_total",
    }

    generated: list[Path] = []
    radial_centers_text = ",".join(f"{float(value):.17g}" for value in radial_centers)
    for file_name, field_name in field_map.items():
        output_path = case_dir / file_name
        lines = [
            "# kind=radial_time_map",
            f"# field={field_name}",
            f"# sample_count={len(samples)}",
            f"# radial_cells={len(radial_centers)}",
            f"# radial_centers={radial_centers_text}",
            "# columns=time_s reference_shock_radius numerical_shock_radius values...",
        ]
        for sample in samples:
            values = sample[field_name]
            row = [
                f"{float(sample['time_s']):.17g}",
                f"{float(sample['reference_shock_radius']):.17g}",
                f"{float(sample['numerical_shock_radius']):.17g}",
            ]
            row.extend(f"{float(value):.17g}" for value in values)
            lines.append(" ".join(row))
        write_text(output_path, "\n".join(lines) + "\n")
        generated.append(output_path)

    shock_lines = [
        "# kind=shock_radius_history",
        "# source=profile_recomputed",
        "# definition=numerical_radius is outer shock front from strongest negative density gradient",
        "# columns=time_s numerical_radius reference_radius relative_error inner_rise_radius density_peak_radius",
    ]
    for sample in samples:
        reference_radius = float(sample["reference_shock_radius"])
        numerical_radius = float(sample["numerical_shock_radius"])
        relative_error = (
            abs(numerical_radius - reference_radius) / reference_radius
            if reference_radius > 0.0
            else 0.0
        )
        shock_lines.append(
            " ".join(
                [
                    f"{float(sample['time_s']):.17g}",
                    f"{numerical_radius:.17g}",
                    f"{reference_radius:.17g}",
                    f"{relative_error:.17g}",
                    f"{float(sample['inner_rise_radius']):.17g}",
                    f"{float(sample['density_peak_radius']):.17g}",
                ]
            )
        )
    recomputed_history_path = case_dir / "shock_radius_vs_time_recomputed.txt"
    write_text(recomputed_history_path, "\n".join(shock_lines) + "\n")
    generated.append(recomputed_history_path)

    shock_feature_path = case_dir / "shock_feature_radii_vs_time.txt"
    write_text(shock_feature_path, "\n".join(shock_lines) + "\n")
    generated.append(shock_feature_path)

    return generated


def plot_polar_slice(path: Path, output_path: Path) -> None:
    rows = read_numeric_rows(path)
    phi = [float(row[1]) for row in rows]
    rho = [float(row[2]) for row in rows]
    te = [float(row[3]) for row in rows]
    velocity = [float(row[4]) for row in rows]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(phi, rho, label="rho")
    ax.plot(phi, te, label="hydro-only Te")
    ax.plot(phi, velocity, label="|v|")
    ax.set_xlabel("phi center")
    ax.set_title(path.stem)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_key_value_bars(path: Path, output_path: Path, title: str) -> None:
    rows = read_numeric_rows(path)
    names: list[str] = []
    values: list[float] = []
    with path.open("r", encoding="utf-8") as handle:
      for line in handle:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or "=" not in stripped:
          continue
        key, value = stripped.split("=", 1)
        try:
          numeric = float(value)
        except ValueError:
          continue
        if math.isfinite(numeric):
          names.append(key)
          values.append(numeric)

    if not names:
      return

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(range(len(names)), values)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=30, ha="right")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_case_directory(case_dir: Path) -> list[Path]:
    generated: list[Path] = []
    synthesize_radial_time_maps(case_dir)
    for path in sorted(case_dir.glob("*.txt")):
        headers = read_headers(path)
        kind = headers.get("kind", "")
        output_path = path.with_suffix(".png")
        if kind == "shell_map":
            plot_shell_map(path, output_path)
        elif kind == "shock_radius_history":
            plot_shock_radius_history(path, output_path)
        elif kind == "radial_profile":
            plot_radial_profile(path, output_path)
        elif kind == "radial_time_map":
            plot_radial_time_map(path, output_path)
        elif kind == "mode_projection":
            plot_mode_projection(path, output_path)
        elif kind == "polar_slice":
            plot_polar_slice(path, output_path)
        elif kind == "seam_diagnostics":
            plot_key_value_bars(path, output_path, "seam / polar diagnostics")
        elif kind == "radial_only_comparison":
            plot_key_value_bars(path, output_path, "radial-only comparison")
        elif kind == "budget_residual":
            plot_key_value_bars(path, output_path, "budget residual")
        elif kind == "budget_residual_comparison":
            plot_key_value_bars(path, output_path, "budget residual comparison")
        elif kind == "mpi_parity_diagnostics":
            plot_key_value_bars(path, output_path, "mpi parity diagnostics")
        elif kind == "moving_mesh_diagnostics":
            plot_key_value_bars(path, output_path, "moving mesh diagnostics")
        else:
            continue
        generated.append(output_path)
    return generated


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", required=True, help="Case output directory under analysis/output")
    args = parser.parse_args()

    case_dir = Path(args.case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    generated = plot_case_directory(case_dir)
    for item in generated:
        print(item)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
