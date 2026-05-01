from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.noh.noh1 import Noh

from plot_hydro_fields import read_headers, read_numeric_rows


EXACTPACK_SOURCE_URL = "https://github.com/lanl/ExactPack/tree/master/exactpack/solvers/noh"


def parse_profile(path: Path) -> dict[str, np.ndarray]:
    rows = read_numeric_rows(path)
    return {
        "radius": np.array([float(row[1]) for row in rows], dtype=float),
        "rho": np.array([float(row[2]) for row in rows], dtype=float),
        "radial_velocity": np.array([float(row[3]) for row in rows], dtype=float),
        "velocity_magnitude": np.array([float(row[4]) for row in rows], dtype=float),
        "mom_r": np.array([float(row[5]) for row in rows], dtype=float),
        "e_fluid_total": np.array([float(row[6]) for row in rows], dtype=float),
    }


def numerical_pressure(gamma: float, profile: dict[str, np.ndarray]) -> np.ndarray:
    kinetic = 0.5 * np.divide(
        profile["mom_r"] * profile["mom_r"],
        profile["rho"],
        out=np.zeros_like(profile["mom_r"]),
        where=profile["rho"] > 0.0,
    )
    return (gamma - 1.0) * (profile["e_fluid_total"] - kinetic)


def relative_l1(numerical: np.ndarray, exact: np.ndarray) -> float:
    denominator = float(np.sum(np.abs(exact)))
    if denominator <= 0.0:
        return float("nan")
    return float(np.sum(np.abs(numerical - exact)) / denominator)


def absolute_l1_mean(numerical: np.ndarray, exact: np.ndarray) -> float:
    if numerical.size == 0:
        return float("nan")
    return float(np.mean(np.abs(numerical - exact)))


def detect_outer_shock_near_exact(radius: np.ndarray, rho: np.ndarray, exact_shock_radius: float) -> float:
    if radius.size < 2 or exact_shock_radius <= 0.0:
        return 0.0
    best_radius = 0.0
    strongest_fall = float("inf")
    for index in range(radius.size - 1):
        interface_radius = 0.5 * (radius[index] + radius[index + 1])
        if interface_radius < 0.4 * exact_shock_radius or interface_radius > 1.6 * exact_shock_radius:
            continue
        gradient = (rho[index + 1] - rho[index]) / (radius[index + 1] - radius[index])
        if gradient < strongest_fall:
            strongest_fall = gradient
            best_radius = interface_radius
    return best_radius


def write_exact_profile(
    path: Path,
    radius: np.ndarray,
    density: np.ndarray,
    pressure: np.ndarray,
    specific_internal_energy: np.ndarray,
    velocity: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=exactpack_noh_radial_profile\n")
        output.write(f"# source={EXACTPACK_SOURCE_URL}\n")
        output.write("# columns=radial_index radius density pressure specific_internal_energy velocity\n")
        for index in range(radius.size):
            output.write(
                f"{index} {radius[index]:.17g} {density[index]:.17g} {pressure[index]:.17g} "
                f"{specific_internal_energy[index]:.17g} {velocity[index]:.17g}\n"
            )


def write_profile_comparison(
    path: Path,
    radius: np.ndarray,
    numerical_rho: np.ndarray,
    exact_rho: np.ndarray,
    numerical_velocity: np.ndarray,
    exact_velocity: np.ndarray,
    numerical_pressure_values: np.ndarray,
    exact_pressure: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=exactpack_noh_profile_comparison\n")
        output.write("# columns=radial_index radius numerical_rho exact_rho numerical_velocity exact_velocity numerical_pressure exact_pressure\n")
        for index in range(radius.size):
            output.write(
                f"{index} {radius[index]:.17g} {numerical_rho[index]:.17g} {exact_rho[index]:.17g} "
                f"{numerical_velocity[index]:.17g} {exact_velocity[index]:.17g} "
                f"{numerical_pressure_values[index]:.17g} {exact_pressure[index]:.17g}\n"
            )


def write_shock_history(path: Path, rows: list[tuple[float, float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=noh_shock_radius_history\n")
        output.write("# source=exactpack_profile_comparison\n")
        output.write("# definition=numerical radius is strongest negative density gradient near ExactPack shock radius\n")
        output.write("# columns=time_s numerical_radius exactpack_radius relative_error\n")
        for time_s, numerical_radius, exact_radius, relative_error in rows:
            output.write(
                f"{time_s:.17g} {numerical_radius:.17g} {exact_radius:.17g} {relative_error:.17g}\n"
            )


def plot_overlay(
    output_path: Path,
    radius: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    y_label: str,
    title: str,
    exact_shock_radius: float,
    x_max: float,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(radius, exact, label="ExactPack", linewidth=1.7, zorder=2)
    ax.plot(
        radius,
        numerical,
        label="DEC3D numerical",
        linewidth=1.7,
        linestyle="--",
        zorder=3,
    )
    ax.axvline(exact_shock_radius, color="black", linestyle="--", linewidth=1.2, label="ExactPack shock")
    ax.set_xlim(0.0, x_max)
    ax.set_xlabel("radius")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_shock_history(path: Path, output_path: Path) -> None:
    rows = read_numeric_rows(path)
    time_s = [float(row[0]) for row in rows]
    numerical = [float(row[1]) for row in rows]
    exact = [float(row[2]) for row in rows]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(time_s, numerical, marker="o", label="DEC3D numerical")
    ax.plot(time_s, exact, marker="s", label="ExactPack")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("shock radius")
    ax.set_title("Spherical Noh shock radius: numerical vs ExactPack")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", required=True)
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--x-max", type=float, default=0.4)
    args = parser.parse_args()

    case_json = Path(args.case_json)
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else (case_dir / "exactpack_compare")
    output_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    vr0 = float(config["init"]["vr0"])
    p0 = float(config["init"]["p0"])
    final_time_s = float(config["runtime"]["t_end_s"])

    profile = parse_profile(case_dir / "radial_profile_t1.txt")
    radius = profile["radius"]
    pressure = numerical_pressure(gamma, profile)

    solver = Noh(geometry=3, gamma=gamma, rho0=rho0, u0=vr0)
    exact = solver(radius, final_time_s)
    exact_rho = np.asarray(exact["density"], dtype=float)
    exact_pressure = np.asarray(exact["pressure"], dtype=float)
    exact_specific_internal_energy = np.asarray(exact["specific_internal_energy"], dtype=float)
    exact_velocity = np.asarray(exact["velocity"], dtype=float)
    exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else final_time_s / 3.0

    numerical_shock_radius = detect_outer_shock_near_exact(radius, profile["rho"], exact_shock_radius)
    shock_relative_error = (
        abs(numerical_shock_radius - exact_shock_radius) / exact_shock_radius
        if exact_shock_radius > 0.0
        else float("nan")
    )

    inside = radius <= exact_shock_radius + 1.0e-14
    density_l1_global = relative_l1(profile["rho"], exact_rho)
    density_l1_inside = relative_l1(profile["rho"][inside], exact_rho[inside])
    velocity_l1_inside = relative_l1(profile["radial_velocity"][inside], exact_velocity[inside])
    velocity_l1_absolute_inside = absolute_l1_mean(profile["radial_velocity"][inside], exact_velocity[inside])
    velocity_max_abs_inside = (
        float(np.max(np.abs(profile["radial_velocity"][inside] - exact_velocity[inside])))
        if np.any(inside)
        else float("nan")
    )
    pressure_l1_inside = relative_l1(pressure[inside], exact_pressure[inside])
    density_plateau_mean = float(np.mean(profile["rho"][inside])) if np.any(inside) else float("nan")
    density_plateau_exact = float(np.mean(exact_rho[inside])) if np.any(inside) else float("nan")

    shock_rows: list[tuple[float, float, float, float]] = []
    for path in sorted(case_dir.glob("radial_profile_step*.txt")):
        headers = read_headers(path)
        rows = read_numeric_rows(path)
        if not rows:
            continue
        time_s = float(headers.get("time_s", "0"))
        profile_radius = np.array([float(row[1]) for row in rows], dtype=float)
        profile_rho = np.array([float(row[2]) for row in rows], dtype=float)
        if time_s <= 0.0:
            exact_radius = 0.0
            numerical_radius = 0.0
            relative_error = 0.0
        else:
            exact_solution = solver(profile_radius, time_s)
            exact_radius = float(exact_solution.jumps[0]) if getattr(exact_solution, "jumps", None) else time_s / 3.0
            numerical_radius = detect_outer_shock_near_exact(profile_radius, profile_rho, exact_radius)
            relative_error = abs(numerical_radius - exact_radius) / exact_radius if exact_radius > 0.0 else float("nan")
        shock_rows.append((time_s, numerical_radius, exact_radius, relative_error))

    write_exact_profile(
        output_dir / "exactpack_radial_profile_t1.txt",
        radius,
        exact_rho,
        exact_pressure,
        exact_specific_internal_energy,
        exact_velocity,
    )
    write_profile_comparison(
        output_dir / "radial_profile_comparison_t1.txt",
        radius,
        profile["rho"],
        exact_rho,
        profile["radial_velocity"],
        exact_velocity,
        pressure,
        exact_pressure,
    )
    write_shock_history(output_dir / "shock_radius_exactpack_compare.txt", shock_rows)

    summary_lines = [
        "case=case_noh_spherical_exactpack_compare",
        f"exactpack_source={EXACTPACK_SOURCE_URL}",
        f"gamma={gamma:.17g}",
        f"rho0={rho0:.17g}",
        f"vr0={vr0:.17g}",
        f"ambient_pressure_floor={p0:.17g}",
        f"t_end_s={final_time_s:.17g}",
        f"exact_shock_radius_t1={exact_shock_radius:.17g}",
        f"numerical_shock_radius_t1={numerical_shock_radius:.17g}",
        f"shock_radius_relative_error_t1={shock_relative_error:.17g}",
        f"density_plateau_mean_inside_exact_shock={density_plateau_mean:.17g}",
        f"density_plateau_exact_inside_shock={density_plateau_exact:.17g}",
        f"density_l1_relative_global={density_l1_global:.17g}",
        f"density_l1_relative_inside_shock={density_l1_inside:.17g}",
        f"velocity_l1_relative_inside_shock={velocity_l1_inside:.17g}",
        f"velocity_l1_absolute_inside_shock={velocity_l1_absolute_inside:.17g}",
        f"velocity_max_abs_inside_shock={velocity_max_abs_inside:.17g}",
        f"pressure_l1_relative_inside_shock={pressure_l1_inside:.17g}",
        "notes=ExactPack Noh uses the ideal zero-pressure self-similar solution; DEC3D run uses a small pressure floor.",
    ]
    (output_dir / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest_lines = [
        "# kind=case_manifest",
        "case=case_noh_spherical_exactpack_compare",
        "summary=summary.txt",
        "exact_profile=exactpack_radial_profile_t1.txt",
        "profile_comparison=radial_profile_comparison_t1.txt",
        "shock_radius_comparison=shock_radius_exactpack_compare.txt",
        "rho_plot=rho_exactpack_compare_t1.png",
        "velocity_plot=velocity_exactpack_compare_t1.png",
        "pressure_plot=pressure_exactpack_compare_t1.png",
        "shock_history_plot=shock_radius_exactpack_compare.png",
    ]
    (output_dir / "case_manifest.txt").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    plot_overlay(
        output_dir / "rho_exactpack_compare_t1.png",
        radius,
        profile["rho"],
        exact_rho,
        "density",
        "Spherical Noh density: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_overlay(
        output_dir / "velocity_exactpack_compare_t1.png",
        radius,
        profile["radial_velocity"],
        exact_velocity,
        "radial velocity",
        "Spherical Noh radial velocity: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_overlay(
        output_dir / "pressure_exactpack_compare_t1.png",
        radius,
        pressure,
        exact_pressure,
        "pressure",
        "Spherical Noh pressure: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_shock_history(
        output_dir / "shock_radius_exactpack_compare.txt",
        output_dir / "shock_radius_exactpack_compare.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
