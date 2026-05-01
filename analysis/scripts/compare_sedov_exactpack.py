from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from exactpack.solvers.sedov.sedov import Sedov

from plot_hydro_fields import detect_shock_features_from_profile, read_headers, read_numeric_rows


EXACTPACK_SOURCE_URL = "https://github.com/lanl/ExactPack/tree/master/exactpack/solvers/sedov"


def read_key_value_file(path: Path) -> dict[str, str]:
    entries: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or "=" not in stripped:
                continue
            key, value = stripped.split("=", 1)
            entries[key] = value
    return entries


def parse_radial_profile(path: Path) -> dict[str, np.ndarray]:
    rows = read_numeric_rows(path)
    return {
        "radius": np.array([float(row[1]) for row in rows], dtype=float),
        "rho": np.array([float(row[2]) for row in rows], dtype=float),
        "hydro_only_te": np.array([float(row[3]) for row in rows], dtype=float),
        "velocity_magnitude": np.array([float(row[4]) for row in rows], dtype=float),
        "mom_r": np.array([float(row[5]) for row in rows], dtype=float),
        "e_fluid_total": np.array([float(row[6]) for row in rows], dtype=float),
    }


def compute_numerical_pressure(gamma: float, rho: np.ndarray, mom_r: np.ndarray, e_total: np.ndarray) -> np.ndarray:
    kinetic = 0.5 * np.divide(mom_r * mom_r, rho, out=np.zeros_like(mom_r), where=rho > 0.0)
    return (gamma - 1.0) * (e_total - kinetic)


def relative_l1(numerical: np.ndarray, exact: np.ndarray) -> float:
    denominator = float(np.sum(np.abs(exact)))
    if denominator <= 0.0:
        return float("nan")
    return float(np.sum(np.abs(numerical - exact)) / denominator)


def mask_inside_shock(radius: np.ndarray, shock_radius: float) -> np.ndarray:
    return radius <= shock_radius + 1.0e-14


def write_exact_profile(
    path: Path,
    radius: np.ndarray,
    density: np.ndarray,
    pressure: np.ndarray,
    specific_internal_energy: np.ndarray,
    velocity: np.ndarray,
    sound_speed: np.ndarray,
    e_fluid_total: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=exactpack_radial_profile\n")
        output.write(f"# source={EXACTPACK_SOURCE_URL}\n")
        output.write("# columns=radial_index radius density pressure specific_internal_energy velocity sound_speed e_fluid_total\n")
        for index in range(radius.size):
            output.write(
                f"{index} {radius[index]:.17g} {density[index]:.17g} {pressure[index]:.17g} "
                f"{specific_internal_energy[index]:.17g} {velocity[index]:.17g} {sound_speed[index]:.17g} {e_fluid_total[index]:.17g}\n"
            )


def write_profile_comparison(
    path: Path,
    radius: np.ndarray,
    numerical_rho: np.ndarray,
    exact_rho: np.ndarray,
    numerical_velocity: np.ndarray,
    exact_velocity: np.ndarray,
    numerical_pressure: np.ndarray,
    exact_pressure: np.ndarray,
    numerical_e_total: np.ndarray,
    exact_e_total: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=exactpack_profile_comparison\n")
        output.write("# columns=radial_index radius numerical_rho exact_rho numerical_velocity exact_velocity numerical_pressure exact_pressure numerical_e_fluid_total exact_e_fluid_total\n")
        for index in range(radius.size):
            output.write(
                f"{index} {radius[index]:.17g} {numerical_rho[index]:.17g} {exact_rho[index]:.17g} "
                f"{numerical_velocity[index]:.17g} {exact_velocity[index]:.17g} "
                f"{numerical_pressure[index]:.17g} {exact_pressure[index]:.17g} "
                f"{numerical_e_total[index]:.17g} {exact_e_total[index]:.17g}\n"
            )


def write_shock_history_comparison(path: Path, rows: list[tuple[float, float, float, float, float, float]]) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=shock_radius_history\n")
        output.write("# source=exactpack_profile_comparison\n")
        output.write("# definition=numerical_radius is outer shock front from strongest negative density gradient\n")
        output.write("# columns=time_s numerical_radius exactpack_radius relative_error inner_rise_radius density_peak_radius\n")
        for time_s, numerical_radius, exact_radius, relative_error, inner_rise_radius, density_peak_radius in rows:
            output.write(
                f"{time_s:.17g} {numerical_radius:.17g} {exact_radius:.17g} {relative_error:.17g} "
                f"{inner_rise_radius:.17g} {density_peak_radius:.17g}\n"
            )


def plot_overlay(
    output_path: Path,
    radius: np.ndarray,
    numerical: np.ndarray,
    exact: np.ndarray,
    y_label: str,
    title: str,
    exact_shock_radius: float,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(radius, numerical, label="numerical", linewidth=1.7)
    ax.plot(radius, exact, label="ExactPack", linewidth=1.7)
    ax.axvline(exact_shock_radius, color="black", linestyle="--", linewidth=1.2, label="ExactPack shock")
    ax.set_xlabel("radius")
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_density_time_slice(
    output_path: Path,
    radius: np.ndarray,
    numerical_density: np.ndarray,
    exact_density: np.ndarray,
    exact_shock_radius: float,
    time_s: float,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(radius, numerical_density, label="numerical", linewidth=1.7)
    ax.plot(radius, exact_density, label="ExactPack", linewidth=1.7)
    ax.axvline(exact_shock_radius, color="black", linestyle="--", linewidth=1.2, label="ExactPack shock")
    ax.set_xlabel("radius")
    ax.set_ylabel("density")
    ax.set_title(f"Sedov density comparison @ t={time_s:.6g}")
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
    ax.plot(time_s, numerical, marker="o", label="numerical")
    ax.plot(time_s, exact, marker="s", label="ExactPack")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("shock radius")
    ax.set_title("Sedov shock radius: numerical vs ExactPack")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", required=True)
    parser.add_argument("--case-dir", required=True)
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    case_json = Path(args.case_json)
    case_dir = Path(args.case_dir)
    output_dir = Path(args.output_dir) if args.output_dir else (case_dir / "exactpack_compare")
    output_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    eblast = float(config["init"]["blast_energy"])
    final_time_s = float(config["runtime"]["t_end_s"])
    blast_radius = float(config["init"]["blast_radius"])
    ambient_pressure = float(config["init"]["p0"])

    numerical_profile = parse_radial_profile(case_dir / "radial_profile_t1.txt")
    radius = numerical_profile["radius"]
    numerical_rho = numerical_profile["rho"]
    numerical_velocity = np.divide(
        numerical_profile["mom_r"],
        numerical_rho,
        out=np.zeros_like(numerical_rho),
        where=numerical_rho > 0.0,
    )
    numerical_pressure = compute_numerical_pressure(
        gamma,
        numerical_rho,
        numerical_profile["mom_r"],
        numerical_profile["e_fluid_total"],
    )

    solver = Sedov(geometry=3, gamma=gamma, rho0=rho0, omega=0.0, eblast=eblast)
    exact = solver(radius, final_time_s)
    exact_rho = np.asarray(exact["density"], dtype=float)
    exact_pressure = np.asarray(exact["pressure"], dtype=float)
    exact_specific_internal_energy = np.asarray(exact["specific_internal_energy"], dtype=float)
    exact_velocity = np.asarray(exact["velocity"], dtype=float)
    exact_sound_speed = np.asarray(exact["sound_speed"], dtype=float)
    exact_e_total = exact_rho * exact_specific_internal_energy + 0.5 * exact_rho * exact_velocity * exact_velocity
    exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")

    inside_mask = mask_inside_shock(radius, exact_shock_radius)
    global_density_l1 = relative_l1(numerical_rho, exact_rho)
    shocked_density_l1 = relative_l1(numerical_rho[inside_mask], exact_rho[inside_mask])
    shocked_velocity_l1 = relative_l1(numerical_velocity[inside_mask], exact_velocity[inside_mask])
    shocked_pressure_l1 = relative_l1(numerical_pressure[inside_mask], exact_pressure[inside_mask])
    shocked_e_total_l1 = relative_l1(numerical_profile["e_fluid_total"][inside_mask], exact_e_total[inside_mask])

    final_features = detect_shock_features_from_profile(radius.tolist(), numerical_rho.tolist())

    shock_rows: list[tuple[float, float, float, float, float, float]] = []
    for path in sorted(case_dir.glob("radial_profile_step*.txt")):
        headers = read_headers(path)
        rows = read_numeric_rows(path)
        if not rows:
            continue
        profile_radius = [float(row[1]) for row in rows]
        profile_rho = [float(row[2]) for row in rows]
        time_s = float(headers.get("time_s", "0"))
        features = detect_shock_features_from_profile(profile_radius, profile_rho)
        numerical_shock = features.shock_radius
        exact_solution = solver(np.array([0.0, profile_radius[-1]], dtype=float), time_s)
        exact_shock = float(exact_solution.jumps[0]) if getattr(exact_solution, "jumps", None) else float("nan")
        relative_error = abs(numerical_shock - exact_shock) / exact_shock if exact_shock > 0.0 else float("nan")
        shock_rows.append(
            (
                time_s,
                numerical_shock,
                exact_shock,
                relative_error,
                features.inner_rise_radius,
                features.density_peak_radius,
            )
        )

        profile_radius_array = np.array(profile_radius, dtype=float)
        exact_profile_solution = solver(profile_radius_array, time_s)
        exact_profile_density = np.asarray(exact_profile_solution["density"], dtype=float)
        step_index = int(headers.get("step_index", "0"))
        plot_density_time_slice(
            output_dir / f"density_exactpack_compare_step{step_index:06d}.png",
            profile_radius_array,
            np.array(profile_rho, dtype=float),
            exact_profile_density,
            exact_shock,
            time_s,
        )

    numerical_final_shock = shock_rows[-1][1] if shock_rows else float("nan")
    final_shock_relative_error = shock_rows[-1][3] if shock_rows else float("nan")

    write_exact_profile(
        output_dir / "exactpack_radial_profile_t1.txt",
        radius,
        exact_rho,
        exact_pressure,
        exact_specific_internal_energy,
        exact_velocity,
        exact_sound_speed,
        exact_e_total,
    )
    write_profile_comparison(
        output_dir / "radial_profile_comparison_t1.txt",
        radius,
        numerical_rho,
        exact_rho,
        numerical_velocity,
        exact_velocity,
        numerical_pressure,
        exact_pressure,
        numerical_profile["e_fluid_total"],
        exact_e_total,
    )
    write_shock_history_comparison(output_dir / "shock_radius_exactpack_compare.txt", shock_rows)

    summary_lines = [
        "case=case_sedov_spherical_exactpack_compare",
        f"exactpack_source={EXACTPACK_SOURCE_URL}",
        f"gamma={gamma:.17g}",
        f"rho0={rho0:.17g}",
        f"eblast={eblast:.17g}",
        f"t_end_s={final_time_s:.17g}",
        f"blast_radius={blast_radius:.17g}",
        f"ambient_pressure={ambient_pressure:.17g}",
        f"exact_shock_radius_t1={exact_shock_radius:.17g}",
        f"numerical_shock_radius_t1={numerical_final_shock:.17g}",
        f"shock_radius_relative_error_t1={final_shock_relative_error:.17g}",
        f"inner_rise_radius_t1={final_features.inner_rise_radius:.17g}",
        f"density_peak_radius_t1={final_features.density_peak_radius:.17g}",
        f"density_l1_relative_global={global_density_l1:.17g}",
        f"density_l1_relative_inside_shock={shocked_density_l1:.17g}",
        f"velocity_l1_relative_inside_shock={shocked_velocity_l1:.17g}",
        f"pressure_l1_relative_inside_shock={shocked_pressure_l1:.17g}",
        f"e_fluid_total_l1_relative_inside_shock={shocked_e_total_l1:.17g}",
        "notes=ExactPack Sedov is a point-blast self-similar solution; DEC3D baseline uses a finite blast_radius and small nonzero ambient pressure.",
    ]
    (output_dir / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest_lines = [
        "# kind=case_manifest",
        "case=case_sedov_spherical_exactpack_compare",
        "summary=summary.txt",
        "exact_profile=exactpack_radial_profile_t1.txt",
        "profile_comparison=radial_profile_comparison_t1.txt",
        "shock_radius_comparison=shock_radius_exactpack_compare.txt",
        "density_time_slices=density_exactpack_compare_step*.png",
        "rho_plot=rho_exactpack_compare_t1.png",
        "velocity_plot=velocity_exactpack_compare_t1.png",
        "pressure_plot=pressure_exactpack_compare_t1.png",
        "e_fluid_total_plot=e_fluid_total_exactpack_compare_t1.png",
        "shock_history_plot=shock_radius_exactpack_compare.png",
    ]
    (output_dir / "case_manifest.txt").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    plot_overlay(
        output_dir / "rho_exactpack_compare_t1.png",
        radius,
        numerical_rho,
        exact_rho,
        "density",
        "Sedov density: numerical vs ExactPack",
        exact_shock_radius,
    )
    plot_overlay(
        output_dir / "velocity_exactpack_compare_t1.png",
        radius,
        numerical_velocity,
        exact_velocity,
        "radial velocity",
        "Sedov velocity: numerical vs ExactPack",
        exact_shock_radius,
    )
    plot_overlay(
        output_dir / "pressure_exactpack_compare_t1.png",
        radius,
        numerical_pressure,
        exact_pressure,
        "pressure",
        "Sedov pressure: numerical vs ExactPack",
        exact_shock_radius,
    )
    plot_overlay(
        output_dir / "e_fluid_total_exactpack_compare_t1.png",
        radius,
        numerical_profile["e_fluid_total"],
        exact_e_total,
        "E_fluid_total",
        "Sedov total energy density: numerical vs ExactPack",
        exact_shock_radius,
    )
    plot_shock_history(
        output_dir / "shock_radius_exactpack_compare.txt",
        output_dir / "shock_radius_exactpack_compare.png",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
