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

from plot_hydro_fields import detect_shock_features_from_profile


EXACTPACK_SOURCE_URL = "https://github.com/lanl/ExactPack/tree/master/exactpack/solvers/sedov"
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
    return profiles


def numerical_pressure(profile: HistoryProfile) -> np.ndarray:
    ion_number_density = profile.rho / MEAN_DT_ION_MASS_G
    electron_number_density = DEFAULT_ZBAR * ion_number_density
    electron_pressure = electron_number_density * profile.te_kev * ERG_PER_KEV
    ion_pressure = ion_number_density * profile.ti_kev * ERG_PER_KEV
    return electron_pressure + ion_pressure


def relative_l1(numerical: np.ndarray, exact: np.ndarray) -> float:
    denominator = float(np.sum(np.abs(exact)))
    if denominator <= 0.0:
        return float("nan")
    return float(np.sum(np.abs(numerical - exact)) / denominator)


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
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            rows.append([float(value) for value in stripped.split()])

    time_s = [row[0] for row in rows]
    numerical = [row[1] for row in rows]
    exact = [row[2] for row in rows]

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(time_s, numerical, marker="o", label="DEC3D numerical")
    ax.plot(time_s, exact, marker="s", label="ExactPack")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("shock radius")
    ax.set_title("Spherical Sedov shock radius: numerical vs ExactPack")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def write_exact_profile(
    path: Path,
    radius: np.ndarray,
    density: np.ndarray,
    pressure: np.ndarray,
    specific_internal_energy: np.ndarray,
    velocity: np.ndarray,
) -> None:
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=exactpack_sedov_radial_profile\n")
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
        output.write("# kind=exactpack_sedov_profile_comparison\n")
        output.write("# columns=radial_index radius numerical_rho exact_rho numerical_velocity exact_velocity numerical_pressure exact_pressure\n")
        for index in range(radius.size):
            output.write(
                f"{index} {radius[index]:.17g} {numerical_rho[index]:.17g} {exact_rho[index]:.17g} "
                f"{numerical_velocity[index]:.17g} {exact_velocity[index]:.17g} "
                f"{numerical_pressure_values[index]:.17g} {exact_pressure[index]:.17g}\n"
            )


def write_shock_history(
    path: Path,
    profiles: list[HistoryProfile],
    solver: Sedov,
) -> tuple[float, float, float, float]:
    final_numerical = float("nan")
    final_exact = float("nan")
    final_relative_error = float("nan")
    final_density_peak = float("nan")
    with path.open("w", encoding="utf-8") as output:
        output.write("# kind=sedov_shock_radius_history\n")
        output.write("# source=exactpack_profile_comparison\n")
        output.write("# definition=numerical radius is outer shock front from strongest negative density gradient\n")
        output.write("# columns=time_s numerical_radius exactpack_radius relative_error density_peak_radius\n")
        for profile in profiles:
            if profile.time_s <= 0.0:
                numerical_radius = 0.0
                exact_radius = 0.0
                relative_error = 0.0
                density_peak_radius = 0.0
            else:
                features = detect_shock_features_from_profile(profile.radius.tolist(), profile.rho.tolist())
                exact_solution = solver(profile.radius, profile.time_s)
                exact_radius = float(exact_solution.jumps[0]) if getattr(exact_solution, "jumps", None) else float("nan")
                numerical_radius = features.shock_radius
                relative_error = (
                    abs(numerical_radius - exact_radius) / exact_radius
                    if exact_radius > 0.0
                    else float("nan")
                )
                density_peak_radius = features.density_peak_radius
            output.write(
                f"{profile.time_s:.17g} {numerical_radius:.17g} {exact_radius:.17g} "
                f"{relative_error:.17g} {density_peak_radius:.17g}\n"
            )
            final_numerical = numerical_radius
            final_exact = exact_radius
            final_relative_error = relative_error
            final_density_peak = density_peak_radius
    return final_numerical, final_exact, final_relative_error, final_density_peak


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-json", required=True)
    parser.add_argument("--history-profile", required=True)
    parser.add_argument("--output-dir")
    parser.add_argument("--x-max", type=float, default=0.12)
    args = parser.parse_args()

    case_json = Path(args.case_json)
    history_profile = Path(args.history_profile)
    output_dir = Path(args.output_dir) if args.output_dir else (history_profile.parent / "exactpack_compare_no_interp_latest")
    output_dir.mkdir(parents=True, exist_ok=True)

    config = json.loads(case_json.read_text(encoding="utf-8"))
    gamma = float(config["hydro"]["gamma"])
    rho0 = float(config["init"]["rho0"])
    eblast = float(config["init"]["blast_energy"])
    p0 = float(config["init"]["p0"])
    blast_radius = float(config["init"]["blast_radius"])

    profiles = read_history_profiles(history_profile)
    final_profile = profiles[-1]
    pressure = numerical_pressure(final_profile)

    solver = Sedov(geometry=3, gamma=gamma, rho0=rho0, omega=0.0, eblast=eblast)
    exact = solver(final_profile.radius, final_profile.time_s)
    exact_rho = np.asarray(exact["density"], dtype=float)
    exact_pressure = np.asarray(exact["pressure"], dtype=float)
    exact_specific_internal_energy = np.asarray(exact["specific_internal_energy"], dtype=float)
    exact_velocity = np.asarray(exact["velocity"], dtype=float)
    exact_shock_radius = float(exact.jumps[0]) if getattr(exact, "jumps", None) else float("nan")

    inside = final_profile.radius <= exact_shock_radius + 1.0e-14
    density_l1_global = relative_l1(final_profile.rho, exact_rho)
    density_l1_inside = relative_l1(final_profile.rho[inside], exact_rho[inside])
    velocity_l1_inside = relative_l1(final_profile.radial_velocity[inside], exact_velocity[inside])
    pressure_l1_inside = relative_l1(pressure[inside], exact_pressure[inside])

    shock_history_path = output_dir / "shock_radius_exactpack_compare.txt"
    numerical_shock, exact_shock, shock_relative_error, density_peak_radius = write_shock_history(
        shock_history_path,
        profiles,
        solver,
    )

    write_exact_profile(
        output_dir / "exactpack_radial_profile_t1.txt",
        final_profile.radius,
        exact_rho,
        exact_pressure,
        exact_specific_internal_energy,
        exact_velocity,
    )
    write_profile_comparison(
        output_dir / "radial_profile_comparison_t1.txt",
        final_profile.radius,
        final_profile.rho,
        exact_rho,
        final_profile.radial_velocity,
        exact_velocity,
        pressure,
        exact_pressure,
    )

    summary_lines = [
        "case=case_sedov_spherical_exactpack_compare",
        f"exactpack_source={EXACTPACK_SOURCE_URL}",
        f"history_profile={history_profile}",
        f"gamma={gamma:.17g}",
        f"rho0={rho0:.17g}",
        f"eblast={eblast:.17g}",
        f"last_history_step={final_profile.step_index}",
        f"last_history_time_s={final_profile.time_s:.17g}",
        f"blast_radius={blast_radius:.17g}",
        f"ambient_pressure={p0:.17g}",
        f"exact_shock_radius_t1={exact_shock:.17g}",
        f"numerical_shock_radius_t1={numerical_shock:.17g}",
        f"shock_radius_relative_error_t1={shock_relative_error:.17g}",
        f"density_peak_radius_t1={density_peak_radius:.17g}",
        f"density_l1_relative_global={density_l1_global:.17g}",
        f"density_l1_relative_inside_shock={density_l1_inside:.17g}",
        f"velocity_l1_relative_inside_shock={velocity_l1_inside:.17g}",
        f"pressure_l1_relative_inside_shock={pressure_l1_inside:.17g}",
        "notes=ExactPack Sedov is a point-blast self-similar solution; DEC3D history profile uses finite blast_radius and angular volume averages.",
    ]
    (output_dir / "summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")

    manifest_lines = [
        "# kind=case_manifest",
        "case=case_sedov_spherical_exactpack_compare",
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
        final_profile.radius,
        final_profile.rho,
        exact_rho,
        "density",
        "Spherical Sedov density: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_overlay(
        output_dir / "velocity_exactpack_compare_t1.png",
        final_profile.radius,
        final_profile.radial_velocity,
        exact_velocity,
        "radial velocity",
        "Spherical Sedov radial velocity: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_overlay(
        output_dir / "pressure_exactpack_compare_t1.png",
        final_profile.radius,
        pressure,
        exact_pressure,
        "pressure",
        "Spherical Sedov pressure: DEC3D vs ExactPack",
        exact_shock_radius,
        args.x_max,
    )
    plot_shock_history(
        shock_history_path,
        output_dir / "shock_radius_exactpack_compare.png",
    )
    print(output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
