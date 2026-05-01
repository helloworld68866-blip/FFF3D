import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def save_profiles(out_dir: Path, serial: pd.DataFrame, distributed: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 7.0), sharex=True)
    for time_s, frame in serial.groupby("time_s"):
        axes[0].plot(
            frame["r_cm"] * 1.0e4,
            frame["epsilon_alpha_erg_cm3"],
            lw=2.0,
            label=f"serial t={time_s:.1e}s",
        )
        axes[1].plot(
            frame["r_cm"] * 1.0e4,
            frame["Te_keV"],
            lw=2.0,
            label=f"serial t={time_s:.1e}s",
        )
    for time_s, frame in distributed.groupby("time_s"):
        axes[0].plot(
            frame["r_cm"] * 1.0e4,
            frame["epsilon_alpha_erg_cm3"],
            "--",
            lw=1.8,
            label=f"distributed t={time_s:.1e}s",
        )
        axes[1].plot(
            frame["r_cm"] * 1.0e4,
            frame["Te_keV"],
            "--",
            lw=1.8,
            label=f"distributed t={time_s:.1e}s",
        )
    axes[0].set_ylabel(r"$\epsilon_\alpha$ (erg/cm$^3$)")
    axes[1].set_ylabel(r"$T_e$ (keV)")
    axes[1].set_xlabel("r (um)")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle("P4-B1 alpha hotspot: serial vs distributed")
    fig.tight_layout()
    fig.savefig(out_dir / "p4_b1_profiles.png", dpi=180)
    plt.close(fig)


def save_budget(out_dir: Path, budget: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    residual_columns = [
        "alpha_equation_budget_residual",
        "alpha_electron_exchange_residual",
        "global_alpha_plus_electron_budget_residual",
    ]
    for run_kind, frame in budget.groupby("run_kind"):
        for column in residual_columns:
            residual = np.maximum(np.abs(frame[column].to_numpy(dtype=float)), 1.0e-45)
            ax.plot(
                frame["time_s"],
                residual,
                marker="o",
                label=f"{run_kind} {column}",
            )
    ax.set_xlabel("time (s)")
    ax.set_ylabel("absolute residual (floored at 1e-45)")
    ax.set_yscale("log")
    ax.set_title("P4-B1 alpha budget closure residuals")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "p4_b1_budget_closure.png", dpi=180)
    plt.close(fig)


def save_error(out_dir: Path, serial: pd.DataFrame, distributed: pd.DataFrame) -> dict:
    merged = serial.merge(
        distributed,
        on=["time_s", "r_cm"],
        suffixes=("_serial", "_distributed"),
    )
    eps_err = np.abs(
        merged["epsilon_alpha_erg_cm3_serial"]
        - merged["epsilon_alpha_erg_cm3_distributed"]
    )
    te_err = np.abs(merged["Te_keV_serial"] - merged["Te_keV_distributed"])
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
    axes[0].plot(merged["r_cm"] * 1.0e4, eps_err, ".", ms=4)
    axes[1].plot(merged["r_cm"] * 1.0e4, te_err, ".", ms=4)
    axes[0].set_ylabel(r"$|\Delta\epsilon_\alpha|$")
    axes[1].set_ylabel(r"$|\Delta T_e|$ (keV)")
    axes[1].set_xlabel("r (um)")
    axes[0].set_title("Serial vs distributed pointwise error")
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "p4_b1_serial_distributed_error.png", dpi=180)
    plt.close(fig)
    return {
        "epsilon_alpha_l1": float(eps_err.mean()) if len(eps_err) else 0.0,
        "epsilon_alpha_linf": float(eps_err.max()) if len(eps_err) else 0.0,
        "Te_l1": float(te_err.mean()) if len(te_err) else 0.0,
        "Te_linf": float(te_err.max()) if len(te_err) else 0.0,
    }


def save_timing(out_dir: Path, timing: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    labels = [f"{int(row.rank_count)} ranks" for row in timing.itertuples()]
    components = [
        "provider_build_wall_s",
        "assembly_wall_s",
        "hypre_solve_wall_s",
        "writeback_wall_s",
    ]
    bottom = np.zeros(len(timing))
    for component in components:
        values = timing[component].to_numpy(dtype=float)
        ax.bar(labels, values, bottom=bottom, label=component.replace("_wall_s", ""))
        bottom += values
    ax.set_ylabel("wall time (s)")
    ax.set_title("P4-B3 distributed alpha timing breakdown")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "p4_b3_timing_breakdown.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    args = parser.parse_args()
    out_dir = Path(args.input_dir)
    serial = pd.read_csv(require(out_dir / "p4_b1_serial_profiles.csv"))
    distributed = pd.read_csv(require(out_dir / "p4_b1_distributed_profiles.csv"))
    timing = pd.read_csv(require(out_dir / "p4_b3_distributed_performance.csv"))
    budget = pd.read_csv(require(out_dir / "p4_b1_budget_summary.csv"))

    save_profiles(out_dir, serial, distributed)
    save_budget(out_dir, budget)
    metrics = save_error(out_dir, serial, distributed)
    save_timing(out_dir, timing)

    summary_path = out_dir / "p4_b1_serial_vs_distributed_summary.json"
    existing = {}
    if summary_path.exists():
        existing = json.loads(summary_path.read_text(encoding="utf-8"))
    existing.update(
        {
            "plot_artifacts_written": True,
            "plot_files": [
                "p4_b1_profiles.png",
                "p4_b1_budget_closure.png",
                "p4_b1_serial_distributed_error.png",
                "p4_b3_timing_breakdown.png",
            ],
            "plot_metrics": metrics,
        }
    )
    summary_path.write_text(json.dumps(existing, indent=2), encoding="utf-8")

    manifest_path = out_dir / "p4_benchmark_manifest.json"
    manifest = {}
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.update(
        {
            "diagnostic_id": "p4.alpha.benchmark_manifest",
            "manifest_scope": "available_artifacts",
            "profile_artifacts_written": (
                (out_dir / "p4_b1_serial_profiles.csv").exists()
                and (out_dir / "p4_b1_distributed_profiles.csv").exists()
            ),
            "budget_artifacts_written": (out_dir / "p4_b1_budget_summary.csv").exists(),
            "performance_artifacts_written": (
                out_dir / "p4_b3_distributed_performance.csv"
            ).exists(),
            "plot_artifacts_written": True,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
