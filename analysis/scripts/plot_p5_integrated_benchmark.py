import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def load_csv(path: Path):
    return np.genfromtxt(path, delimiter=",", names=True, encoding="utf-8")


def save_radial_profiles(root: Path) -> None:
    data = load_csv(root / "p5_radial_profiles.csv")
    steps = np.unique(data["step"])
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    fields = [
        ("Te_eV", "Te (eV)"),
        ("Ti_eV", "Ti (eV)"),
        ("radiation_energy_erg_cm3", "sum Ug (erg/cm^3)"),
        ("epsilon_alpha_erg_cm3", "epsilon alpha (erg/cm^3)"),
    ]
    for ax, (field, label) in zip(axes.ravel(), fields):
        for step in steps:
            mask = data["step"] == step
            ax.plot(data["r_um"][mask], data[field][mask], label=f"step {int(step)}")
        ax.set_xlabel("r (um)")
        ax.set_ylabel(label)
        ax.grid(True, alpha=0.3)
    axes[0, 0].legend(loc="best", fontsize=8)
    fig.suptitle("P5 integrated radial profiles")
    fig.savefig(root / "p5_radial_profiles.png", dpi=160)
    plt.close(fig)


def save_mode_amplitudes(root: Path) -> None:
    data = load_csv(root / "p5_mode_amplitudes.csv")
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.plot(data["time_s"], data["amplitude"], marker="o", label="target mode")
    ax.plot(data["time_s"], data["phi_leakage"], marker="s", label="phi leakage")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("amplitude")
    ax.set_title("P5 mode diagnostics")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(root / "p5_mode_amplitudes.png", dpi=160)
    plt.close(fig)


def save_budget_residuals(root: Path) -> None:
    data = load_csv(root / "p5_budget_history.csv")
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for field in [
        "hydro_budget_residual",
        "thermal_budget_residual",
        "radiation_budget_residual",
        "alpha_budget_residual",
        "max_stage_budget_residual",
    ]:
        ax.plot(data["time_s"], data[field], marker="o", label=field)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("residual")
    ax.set_title("P5 stage budget residuals")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.savefig(root / "p5_budget_residuals.png", dpi=160)
    plt.close(fig)


def save_rz_gif(root: Path) -> None:
    frame_path = root / "p5_rz_frame_data.csv"
    if not frame_path.exists():
        return
    data = load_csv(frame_path)
    steps = np.unique(data["step"])
    frame_dir = root / "_p5_rz_frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    images = []
    fields = [
        ("rho_g_cm3", "rho (g/cm^3)"),
        ("Te_eV", "Te (eV)"),
        ("radiation_energy_erg_cm3", "sum Ug"),
        ("epsilon_alpha_erg_cm3", "epsilon alpha"),
    ]
    for step in steps:
        mask = data["step"] == step
        fig, axes = plt.subplots(2, 2, figsize=(9, 8), constrained_layout=True)
        for ax, (field, label) in zip(axes.ravel(), fields):
            sc = ax.scatter(
                data["x_um"][mask],
                data["z_um"][mask],
                c=data[field][mask],
                s=8,
                cmap="viridis",
            )
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x (um)")
            ax.set_ylabel("z (um)")
            ax.set_title(label)
            fig.colorbar(sc, ax=ax, shrink=0.85)
        fig.suptitle(f"P5 r-z slice, step {int(step)}")
        png = frame_dir / f"frame_{int(step):04d}.png"
        fig.savefig(png, dpi=130)
        plt.close(fig)
        with Image.open(png) as im:
            images.append(im.convert("P").copy())
    if images:
        images[0].save(
            root / "p5_p5_legendre_l2_m0_primary_rz_evolution.gif",
            save_all=True,
            append_images=images[1:],
            duration=350,
            loop=0,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    args = parser.parse_args()
    root = Path(args.input_dir)
    save_radial_profiles(root)
    save_mode_amplitudes(root)
    save_budget_residuals(root)
    save_rz_gif(root)


if __name__ == "__main__":
    main()
