from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RuntimeStatus:
    step: int
    time_s: float
    dt_s: float
    wall_s: float
    status: str


def latest_runtime_status(dec3d_out: Path) -> RuntimeStatus | None:
    if not dec3d_out.exists():
        return None
    latest: RuntimeStatus | None = None
    with dec3d_out.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("step "):
                continue
            parts = stripped.split()
            if len(parts) < 13:
                continue
            try:
                latest = RuntimeStatus(
                    step=int(parts[0]),
                    time_s=float(parts[1]),
                    dt_s=float(parts[2]),
                    wall_s=float(parts[3]),
                    status=parts[12],
                )
            except ValueError:
                continue
    return latest


def plot_new_snaps(args: argparse.Namespace) -> list[Path]:
    output_dir = Path(args.output_dir)
    plot_dir = output_dir / "rz_snap_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    plot_script = Path(args.plot_script)
    streamline_script = (
        Path(args.streamline_script)
        if args.streamline_script
        else plot_script.with_name("plot_snap_rz_velocity_streamlines.py")
    )
    generated: list[Path] = []

    for snap in sorted(output_dir.glob("fields_*.snap")):
        plot_jobs = [
            (plot_script, plot_dir / f"{snap.stem}_rz6.png"),
            (streamline_script, plot_dir / f"{snap.stem}_velocity_arrows.png"),
        ]
        for script, png in plot_jobs:
            if png.exists() and png.stat().st_mtime >= snap.stat().st_mtime:
                continue
            command = [
                sys.executable,
                str(script),
                "--snap",
                str(snap),
                "--output",
                str(png),
                "--radial-max-cm",
                str(args.radial_max_cm),
                "--radial-min-cm",
                str(args.radial_min_cm),
                "--dpi",
                str(args.dpi),
            ]
            subprocess.run(command, check=True)
            generated.append(png)
    return generated


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plot-script", required=True)
    parser.add_argument("--streamline-script", default="")
    parser.add_argument("--radial-min-cm", type=float, default=0.0)
    parser.add_argument("--radial-max-cm", type=float, required=True)
    parser.add_argument("--target-time-s", type=float, default=1.0e-10)
    parser.add_argument("--dpi", type=int, default=170)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    generated = plot_new_snaps(args)
    status = latest_runtime_status(output_dir / "dec3d.out")

    lines: list[str] = []
    lines.append(f"output_dir={output_dir}")
    if status is None:
        lines.append("latest_status=unavailable")
    else:
        lines.append(
            "latest_status="
            f"step={status.step}; time_ps={status.time_s * 1.0e12:.9g}; "
            f"dt_ps={status.dt_s * 1.0e12:.9g}; wall_s={status.wall_s:.3f}; "
            f"status={status.status}; target_reached={status.time_s >= args.target_time_s}"
        )
    lines.append(f"new_plots={len(generated)}")
    for png in generated:
        lines.append(f"plot={png}")

    summary = "\n".join(lines) + "\n"
    (output_dir / "monitor_summary.txt").write_text(summary, encoding="utf-8")
    print(summary, end="")


if __name__ == "__main__":
    main()
