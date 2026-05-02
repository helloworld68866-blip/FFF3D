from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis" / "scripts"))

from animate_hydro_history_rho_pressure_velocity import (  # noqa: E402
    ERG_PER_KEV,
    HistoryProfile,
    ReferenceProfile,
    compute_plot_limits,
    compute_total_pressure,
    render_frame,
)


def test_total_pressure_uses_both_temperatures() -> None:
    rho = [2.0, 4.0]
    te_kev = [3.0, 5.0]
    ti_kev = [7.0, 11.0]
    ion_mass_g = 2.0
    zbar = 1.5

    pressure = compute_total_pressure(rho, te_kev, ti_kev, ion_mass_g=ion_mass_g, zbar=zbar)

    expected0 = (rho[0] / ion_mass_g) * (zbar * te_kev[0] + ti_kev[0]) * ERG_PER_KEV
    expected1 = (rho[1] / ion_mass_g) * (zbar * te_kev[1] + ti_kev[1]) * ERG_PER_KEV
    assert abs(float(pressure[0]) - expected0) < 1.0e-30
    assert abs(float(pressure[1]) - expected1) < 1.0e-30


def test_render_frame_writes_wide_three_panel_png(tmp_path: Path) -> None:
    radius = [0.0, 0.5, 1.0]
    current = HistoryProfile(
        step_index=4,
        time_s=0.25,
        radius=radius,
        rho=[1.0, 3.0, 1.0],
        te_kev=[1.0e-6, 2.0e-6, 1.0e-6],
        ti_kev=[1.0e-6, 2.0e-6, 1.0e-6],
        radial_velocity=[-1.0, -0.5, 0.0],
    )
    reference = ReferenceProfile(
        radius=radius,
        rho=[1.0, 4.0, 1.0],
        pressure=[1.0e-6, 4.0e-6, 1.0e-6],
        radial_velocity=[-1.0, 0.0, 0.0],
        shock_radius=0.5,
    )
    output_path = tmp_path / "three_panel.png"

    render_frame(
        current=current,
        reference=reference,
        output_path=output_path,
        title_prefix="Synthetic",
        reference_label="Exact reference",
        x_limit=(0.0, 1.0),
        rho_limit=(0.0, 5.0),
        pressure_limit=(0.0, 5.0e-6),
        velocity_limit=(-1.2, 0.2),
    )

    with Image.open(output_path) as image:
        width, height = image.size
    assert width > 2.5 * height


def test_plot_limits_do_not_let_initial_blast_pressure_hide_later_profiles() -> None:
    initial = HistoryProfile(
        step_index=0,
        time_s=0.0,
        radius=[0.0, 0.5, 1.0],
        rho=[1.0e-24, 1.0e-24, 1.0e-24],
        te_kev=[1000.0, 1.0, 1.0],
        ti_kev=[1000.0, 1.0, 1.0],
        radial_velocity=[0.0, 0.0, 0.0],
    )
    evolved = HistoryProfile(
        step_index=10,
        time_s=0.1,
        radius=[0.0, 0.5, 1.0],
        rho=[1.0e-24, 1.0e-24, 1.0e-24],
        te_kev=[1.0, 2.0, 3.0],
        ti_kev=[1.0, 2.0, 3.0],
        radial_velocity=[0.0, 0.5, 1.0],
    )
    references = [
        ReferenceProfile([0.0, 0.5, 1.0], [1.0, 1.0, 1.0], [1000.0, 1.0, 1.0], [0.0, 0.0, 0.0], 0.0),
        ReferenceProfile([0.0, 0.5, 1.0], [1.0, 1.0, 1.0], [1.0, 2.0, 3.0], [0.0, 0.5, 1.0], 0.5),
    ]

    _, pressure_limit, _ = compute_plot_limits([initial, evolved], references, (0.0, 1.0))

    assert pressure_limit[1] < 10.0


def main() -> int:
    import tempfile

    test_total_pressure_uses_both_temperatures()
    with tempfile.TemporaryDirectory() as temp_dir:
        test_render_frame_writes_wide_three_panel_png(Path(temp_dir))
    test_plot_limits_do_not_let_initial_blast_pressure_hide_later_profiles()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
