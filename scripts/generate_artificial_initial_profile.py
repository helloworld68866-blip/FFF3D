#!/usr/bin/env python3
"""Generate artificial DEC3D radial initial profiles in .pro format.

The default profile is an artificial ICF deceleration-stage initial profile:
it follows docs/artificial_initial_setup_2.md for a hot-spot/shell/buffer
state, while keeping the Woo Eq. 2.2 perturbation shape as metadata for
downstream 3D velocity perturbations. It is not an exact LILAC profile.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, NamedTuple


ERG_PER_KEV = 1.602176634e-9
HBAR_ERG_S = 1.054571817e-27
ELECTRON_MASS_G = 9.1093837015e-28
DEUTERON_MASS_G = 3.3435837724e-24
TRITON_MASS_G = 5.0073567446e-24
MEAN_DT_ION_MASS_G = 0.5 * (DEUTERON_MASS_G + TRITON_MASS_G)
GBAR_TO_DYN_PER_CM2 = 1.0e15

HEADER = (
    "r_um",
    "rho_g_cm3",
    "Te_keV",
    "Ti_keV",
    "vr_cm_s",
    "vt_cm_s",
    "vp_cm_s",
    "epsilon_alpha_erg_cm3",
    "radiation_scale",
)


class ProfileConfig:
    def __init__(
        self,
        *,
        points: int = 1025,
        r_max_um: float = 115.0,
        r0_um: float = 68.0,
        shell_outer_um: float = 86.0,
        ell: int = 2,
        smooth_width_cells: float = 3.0,
        smooth_width_um: float = 4.0,
        hotspot_center_density_g_cm3: float = 1.5,
        hotspot_density_slope: float = 1.0,
        shell_inner_density_g_cm3: float = 35.0,
        shell_density_gradient: float = 0.25,
        outer_density_g_cm3: float = 0.03,
        hotspot_pressure_gbar: float = 3.0,
        shell_adiabat: float = 3.2,
        outer_pressure_gbar: float = 0.003,
        implosion_speed_cm_s: float = 3.0e7,
        hotspot_velocity_fraction: float = 0.45,
        shell_velocity_slope: float = 0.1,
        velocity_model: str = "decelerating-shell",
        shell_inner_velocity_fraction: float = 0.65,
        shell_outer_velocity_fraction: float = 1.0,
    ) -> None:
        self.points = points
        self.r_max_um = r_max_um
        self.r0_um = r0_um
        self.shell_outer_um = shell_outer_um
        self.ell = ell
        self.smooth_width_cells = smooth_width_cells
        self.smooth_width_um = smooth_width_um
        self.hotspot_center_density_g_cm3 = hotspot_center_density_g_cm3
        self.hotspot_density_slope = hotspot_density_slope
        self.shell_inner_density_g_cm3 = shell_inner_density_g_cm3
        self.shell_density_gradient = shell_density_gradient
        self.outer_density_g_cm3 = outer_density_g_cm3
        self.hotspot_pressure_gbar = hotspot_pressure_gbar
        self.shell_adiabat = shell_adiabat
        self.outer_pressure_gbar = outer_pressure_gbar
        self.implosion_speed_cm_s = implosion_speed_cm_s
        self.hotspot_velocity_fraction = hotspot_velocity_fraction
        self.shell_velocity_slope = shell_velocity_slope
        self.velocity_model = velocity_model
        self.shell_inner_velocity_fraction = shell_inner_velocity_fraction
        self.shell_outer_velocity_fraction = shell_outer_velocity_fraction


class ProfileRow(NamedTuple):
    r_um: float
    rho_g_cm3: float
    Te_keV: float
    Ti_keV: float
    vr_cm_s: float
    vt_cm_s: float
    vp_cm_s: float
    epsilon_alpha_erg_cm3: float
    radiation_scale: float


class DecelerationDiagnostics(NamedTuple):
    passed: bool
    checks: dict[str, bool]
    metrics: dict[str, float | str]


def _smooth_step(r_um: float, center_um: float, width_um: float) -> float:
    return 0.5 * (1.0 + math.tanh((r_um - center_um) / width_um))


def _smooth_unit_interval(value: float) -> float:
    x = min(1.0, max(0.0, value))
    return x * x * (3.0 - 2.0 * x)


def _transition_width_um(config: ProfileConfig) -> float:
    dr_um = config.r_max_um / float(config.points - 1)
    return max(config.smooth_width_cells * dr_um, config.smooth_width_um)


def region_weights(r_um: float, config: ProfileConfig) -> tuple[float, float, float]:
    """Return C/S/O weights from artificial_initial_setup.md."""

    width_um = _transition_width_um(config)
    h_if = _smooth_step(r_um, config.r0_um, width_um)
    h_out = _smooth_step(r_um, config.shell_outer_um, width_um)
    center = 1.0 - h_if
    shell = h_if * (1.0 - h_out)
    outer = h_out
    return center, shell, outer


def woo_shape_function(r_um: float, r0_um: float, ell: int) -> float:
    """Woo thesis Eq. 2.2 radial shape function f(r)."""

    if r0_um <= 0.0:
        raise ValueError("r0_um must be positive")
    if ell <= 0:
        raise ValueError("ell must be positive")
    if r_um <= 0.0:
        return 0.0

    width_um = 0.015 * r0_um
    transition = math.tanh((r_um - r0_um) / width_um)
    inner = 0.5 * (r_um / r0_um) ** ell * (1.0 - transition)
    outer = 0.5 * (r0_um / r_um) ** ell * (1.0 + transition)
    return inner + outer


def density_g_cm3(r_um: float, config: ProfileConfig) -> float:
    """Density from docs/artificial_initial_setup.md section 3.1."""

    center_weight, shell_weight, outer_weight = region_weights(r_um, config)
    hotspot = config.hotspot_center_density_g_cm3 * (
        1.0 + config.hotspot_density_slope * (r_um / config.r0_um) ** 2
    )
    shell_x = (r_um - config.r0_um) / (config.shell_outer_um - config.r0_um)
    shell = config.shell_inner_density_g_cm3 * (
        1.0 - config.shell_density_gradient * shell_x
    )
    shell = max(shell, config.outer_density_g_cm3)
    return (
        center_weight * hotspot
        + shell_weight * shell
        + outer_weight * config.outer_density_g_cm3
    )


def fermi_pressure_dyn_cm2(rho_g_cm3: float) -> float:
    """Electron Fermi pressure for equal-molar fully ionized DT."""

    electron_number_density = rho_g_cm3 / MEAN_DT_ION_MASS_G
    return (
        (3.0 * math.pi**2) ** (2.0 / 3.0)
        / 5.0
        * HBAR_ERG_S**2
        / ELECTRON_MASS_G
        * electron_number_density ** (5.0 / 3.0)
    )


def pressure_dyn_cm2(r_um: float, config: ProfileConfig) -> float:
    """Total pressure from docs/artificial_initial_setup.md sections 3.2-3.3."""

    center_weight, shell_weight, outer_weight = region_weights(r_um, config)
    shell_x = (r_um - config.r0_um) / (config.shell_outer_um - config.r0_um)
    shell_density = config.shell_inner_density_g_cm3 * (
        1.0 - config.shell_density_gradient * shell_x
    )
    shell_density = max(shell_density, config.outer_density_g_cm3)
    shell_pressure = config.shell_adiabat * fermi_pressure_dyn_cm2(shell_density)
    hotspot_pressure = config.hotspot_pressure_gbar * GBAR_TO_DYN_PER_CM2
    outer_pressure = config.outer_pressure_gbar * GBAR_TO_DYN_PER_CM2
    return (
        center_weight * hotspot_pressure
        + shell_weight * shell_pressure
        + outer_weight * outer_pressure
    )


def temperature_keV(r_um: float, config: ProfileConfig) -> float:
    """Recover Te=Ti from P_e=P_i=P/2."""

    density = density_g_cm3(r_um, config)
    pressure = pressure_dyn_cm2(r_um, config)
    ion_number_density = density / MEAN_DT_ION_MASS_G
    return 0.5 * pressure / ion_number_density / ERG_PER_KEV


def radial_velocity_cm_s(r_um: float, config: ProfileConfig) -> float:
    """Velocity from docs/artificial_initial_setup_2.md section 3.5.

    The default model is a deceleration-stage interpretation: the inner
    shell is already being braked by hot-spot pressure, while the outer shell
    remains closer to the peak implosion speed. The setup2-original model is
    kept for comparison with the literal formula in the design note.
    """

    center_weight, shell_weight, outer_weight = region_weights(r_um, config)
    shell_x = (r_um - config.r0_um) / (config.shell_outer_um - config.r0_um)
    hotspot_velocity = (
        -config.hotspot_velocity_fraction
        * config.implosion_speed_cm_s
        * r_um
        / config.r0_um
    )
    if config.velocity_model == "setup2-original":
        shell_velocity = -config.implosion_speed_cm_s * (
            1.0 - config.shell_velocity_slope * shell_x
        )
        outer_velocity = -0.5 * config.implosion_speed_cm_s
    elif config.velocity_model == "decelerating-shell":
        shell_fraction = (
            config.shell_inner_velocity_fraction
            + (
                config.shell_outer_velocity_fraction
                - config.shell_inner_velocity_fraction
            )
            * _smooth_unit_interval(shell_x)
        )
        shell_velocity = -config.implosion_speed_cm_s * shell_fraction
        outer_x = (r_um - config.shell_outer_um) / (
            config.r_max_um - config.shell_outer_um
        )
        outer_velocity = (
            -config.implosion_speed_cm_s
            * config.shell_outer_velocity_fraction
            * (1.0 - _smooth_unit_interval(outer_x))
        )
    else:
        raise ValueError(f"unknown velocity_model: {config.velocity_model}")
    return (
        center_weight * hotspot_velocity
        + shell_weight * shell_velocity
        + outer_weight * outer_velocity
    )


def validate_config(config: ProfileConfig) -> None:
    if config.points < 2:
        raise ValueError("points must be at least 2")
    if config.r_max_um <= 0.0:
        raise ValueError("r_max_um must be positive")
    if not (0.0 < config.r0_um < config.shell_outer_um < config.r_max_um):
        raise ValueError("r0_um, shell_outer_um, and r_max_um must be ordered")
    if config.ell <= 0:
        raise ValueError("ell must be positive")
    if config.smooth_width_cells <= 0.0:
        raise ValueError("smooth_width_cells must be positive")
    if config.smooth_width_um < 0.0:
        raise ValueError("smooth_width_um must be non-negative")
    if (
        config.hotspot_center_density_g_cm3 <= 0.0
        or config.shell_inner_density_g_cm3 <= 0.0
        or config.outer_density_g_cm3 <= 0.0
    ):
        raise ValueError("density controls must be positive")
    if (
        config.hotspot_pressure_gbar <= 0.0
        or config.shell_adiabat <= 0.0
        or config.outer_pressure_gbar <= 0.0
    ):
        raise ValueError("pressure controls must be positive")
    if config.implosion_speed_cm_s <= 0.0:
        raise ValueError("implosion_speed_cm_s must be positive")
    if config.hotspot_velocity_fraction < 0.0 or config.shell_velocity_slope < 0.0:
        raise ValueError("velocity shape controls must be non-negative")
    if config.velocity_model not in {"decelerating-shell", "setup2-original"}:
        raise ValueError("velocity_model must be decelerating-shell or setup2-original")
    if (
        config.hotspot_velocity_fraction < 0.0
        or config.shell_inner_velocity_fraction < 0.0
        or config.shell_outer_velocity_fraction < 0.0
    ):
        raise ValueError("velocity fractions must be non-negative")
    if config.shell_inner_velocity_fraction >= config.shell_outer_velocity_fraction:
        raise ValueError("shell_inner_velocity_fraction must be less than shell_outer_velocity_fraction")
    if config.shell_outer_velocity_fraction > 1.2:
        raise ValueError("shell_outer_velocity_fraction is unexpectedly large")


def generate_profile(config: ProfileConfig) -> list[ProfileRow]:
    validate_config(config)
    dr_um = config.r_max_um / float(config.points - 1)
    rows: list[ProfileRow] = []
    for i in range(config.points):
        r_um = dr_um * i
        temp = temperature_keV(r_um, config)
        rows.append(
            ProfileRow(
                r_um=r_um,
                rho_g_cm3=density_g_cm3(r_um, config),
                Te_keV=temp,
                Ti_keV=temp,
                vr_cm_s=radial_velocity_cm_s(r_um, config),
                vt_cm_s=0.0,
                vp_cm_s=0.0,
                epsilon_alpha_erg_cm3=0.0,
                radiation_scale=1.0,
            )
        )
    return rows


def _row_pressure_dyn_cm2(row: ProfileRow) -> float:
    ion_number_density = row.rho_g_cm3 / MEAN_DT_ION_MASS_G
    return ion_number_density * (row.Te_keV + row.Ti_keV) * ERG_PER_KEV


def _gradient(values: list[float], coordinates: list[float]) -> list[float]:
    if len(values) != len(coordinates):
        raise ValueError("values and coordinates must have the same length")
    if len(values) < 2:
        raise ValueError("at least two points are required")

    gradient: list[float] = []
    for index in range(len(values)):
        if index == 0:
            numerator = values[1] - values[0]
            denominator = coordinates[1] - coordinates[0]
        elif index == len(values) - 1:
            numerator = values[-1] - values[-2]
            denominator = coordinates[-1] - coordinates[-2]
        else:
            numerator = values[index + 1] - values[index - 1]
            denominator = coordinates[index + 1] - coordinates[index - 1]
        gradient.append(numerator / denominator)
    return gradient


def _nearest_index(rows: list[ProfileRow], r_um: float) -> int:
    return min(range(len(rows)), key=lambda index: abs(rows[index].r_um - r_um))


def _window_indices(rows: list[ProfileRow], center_um: float, half_width_um: float) -> list[int]:
    selected = [
        index
        for index, row in enumerate(rows)
        if abs(row.r_um - center_um) <= half_width_um
    ]
    if selected:
        return selected
    return [_nearest_index(rows, center_um)]


def _mean(values: list[float]) -> float:
    return sum(values) / float(len(values))


def deceleration_diagnostics(
    rows: list[ProfileRow], config: ProfileConfig
) -> DecelerationDiagnostics:
    """Check whether the generated profile is a deceleration-stage state."""

    if len(rows) < 5:
        raise ValueError("at least five rows are required for diagnostics")

    r_cm = [row.r_um * 1.0e-4 for row in rows]
    rho = [row.rho_g_cm3 for row in rows]
    pressure = [_row_pressure_dyn_cm2(row) for row in rows]
    vr = [row.vr_cm_s for row in rows]
    pressure_gradient = _gradient(pressure, r_cm)
    pressure_acceleration = [
        -pressure_gradient[index] / rho[index] for index in range(len(rows))
    ]
    vr_gradient = _gradient(vr, r_cm)
    r2v = [r_cm[index] * r_cm[index] * vr[index] for index in range(len(rows))]
    r2v_gradient = _gradient(r2v, r_cm)
    divergence = [
        (3.0 * vr_gradient[0])
        if r_cm[index] == 0.0
        else r2v_gradient[index] / (r_cm[index] * r_cm[index])
        for index in range(len(rows))
    ]

    width_um = _transition_width_um(config)
    interface_indices = _window_indices(rows, config.r0_um, width_um)
    hotspot_probe = max(0.0, config.r0_um - 2.0 * width_um)
    inner_shell_probe = min(
        config.shell_outer_um - width_um, config.r0_um + 2.0 * width_um
    )
    outer_shell_probe = max(
        config.r0_um + width_um, config.shell_outer_um - 2.0 * width_um
    )

    hotspot_index = _nearest_index(rows, hotspot_probe)
    inner_shell_index = _nearest_index(rows, inner_shell_probe)
    outer_shell_index = _nearest_index(rows, outer_shell_probe)

    interface_velocity = _mean([vr[index] for index in interface_indices])
    interface_acceleration = _mean(
        [pressure_acceleration[index] for index in interface_indices]
    )
    interface_convergence = max(
        0.0, -min(divergence[index] for index in interface_indices)
    )
    pressure_ratio = pressure[hotspot_index] / pressure[inner_shell_index]
    shell_density_ratio = rho[inner_shell_index] / rho[hotspot_index]
    shell_inner_speed_fraction = abs(vr[inner_shell_index]) / config.implosion_speed_cm_s
    shell_outer_speed_fraction = abs(vr[outer_shell_index]) / config.implosion_speed_cm_s

    checks = {
        "interface_velocity_is_inward": interface_velocity < 0.0,
        "pressure_force_is_outward_at_interface": interface_acceleration > 0.0,
        "pressure_force_opposes_interface_motion": interface_velocity
        * interface_acceleration
        < 0.0,
        "hotspot_pressure_exceeds_inner_shell_pressure": 1.05 <= pressure_ratio <= 1.8,
        "shell_is_much_denser_than_hotspot_edge": shell_density_ratio >= 6.0,
        "hotspot_temperature_is_reasonable": 0.5
        <= rows[hotspot_index].Te_keV
        <= 5.0,
        "shell_temperature_is_cold": 0.02 <= rows[inner_shell_index].Te_keV <= 0.25,
        "outer_shell_is_faster_than_inner_shell": vr[outer_shell_index]
        < vr[inner_shell_index],
        "outer_shell_speed_is_omega_like": 0.75 <= shell_outer_speed_fraction <= 1.10,
        "interface_convergence_is_limited": interface_convergence < 2.5e10,
        "transition_width_is_not_too_sharp": width_um >= 3.0,
        "outer_boundary_velocity_matches_static_boundary": abs(vr[-1])
        <= 1.0e-3 * config.implosion_speed_cm_s,
    }
    metrics: dict[str, float | str] = {
        "diagnostic_id": "artificial_profile.deceleration_check",
        "velocity_model": config.velocity_model,
        "r0_um": config.r0_um,
        "shell_outer_um": config.shell_outer_um,
        "effective_smooth_width_um": width_um,
        "hotspot_probe_um": rows[hotspot_index].r_um,
        "inner_shell_probe_um": rows[inner_shell_index].r_um,
        "outer_shell_probe_um": rows[outer_shell_index].r_um,
        "interface_velocity_cm_s": interface_velocity,
        "interface_pressure_acceleration_cm_s2": interface_acceleration,
        "interface_convergence_1_s": interface_convergence,
        "hotspot_to_inner_shell_pressure_ratio": pressure_ratio,
        "inner_shell_to_hotspot_edge_density_ratio": shell_density_ratio,
        "hotspot_probe_te_kev": rows[hotspot_index].Te_keV,
        "inner_shell_probe_te_kev": rows[inner_shell_index].Te_keV,
        "inner_shell_speed_fraction_of_v0": shell_inner_speed_fraction,
        "outer_shell_speed_fraction_of_v0": shell_outer_speed_fraction,
        "pressure_max_gbar": max(pressure) / GBAR_TO_DYN_PER_CM2,
        "rho_peak_g_cm3": max(rho),
        "outer_boundary_velocity_cm_s": vr[-1],
    }
    return DecelerationDiagnostics(
        passed=all(checks.values()), checks=checks, metrics=metrics
    )


def write_diagnostics_json(path: Path, diagnostics: DecelerationDiagnostics) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "passed": diagnostics.passed,
        "checks": diagnostics.checks,
        "metrics": diagnostics.metrics,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def plot_profile_diagnostics(
    path: Path, rows: list[ProfileRow], config: ProfileConfig, diagnostics: DecelerationDiagnostics
) -> None:
    import matplotlib.pyplot as plt

    r_um = [row.r_um for row in rows]
    rho = [row.rho_g_cm3 for row in rows]
    te = [row.Te_keV for row in rows]
    ti = [row.Ti_keV for row in rows]
    vr = [row.vr_cm_s / 1.0e7 for row in rows]
    pressure = [_row_pressure_dyn_cm2(row) / GBAR_TO_DYN_PER_CM2 for row in rows]
    r_cm = [row.r_um * 1.0e-4 for row in rows]
    pressure_gradient = _gradient([p * GBAR_TO_DYN_PER_CM2 for p in pressure], r_cm)
    acceleration = [
        -pressure_gradient[index] / rho[index] / 1.0e18 for index in range(len(rows))
    ]
    r2v = [
        r_cm[index] * r_cm[index] * rows[index].vr_cm_s
        for index in range(len(rows))
    ]
    vr_gradient = _gradient([row.vr_cm_s for row in rows], r_cm)
    divergence = _gradient(r2v, r_cm)
    divergence = [
        (3.0 * vr_gradient[0])
        if r_cm[index] == 0.0
        else divergence[index] / (r_cm[index] * r_cm[index])
        for index in range(len(rows))
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.2), constrained_layout=True)
    axes_flat = list(axes.ravel())

    def mark_regions(ax) -> None:
        width_um = _transition_width_um(config)
        ax.axvline(config.r0_um, color="0.25", lw=1.0, ls="--")
        ax.axvline(config.shell_outer_um, color="0.35", lw=1.0, ls=":")
        ax.axvspan(config.r0_um - width_um, config.r0_um + width_um, color="0.9", zorder=-10)
        ax.axvspan(
            config.shell_outer_um - width_um,
            config.shell_outer_um + width_um,
            color="0.94",
            zorder=-10,
        )
        ax.set_xlim(0.0, config.r_max_um)
        ax.set_xlabel(r"$r$ ($\mu$m)")

    axes_flat[0].plot(r_um, rho, color="#1f77b4")
    axes_flat[0].set_ylabel(r"$\rho$ (g/cm$^3$)")
    axes_flat[0].set_title("Density")

    axes_flat[1].plot(r_um, pressure, color="#d62728")
    axes_flat[1].set_ylabel("Pressure (Gbar)")
    axes_flat[1].set_title("Total pressure")

    axes_flat[2].plot(r_um, te, color="#ff7f0e", label=r"$T_e$")
    axes_flat[2].plot(r_um, ti, color="#2ca02c", ls="--", label=r"$T_i$")
    axes_flat[2].set_ylabel("Temperature (keV)")
    axes_flat[2].set_title("Temperatures")
    axes_flat[2].legend(frameon=False)

    axes_flat[3].plot(r_um, vr, color="#9467bd")
    axes_flat[3].axhline(0.0, color="0.4", lw=0.8)
    axes_flat[3].set_ylabel(r"$v_r$ ($10^7$ cm/s)")
    axes_flat[3].set_title("Radial velocity")

    axes_flat[4].plot(r_um, acceleration, color="#8c564b")
    axes_flat[4].axhline(0.0, color="0.4", lw=0.8)
    axes_flat[4].set_ylabel(r"$a_P$ ($10^{18}$ cm/s$^2$)")
    axes_flat[4].set_title(r"Pressure force $a_P=-\nabla P/\rho$")

    axes_flat[5].plot(r_um, [value / 1.0e10 for value in divergence], color="#17becf")
    axes_flat[5].axhline(0.0, color="0.4", lw=0.8)
    axes_flat[5].set_ylabel(r"$\nabla\cdot v$ ($10^{10}$ s$^{-1}$)")
    axes_flat[5].set_title("Spherical velocity divergence")

    for ax in axes_flat:
        mark_regions(ax)
        ax.grid(True, alpha=0.25)

    status = "PASS" if diagnostics.passed else "FAIL"
    fig.suptitle(
        f"Artificial deceleration profile check: {status}  "
        f"r0={config.r0_um:g} um, w={_transition_width_um(config):g} um",
        fontsize=13,
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _format_float(value: float) -> str:
    return f"{value:.17e}"


def write_profile(path: Path, rows: Iterable[ProfileRow], config: ProfileConfig) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    materialized = list(rows)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# artificial DEC3D radial profile surrogate\n")
        out.write("# source_density=docs/artificial_initial_setup_2.md section 3.1\n")
        out.write("# source_pressure=docs/artificial_initial_setup_2.md sections 3.2-3.4\n")
        out.write("# source_velocity=docs/artificial_initial_setup_2.md section 3.5\n")
        out.write("# perturbation_shape=Woo_thesis_Eq_2_2_not_applied_to_1d_profile\n")
        out.write("# closure=Pe_equals_Pi_from_total_pressure_not_LILAC_exact\n")
        out.write(
            "# "
            f"r0_um={config.r0_um:.17g} ell={config.ell} "
            f"shell_outer_um={config.shell_outer_um:.17g} "
            f"velocity_model={config.velocity_model} "
            f"smooth_width_um={config.smooth_width_um:.17g} "
            f"effective_smooth_width_um={_transition_width_um(config):.17g} "
            f"hotspot_pressure_gbar={config.hotspot_pressure_gbar:.17g} "
            f"shell_adiabat={config.shell_adiabat:.17g} "
            f"implosion_speed_cm_s={config.implosion_speed_cm_s:.17g} "
            f"shell_inner_velocity_fraction={config.shell_inner_velocity_fraction:.17g} "
            f"shell_outer_velocity_fraction={config.shell_outer_velocity_fraction:.17g}\n"
        )
        out.write(" ".join(HEADER) + "\n")
        for row in materialized:
            out.write(
                " ".join(
                    _format_float(value)
                    for value in (
                        row.r_um,
                        row.rho_g_cm3,
                        row.Te_keV,
                        row.Ti_keV,
                        row.vr_cm_s,
                        row.vt_cm_s,
                        row.vp_cm_s,
                        row.epsilon_alpha_erg_cm3,
                        row.radiation_scale,
                    )
                )
                + "\n"
            )


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate an artificial DEC3D .pro radial initial profile."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("cases/artificial_profile_1024.pro"),
        help="Output .pro path.",
    )
    parser.add_argument("--points", type=int, default=1025, help="Number of radial rows.")
    parser.add_argument("--r-max-um", type=float, default=115.0, help="Outer profile radius.")
    parser.add_argument("--r0-um", type=float, default=68.0, help="Perturbation center radius.")
    parser.add_argument("--shell-outer-um", type=float, default=86.0, help="Compressed shell outer radius.")
    parser.add_argument("--ell", type=int, default=2, help="Legendre mode for Eq. 2.2 metadata.")
    parser.add_argument("--smooth-width-cells", type=float, default=3.0, help="Interface width in grid cells.")
    parser.add_argument(
        "--smooth-width-um",
        type=float,
        default=4.0,
        help="Minimum physical interface width in microns.",
    )
    parser.add_argument("--hotspot-center-density", type=float, default=1.5, help="Central hot-spot density.")
    parser.add_argument("--hotspot-density-slope", type=float, default=1.0, help="Quadratic hot-spot density slope.")
    parser.add_argument("--shell-inner-density", type=float, default=35.0, help="Compressed shell density at R_if.")
    parser.add_argument("--shell-density-gradient", type=float, default=0.25, help="Fractional shell density drop.")
    parser.add_argument("--outer-density", type=float, default=0.03, help="Outer buffer density.")
    parser.add_argument("--hotspot-pressure-gbar", type=float, default=3.0, help="Uniform hot-spot pressure.")
    parser.add_argument("--shell-adiabat", type=float, default=3.2, help="Shell adiabat multiplying Fermi pressure.")
    parser.add_argument("--outer-pressure-gbar", type=float, default=0.003, help="Outer buffer pressure.")
    parser.add_argument(
        "--implosion-speed",
        type=float,
        default=3.0e7,
        help="Characteristic inward shell speed in cm/s.",
    )
    parser.add_argument("--hotspot-velocity-fraction", type=float, default=0.45, help="Hot-spot velocity fraction.")
    parser.add_argument("--shell-velocity-slope", type=float, default=0.1, help="Literal setup2 shell velocity fractional slope.")
    parser.add_argument(
        "--velocity-model",
        choices=("decelerating-shell", "setup2-original"),
        default="decelerating-shell",
        help="Velocity model for the 1D base profile.",
    )
    parser.add_argument(
        "--shell-inner-velocity-fraction",
        type=float,
        default=0.65,
        help="Inner shell inward speed fraction for decelerating-shell mode.",
    )
    parser.add_argument(
        "--shell-outer-velocity-fraction",
        type=float,
        default=1.0,
        help="Outer shell inward speed fraction for decelerating-shell mode.",
    )
    parser.add_argument(
        "--diagnostics-json",
        type=Path,
        default=None,
        help="Optional path for deceleration check diagnostics.",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Optional path for the profile diagnostic plot.",
    )
    parser.add_argument(
        "--require-deceleration",
        action="store_true",
        help="Exit non-zero if the deceleration check fails.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output file.")
    return parser


def config_from_args(args: argparse.Namespace) -> ProfileConfig:
    return ProfileConfig(
        points=args.points,
        r_max_um=args.r_max_um,
        r0_um=args.r0_um,
        shell_outer_um=args.shell_outer_um,
        ell=args.ell,
        smooth_width_cells=args.smooth_width_cells,
        smooth_width_um=args.smooth_width_um,
        hotspot_center_density_g_cm3=args.hotspot_center_density,
        hotspot_density_slope=args.hotspot_density_slope,
        shell_inner_density_g_cm3=args.shell_inner_density,
        shell_density_gradient=args.shell_density_gradient,
        outer_density_g_cm3=args.outer_density,
        hotspot_pressure_gbar=args.hotspot_pressure_gbar,
        shell_adiabat=args.shell_adiabat,
        outer_pressure_gbar=args.outer_pressure_gbar,
        implosion_speed_cm_s=args.implosion_speed,
        hotspot_velocity_fraction=args.hotspot_velocity_fraction,
        shell_velocity_slope=args.shell_velocity_slope,
        velocity_model=args.velocity_model,
        shell_inner_velocity_fraction=args.shell_inner_velocity_fraction,
        shell_outer_velocity_fraction=args.shell_outer_velocity_fraction,
    )


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    output = args.output
    if output.exists() and not args.force:
        parser.error(f"output exists: {output}; pass --force to overwrite")

    config = config_from_args(args)
    rows = generate_profile(config)
    diagnostics = deceleration_diagnostics(rows, config)
    write_profile(output, rows, config)
    if args.diagnostics_json is not None:
        write_diagnostics_json(args.diagnostics_json, diagnostics)
    if args.plot is not None:
        plot_profile_diagnostics(args.plot, rows, config, diagnostics)
    peak = max(rows, key=lambda row: row.rho_g_cm3)
    center_pressure_gbar = pressure_dyn_cm2(rows[0].r_um, config) / GBAR_TO_DYN_PER_CM2
    print(f"output={output}")
    print(f"rows={len(rows)}")
    print(f"r_range_um={rows[0].r_um:.8g},{rows[-1].r_um:.8g}")
    print(f"rho_peak_g_cm3={peak.rho_g_cm3:.8g}")
    print(f"rho_peak_r_um={peak.r_um:.8g}")
    print(f"pressure_center_gbar={center_pressure_gbar:.8g}")
    print(f"effective_smooth_width_um={_transition_width_um(config):.8g}")
    print(f"r0_um={config.r0_um:.8g}")
    print(f"f_r0={woo_shape_function(config.r0_um, config.r0_um, config.ell):.8g}")
    print(f"deceleration_check_passed={str(diagnostics.passed).lower()}")
    for name, passed in diagnostics.checks.items():
        print(f"check.{name}={str(passed).lower()}")
    for name, value in diagnostics.metrics.items():
        if isinstance(value, float):
            print(f"metric.{name}={value:.8g}")
        else:
            print(f"metric.{name}={value}")
    if args.require_deceleration and not diagnostics.passed:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
