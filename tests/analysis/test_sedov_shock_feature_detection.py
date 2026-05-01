from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis" / "scripts"))

from plot_hydro_fields import detect_shock_features_from_profile, recompute_shock_radius_from_profile


def assert_close(actual: float, expected: float, tolerance: float = 1.0e-12) -> None:
    if abs(actual - expected) > tolerance:
        raise AssertionError(f"expected {expected:.17g}, got {actual:.17g}")


def assert_between(value: float, lower: float, upper: float) -> None:
    if not lower <= value <= upper:
        raise AssertionError(f"expected {lower:.17g} <= {value:.17g} <= {upper:.17g}")


def main() -> int:
    radial_centers = [0.07, 0.08, 0.09, 0.10, 0.11, 0.12]
    rho_values = [0.20, 1.00, 2.50, 4.00, 1.00, 1.00]

    features = detect_shock_features_from_profile(radial_centers, rho_values)

    assert_between(features.inner_rise_radius, 0.09, 0.10)
    assert_between(features.shock_radius, 0.10, 0.11)
    assert_close(features.density_peak_radius, 0.10)
    if not features.inner_rise_radius < features.density_peak_radius < features.shock_radius:
        raise AssertionError("expected inner rise, density peak, and outer shock to be ordered radially")
    assert_close(recompute_shock_radius_from_profile(radial_centers, rho_values), features.shock_radius)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
