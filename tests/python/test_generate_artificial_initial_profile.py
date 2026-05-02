import importlib.util
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "generate_artificial_initial_profile.py"
PYTHON = r"C:\Users\Administrator\anaconda3\envs\spyder55-pip\python.exe"
ERG_PER_KEV = 1.602176634e-9
MEAN_DT_ION_MASS_G = 0.5 * (3.3435837724e-24 + 5.0073567446e-24)


def load_generator_module():
    spec = importlib.util.spec_from_file_location("generate_artificial_initial_profile", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def total_pressure_gbar(row):
    ion_number_density = row.rho_g_cm3 / MEAN_DT_ION_MASS_G
    electron_number_density = ion_number_density
    pressure = (
        electron_number_density * row.Te_keV * ERG_PER_KEV
        + ion_number_density * row.Ti_keV * ERG_PER_KEV
    )
    return pressure / 1.0e15


class ArtificialInitialProfileGeneratorTest(unittest.TestCase):
    def test_woo_shape_function_is_centered_and_decays(self):
        generator = load_generator_module()

        self.assertAlmostEqual(generator.woo_shape_function(68.0, 68.0, 2), 1.0, places=12)
        self.assertLess(generator.woo_shape_function(10.0, 68.0, 2), 0.03)
        self.assertLess(generator.woo_shape_function(120.0, 68.0, 2), 0.35)

    def test_default_surrogate_profile_matches_artificial_setup_document(self):
        generator = load_generator_module()
        config = generator.ProfileConfig(points=121)
        rows = generator.generate_profile(config)

        self.assertEqual(len(rows), 121)
        self.assertAlmostEqual(rows[0].r_um, 0.0)
        self.assertAlmostEqual(rows[-1].r_um, 115.0)
        self.assertTrue(all(right.r_um > left.r_um for left, right in zip(rows, rows[1:])))
        self.assertTrue(all(row.rho_g_cm3 > 0.0 for row in rows))
        self.assertTrue(all(row.Te_keV > 0.0 and row.Ti_keV > 0.0 for row in rows))
        self.assertTrue(all(math.isfinite(row.vr_cm_s) for row in rows))
        self.assertAlmostEqual(rows[0].vr_cm_s, 0.0, delta=1.0e5)

        self.assertGreater(rows[0].rho_g_cm3, 1.0)
        self.assertLess(rows[0].rho_g_cm3, 2.0)
        self.assertGreater(total_pressure_gbar(rows[0]), 2.5)
        self.assertLess(total_pressure_gbar(rows[0]), 4.0)

        peak = max(rows, key=lambda row: row.rho_g_cm3)
        self.assertGreater(peak.rho_g_cm3, 25.0)
        self.assertLess(peak.rho_g_cm3, 40.0)
        self.assertGreater(peak.r_um, 68.0)
        self.assertLess(peak.r_um, 88.0)

    def test_default_interfaces_have_physical_minimum_smoothing(self):
        generator = load_generator_module()
        config = generator.ProfileConfig(points=1025)

        self.assertAlmostEqual(generator._transition_width_um(config), 4.0)

    def test_default_profile_satisfies_deceleration_checks(self):
        generator = load_generator_module()
        config = generator.ProfileConfig(points=1025)
        rows = generator.generate_profile(config)
        diagnostics = generator.deceleration_diagnostics(rows, config)

        self.assertTrue(diagnostics.passed, msg=diagnostics)
        self.assertGreater(
            diagnostics.metrics["hotspot_to_inner_shell_pressure_ratio"], 1.05
        )
        self.assertLess(
            diagnostics.metrics["interface_convergence_1_s"], 2.5e10
        )
        self.assertAlmostEqual(
            diagnostics.metrics["outer_boundary_velocity_cm_s"], 0.0, delta=3.0e4
        )
        self.assertGreater(
            diagnostics.metrics["outer_shell_speed_fraction_of_v0"],
            diagnostics.metrics["inner_shell_speed_fraction_of_v0"],
        )

    def test_write_profile_uses_dec3d_pro_columns(self):
        generator = load_generator_module()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "artificial.pro"
            rows = generator.generate_profile(generator.ProfileConfig(points=9))
            generator.write_profile(path, rows, generator.ProfileConfig(points=9))

            lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
            header = next(line for line in lines if line and not line.startswith("#"))
            self.assertEqual(
                header,
                "r_um rho_g_cm3 Te_keV Ti_keV vr_cm_s vt_cm_s vp_cm_s "
                "epsilon_alpha_erg_cm3 radiation_scale",
            )
            data_lines = [line for line in lines if line and not line.startswith("#")][1:]
            self.assertEqual(len(data_lines), 9)
            self.assertEqual(len(data_lines[0].split()), 9)

    def test_cli_generates_profile_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "cli.pro"
            result = subprocess.run(
                [
                    PYTHON,
                    str(SCRIPT_PATH),
                    "--output",
                    str(output),
                    "--points",
                    "33",
                    "--force",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, msg=result.stdout + result.stderr)
            self.assertTrue(output.exists())
            self.assertIn("rows=33", result.stdout)


if __name__ == "__main__":
    unittest.main()
