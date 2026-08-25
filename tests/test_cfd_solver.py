import math
import unittest

import numpy as np

from cfd_mesh import CFDMeshGenerator
from cfd_model import CFDObstacle, EnclosureModel
from cfd_solver import EnclosureCFDSolver
from models import CFDBoundaryPatch, EnclosureCFDSettings


class TestEnclosureCFDSolver(unittest.TestCase):
    def _settings(self):
        settings = EnclosureCFDSettings()
        settings.geometry.width_mm = 30.0
        settings.geometry.depth_mm = 30.0
        settings.geometry.height_mm = 30.0
        settings.geometry.wall_heat_transfer_w_m2k = 8.0
        settings.solver.cell_size_mm = 6.0
        settings.solver.max_iterations = 12
        settings.solver.pressure_iterations = 12
        settings.solver.tolerance = 1e-8
        settings.solver.pseudo_time_step_s = 0.002
        return settings

    def test_heated_solid_raises_solid_and_air_temperature(self):
        settings = self._settings()
        model = EnclosureModel(
            dimensions_mm=(30.0, 30.0, 30.0),
            obstacles=[CFDObstacle("HOT", (12, 12, 12, 18, 18, 18), 5.0, 0.5)],
        )
        mesh = CFDMeshGenerator().generate_mesh(model, settings)

        result = EnclosureCFDSolver().solve(mesh, settings)

        self.assertEqual(len(result.pressure_pa), mesh.cell_count)
        self.assertTrue(all(math.isfinite(value) for value in result.pressure_pa))
        self.assertGreater(result.maximum_solid_temperature_c, settings.ambient_c)
        self.assertGreater(result.maximum_air_temperature_c, settings.ambient_c)
        self.assertAlmostEqual(result.total_heat_w, 0.5)
        self.assertEqual(len(result.residuals.continuity), result.iterations)
        self.assertTrue(math.isfinite(result.energy_balance_error_pct))

    def test_forced_inlet_produces_a_finite_velocity_field(self):
        settings = self._settings()
        settings.solver.include_buoyancy = False
        inlet = CFDBoundaryPatch("Fan", "FAN", "XMIN", 0.5, 0.5, 0.5, 0.5, 0.8, 25.0)
        outlet = CFDBoundaryPatch("Outlet", "OUTLET", "XMAX", 0.5, 0.5, 0.5, 0.5)
        model = EnclosureModel((30.0, 30.0, 30.0), patches=[inlet, outlet])
        mesh = CFDMeshGenerator().generate_mesh(model, settings)

        result = EnclosureCFDSolver().solve(mesh, settings)

        self.assertGreater(result.maximum_velocity_m_s, 0.0)
        self.assertTrue(math.isfinite(result.mass_balance_error_pct))
        self.assertTrue(np.isfinite(np.asarray(result.velocity_u_m_s)).all())

    def test_cancellation_stops_before_iteration(self):
        settings = self._settings()
        model = EnclosureModel((30.0, 30.0, 30.0))
        mesh = CFDMeshGenerator().generate_mesh(model, settings)

        with self.assertRaisesRegex(RuntimeError, "cancelled"):
            EnclosureCFDSolver().solve(mesh, settings, cancel_callback=lambda: True)


if __name__ == "__main__":
    unittest.main()
