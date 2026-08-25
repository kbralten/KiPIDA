import unittest

import numpy as np

from cfd_mesh import CFDMeshGenerator
from cfd_model import CFDObstacle, EnclosureModel
from models import CFDBoundaryPatch, EnclosureCFDSettings


class TestCFDMeshGenerator(unittest.TestCase):
    def _settings(self, cell=5.0):
        settings = EnclosureCFDSettings()
        settings.solver.cell_size_mm = cell
        return settings

    def test_structured_mesh_maps_solids_heat_and_boundary_patches(self):
        model = EnclosureModel(
            dimensions_mm=(30.0, 20.0, 10.0),
            obstacles=[CFDObstacle("HOT", (10, 5, 2, 20, 15, 8), 8.0, 2.5)],
            patches=[CFDBoundaryPatch(
                "Inlet", "INLET", "XMIN", 0.5, 0.5, 0.5, 0.5, 1.0, 25.0,
            ), CFDBoundaryPatch(
                "Outlet", "OUTLET", "XMAX", 0.5, 0.5, 0.5, 0.5,
            )],
        )

        mesh = CFDMeshGenerator().generate_mesh(model, self._settings())

        self.assertEqual(mesh.shape, (6, 4, 3))
        self.assertGreater(np.count_nonzero(mesh.solid_mask), 0)
        self.assertAlmostEqual(float(np.sum(mesh.heat_sources_w)), 2.5)
        self.assertTrue(mesh.patch_cells[0].cells)
        self.assertTrue(all(cell[0] == 0 for cell in mesh.patch_cells[0].cells))
        self.assertTrue(all(mesh.fluid_mask[cell] for cell in mesh.patch_cells[0].cells))

    def test_subcell_obstacle_is_stamped_into_nearest_cell(self):
        model = EnclosureModel(
            dimensions_mm=(30.0, 30.0, 30.0),
            obstacles=[CFDObstacle("TINY", (14.9, 14.9, 14.9, 15.1, 15.1, 15.1), 1.0, 1.0)],
        )
        mesh = CFDMeshGenerator().generate_mesh(model, self._settings(10.0))

        self.assertEqual(np.count_nonzero(mesh.solid_mask), 1)
        self.assertAlmostEqual(float(np.sum(mesh.heat_sources_w)), 1.0)

    def test_overlapping_boundary_patches_are_rejected(self):
        patch = CFDBoundaryPatch("A", "VENT", "XMIN", 0.5, 0.5, 0.5, 0.5)
        model = EnclosureModel(
            dimensions_mm=(30.0, 30.0, 30.0),
            patches=[patch, CFDBoundaryPatch("B", "OUTLET", "XMIN", 0.5, 0.5, 0.5, 0.5)],
        )

        with self.assertRaisesRegex(ValueError, "overlaps"):
            CFDMeshGenerator().generate_mesh(model, self._settings())

    def test_maximum_cell_guard_prevents_oversized_run(self):
        model = EnclosureModel(dimensions_mm=(100.0, 100.0, 100.0))
        settings = self._settings(1.0)
        settings.solver.max_cells = 10_000

        with self.assertRaisesRegex(ValueError, "Increase the CFD cell size"):
            CFDMeshGenerator().generate_mesh(model, settings)

    def test_invalid_inlet_velocity_is_rejected(self):
        model = EnclosureModel(
            dimensions_mm=(30.0, 30.0, 30.0),
            patches=[CFDBoundaryPatch("Bad fan", "FAN", "XMIN", velocity_m_s=0.0)],
        )

        with self.assertRaisesRegex(ValueError, "positive velocity"):
            CFDMeshGenerator().generate_mesh(model, self._settings())

    def test_forced_inlet_requires_an_exhaust_patch(self):
        model = EnclosureModel(
            dimensions_mm=(30.0, 30.0, 30.0),
            patches=[CFDBoundaryPatch(
                "Fan", "FAN", "XMIN", velocity_m_s=1.0,
            )],
        )

        with self.assertRaisesRegex(ValueError, "OUTLET or VENT"):
            CFDMeshGenerator().generate_mesh(model, self._settings())


if __name__ == "__main__":
    unittest.main()
