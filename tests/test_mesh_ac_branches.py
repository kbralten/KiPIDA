import os
import sys
import unittest

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from mesh import Mesh


class TestMeshACBranches(unittest.TestCase):
    def test_edge_retains_r_and_l_without_changing_dc_stamp(self):
        mesh = Mesh()
        mesh.add_edge_direct(0, 1, 4.0, inductance_h=0.5e-9, kind="via")

        self.assertEqual(mesh.G_coo_row, [0, 1, 0, 1])
        self.assertEqual(mesh.G_coo_col, [0, 1, 1, 0])
        self.assertEqual(mesh.G_coo_data, [4.0, 4.0, -4.0, -4.0])
        self.assertEqual(len(mesh.branches), 1)
        self.assertAlmostEqual(mesh.branches[0].resistance_ohm, 0.25)
        self.assertAlmostEqual(mesh.branches[0].inductance_h, 0.5e-9)
        self.assertEqual(mesh.branches[0].kind, "via")


if __name__ == "__main__":
    unittest.main()
