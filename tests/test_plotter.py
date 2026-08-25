
import unittest
import sys
import os

# Add plugin root to path
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

# Mock wx if strictly necessary, but let's try assuming it might be present or we can mock it
# For CI/headless correctness without wx installed, we often mock wx.
# Let's simple-mock wx to ensure logic runs even if system python doesn't have wx (though user python likely does)
import types
from types import SimpleNamespace
if 'wx' not in sys.modules:
    wx_mock = types.ModuleType('wx')
    wx_mock.Bitmap = lambda *args: "BITMAP_OBJECT"
    wx_mock.Image = lambda *args: "IMAGE_OBJECT"
    wx_mock.BITMAP_TYPE_PNG = 1
    sys.modules['wx'] = wx_mock

from plotter import Plotter
from mesh import Mesh

class TestPlotter(unittest.TestCase):
    def setUp(self):
        self.plotter = Plotter(debug=True)
        self.mesh = Mesh()
        # Create a simple dummy mesh
        # Nodes 0,1 on layer 0
        # Nodes 2,3 on layer 1
        self.mesh.nodes = [0, 1, 2, 3]
        self.mesh.node_coords = {
            0: (0.0, 0.0, 0),
            1: (1.0, 1.0, 0),
            2: (0.0, 0.0, 1),
            3: (1.0, 1.0, 1)
        }
        self.mesh.results = {
            0: 3.3,
            1: 3.2,
            2: 3.3,
            3: 3.1
        }
        self.stackup = {
            'copper': {0: {}, 1: {}},
            'layer_order': [0, 1]
        }

    def test_plot_3d(self):
        # Should return a bitmap (or mock string)
        bmp = self.plotter.plot_3d_mesh(self.mesh, self.stackup)
        print(f"3D Plot result: {bmp}")
        self.assertIsNotNone(bmp)
        
    def test_plot_2d_layer(self):
        # Should return a bitmap for layer 0
        bmp = self.plotter.plot_layer_2d(self.mesh, 0, self.stackup, vmin=3.0, vmax=3.5, layer_name="F.Cu (Test)")
        print(f"2D Plot Layer 0 result: {bmp}")
        self.assertIsNotNone(bmp)
        
    def test_plot_2d_empty_layer(self):
        # Layer 99 empty
        bmp = self.plotter.plot_layer_2d(self.mesh, 99, self.stackup)
        self.assertIsNone(bmp)

    def test_plot_impedance_sweep(self):
        baseline = SimpleNamespace(
            frequencies_hz=[1e3, 1e4, 1e5],
            impedance_ohm=[0.05 + 0j, 0.04 - 0.01j, 0.08 + 0.02j],
            target_impedance_ohm=0.06,
        )
        optimized = SimpleNamespace(
            frequencies_hz=baseline.frequencies_hz,
            impedance_ohm=[0.04 + 0j, 0.03 - 0.01j, 0.05 + 0.01j],
        )

        bitmap = self.plotter.plot_impedance_sweep(baseline, optimized)

        self.assertIsNotNone(bitmap)

    def test_plot_thermal_views(self):
        thermal_mesh = SimpleNamespace(
            nodes=[0, 1, 2, 3],
            node_coords={
                0: (0.0, 0.0, 0.0),
                1: (1.0, 0.0, 0.0),
                2: (0.0, 1.0, 1.6),
                3: (1.0, 1.0, 1.6),
            },
            layer_specs=[
                SimpleNamespace(layer_id=0),
                SimpleNamespace(layer_id=31),
            ],
            node_layers={0: 0, 1: 0, 2: 31, 3: 31},
            node_map={
                (0, 0, 0): 0,
                (1, 0, 0): 1,
                (0, 1, 1): 2,
                (1, 1, 1): 3,
            },
        )
        result = SimpleNamespace(
            temperatures_c={0: 28.0, 1: 30.0, 2: 35.0, 3: 32.0},
        )

        self.assertIsNotNone(self.plotter.plot_thermal_3d(thermal_mesh, result))
        self.assertIsNotNone(self.plotter.plot_thermal_surface(thermal_mesh, result, side="TOP"))
        self.assertIsNotNone(self.plotter.plot_thermal_surface(thermal_mesh, result, side="BOTTOM"))

        thermal_3d_png = self.plotter.plot_thermal_3d(thermal_mesh, result, as_png=True)
        top_png = self.plotter.plot_thermal_surface(
            thermal_mesh, result, side="TOP", as_png=True
        )
        self.assertIsInstance(thermal_3d_png, bytes)
        self.assertIsInstance(top_png, bytes)
        self.assertTrue(thermal_3d_png.startswith(b"\x89PNG"))
        self.assertTrue(top_png.startswith(b"\x89PNG"))

    def test_plot_cfd_views(self):
        shape = (3, 3, 3)
        count = 27
        cfd_mesh = SimpleNamespace(
            shape=shape,
            spacing_m=(0.01, 0.01, 0.01),
            dimensions_m=(0.03, 0.03, 0.03),
        )
        air = [25.0 + index / 100.0 for index in range(count)]
        solid = [float('nan')] * count
        solid[13] = 31.0
        result = SimpleNamespace(
            pressure_pa=[0.1 * index for index in range(count)],
            velocity_u_m_s=[0.2] * count,
            velocity_v_m_s=[0.0] * count,
            velocity_w_m_s=[0.01] * count,
            air_temperature_c=air,
            solid_temperature_c=solid,
            residuals=SimpleNamespace(
                continuity=[1.0, 0.1, 0.01],
                momentum=[0.5, 0.05, 0.005],
                energy=[2.0, 0.2, 0.02],
            ),
        )

        self.assertIsNotNone(self.plotter.plot_cfd_3d(cfd_mesh, result))
        self.assertIsNotNone(self.plotter.plot_cfd_slice(cfd_mesh, result, "TEMPERATURE", "XY"))
        self.assertIsNotNone(self.plotter.plot_cfd_slice(cfd_mesh, result, "VELOCITY", "XZ"))
        self.assertIsNotNone(self.plotter.plot_cfd_slice(cfd_mesh, result, "PRESSURE", "YZ"))
        self.assertIsNotNone(self.plotter.plot_cfd_residuals(result))

if __name__ == '__main__':
    unittest.main()
