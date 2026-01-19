
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

if __name__ == '__main__':
    unittest.main()
