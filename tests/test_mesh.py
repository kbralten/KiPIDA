
import unittest
import sys
import os

try:
    import numpy as np
    from shapely.geometry import Polygon
except ImportError:
    pass

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

# 1. Pcbnew Mock (Base for others)
class MockPcbnew:
    PAD_ATTRIB_PTH = 0
    class ActionPlugin: pass
    class PCB_VIA: pass
    class PCB_TRACK: pass
    @staticmethod
    def ToMM(v): return v / 1e6

if 'pcbnew' not in sys.modules:
    sys.modules['pcbnew'] = MockPcbnew

import pcbnew # Now we can import it

# 2. Other Mocks
class MockBoard:
    def __init__(self):
        self.tracks = []
        self.footprints = []
        self.layers_copper = {0: "F_Cu", 31: "B_Cu"}
    def FindNet(self, name): return MockNet(1)
    def GetNetCode(self): return 1
    def GetTracks(self): return self.tracks
    def GetFootprints(self): return self.footprints
    def GetEnabledLayers(self): return MockSeq([0, 31])
    def IsLayerCopper(self, id): return id in self.layers_copper

class MockNet:
    def __init__(self, c): self.c = c
    def GetNetCode(self): return self.c

class MockSeq:
    def __init__(self, s): self.s = s
    def Seq(self): return self.s
    
class MockTrack(MockPcbnew.PCB_VIA):
    def __init__(self, net, layer, pos, w, cls):
        self.net = net
        self.pos = pos
        self.w = w
        self.cls = cls
        self.top = 0
        self.bot = 31
    def GetNetCode(self): return self.net
    def GetClass(self): return self.cls
    def GetPosition(self): return self.pos
    def GetWidth(self): return self.w
    def TopLayer(self): return self.top
    def BottomLayer(self): return self.bot
    def GetLayerSet(self):
        class MockLayerSet:
            def Seq(self): return [0, 31]
        return MockLayerSet()

class MockPad:
    def __init__(self, net, pos, size, attrib):
        self.net = net
        self.pos = pos
        self.size = size
        self.attrib = attrib
    def GetNetCode(self): return self.net
    def GetPosition(self): return self.pos
    def GetSize(self): return self.size
    def GetAttribute(self): return self.attrib

class MockFootprint:
    def __init__(self, pads): self.pads = pads
    def Pads(self): return self.pads

class MockPoint:
    def __init__(self, x, y): self.x, self.y = x, y

from mesh import Mesher, Mesh

class TestMesher(unittest.TestCase):
    def setUp(self):
        self.board = MockBoard()
        self.mesher = Mesher(self.board)
        
    def test_simple_rect_mesh(self):
        # 20x10 mm rect on layer 0
        poly = Polygon([(0,0), (20,0), (20,10), (0,10)])
        geo = {0: poly}
        stackup = {0: {'thickness_mm': 0.035, 'resistivity': 1}}
        
        mesh = self.mesher.generate_mesh("Test", geo, stackup, grid_size_mm=5.0)
        
        self.assertEqual(len(mesh.nodes), 15)
        self.assertTrue(len(mesh.edges) > 0)
        
    def test_vertical_link_via(self):
        # 2 layers, node at 0,0 on both
        poly = Polygon([(-1,-1), (1,-1), (1,1), (-1,1)]) # Small square at origin
        geo = {0: poly, 31: poly}
        stackup = {0: {}, 31: {}}
        
        # Add Via at 0,0
        via = MockTrack(1, 0, MockPoint(0,0), 300000, "PCB_VIA")
        self.board.tracks.append(via)
        
        mesh = self.mesher.generate_mesh("Test", geo, stackup, grid_size_mm=0.5)
        
        # Just check edges for vertical type
        vertical_edges = [e for e in mesh.edges if e[2] > 500] # G > 500
        self.assertTrue(len(vertical_edges) >= 1)

if __name__ == '__main__':
    unittest.main()
