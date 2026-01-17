
import unittest
import sys
import os

# Inject Mock modules to run this standalone (without KiCad) if needed
# But for Solver, we mostly need numpy/scipy, which we assume are available in the env we run.
# The solver doesn't depend on pcbnew directly, only the Mesh object structure.

class MockMesh:
    def __init__(self):
        self.nodes = []
        self.edges = [] 

class TestSolver(unittest.TestCase):
    def setUp(self):
        # Ensure we can import solver
        plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        if plugin_dir not in sys.path:
            sys.path.insert(0, plugin_dir)
            
        from solver import Solver
        self.solver_class = Solver

    def test_series_resistors(self):
        """
        Test 3 nodes in series:
        (Src=10V) --[R1=1ohm]-- (Node1) --[R2=1ohm]-- (Load=1A)
        
        Edges: (0, 1, G=1), (1, 2, G=1)
        Src at 0: 10V
        Load at 2: 1.0A
        
        Expected:
        Current I = 1A flows through entire string.
        V2 = ? 
        Voltage drop across R2 = I*R = 1*1 = 1V.
        Voltage drop across R1 = I*R = 1*1 = 1V.
        
        Wait, load at Node 2 sinks 1A.
        KCL Node 2: (V1-V2)*G = 1A -> V1-V2 = 1 -> V2 = V1-1
        KCL Node 1: (V0-V1)*G + (V2-V1)*G = 0 
                     (10-V1) + (V1-1-V1) = 0 ?? No
                     
        Let's solve by hand:
        V0 = 10
        I_path = 1A (Since load sinks 1A at end and no other path)
        V_drop_R1 = 1A * 1ohm = 1V => V1 = 9V
        V_drop_R2 = 1A * 1ohm = 1V => V2 = 8V
        """
        solver = self.solver_class()
        mesh = MockMesh()
        mesh.nodes = [0, 1, 2]
        # Conductance 1.0 = Resistance 1.0
        mesh.edges = [
            (0, 1, 1.0),
            (1, 2, 1.0)
        ]
        
        sources = [{'node_id': 0, 'voltage': 10.0}]
        loads = [{'node_id': 2, 'current': 1.0}]
        
        result = solver.solve(mesh, sources, loads)
        
        print(f"Series Result: {result}")
        
        self.assertAlmostEqual(result[0], 10.0, places=5)
        self.assertAlmostEqual(result[1], 9.0, places=5)
        self.assertAlmostEqual(result[2], 8.0, places=5)

    def test_parallel_resistors(self):
        """
        (Src=10V) ----+----[R1=1]---- (Node1)
                      |
                      +----[R2=1]---- (Node1)
                      
        Wait, parallel resistors between same two nodes?
        (Src=10V @ 0) --[R1=1]-- (Node1)
                      --[R2=1]--
                      
        Edges: (0, 1, 1.0), (0, 1, 1.0). Total G = 2.0.
        Combined R = 0.5.
        
        Load 2A at Node 1.
        
        V_drop = I * R = 2A * 0.5ohm = 1V.
        V1 should be 9V.
        """
        solver = self.solver_class()
        mesh = MockMesh()
        mesh.nodes = [0, 1]
        mesh.edges = [
            (0, 1, 1.0),
            (0, 1, 1.0)
        ]
        
        sources = [{'node_id': 0, 'voltage': 10.0}]
        loads = [{'node_id': 1, 'current': 2.0}]
        
        result = solver.solve(mesh, sources, loads)
        
        print(f"Parallel Result: {result}")
        
        self.assertAlmostEqual(result[0], 10.0, places=5)
        self.assertAlmostEqual(result[1], 9.0, places=5)

if __name__ == '__main__':
    unittest.main()
