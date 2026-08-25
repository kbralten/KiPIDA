import os
import sys
import unittest

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

try:
    from shapely.geometry import box
    from electrothermal import ElectroThermalSolver
    from mesh import Mesh
    from models import AirflowSettings, ThermalAnalysisSettings
    from thermal_mesh import ThermalMesher
    from thermal_model import ThermalBoardModel
    ELECTROTHERMAL_AVAILABLE = True
except ImportError:
    ELECTROTHERMAL_AVAILABLE = False


@unittest.skipUnless(ELECTROTHERMAL_AVAILABLE, "thermal numerical dependencies unavailable")
class TestElectroThermalSolver(unittest.TestCase):
    def test_dc_loss_heats_board_and_updates_resistance(self):
        settings = ThermalAnalysisSettings(
            ambient_c=25.0,
            grid_size_mm=2.0,
            airflow=AirflowSettings(
                mode="CUSTOM", custom_h_w_m2k=500.0,
                expose_top=True, expose_bottom=True, expose_edges=True,
            ),
            include_radiation=False,
            coupled_iterations=10,
            convergence_c=0.01,
            relaxation=1.0,
        )
        board_model = ThermalBoardModel(
            bounds_mm=(0.0, 0.0, 10.0, 10.0),
            outline=box(0.0, 0.0, 10.0, 10.0),
            stackup={
                "copper": {
                    0: {"name": "F.Cu", "thickness_mm": 0.035},
                    31: {"name": "B.Cu", "thickness_mm": 0.035},
                },
                "layer_order": [0, 31],
                "substrate": [{"between": [0, 31], "thickness_mm": 1.53}],
            },
            copper_by_layer={
                0: box(0.0, 0.0, 10.0, 10.0),
                31: box(0.0, 0.0, 10.0, 10.0),
            },
        )
        thermal_mesh = ThermalMesher().generate_mesh(board_model, settings)

        electrical_mesh = Mesh()
        electrical_mesh.nodes = [0, 1]
        electrical_mesh.node_coords = {0: (2.0, 5.0, 0), 1: (8.0, 5.0, 0)}
        electrical_mesh.add_edge_direct(0, 1, g=1.0)
        contexts = {
            "5V": {
                "mesh": electrical_mesh,
                "sources": [{"node_id": 0, "voltage": 5.0}],
                "loads": [{"node_id": 1, "current": 1.0}],
            }
        }

        result = ElectroThermalSolver().solve(thermal_mesh, settings, contexts)

        self.assertGreater(result.thermal.hotspot.temperature_c, settings.ambient_c)
        self.assertGreater(result.dc_results["5V"].total_loss_w, 1.0)
        self.assertLessEqual(result.iterations, settings.coupled_iterations)
        self.assertTrue(result.converged)


if __name__ == "__main__":
    unittest.main()
