import os
import sys
import unittest

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

try:
    from shapely.geometry import box
    from thermal_mesh import ThermalMesher
    from thermal_model import ThermalBoardModel, ThermalPlacement, ThermalVia
    from thermal_solver import ThermalSolver
    from models import AirflowSettings, ThermalAnalysisSettings, ThermalComponentModel
    THERMAL_AVAILABLE = True
except ImportError:
    THERMAL_AVAILABLE = False


@unittest.skipUnless(THERMAL_AVAILABLE, "thermal numerical dependencies unavailable")
class TestThermalSolver(unittest.TestCase):
    def _model(self, power_w=1.0, vias=None):
        component = ThermalComponentModel(
            ref_des="U1", power_w=power_w, width_mm=2.0, depth_mm=2.0,
            theta_jb_c_per_w=5.0, max_junction_c=125.0,
        )
        return ThermalBoardModel(
            bounds_mm=(0.0, 0.0, 10.0, 10.0),
            outline=box(0.0, 0.0, 10.0, 10.0),
            stackup={
                "copper": {
                    0: {"name": "F.Cu", "thickness_mm": 0.035},
                    31: {"name": "B.Cu", "thickness_mm": 0.035},
                },
                "layer_order": [0, 31],
                "substrate": [{
                    "between": [0, 31], "material": "FR4", "thickness_mm": 1.53,
                }],
            },
            copper_by_layer={
                0: box(0.0, 0.0, 10.0, 10.0),
                31: box(0.0, 0.0, 10.0, 10.0),
            },
            vias=list(vias or []),
            placements={"U1": ThermalPlacement("U1", 5.0, 5.0, 2.0, 2.0, "TOP")},
            components=[component],
        )

    def _settings(self, power=True, velocity=0.0):
        return ThermalAnalysisSettings(
            ambient_c=25.0,
            grid_size_mm=2.0,
            airflow=AirflowSettings(
                mode="FORCED" if velocity else "NATURAL",
                velocity_m_s=velocity,
                expose_top=True,
                expose_bottom=True,
                expose_edges=True,
            ),
        )

    def test_zero_heat_settles_at_ambient(self):
        mesh = ThermalMesher().generate_mesh(self._model(power_w=0.0), self._settings())
        result = ThermalSolver().solve(mesh, ambient_c=25.0)

        self.assertAlmostEqual(result.hotspot.temperature_c, 25.0, places=6)
        self.assertAlmostEqual(result.total_input_power_w, 0.0, places=9)

    def test_energy_balance_and_junction_estimate(self):
        mesh = ThermalMesher().generate_mesh(self._model(), self._settings())
        result = ThermalSolver().solve(mesh, ambient_c=25.0)

        self.assertGreater(result.hotspot.temperature_c, 25.0)
        self.assertLess(result.energy_balance_error_pct, 1e-5)
        self.assertEqual(result.component_results[0].ref_des, "U1")
        self.assertAlmostEqual(
            result.component_results[0].junction_temperature_c,
            result.component_results[0].board_temperature_c + 5.0,
        )

    def test_forced_air_reduces_hotspot(self):
        model = self._model()
        natural_mesh = ThermalMesher().generate_mesh(model, self._settings())
        forced_mesh = ThermalMesher().generate_mesh(model, self._settings(velocity=3.0))

        natural = ThermalSolver().solve(natural_mesh, ambient_c=25.0)
        forced = ThermalSolver().solve(forced_mesh, ambient_c=25.0)

        self.assertLess(forced.hotspot.temperature_c, natural.hotspot.temperature_c)

    def test_thermal_via_adds_vertical_branch(self):
        mesh = ThermalMesher().generate_mesh(
            self._model(vias=[ThermalVia(5.0, 5.0, 0.6)]), self._settings()
        )

        self.assertTrue(any(branch.kind == "via" for branch in mesh.branches))

    def test_forced_air_direction_sets_leading_edge_transfer(self):
        settings = self._settings(velocity=2.0)
        settings.airflow.direction_deg = 0.0
        mesh = ThermalMesher().generate_mesh(self._model(), settings)
        top = [boundary for boundary in mesh.boundaries if boundary.kind == "top"]
        upstream = min(top, key=lambda item: mesh.node_coords[item.node_id][0])
        downstream = max(top, key=lambda item: mesh.node_coords[item.node_id][0])

        self.assertGreater(upstream.conductance_w_k, downstream.conductance_w_k)


if __name__ == "__main__":
    unittest.main()
