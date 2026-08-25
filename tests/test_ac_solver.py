import os
import sys
import unittest

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

import ac_solver
from ac_model import ACNetwork, ACNodeConnection
from ac_solver import ACSolver
from decoupling_optimizer import DecouplingOptimizer
from models import ACAnalysisSettings, ACSourceModel, CapacitorModel


HAS_NUMERICAL_DEPS = ac_solver.np is not None and ac_solver.scipy is not None


@unittest.skipUnless(HAS_NUMERICAL_DEPS, "NumPy/SciPy are not installed in this test interpreter")
class TestACSolver(unittest.TestCase):
    @staticmethod
    def network(node_count=2, measurement=None, capacitor_nodes=None):
        connection = ACNodeConnection([0], [1])
        return ACNetwork(
            node_count=node_count,
            branches=[],
            source=connection,
            measurement=measurement or connection,
            capacitor_nodes=capacitor_nodes or {},
        )

    @staticmethod
    def settings(**overrides):
        values = dict(
            rail_name="3V3",
            frequency_start_hz=1e3,
            frequency_stop_hz=1e5,
            frequency_points=3,
            target_impedance_ohm=0.11,
            source=ACSourceModel(resistance_ohm=0.1, inductance_h=0.0),
        )
        values.update(overrides)
        return ACAnalysisSettings(**values)

    def test_source_impedance_is_recovered_by_one_ampere_injection(self):
        result = ACSolver().solve_sweep(self.network(), self.settings())

        for impedance in result.impedance_ohm:
            self.assertAlmostEqual(abs(impedance), 0.1, places=9)
        self.assertTrue(result.meets_target)

    def test_rlc_capacitor_resonance_approaches_esr(self):
        capacitor = CapacitorModel(
            ref_des="C1", capacitance_f=1e-6, esr_ohm=0.02, esl_h=1e-9,
        )
        network = self.network(capacitor_nodes={"C1": ACNodeConnection([0], [1])})
        resonance_hz = 1.0 / (2.0 * ac_solver.np.pi * (1e-6 * 1e-9) ** 0.5)
        settings = self.settings(
            frequency_start_hz=resonance_hz * 0.999,
            frequency_stop_hz=resonance_hz * 1.001,
            frequency_points=3,
            source=ACSourceModel(resistance_ohm=1000.0, inductance_h=0.0),
            capacitors=[capacitor],
        )

        result = ACSolver().solve_sweep(network, settings)

        self.assertAlmostEqual(abs(result.impedance_ohm[1]), 0.02, delta=5e-5)

    def test_disconnected_measurement_port_is_rejected(self):
        network = self.network(
            node_count=4,
            measurement=ACNodeConnection([2], [3]),
        )
        with self.assertRaisesRegex(ValueError, "not electrically connected"):
            ACSolver().solve_sweep(network, self.settings())

    def test_optimizer_populates_existing_candidate(self):
        candidate = CapacitorModel(
            ref_des="C_DNP", capacitance_f=100e-9, esr_ohm=0.01,
            esl_h=0.0, enabled=False, candidate=True,
        )
        network = self.network(capacitor_nodes={"C_DNP": ACNodeConnection([0], [1])})
        settings = self.settings(
            frequency_start_hz=100.0,
            frequency_stop_hz=200.0,
            frequency_points=5,
            target_impedance_ohm=0.1,
            source=ACSourceModel(resistance_ohm=0.5, inductance_h=0.0),
            capacitors=[candidate],
            optimizer_values_f=[0.1],
            optimizer_max_additions=1,
        )

        result = DecouplingOptimizer(ACSolver()).optimize(network, settings)

        self.assertEqual([item.ref_des for item in result.recommendations], ["C_DNP"])
        self.assertLess(result.optimized.worst_impedance_ohm, result.baseline.worst_impedance_ohm)
        self.assertTrue(result.reached_target)


if __name__ == "__main__":
    unittest.main()
