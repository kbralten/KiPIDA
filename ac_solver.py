"""Sparse frequency-domain solver for Ki-PIDA power distribution networks."""

import warnings

try:
    import numpy as np
    import scipy.sparse
    import scipy.sparse.linalg
except ImportError:
    np = None
    scipy = None

try:
    from .models import ACAnalysisSettings, ImpedanceSweepResult
except (ImportError, ValueError):
    from models import ACAnalysisSettings, ImpedanceSweepResult


class ACSolver:
    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback
        if np is None or scipy is None:
            raise ImportError("NumPy and SciPy are required for AC analysis.")

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[AC SOLVER] {message}")

    @staticmethod
    def capacitor_impedance(capacitor, frequency_hz):
        omega = 2.0 * np.pi * frequency_hz
        capacitance = max(float(capacitor.capacitance_f), 1e-18)
        return complex(capacitor.esr_ohm, omega * capacitor.esl_h - 1.0 / (omega * capacitance))

    @staticmethod
    def _stamp_branch(matrix, node_a, node_b, admittance):
        matrix[node_a, node_a] += admittance
        matrix[node_b, node_b] += admittance
        matrix[node_a, node_b] -= admittance
        matrix[node_b, node_a] -= admittance

    def _stamp_group_branch(self, matrix, connection, admittance):
        if not connection.rail_nodes or not connection.ground_nodes:
            return
        pairs = len(connection.rail_nodes) * len(connection.ground_nodes)
        distributed = admittance / pairs
        for rail_node in connection.rail_nodes:
            for ground_node in connection.ground_nodes:
                self._stamp_branch(matrix, rail_node, ground_node, distributed)

    @staticmethod
    def _topology_anchors(network, capacitors):
        """Find floating islands and verify that the measurement reaches the source."""
        adjacency = [set() for _ in range(network.node_count)]

        def connect(node_a, node_b):
            if not (0 <= node_a < network.node_count and 0 <= node_b < network.node_count):
                raise ValueError("The AC network contains an invalid node reference.")
            adjacency[node_a].add(node_b)
            adjacency[node_b].add(node_a)

        for branch in network.branches:
            connect(branch.node_a, branch.node_b)
        for rail_node in network.source.rail_nodes:
            for ground_node in network.source.ground_nodes:
                connect(rail_node, ground_node)
        for capacitor in capacitors:
            if not capacitor.enabled or capacitor.capacitance_f <= 0:
                continue
            connection = network.capacitor_nodes.get(capacitor.ref_des)
            if connection is None:
                continue
            for rail_node in connection.rail_nodes:
                for ground_node in connection.ground_nodes:
                    connect(rail_node, ground_node)

        components = []
        component_by_node = {}
        remaining = set(range(network.node_count))
        while remaining:
            start = next(iter(remaining))
            stack = [start]
            component = []
            remaining.remove(start)
            while stack:
                node = stack.pop()
                component_by_node[node] = len(components)
                component.append(node)
                for neighbor in adjacency[node]:
                    if neighbor in remaining:
                        remaining.remove(neighbor)
                        stack.append(neighbor)
            components.append(component)

        source_anchor = network.source.ground_nodes[0]
        source_component = component_by_node[source_anchor]
        measurement_nodes = network.measurement.rail_nodes + network.measurement.ground_nodes
        if any(node < 0 or node >= network.node_count for node in measurement_nodes):
            raise ValueError("The measurement port contains an invalid node reference.")
        if any(component_by_node[node] != source_component for node in measurement_nodes):
            raise ValueError("The measurement port is not electrically connected to the AC source.")

        anchors = [source_anchor]
        anchors.extend(component[0] for index, component in enumerate(components)
                       if index != source_component)
        return anchors

    def solve_sweep(self, network, settings: ACAnalysisSettings, capacitors=None, progress_callback=None):
        if settings.frequency_start_hz <= 0:
            raise ValueError("Start frequency must be greater than zero.")
        if settings.frequency_stop_hz <= settings.frequency_start_hz:
            raise ValueError("Stop frequency must be greater than start frequency.")
        if settings.frequency_points < 2:
            raise ValueError("At least two frequency points are required.")
        if not network.measurement.rail_nodes or not network.measurement.ground_nodes:
            raise ValueError("Measurement port is not connected to the AC network.")

        capacitors = list(settings.capacitors if capacitors is None else capacitors)
        if not network.source.rail_nodes or not network.source.ground_nodes:
            raise ValueError("AC source is not connected to the network.")
        anchors = self._topology_anchors(network, capacitors)
        frequencies = np.logspace(
            np.log10(settings.frequency_start_hz),
            np.log10(settings.frequency_stop_hz),
            int(settings.frequency_points),
        )
        impedances = []

        for index, frequency in enumerate(frequencies):
            omega = 2.0 * np.pi * frequency
            matrix = scipy.sparse.lil_matrix((network.node_count, network.node_count), dtype=np.complex128)

            for branch in network.branches:
                resistance = max(float(branch.resistance_ohm), 1e-15)
                impedance = complex(resistance, omega * max(0.0, float(branch.inductance_h)))
                self._stamp_branch(matrix, branch.node_a, branch.node_b, 1.0 / impedance)

            source_impedance = complex(
                max(float(settings.source.resistance_ohm), 1e-9),
                omega * max(0.0, float(settings.source.inductance_h)),
            )
            self._stamp_group_branch(matrix, network.source, 1.0 / source_impedance)

            for capacitor in capacitors:
                if not capacitor.enabled or capacitor.capacitance_f <= 0:
                    continue
                connection = network.capacitor_nodes.get(capacitor.ref_des)
                if connection is None:
                    continue
                impedance = self.capacitor_impedance(capacitor, frequency)
                self._stamp_group_branch(matrix, connection, 1.0 / impedance)

            current = np.zeros(network.node_count, dtype=np.complex128)
            for node_id in network.measurement.rail_nodes:
                current[node_id] += 1.0 / len(network.measurement.rail_nodes)
            for node_id in network.measurement.ground_nodes:
                current[node_id] -= 1.0 / len(network.measurement.ground_nodes)

            # Fix one return node to zero volts to remove the arbitrary common-mode potential.
            for anchor in anchors:
                matrix.rows[anchor] = [anchor]
                matrix.data[anchor] = [1.0 + 0.0j]
                current[anchor] = 0.0

            with warnings.catch_warnings():
                warnings.simplefilter("error", scipy.sparse.linalg.MatrixRankWarning)
                try:
                    voltage = scipy.sparse.linalg.spsolve(matrix.tocsr(), current)
                except Exception as exc:
                    raise ValueError(f"AC solve failed at {frequency:g} Hz: {exc}") from exc

            rail_voltage = np.mean([voltage[node] for node in network.measurement.rail_nodes])
            ground_voltage = np.mean([voltage[node] for node in network.measurement.ground_nodes])
            impedances.append(complex(rail_voltage - ground_voltage))

            if progress_callback:
                progress_callback(index + 1, len(frequencies), float(frequency))

        magnitudes = np.abs(np.asarray(impedances))
        worst_index = int(np.argmax(magnitudes))
        target = max(0.0, float(settings.target_impedance_ohm))
        meets_target = bool(target > 0 and np.all(magnitudes <= target))
        if self.debug:
            self._log(
                f"Solved {len(frequencies)} points; worst |Z|={magnitudes[worst_index]:.6g} ohm "
                f"at {frequencies[worst_index]:.6g} Hz."
            )

        return ImpedanceSweepResult(
            frequencies_hz=[float(value) for value in frequencies],
            impedance_ohm=impedances,
            target_impedance_ohm=target,
            worst_frequency_hz=float(frequencies[worst_index]),
            worst_impedance_ohm=float(magnitudes[worst_index]),
            meets_target=meets_target,
        )
