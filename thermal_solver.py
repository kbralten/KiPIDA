"""Sparse steady-state thermal solver with convective boundaries."""

try:
    import numpy as np
    import scipy.sparse
    import scipy.sparse.linalg
except ImportError:
    np = None
    scipy = None

try:
    import pypardiso
except ImportError:
    pypardiso = None

try:
    from .models import ComponentThermalResult, ThermalHotspot, ThermalResult
except (ImportError, ValueError):
    from models import ComponentThermalResult, ThermalHotspot, ThermalResult


class ThermalSolver:
    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback
        if np is None or scipy is None:
            raise ImportError("NumPy and SciPy are required for 3D thermal analysis.")

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[THERMAL SOLVER] {message}")

    def solve(self, mesh, ambient_c=25.0, progress_callback=None):
        if not mesh.nodes:
            raise ValueError("The thermal mesh is empty.")
        if not mesh.boundaries:
            raise ValueError("At least one exposed convective surface is required.")

        node_to_index = {node: index for index, node in enumerate(mesh.nodes)}
        count = len(mesh.nodes)
        matrix = scipy.sparse.lil_matrix((count, count), dtype=float)
        rhs = np.zeros(count, dtype=float)

        for branch in mesh.branches:
            a = node_to_index.get(branch.node_a)
            b = node_to_index.get(branch.node_b)
            conductance = float(branch.conductance_w_k)
            if a is None or b is None or conductance <= 0:
                continue
            matrix[a, a] += conductance
            matrix[b, b] += conductance
            matrix[a, b] -= conductance
            matrix[b, a] -= conductance

        for boundary in mesh.boundaries:
            index = node_to_index.get(boundary.node_id)
            conductance = float(boundary.conductance_w_k)
            if index is None or conductance <= 0:
                continue
            matrix[index, index] += conductance
            rhs[index] += conductance * float(ambient_c)

        for node, power in mesh.heat_sources_w.items():
            index = node_to_index.get(node)
            if index is not None:
                rhs[index] += float(power)

        if progress_callback:
            progress_callback(1, 3, "matrix")
        csr = matrix.tocsr()
        if pypardiso is not None:
            temperatures = pypardiso.spsolve(csr, rhs)
        else:
            temperatures = scipy.sparse.linalg.spsolve(csr, rhs)
        if progress_callback:
            progress_callback(2, 3, "solve")

        if np.any(~np.isfinite(temperatures)):
            raise ValueError("Thermal solution contains non-finite temperatures.")
        result_temperatures = {
            node: float(temperatures[index]) for node, index in node_to_index.items()
        }
        hotspot_node = max(result_temperatures, key=result_temperatures.get)
        hotspot_x, hotspot_y, hotspot_z = mesh.node_coords[hotspot_node]
        hotspot = ThermalHotspot(
            node_id=hotspot_node,
            x_mm=float(hotspot_x),
            y_mm=float(hotspot_y),
            z_mm=float(hotspot_z),
            temperature_c=result_temperatures[hotspot_node],
        )

        component_results = []
        for ref_des, nodes in mesh.component_nodes.items():
            component = mesh.component_models[ref_des]
            valid = [result_temperatures[node] for node in nodes if node in result_temperatures]
            if not valid:
                continue
            board_temperature = float(sum(valid) / len(valid))
            junction = board_temperature + component.power_w * max(0.0, component.theta_jb_c_per_w)
            component_results.append(ComponentThermalResult(
                ref_des=ref_des,
                board_temperature_c=board_temperature,
                junction_temperature_c=junction,
                power_w=component.power_w,
                max_junction_c=component.max_junction_c,
                margin_c=component.max_junction_c - junction,
                model_source=component.model_source,
            ))
        component_results.sort(key=lambda item: item.junction_temperature_c, reverse=True)

        boundary_power = sum(
            boundary.conductance_w_k *
            (result_temperatures[boundary.node_id] - float(ambient_c))
            for boundary in mesh.boundaries
        )
        input_power = float(sum(mesh.heat_sources_w.values()))
        denominator = max(abs(input_power), 1e-12)
        balance_error = abs(boundary_power - input_power) / denominator * 100.0
        if progress_callback:
            progress_callback(3, 3, "results")
        self._log(
            f"Solved {count:,} nodes; hotspot {hotspot.temperature_c:.2f} C, "
            f"energy error {balance_error:.3g}%."
        )
        return ThermalResult(
            temperatures_c=result_temperatures,
            hotspot=hotspot,
            component_results=component_results,
            total_input_power_w=input_power,
            total_boundary_power_w=float(boundary_power),
            energy_balance_error_pct=float(balance_error),
            convection_coefficient_w_m2k=float(mesh.convection_coefficient_w_m2k),
        )
