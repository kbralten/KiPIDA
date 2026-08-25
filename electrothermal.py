"""Iterative coupling between Ki-PIDA DC meshes and the 3D thermal model."""

try:
    from .models import ElectroThermalResult
    from .solver import Solver
    from .thermal_solver import ThermalSolver
except (ImportError, ValueError):
    from models import ElectroThermalResult
    from solver import Solver
    from thermal_solver import ThermalSolver


class ElectroThermalSolver:
    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[ELECTRO-THERMAL] {message}")

    @staticmethod
    def _branch_temperature(thermal_mesh, temperatures, electrical_mesh, branch, ambient_c):
        coord_a = electrical_mesh.node_coords.get(branch.node_a)
        coord_b = electrical_mesh.node_coords.get(branch.node_b)
        if coord_a is None or coord_b is None:
            return float(ambient_c)
        x = (coord_a[0] + coord_b[0]) / 2.0
        y = (coord_a[1] + coord_b[1]) / 2.0
        layer = coord_a[2]
        node = thermal_mesh.nearest_node(x, y, layer)
        return float(temperatures.get(node, ambient_c)) if node is not None else float(ambient_c)

    def solve(self, thermal_mesh, settings, rail_contexts, progress_callback=None):
        if not rail_contexts:
            raise ValueError("Run a DC analysis before coupled electro-thermal analysis.")
        thermal_solver = ThermalSolver(debug=self.debug, log_callback=self.log_callback)
        dc_solver = Solver(debug=self.debug, log_callback=self.log_callback)
        base_heat = dict(thermal_mesh.heat_sources_w)
        previous_temperatures = {
            node: float(settings.ambient_c) for node in thermal_mesh.nodes
        }
        dc_results = {}
        converged = False
        thermal_result = None
        iterations = max(1, int(settings.coupled_iterations))

        for iteration in range(iterations):
            thermal_mesh.heat_sources_w = dict(base_heat)
            dc_results = {}
            for rail_name, context in rail_contexts.items():
                electrical_mesh = context["mesh"]
                scales = []
                for branch in electrical_mesh.branches:
                    temperature = self._branch_temperature(
                        thermal_mesh, previous_temperatures, electrical_mesh, branch, settings.ambient_c
                    )
                    scale = 1.0 + settings.copper_temp_coefficient_per_c * (temperature - 20.0)
                    scales.append(max(0.2, min(3.0, scale)))
                detailed = dc_solver.solve_detailed(
                    electrical_mesh,
                    context.get("sources", []),
                    context.get("loads", []),
                    branch_resistance_scales=scales,
                )
                dc_results[rail_name] = detailed
                for branch, power in zip(electrical_mesh.branches, detailed.branch_losses_w):
                    if power <= 0:
                        continue
                    coord_a = electrical_mesh.node_coords.get(branch.node_a)
                    coord_b = electrical_mesh.node_coords.get(branch.node_b)
                    if coord_a is None or coord_b is None:
                        continue
                    x = (coord_a[0] + coord_b[0]) / 2.0
                    y = (coord_a[1] + coord_b[1]) / 2.0
                    node = thermal_mesh.nearest_node(x, y, coord_a[2])
                    if node is not None:
                        thermal_mesh.add_heat(node, power)

            thermal_result = thermal_solver.solve(thermal_mesh, settings.ambient_c)
            max_delta = max(
                abs(thermal_result.temperatures_c[node] - previous_temperatures.get(node, settings.ambient_c))
                for node in thermal_mesh.nodes
            )
            relaxation = max(0.05, min(1.0, float(settings.relaxation)))
            previous_temperatures = {
                node: previous_temperatures.get(node, settings.ambient_c) + relaxation * (
                    thermal_result.temperatures_c[node] - previous_temperatures.get(node, settings.ambient_c)
                ) for node in thermal_mesh.nodes
            }
            self._log(
                f"Iteration {iteration + 1}/{iterations}: max delta {max_delta:.4g} C, "
                f"copper loss {sum(result.total_loss_w for result in dc_results.values()):.4g} W."
            )
            if progress_callback:
                progress_callback(iteration + 1, iterations, f"delta={max_delta:.3g} C")
            if max_delta <= max(1e-6, float(settings.convergence_c)):
                converged = True
                break

        thermal_result.iterations = iteration + 1
        thermal_result.converged = converged
        return ElectroThermalResult(
            thermal=thermal_result,
            dc_results=dc_results,
            iterations=iteration + 1,
            converged=converged,
        )
