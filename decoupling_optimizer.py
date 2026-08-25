"""Deterministic target-impedance decoupling optimization."""

from dataclasses import replace

try:
    import numpy as np
except ImportError:
    np = None

try:
    from .models import (
        DecouplingOptimizationResult,
        OptimizationRecommendation,
    )
except (ImportError, ValueError):
    from models import DecouplingOptimizationResult, OptimizationRecommendation


class DecouplingOptimizer:
    """Greedily populate existing/DNP capacitor footprints.

    The optimizer never invents placement locations and never writes the PCB.
    Each recommendation is tied to a capacitor footprint already present in the
    AC network.
    """

    def __init__(self, solver, debug=False, log_callback=None):
        self.solver = solver
        self.debug = debug
        self.log_callback = log_callback

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[DECAP OPT] {message}")

    @staticmethod
    def _score(result, target):
        magnitudes = np.abs(np.asarray(result.impedance_ohm))
        if target > 0:
            ratios = magnitudes / target
            return float(np.max(ratios) + 0.05 * np.mean(np.maximum(ratios - 1.0, 0.0)))
        return float(np.max(magnitudes))

    @staticmethod
    def _clone_capacitors(capacitors):
        return [replace(capacitor) for capacitor in capacitors]

    def optimize(self, network, settings, progress_callback=None):
        if np is None:
            raise ImportError("NumPy is required for decoupling optimization.")

        baseline_caps = self._clone_capacitors(settings.capacitors)
        baseline = self.solver.solve_sweep(network, settings, baseline_caps)
        target = max(0.0, float(settings.target_impedance_ohm))
        if baseline.meets_target:
            return DecouplingOptimizationResult(
                baseline=baseline,
                optimized=baseline,
                recommendations=[],
                reached_target=True,
            )

        working = self._clone_capacitors(baseline_caps)
        available_refs = sorted(
            capacitor.ref_des for capacitor in working
            if not capacitor.enabled and capacitor.ref_des in network.capacitor_nodes
        )
        values = sorted(set(float(value) for value in settings.optimizer_values_f if value > 0))
        max_additions = max(0, int(settings.optimizer_max_additions))
        recommendations = []

        coarse_settings = replace(settings, frequency_points=min(41, settings.frequency_points))
        current_result = self.solver.solve_sweep(network, coarse_settings, working)
        current_score = self._score(current_result, target)
        total_trials = max(1, max_additions * max(1, len(available_refs)) * max(1, len(values)))
        completed_trials = 0

        for _ in range(max_additions):
            best = None
            best_result = None
            best_score = current_score

            for ref_des in list(available_refs):
                for capacitance in values:
                    trial = self._clone_capacitors(working)
                    selected = next(cap for cap in trial if cap.ref_des == ref_des)
                    selected.enabled = True
                    selected.candidate = False
                    selected.capacitance_f = capacitance
                    result = self.solver.solve_sweep(network, coarse_settings, trial)
                    score = self._score(result, target)
                    completed_trials += 1
                    if progress_callback:
                        progress_callback(completed_trials, total_trials, ref_des)

                    key = (score, ref_des, capacitance)
                    if best is None or key < best[0]:
                        best = (key, ref_des, capacitance, trial)
                        best_result = result
                        best_score = score

            if best is None or best_score >= current_score * (1.0 - 1e-9):
                self._log("No remaining candidate improves the target-impedance score.")
                break

            _, ref_des, capacitance, working = best
            available_refs.remove(ref_des)
            current_result = best_result
            current_score = best_score
            recommendations.append(OptimizationRecommendation(
                ref_des=ref_des,
                capacitance_f=capacitance,
                action="populate",
            ))
            self._log(f"Selected {ref_des} = {capacitance:g} F (score {current_score:.4g}).")

            if target > 0 and current_result.meets_target:
                break
            if not available_refs:
                break

        optimized = self.solver.solve_sweep(network, settings, working)
        return DecouplingOptimizationResult(
            baseline=baseline,
            optimized=optimized,
            recommendations=recommendations,
            reached_target=optimized.meets_target,
        )
