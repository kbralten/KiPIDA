"""Build board-level thermal models from KiCad geometry and power-tree data."""

from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Tuple

try:
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
except ImportError:
    Polygon = box = unary_union = None

try:
    from .extractor import GeometryExtractor
    from .models import PowerRail, ThermalAnalysisSettings, ThermalComponentModel
except (ImportError, ValueError):
    from extractor import GeometryExtractor
    from models import PowerRail, ThermalAnalysisSettings, ThermalComponentModel


@dataclass
class ThermalVia:
    x_mm: float
    y_mm: float
    diameter_mm: float = 0.6


@dataclass
class ThermalPlacement:
    ref_des: str
    x_mm: float
    y_mm: float
    width_mm: float
    depth_mm: float
    side: str = "TOP"


@dataclass
class CopperLossPoint:
    x_mm: float
    y_mm: float
    layer_id: int
    power_w: float


@dataclass
class ThermalBoardModel:
    bounds_mm: Tuple[float, float, float, float]
    outline: object
    stackup: dict
    copper_by_layer: Dict[int, object] = field(default_factory=dict)
    vias: List[ThermalVia] = field(default_factory=list)
    placements: Dict[str, ThermalPlacement] = field(default_factory=dict)
    components: List[ThermalComponentModel] = field(default_factory=list)
    copper_losses: List[CopperLossPoint] = field(default_factory=list)


class PowerLossEstimator:
    """Estimate component heat from the existing power-tree semantics."""

    @staticmethod
    def estimate(rails: List[PowerRail]) -> List[ThermalComponentModel]:
        rail_by_name = {rail.net_name: rail for rail in rails}
        component_power = {}
        model_source = {}

        for rail in rails:
            voltage = max(0.0, float(rail.nominal_voltage))
            for load in rail.loads:
                ref_des = load.component_ref.ref_des
                component_power[ref_des] = component_power.get(ref_des, 0.0) + voltage * max(
                    0.0, float(load.total_current)
                )
                model_source[ref_des] = "power-tree-load"

        memo = {}
        visiting = set()

        def rail_current(rail_name):
            if rail_name in memo:
                return memo[rail_name]
            if rail_name in visiting:
                return 0.0
            visiting.add(rail_name)
            rail = rail_by_name.get(rail_name)
            if rail is None:
                visiting.remove(rail_name)
                return 0.0
            total = sum(max(0.0, float(load.total_current)) for load in rail.loads)
            for regulator in rail.child_regulators:
                output_rail = rail_by_name.get(regulator.output_rail_name)
                output_current = rail_current(regulator.output_rail_name)
                output_voltage = max(0.0, float(output_rail.nominal_voltage)) if output_rail else 0.0
                input_voltage = max(0.0, float(rail.nominal_voltage))
                output_power = output_voltage * output_current
                if regulator.reg_type == "SWITCHING":
                    efficiency = max(1e-6, min(1.0, float(regulator.efficiency)))
                    loss = output_power * (1.0 / efficiency - 1.0)
                    input_current = output_power / efficiency / input_voltage if input_voltage > 0 else 0.0
                else:
                    loss = max(0.0, input_voltage - output_voltage) * output_current
                    input_current = output_current
                ref_des = regulator.output_ref_des or regulator.input_ref_des
                if ref_des:
                    component_power[ref_des] = component_power.get(ref_des, 0.0) + loss
                    model_source[ref_des] = "regulator-loss"
                total += input_current
            visiting.remove(rail_name)
            memo[rail_name] = total
            return total

        for rail in rails:
            rail_current(rail.net_name)

        return [ThermalComponentModel(
            ref_des=ref_des,
            power_w=power,
            model_source=model_source.get(ref_des, "estimated"),
        ) for ref_des, power in sorted(component_power.items()) if power > 0]


class ThermalModelBuilder:
    def __init__(self, board, debug=False, log_callback=None):
        self.board = board
        self.debug = debug
        self.log_callback = log_callback

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[THERMAL MODEL] {message}")

    def _get_val(self, obj, attr_name, default=None):
        if obj is None:
            return default
        if hasattr(obj, attr_name):
            value = getattr(obj, attr_name)
            if value is not None and not callable(value):
                return value
        for method_name in (f"get_{attr_name}", attr_name):
            if hasattr(obj, method_name):
                try:
                    value = getattr(obj, method_name)()
                    if value is not None:
                        return value
                except Exception:
                    pass
        return default

    def _items(self, name):
        return self._get_val(self.board, name, []) or []

    @staticmethod
    def _to_mm(value):
        value = float(value or 0.0)
        return value / 1e6 if abs(value) > 10000 else value

    def _position(self, item):
        position = self._get_val(item, "position", self._get_val(item, "start"))
        if position is None:
            return None
        return self._to_mm(self._get_val(position, "x", 0.0)), self._to_mm(
            self._get_val(position, "y", 0.0)
        )

    def _reference(self, footprint):
        value = self._get_val(footprint, "reference", self._get_val(footprint, "ref_des", ""))
        if value:
            return str(value)
        field = self._get_val(footprint, "reference_field")
        text = self._get_val(field, "text")
        return str(self._get_val(text, "value", "") or "")

    def _pads(self, footprint):
        pads = self._get_val(footprint, "pads")
        if pads is None:
            pads = self._get_val(self._get_val(footprint, "definition"), "pads", [])
        return pads or []

    def _footprint_size(self, footprint, default_width, default_depth):
        points = [self._position(pad) for pad in self._pads(footprint)]
        points = [point for point in points if point is not None]
        if len(points) >= 2:
            width = max(point[0] for point in points) - min(point[0] for point in points)
            depth = max(point[1] for point in points) - min(point[1] for point in points)
            return max(default_width, width + 0.5), max(default_depth, depth + 0.5)
        return default_width, default_depth

    def _discover_net_names(self):
        names = set()
        for collection in ("tracks", "vias", "zones", "footprints"):
            for item in self._items(collection):
                candidates = self._pads(item) if collection == "footprints" else [item]
                for candidate in candidates:
                    net = self._get_val(candidate, "net")
                    name = self._get_val(net, "name", "")
                    if name:
                        names.add(str(name))
        return sorted(names)

    def _extract_copper(self, extractor):
        by_layer = {}
        names = self._discover_net_names()
        for index, net_name in enumerate(names):
            geometry = extractor.get_net_geometry(net_name)
            for layer, polygon in geometry.items():
                if polygon is None or polygon.is_empty:
                    continue
                if layer in by_layer:
                    by_layer[layer] = unary_union([by_layer[layer], polygon])
                else:
                    by_layer[layer] = polygon
            if self.debug and (index + 1) % 20 == 0:
                self._log(f"Aggregated copper for {index + 1}/{len(names)} nets.")
        return by_layer

    def _outline_and_bounds(self, copper_by_layer, placements):
        polygons = [
            polygon for polygon in copper_by_layer.values()
            if polygon is not None and not polygon.is_empty
        ]
        extents = list(polygons)
        if box is not None:
            extents.extend(box(
                placement.x_mm - placement.width_mm / 2.0,
                placement.y_mm - placement.depth_mm / 2.0,
                placement.x_mm + placement.width_mm / 2.0,
                placement.y_mm + placement.depth_mm / 2.0,
            ) for placement in placements.values())
        if extents and unary_union is not None:
            bounds = unary_union(extents).bounds
        elif placements:
            xs = [placement.x_mm for placement in placements.values()]
            ys = [placement.y_mm for placement in placements.values()]
            bounds = min(xs), min(ys), max(xs), max(ys)
        else:
            bounds = (0.0, 0.0, 50.0, 50.0)
        min_x, min_y, max_x, max_y = bounds
        if max_x - min_x < 1.0:
            max_x = min_x + 10.0
        if max_y - min_y < 1.0:
            max_y = min_y + 10.0
        margin = 1.0
        bounds = min_x - margin, min_y - margin, max_x + margin, max_y + margin
        return (box(*bounds) if box is not None else None), bounds

    def build(self, settings: ThermalAnalysisSettings, rails=None, copper_losses=None):
        if Polygon is None or unary_union is None:
            raise ImportError("Shapely is required for 3D thermal model extraction.")
        extractor = GeometryExtractor(self.board, debug=self.debug, log_callback=self.log_callback)
        stackup = extractor.get_board_stackup()

        saved_components = {component.ref_des: replace(component) for component in settings.components}
        if not saved_components and rails:
            saved_components = {
                component.ref_des: component for component in PowerLossEstimator.estimate(list(rails))
            }

        placements = {}
        for footprint in self._items("footprints"):
            ref_des = self._reference(footprint)
            position = self._position(footprint)
            if not ref_des or position is None:
                continue
            component = saved_components.get(ref_des, ThermalComponentModel(ref_des=ref_des))
            width, depth = self._footprint_size(footprint, component.width_mm, component.depth_mm)
            layer = self._get_val(footprint, "layer", 0)
            placements[ref_des] = ThermalPlacement(
                ref_des=ref_des,
                x_mm=position[0],
                y_mm=position[1],
                width_mm=width,
                depth_mm=depth,
                side="BOTTOM" if str(layer) in {"31", "B_Cu", "B.Cu"} else "TOP",
            )
            if ref_des in saved_components:
                saved_components[ref_des].width_mm = width
                saved_components[ref_des].depth_mm = depth

        copper_by_layer = self._extract_copper(extractor)
        outline, bounds = self._outline_and_bounds(copper_by_layer, placements)

        vias = []
        for via in self._items("vias"):
            position = self._position(via)
            if position is None:
                continue
            width = self._to_mm(self._get_val(via, "width", 0.6e6))
            vias.append(ThermalVia(position[0], position[1], max(width, 0.1)))

        return ThermalBoardModel(
            bounds_mm=bounds,
            outline=outline,
            stackup=stackup,
            copper_by_layer=copper_by_layer,
            vias=vias,
            placements=placements,
            components=list(saved_components.values()),
            copper_losses=list(copper_losses or []),
        )
