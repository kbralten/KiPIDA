"""Build frequency-domain PDN models from the existing Ki-PIDA mesh pipeline."""

from dataclasses import dataclass, field
import re
from typing import Dict, List, Optional, Tuple

try:
    from .extractor import GeometryExtractor
    from .mesh import Mesher
    from .models import ACAnalysisSettings, CapacitorModel, MeshBranch, PowerRail
except (ImportError, ValueError):
    from extractor import GeometryExtractor
    from mesh import Mesher
    from models import ACAnalysisSettings, CapacitorModel, MeshBranch, PowerRail


@dataclass
class ACNodeConnection:
    rail_nodes: List[int] = field(default_factory=list)
    ground_nodes: List[int] = field(default_factory=list)


@dataclass
class ACNetwork:
    node_count: int
    branches: List[MeshBranch]
    source: ACNodeConnection
    measurement: ACNodeConnection
    capacitor_nodes: Dict[str, ACNodeConnection]
    node_coords: Dict[int, Tuple[float, float, int]] = field(default_factory=dict)


def parse_capacitance(value) -> Optional[float]:
    """Parse common KiCad capacitor values such as 100n, 4u7 and 10uF."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value) if value > 0 else None

    text = str(value).strip().replace("µ", "u").replace("μ", "u")
    text = text.replace(" ", "")
    match = re.search(r"(?i)(\d+(?:\.\d+)?)([pnum]?)(\d*)\s*f?", text)
    if not match:
        return None

    major, prefix, fractional = match.groups()
    if fractional and "." not in major:
        major = f"{major}.{fractional}"

    scale = {"": 1.0, "p": 1e-12, "n": 1e-9, "u": 1e-6, "m": 1e-3}
    try:
        capacitance = float(major) * scale[prefix.lower()]
    except (ValueError, KeyError):
        return None
    return capacitance if capacitance > 0 else None


def format_capacitance(value_f: float) -> str:
    for scale, suffix in ((1e-3, "mF"), (1e-6, "uF"), (1e-9, "nF"), (1e-12, "pF")):
        if value_f >= scale:
            return f"{value_f / scale:g}{suffix}"
    return f"{value_f:g}F"


class ACModelBuilder:
    """Adapter between KiCad board objects and the pure numerical AC solver."""

    def __init__(self, board, debug=False, log_callback=None):
        self.board = board
        self.debug = debug
        self.log_callback = log_callback

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[AC MODEL] {message}")

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

    def _get_board_items(self, attr_name):
        value = self._get_val(self.board, attr_name, [])
        return value if value is not None else []

    def _get_footprint_reference(self, footprint):
        reference = self._get_val(footprint, "reference", self._get_val(footprint, "ref_des", ""))
        if reference:
            return str(reference)
        field = self._get_val(footprint, "reference_field")
        text = self._get_val(field, "text")
        return str(self._get_val(text, "value", "") or "")

    def _get_footprint_value(self, footprint):
        value = self._get_val(footprint, "value", "")
        if value:
            return str(value)
        field = self._get_val(footprint, "value_field")
        text = self._get_val(field, "text")
        return str(self._get_val(text, "value", "") or "")

    def _get_footprint_name(self, footprint):
        for attr in ("footprint_name", "library_link", "lib_id", "name"):
            value = self._get_val(footprint, attr, "")
            if value:
                return str(value)
        definition = self._get_val(footprint, "definition")
        return str(self._get_val(definition, "name", "") or "")

    def _get_pads(self, footprint):
        pads = self._get_val(footprint, "pads")
        if pads is None:
            definition = self._get_val(footprint, "definition")
            pads = self._get_val(definition, "pads", [])
        return pads or []

    def _pad_number(self, pad):
        return str(self._get_val(pad, "number", self._get_val(pad, "name", "")) or "")

    def _pad_net_name(self, pad):
        return str(self._get_val(self._get_val(pad, "net"), "name", "") or "")

    def _is_dnp(self, footprint):
        for attr in ("dnp", "do_not_populate", "exclude_from_bom", "not_in_bom"):
            if bool(self._get_val(footprint, attr, False)):
                return True
        attrs = str(self._get_val(footprint, "attributes", "")).lower()
        return "dnp" in attrs or "do_not_populate" in attrs

    def discover_ground_nets(self):
        names = set()
        for footprint in self._get_board_items("footprints"):
            for pad in self._get_pads(footprint):
                name = self._pad_net_name(pad)
                if re.search(r"(?i)(^|[/_+\-])(gnd|vss|ground)(\d*|$)", name) or name.upper() in {
                    "GND", "AGND", "DGND", "PGND", "VSS"
                }:
                    names.add(name)
        return sorted(names)

    def _default_esl(self, footprint_name):
        name = footprint_name.lower()
        for package, esl in (("0201", 0.25e-9), ("0402", 0.4e-9), ("0603", 0.6e-9),
                             ("0805", 0.8e-9), ("1206", 1.0e-9), ("1210", 1.2e-9)):
            if package in name:
                return esl
        return 0.8e-9

    def discover_capacitors(self, rail_name, ground_name):
        capacitors = []
        for footprint in self._get_board_items("footprints"):
            reference = self._get_footprint_reference(footprint)
            value_text = self._get_footprint_value(footprint)
            capacitance = parse_capacitance(value_text)
            if not reference.upper().startswith("C") and capacitance is None:
                continue

            rail_pads = []
            ground_pads = []
            for pad in self._get_pads(footprint):
                net_name = self._pad_net_name(pad)
                if net_name == rail_name:
                    rail_pads.append(self._pad_number(pad))
                elif net_name == ground_name:
                    ground_pads.append(self._pad_number(pad))

            if not rail_pads or not ground_pads:
                continue

            candidate = self._is_dnp(footprint) or capacitance is None
            capacitors.append(CapacitorModel(
                ref_des=reference,
                rail_pad_names=rail_pads,
                ground_pad_names=ground_pads,
                capacitance_f=capacitance or 100e-9,
                esr_ohm=0.01,
                esl_h=self._default_esl(self._get_footprint_name(footprint)),
                enabled=not candidate,
                candidate=candidate,
                model_source="value-and-package-estimate",
            ))

        capacitors.sort(key=lambda cap: cap.ref_des)
        return capacitors

    def _find_footprint(self, ref_des):
        for footprint in self._get_board_items("footprints"):
            if self._get_footprint_reference(footprint) == ref_des:
                return footprint
        return None

    def pad_names_for_net(self, ref_des, net_name):
        footprint = self._find_footprint(ref_des)
        if not footprint:
            return []
        return [self._pad_number(pad) for pad in self._get_pads(footprint)
                if self._pad_net_name(pad) == net_name]

    def find_pad_nodes(self, mesh, ref_des, pad_names, net_name):
        footprint = self._find_footprint(ref_des)
        if not footprint:
            return []

        wanted = set(str(name) for name in (pad_names or []))
        found = []
        for pad in self._get_pads(footprint):
            if self._pad_net_name(pad) != net_name:
                continue
            if wanted and self._pad_number(pad) not in wanted:
                continue

            position = self._get_val(pad, "position")
            if position is None:
                continue
            x = float(self._get_val(position, "x", 0.0) or 0.0)
            y = float(self._get_val(position, "y", 0.0) or 0.0)
            if abs(x) > 10000 or abs(y) > 10000:
                x /= 1e6
                y /= 1e6

            ix = int(round((x - mesh.grid_origin[0]) / mesh.grid_step))
            iy = int(round((y - mesh.grid_origin[1]) / mesh.grid_step))
            layers = self._get_val(pad, "layers") or []
            if not layers:
                padstack = self._get_val(pad, "padstack")
                layers = self._get_val(padstack, "layers", []) or []
            if not layers:
                layers = sorted({coords[2] for coords in mesh.node_coords.values()})

            pad_nodes = []
            for radius in range(0, 3):
                for dx in range(-radius, radius + 1):
                    for dy in range(-radius, radius + 1):
                        if radius and abs(dx) < radius and abs(dy) < radius:
                            continue
                        for layer in layers:
                            node_id = mesh.node_map.get((ix + dx, iy + dy, layer))
                            if node_id is not None:
                                pad_nodes.append(int(node_id))
                if pad_nodes:
                    break
            found.extend(pad_nodes)

        return sorted(set(found))

    def _translate_mesh(self, mesh, offset):
        branches = [MeshBranch(
            node_a=branch.node_a + offset,
            node_b=branch.node_b + offset,
            resistance_ohm=branch.resistance_ohm,
            inductance_h=branch.inductance_h,
            kind=branch.kind,
        ) for branch in mesh.branches]
        coords = {node_id + offset: coord for node_id, coord in mesh.node_coords.items()}
        return branches, coords

    def build(self, rail: PowerRail, settings: ACAnalysisSettings, grid_size_mm=0.5):
        extractor = GeometryExtractor(self.board, debug=self.debug, log_callback=self.log_callback)
        stackup = extractor.get_board_stackup()
        rail_geometry = extractor.get_net_geometry(settings.rail_name)
        ground_geometry = extractor.get_net_geometry(settings.ground_net_name)
        if not rail_geometry:
            raise ValueError(f"No copper geometry found for rail '{settings.rail_name}'.")
        if not ground_geometry:
            raise ValueError(f"No copper geometry found for return net '{settings.ground_net_name}'.")

        mesher = Mesher(self.board, debug=self.debug, log_callback=self.log_callback)
        rail_mesh = mesher.generate_mesh(settings.rail_name, rail_geometry, stackup, grid_size_mm)
        ground_mesh = mesher.generate_mesh(settings.ground_net_name, ground_geometry, stackup, grid_size_mm)
        if not rail_mesh.nodes or not ground_mesh.nodes:
            raise ValueError("The AC rail-to-ground mesh is empty.")

        offset = len(rail_mesh.nodes)
        rail_branches, rail_coords = self._translate_mesh(rail_mesh, 0)
        ground_branches, ground_coords = self._translate_mesh(ground_mesh, offset)

        source_ref = settings.source.ref_des
        source_rail_pads = settings.source.rail_pad_names
        if not source_ref and rail.sources:
            source_ref = rail.sources[0].component_ref.ref_des
            source_rail_pads = rail.sources[0].pad_names
        source_ground_pads = settings.source.ground_pad_names or self.pad_names_for_net(
            source_ref, settings.ground_net_name
        )

        port_ref = settings.measurement_port.ref_des
        port_rail_pads = settings.measurement_port.rail_pad_names
        if not port_ref and rail.loads:
            port_ref = rail.loads[0].component_ref.ref_des
            port_rail_pads = rail.loads[0].pad_names
        port_ground_pads = settings.measurement_port.ground_pad_names or self.pad_names_for_net(
            port_ref, settings.ground_net_name
        )

        source = ACNodeConnection(
            rail_nodes=self.find_pad_nodes(rail_mesh, source_ref, source_rail_pads, settings.rail_name),
            ground_nodes=[node + offset for node in self.find_pad_nodes(
                ground_mesh, source_ref, source_ground_pads, settings.ground_net_name
            )],
        )
        measurement = ACNodeConnection(
            rail_nodes=self.find_pad_nodes(rail_mesh, port_ref, port_rail_pads, settings.rail_name),
            ground_nodes=[node + offset for node in self.find_pad_nodes(
                ground_mesh, port_ref, port_ground_pads, settings.ground_net_name
            )],
        )
        if not source.rail_nodes or not source.ground_nodes:
            raise ValueError("The AC source must map to pads on both the rail and the return net.")
        if not measurement.rail_nodes or not measurement.ground_nodes:
            raise ValueError("The measurement port must map to pads on both the rail and the return net.")

        capacitor_nodes = {}
        for capacitor in settings.capacitors:
            rail_nodes = self.find_pad_nodes(
                rail_mesh, capacitor.ref_des, capacitor.rail_pad_names, settings.rail_name
            )
            ground_nodes = [node + offset for node in self.find_pad_nodes(
                ground_mesh, capacitor.ref_des, capacitor.ground_pad_names, settings.ground_net_name
            )]
            if rail_nodes and ground_nodes:
                capacitor_nodes[capacitor.ref_des] = ACNodeConnection(rail_nodes, ground_nodes)
            else:
                self._log(f"Skipping {capacitor.ref_des}: pads could not be mapped to both meshes.")

        return ACNetwork(
            node_count=len(rail_mesh.nodes) + len(ground_mesh.nodes),
            branches=rail_branches + ground_branches,
            source=source,
            measurement=measurement,
            capacitor_nodes=capacitor_nodes,
            node_coords={**rail_coords, **ground_coords},
        )
