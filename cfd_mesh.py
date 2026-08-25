"""Structured volumetric mesh for enclosure airflow and conjugate heat transfer."""

from dataclasses import dataclass, field
import math

try:
    import numpy as np
except ImportError:
    np = None


@dataclass
class CFDPatchCells:
    patch: object
    cells: list = field(default_factory=list)


@dataclass
class CFDMesh:
    shape: tuple
    spacing_m: tuple
    dimensions_m: tuple
    fluid_mask: object
    solid_mask: object
    solid_conductivity_w_mk: object
    heat_sources_w: object
    patch_cells: list = field(default_factory=list)
    obstacle_names: object = None

    @property
    def cell_count(self):
        return int(self.shape[0] * self.shape[1] * self.shape[2])

    @property
    def cell_volume_m3(self):
        return self.spacing_m[0] * self.spacing_m[1] * self.spacing_m[2]

    def cell_center_mm(self, i, j, k):
        dx, dy, dz = self.spacing_m
        return ((i + 0.5) * dx * 1000.0,
                (j + 0.5) * dy * 1000.0,
                (k + 0.5) * dz * 1000.0)


class CFDMeshGenerator:
    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback

    def _log(self, message):
        if self.log_callback:
            self.log_callback(f"[CFD MESH] {message}")

    @staticmethod
    def _patch_cells(patch, shape):
        nx, ny, nz = shape
        face = str(patch.face).upper()
        if face in {"XMIN", "XMAX"}:
            fixed = 0 if face == "XMIN" else nx - 1
            nu, nv = ny, nz
            make = lambda u, v: (fixed, u, v)
        elif face in {"YMIN", "YMAX"}:
            fixed = 0 if face == "YMIN" else ny - 1
            nu, nv = nx, nz
            make = lambda u, v: (u, fixed, v)
        elif face in {"ZMIN", "ZMAX"}:
            fixed = 0 if face == "ZMIN" else nz - 1
            nu, nv = nx, ny
            make = lambda u, v: (u, v, fixed)
        else:
            raise ValueError(f"Unsupported CFD patch face: {patch.face}")
        u0 = max(0, int(math.floor((patch.center_u - patch.size_u / 2.0) * nu)))
        u1 = min(nu, int(math.ceil((patch.center_u + patch.size_u / 2.0) * nu)))
        v0 = max(0, int(math.floor((patch.center_v - patch.size_v / 2.0) * nv)))
        v1 = min(nv, int(math.ceil((patch.center_v + patch.size_v / 2.0) * nv)))
        return [make(u, v) for u in range(u0, u1) for v in range(v0, v1)]

    def generate_mesh(self, model, settings):
        if np is None:
            raise ImportError("NumPy is required for enclosure CFD meshing.")
        requested = max(0.5, float(settings.solver.cell_size_mm))
        dims_mm = model.dimensions_mm
        shape = tuple(max(3, int(math.ceil(size / requested))) for size in dims_mm)
        count = shape[0] * shape[1] * shape[2]
        if count > int(settings.solver.max_cells):
            raise ValueError(
                f"CFD mesh would contain {count:,} cells. Increase the CFD cell size."
            )
        dims_m = tuple(value * 1e-3 for value in dims_mm)
        spacing = tuple(dims_m[index] / shape[index] for index in range(3))
        solid = np.zeros(shape, dtype=bool)
        conductivity = np.zeros(shape, dtype=float)
        heat = np.zeros(shape, dtype=float)
        names = np.empty(shape, dtype=object)
        names[:] = "AIR"

        xs = (np.arange(shape[0]) + 0.5) * spacing[0] * 1000.0
        ys = (np.arange(shape[1]) + 0.5) * spacing[1] * 1000.0
        zs = (np.arange(shape[2]) + 0.5) * spacing[2] * 1000.0
        for obstacle in model.obstacles:
            x0, y0, z0, x1, y1, z1 = obstacle.bounds_mm
            mask = (
                (xs[:, None, None] >= x0) & (xs[:, None, None] <= x1) &
                (ys[None, :, None] >= y0) & (ys[None, :, None] <= y1) &
                (zs[None, None, :] >= z0) & (zs[None, None, :] <= z1)
            )
            if not np.any(mask):
                center = ((x0 + x1) / 2.0, (y0 + y1) / 2.0, (z0 + z1) / 2.0)
                index = tuple(min(shape[a] - 1, max(0, int(center[a] / dims_mm[a] * shape[a])))
                              for a in range(3))
                mask[index] = True
            solid |= mask
            conductivity[mask] = max(0.05, float(obstacle.conductivity_w_mk))
            names[mask] = obstacle.name
            cells = int(np.count_nonzero(mask))
            if obstacle.heat_w > 0 and cells:
                heat[mask] += float(obstacle.heat_w) / cells

        fluid = ~solid
        patches = []
        occupied = set()
        for patch in model.patches:
            kind = str(patch.kind).upper()
            if kind not in {"WALL", "INLET", "OUTLET", "VENT", "FAN"}:
                raise ValueError(f"Unsupported CFD patch type: {patch.kind}")
            normalized = (patch.center_u, patch.center_v, patch.size_u, patch.size_v)
            if any(float(value) < 0.0 or float(value) > 1.0 for value in normalized):
                raise ValueError(
                    f"CFD boundary patch '{patch.name}' coordinates must be between 0 and 1."
                )
            if patch.size_u <= 0.0 or patch.size_v <= 0.0:
                raise ValueError(f"CFD boundary patch '{patch.name}' must have a positive size.")
            if kind in {"INLET", "FAN"} and patch.velocity_m_s <= 0.0:
                raise ValueError(
                    f"CFD boundary patch '{patch.name}' requires a positive velocity."
                )
            cells = [cell for cell in self._patch_cells(patch, shape) if fluid[cell]]
            if not cells:
                raise ValueError(
                    f"CFD boundary patch '{patch.name}' has no fluid cells on its face."
                )
            duplicate = occupied.intersection(cells)
            if duplicate:
                raise ValueError(f"CFD boundary patch '{patch.name}' overlaps another patch.")
            occupied.update(cells)
            patches.append(CFDPatchCells(patch=patch, cells=cells))
        inlet_count = sum(
            1 for mapped in patches if str(mapped.patch.kind).upper() in {"INLET", "FAN"}
        )
        outlet_count = sum(
            1 for mapped in patches if str(mapped.patch.kind).upper() in {"OUTLET", "VENT"}
        )
        if inlet_count and not outlet_count:
            raise ValueError("Forced enclosure flow requires at least one OUTLET or VENT patch.")
        mesh = CFDMesh(
            shape=shape,
            spacing_m=spacing,
            dimensions_m=dims_m,
            fluid_mask=fluid,
            solid_mask=solid,
            solid_conductivity_w_mk=conductivity,
            heat_sources_w=heat,
            patch_cells=patches,
            obstacle_names=names,
        )
        self._log(
            f"Generated {shape[0]} x {shape[1]} x {shape[2]} = {count:,} cells "
            f"({np.count_nonzero(fluid):,} fluid)."
        )
        return mesh
