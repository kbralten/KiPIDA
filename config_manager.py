import json
from typing import List, Optional, Union
from pathlib import Path

try:
    from .models import (
        ACAnalysisSettings, ACMeasurementPort, ACSourceModel, CapacitorModel,
        AirflowSettings, CFDBoundaryPatch, CFDSolverSettings, ComponentRef,
        EnclosureCFDSettings, EnclosureGeometrySettings, FluidProperties,
        PowerRail, ProjectConfig,
        ThermalAnalysisSettings, ThermalComponentModel, UnifiedLoad,
        UnifiedSource, VoltageRegulator,
    )
except (ImportError, ValueError):
    from models import (
        ACAnalysisSettings, ACMeasurementPort, ACSourceModel, CapacitorModel,
        AirflowSettings, CFDBoundaryPatch, CFDSolverSettings, ComponentRef,
        EnclosureCFDSettings, EnclosureGeometrySettings, FluidProperties,
        PowerRail, ProjectConfig,
        ThermalAnalysisSettings, ThermalComponentModel, UnifiedLoad,
        UnifiedSource, VoltageRegulator,
    )

CONFIG_VERSION = "1.4"
SUPPORTED_CONFIG_VERSIONS = {"1.0", "1.1", "1.2", "1.3", CONFIG_VERSION}


def get_project_config_path(
    project_path: Union[str, Path], project_name: Optional[str] = None
) -> Path:
    """Return the sidecar config path for a KiCad project or board path.

    KiCad 10's IPC API exposes ``Project.path`` as the ``.kicad_pro`` file,
    whereas older integrations can expose a project directory. Supporting both
    shapes keeps the config file beside the project rather than trying to
    create a directory underneath the project file.
    """
    path = Path(project_path)
    is_project_file = path.suffix.lower() in {".kicad_pro", ".pro"}

    if is_project_file:
        directory = path.parent
        config_stem = path.stem
    else:
        directory = path
        config_stem = project_name or path.name

    return directory / f"{config_stem}.kipida.json"

def save_config(
    rails: List[PowerRail], filepath: str, ac_profiles=None,
    thermal_profile=None, cfd_profile=None,
):
    """Save the power tree and optional AC, thermal, and CFD profiles."""
    config = {
        "version": CONFIG_VERSION,
        "rails": [_rail_to_dict(rail) for rail in rails],
        "ac_profiles": {
            name: _ac_settings_to_dict(settings)
            for name, settings in (ac_profiles or {}).items()
        },
        "thermal_profile": (
            _thermal_settings_to_dict(thermal_profile) if thermal_profile is not None else None
        ),
        "cfd_profile": _cfd_settings_to_dict(cfd_profile) if cfd_profile is not None else None,
    }
    
    with open(filepath, 'w') as f:
        json.dump(config, f, indent=2)

def load_project_config(filepath: str) -> ProjectConfig:
    """Load the complete project configuration, including legacy v1.0 files."""
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Config file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        config = json.load(f)
    
    # Validate version
    version = config.get("version", "unknown")
    if version not in SUPPORTED_CONFIG_VERSIONS:
        raise ValueError(f"Unsupported config version: {version}")
    
    rails = [_dict_to_rail(rail_dict) for rail_dict in config.get("rails", [])]
    ac_profiles = {
        name: _dict_to_ac_settings(settings)
        for name, settings in config.get("ac_profiles", {}).items()
    }
    thermal_data = config.get("thermal_profile")
    thermal_profile = _dict_to_thermal_settings(thermal_data) if thermal_data else None
    cfd_data = config.get("cfd_profile")
    cfd_profile = _dict_to_cfd_settings(cfd_data) if cfd_data else None
    return ProjectConfig(
        rails=rails,
        ac_profiles=ac_profiles,
        thermal_profile=thermal_profile,
        cfd_profile=cfd_profile,
    )


def load_config(filepath: str) -> List[PowerRail]:
    """Backward-compatible rail-only loader."""
    return load_project_config(filepath).rails


def load_ac_profiles(filepath: str):
    return load_project_config(filepath).ac_profiles


def load_thermal_profile(filepath: str):
    return load_project_config(filepath).thermal_profile


def load_cfd_profile(filepath: str):
    return load_project_config(filepath).cfd_profile

def _rail_to_dict(rail: PowerRail) -> dict:
    """Convert PowerRail to dictionary."""
    return {
        "net_name": rail.net_name,
        "nominal_voltage": rail.nominal_voltage,
        "sources": [_source_to_dict(src) for src in rail.sources],
        "loads": [_load_to_dict(load) for load in rail.loads],
        "child_regulators": [_regulator_to_dict(reg) for reg in rail.child_regulators]
    }

def _dict_to_rail(data: dict) -> PowerRail:
    """Convert dictionary to PowerRail."""
    rail = PowerRail(
        net_name=data["net_name"],
        nominal_voltage=data.get("nominal_voltage", 0.0)
    )
    
    # Deserialize sources
    for src_dict in data.get("sources", []):
        rail.sources.append(_dict_to_source(src_dict))
    
    # Deserialize loads
    for load_dict in data.get("loads", []):
        rail.loads.append(_dict_to_load(load_dict))
    
    # Deserialize regulators
    for reg_dict in data.get("child_regulators", []):
        rail.child_regulators.append(_dict_to_regulator(reg_dict))
    
    return rail

def _source_to_dict(source: UnifiedSource) -> dict:
    """Convert UnifiedSource to dictionary."""
    return {
        "ref_des": source.component_ref.ref_des,
        "pad_names": source.pad_names
    }

def _dict_to_source(data: dict) -> UnifiedSource:
    """Convert dictionary to UnifiedSource."""
    return UnifiedSource(
        component_ref=ComponentRef(ref_des=data["ref_des"]),
        pad_names=data.get("pad_names", [])
    )

def _load_to_dict(load: UnifiedLoad) -> dict:
    """Convert UnifiedLoad to dictionary."""
    return {
        "ref_des": load.component_ref.ref_des,
        "total_current": load.total_current,
        "pad_names": load.pad_names,
        "distribution_mode": load.distribution_mode,
        "thermal_mode": load.thermal_mode,
    }

def _dict_to_load(data: dict) -> UnifiedLoad:
    """Convert dictionary to UnifiedLoad."""
    return UnifiedLoad(
        component_ref=ComponentRef(ref_des=data["ref_des"]),
        total_current=data.get("total_current", 0.0),
        pad_names=data.get("pad_names", []),
        distribution_mode=data.get("distribution_mode", "UNIFORM"),
        thermal_mode=data.get("thermal_mode", "AUTO"),
    )

def _regulator_to_dict(reg: VoltageRegulator) -> dict:
    """Convert VoltageRegulator to dictionary."""
    return {
        "name": reg.name,
        "input_rail_name": reg.input_rail_name,
        "input_ref_des": reg.input_ref_des,
        "input_pad_names": reg.input_pad_names,
        "output_rail_name": reg.output_rail_name,
        "output_ref_des": reg.output_ref_des,
        "output_pad_names": reg.output_pad_names,
        "reg_type": reg.reg_type,
        "efficiency": reg.efficiency,
        "thermal_ref_des": reg.thermal_ref_des,
    }

def _dict_to_regulator(data: dict) -> VoltageRegulator:
    """Convert dictionary to VoltageRegulator."""
    return VoltageRegulator(
        name=data["name"],
        input_rail_name=data["input_rail_name"],
        input_ref_des=data["input_ref_des"],
        input_pad_names=data.get("input_pad_names", []),
        output_rail_name=data["output_rail_name"],
        output_ref_des=data["output_ref_des"],
        output_pad_names=data.get("output_pad_names", []),
        reg_type=data.get("reg_type", "LINEAR"),
        efficiency=data.get("efficiency", 0.85),
        thermal_ref_des=data.get("thermal_ref_des", ""),
    )


def _ac_settings_to_dict(settings: ACAnalysisSettings) -> dict:
    return {
        "rail_name": settings.rail_name,
        "ground_net_name": settings.ground_net_name,
        "frequency_start_hz": settings.frequency_start_hz,
        "frequency_stop_hz": settings.frequency_stop_hz,
        "frequency_points": settings.frequency_points,
        "target_impedance_ohm": settings.target_impedance_ohm,
        "source": {
            "ref_des": settings.source.ref_des,
            "rail_pad_names": settings.source.rail_pad_names,
            "ground_pad_names": settings.source.ground_pad_names,
            "resistance_ohm": settings.source.resistance_ohm,
            "inductance_h": settings.source.inductance_h,
        },
        "measurement_port": {
            "ref_des": settings.measurement_port.ref_des,
            "rail_pad_names": settings.measurement_port.rail_pad_names,
            "ground_pad_names": settings.measurement_port.ground_pad_names,
        },
        "capacitors": [{
            "ref_des": cap.ref_des,
            "rail_pad_names": cap.rail_pad_names,
            "ground_pad_names": cap.ground_pad_names,
            "capacitance_f": cap.capacitance_f,
            "esr_ohm": cap.esr_ohm,
            "esl_h": cap.esl_h,
            "enabled": cap.enabled,
            "candidate": cap.candidate,
            "model_source": cap.model_source,
        } for cap in settings.capacitors],
        "optimizer_values_f": settings.optimizer_values_f,
        "optimizer_max_additions": settings.optimizer_max_additions,
    }


def _dict_to_ac_settings(data: dict) -> ACAnalysisSettings:
    source_data = data.get("source", {})
    port_data = data.get("measurement_port", {})
    return ACAnalysisSettings(
        rail_name=data.get("rail_name", ""),
        ground_net_name=data.get("ground_net_name", "GND"),
        frequency_start_hz=float(data.get("frequency_start_hz", 1e3)),
        frequency_stop_hz=float(data.get("frequency_stop_hz", 1e8)),
        frequency_points=int(data.get("frequency_points", 121)),
        target_impedance_ohm=float(data.get("target_impedance_ohm", 0.05)),
        source=ACSourceModel(
            ref_des=source_data.get("ref_des", ""),
            rail_pad_names=source_data.get("rail_pad_names", []),
            ground_pad_names=source_data.get("ground_pad_names", []),
            resistance_ohm=float(source_data.get("resistance_ohm", 0.01)),
            inductance_h=float(source_data.get("inductance_h", 1e-9)),
        ),
        measurement_port=ACMeasurementPort(
            ref_des=port_data.get("ref_des", ""),
            rail_pad_names=port_data.get("rail_pad_names", []),
            ground_pad_names=port_data.get("ground_pad_names", []),
        ),
        capacitors=[CapacitorModel(
            ref_des=cap["ref_des"],
            rail_pad_names=cap.get("rail_pad_names", []),
            ground_pad_names=cap.get("ground_pad_names", []),
            capacitance_f=float(cap.get("capacitance_f", 0.0)),
            esr_ohm=float(cap.get("esr_ohm", 0.01)),
            esl_h=float(cap.get("esl_h", 0.8e-9)),
            enabled=bool(cap.get("enabled", True)),
            candidate=bool(cap.get("candidate", False)),
            model_source=cap.get("model_source", "estimated"),
        ) for cap in data.get("capacitors", [])],
        optimizer_values_f=[float(value) for value in data.get(
            "optimizer_values_f", [10e-9, 47e-9, 100e-9, 470e-9, 1e-6, 4.7e-6, 10e-6]
        )],
        optimizer_max_additions=int(data.get("optimizer_max_additions", 8)),
    )


def _thermal_settings_to_dict(settings: ThermalAnalysisSettings) -> dict:
    return {
        "ambient_c": settings.ambient_c,
        "grid_size_mm": settings.grid_size_mm,
        "airflow": {
            "mode": settings.airflow.mode,
            "velocity_m_s": settings.airflow.velocity_m_s,
            "direction_deg": settings.airflow.direction_deg,
            "custom_h_w_m2k": settings.airflow.custom_h_w_m2k,
            "expose_top": settings.airflow.expose_top,
            "expose_bottom": settings.airflow.expose_bottom,
            "expose_edges": settings.airflow.expose_edges,
        },
        "include_radiation": settings.include_radiation,
        "emissivity": settings.emissivity,
        "include_dc_copper_losses": settings.include_dc_copper_losses,
        "coupled_iterations": settings.coupled_iterations,
        "convergence_c": settings.convergence_c,
        "relaxation": settings.relaxation,
        "copper_temp_coefficient_per_c": settings.copper_temp_coefficient_per_c,
        "components": [{
            "ref_des": component.ref_des,
            "power_w": component.power_w,
            "width_mm": component.width_mm,
            "depth_mm": component.depth_mm,
            "height_mm": component.height_mm,
            "theta_jb_c_per_w": component.theta_jb_c_per_w,
            "max_junction_c": component.max_junction_c,
            "enabled": component.enabled,
            "model_source": component.model_source,
        } for component in settings.components],
    }


def _dict_to_thermal_settings(data: dict) -> ThermalAnalysisSettings:
    airflow = data.get("airflow", {})
    return ThermalAnalysisSettings(
        ambient_c=float(data.get("ambient_c", 25.0)),
        grid_size_mm=float(data.get("grid_size_mm", 1.0)),
        airflow=AirflowSettings(
            mode=airflow.get("mode", "NATURAL"),
            velocity_m_s=float(airflow.get("velocity_m_s", 0.0)),
            direction_deg=float(airflow.get("direction_deg", 0.0)),
            custom_h_w_m2k=float(airflow.get("custom_h_w_m2k", 10.0)),
            expose_top=bool(airflow.get("expose_top", True)),
            expose_bottom=bool(airflow.get("expose_bottom", True)),
            expose_edges=bool(airflow.get("expose_edges", True)),
        ),
        include_radiation=bool(data.get("include_radiation", True)),
        emissivity=float(data.get("emissivity", 0.9)),
        include_dc_copper_losses=bool(data.get("include_dc_copper_losses", True)),
        coupled_iterations=int(data.get("coupled_iterations", 10)),
        convergence_c=float(data.get("convergence_c", 0.1)),
        relaxation=float(data.get("relaxation", 0.6)),
        copper_temp_coefficient_per_c=float(data.get("copper_temp_coefficient_per_c", 0.00393)),
        components=[ThermalComponentModel(
            ref_des=component["ref_des"],
            power_w=float(component.get("power_w", 0.0)),
            width_mm=float(component.get("width_mm", 3.0)),
            depth_mm=float(component.get("depth_mm", 3.0)),
            height_mm=float(component.get("height_mm", 1.0)),
            theta_jb_c_per_w=float(component.get("theta_jb_c_per_w", 20.0)),
            max_junction_c=float(component.get("max_junction_c", 125.0)),
            enabled=bool(component.get("enabled", True)),
            model_source=component.get("model_source", "estimated"),
        ) for component in data.get("components", [])],
    )


def _cfd_settings_to_dict(settings: EnclosureCFDSettings) -> dict:
    geometry = settings.geometry
    fluid = settings.fluid
    solver = settings.solver
    return {
        "ambient_c": settings.ambient_c,
        "geometry": {
            "width_mm": geometry.width_mm,
            "depth_mm": geometry.depth_mm,
            "height_mm": geometry.height_mm,
            "board_orientation": geometry.board_orientation,
            "board_offset_x_mm": geometry.board_offset_x_mm,
            "board_offset_y_mm": geometry.board_offset_y_mm,
            "board_offset_z_mm": geometry.board_offset_z_mm,
            "wall_heat_transfer_w_m2k": geometry.wall_heat_transfer_w_m2k,
        },
        "fluid": {
            "density_kg_m3": fluid.density_kg_m3,
            "dynamic_viscosity_pa_s": fluid.dynamic_viscosity_pa_s,
            "heat_capacity_j_kgk": fluid.heat_capacity_j_kgk,
            "conductivity_w_mk": fluid.conductivity_w_mk,
            "thermal_expansion_per_k": fluid.thermal_expansion_per_k,
        },
        "solver": {
            "cell_size_mm": solver.cell_size_mm,
            "max_iterations": solver.max_iterations,
            "tolerance": solver.tolerance,
            "relaxation": solver.relaxation,
            "pseudo_time_step_s": solver.pseudo_time_step_s,
            "pressure_iterations": solver.pressure_iterations,
            "include_buoyancy": solver.include_buoyancy,
            "gravity_x_m_s2": solver.gravity_x_m_s2,
            "gravity_y_m_s2": solver.gravity_y_m_s2,
            "gravity_z_m_s2": solver.gravity_z_m_s2,
            "max_cells": solver.max_cells,
        },
        "patches": [{
            "name": patch.name,
            "kind": patch.kind,
            "face": patch.face,
            "center_u": patch.center_u,
            "center_v": patch.center_v,
            "size_u": patch.size_u,
            "size_v": patch.size_v,
            "velocity_m_s": patch.velocity_m_s,
            "temperature_c": patch.temperature_c,
            "pressure_pa": patch.pressure_pa,
        } for patch in settings.patches],
        "use_phase3_heat_sources": settings.use_phase3_heat_sources,
        "include_dc_copper_losses": settings.include_dc_copper_losses,
    }


def _dict_to_cfd_settings(data: dict) -> EnclosureCFDSettings:
    geometry = data.get("geometry", {})
    fluid = data.get("fluid", {})
    solver = data.get("solver", {})
    return EnclosureCFDSettings(
        ambient_c=float(data.get("ambient_c", 25.0)),
        geometry=EnclosureGeometrySettings(
            width_mm=float(geometry.get("width_mm", 120.0)),
            depth_mm=float(geometry.get("depth_mm", 100.0)),
            height_mm=float(geometry.get("height_mm", 50.0)),
            board_orientation=geometry.get("board_orientation", "XY"),
            board_offset_x_mm=float(geometry.get("board_offset_x_mm", 0.0)),
            board_offset_y_mm=float(geometry.get("board_offset_y_mm", 0.0)),
            board_offset_z_mm=float(geometry.get("board_offset_z_mm", 15.0)),
            wall_heat_transfer_w_m2k=float(
                geometry.get("wall_heat_transfer_w_m2k", 5.0)
            ),
        ),
        fluid=FluidProperties(
            density_kg_m3=float(fluid.get("density_kg_m3", 1.184)),
            dynamic_viscosity_pa_s=float(fluid.get("dynamic_viscosity_pa_s", 1.85e-5)),
            heat_capacity_j_kgk=float(fluid.get("heat_capacity_j_kgk", 1007.0)),
            conductivity_w_mk=float(fluid.get("conductivity_w_mk", 0.0262)),
            thermal_expansion_per_k=float(fluid.get("thermal_expansion_per_k", 0.00335)),
        ),
        solver=CFDSolverSettings(
            cell_size_mm=float(solver.get("cell_size_mm", 5.0)),
            max_iterations=int(solver.get("max_iterations", 250)),
            tolerance=float(solver.get("tolerance", 1e-4)),
            relaxation=float(solver.get("relaxation", 0.45)),
            pseudo_time_step_s=float(solver.get("pseudo_time_step_s", 0.02)),
            pressure_iterations=int(solver.get("pressure_iterations", 60)),
            include_buoyancy=bool(solver.get("include_buoyancy", True)),
            gravity_x_m_s2=float(solver.get("gravity_x_m_s2", 0.0)),
            gravity_y_m_s2=float(solver.get("gravity_y_m_s2", 0.0)),
            gravity_z_m_s2=float(solver.get("gravity_z_m_s2", -9.81)),
            max_cells=int(solver.get("max_cells", 250000)),
        ),
        patches=[CFDBoundaryPatch(
            name=patch.get("name", "Patch"),
            kind=patch.get("kind", "VENT"),
            face=patch.get("face", "XMIN"),
            center_u=float(patch.get("center_u", 0.5)),
            center_v=float(patch.get("center_v", 0.5)),
            size_u=float(patch.get("size_u", 0.25)),
            size_v=float(patch.get("size_v", 0.25)),
            velocity_m_s=float(patch.get("velocity_m_s", 0.0)),
            temperature_c=float(patch.get("temperature_c", 25.0)),
            pressure_pa=float(patch.get("pressure_pa", 0.0)),
        ) for patch in data.get("patches", [])],
        use_phase3_heat_sources=bool(data.get("use_phase3_heat_sources", True)),
        include_dc_copper_losses=bool(data.get("include_dc_copper_losses", True)),
    )
