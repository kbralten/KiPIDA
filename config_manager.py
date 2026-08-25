import json
from typing import List, Optional, Union
from pathlib import Path

try:
    from .models import (
        ACAnalysisSettings, ACMeasurementPort, ACSourceModel, CapacitorModel,
        ComponentRef, PowerRail, ProjectConfig, UnifiedLoad, UnifiedSource,
        VoltageRegulator,
    )
except (ImportError, ValueError):
    from models import (
        ACAnalysisSettings, ACMeasurementPort, ACSourceModel, CapacitorModel,
        ComponentRef, PowerRail, ProjectConfig, UnifiedLoad, UnifiedSource,
        VoltageRegulator,
    )

CONFIG_VERSION = "1.1"
SUPPORTED_CONFIG_VERSIONS = {"1.0", CONFIG_VERSION}


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

def save_config(rails: List[PowerRail], filepath: str, ac_profiles=None):
    """Save the power tree and optional AC profiles to a JSON sidecar."""
    config = {
        "version": CONFIG_VERSION,
        "rails": [_rail_to_dict(rail) for rail in rails],
        "ac_profiles": {
            name: _ac_settings_to_dict(settings)
            for name, settings in (ac_profiles or {}).items()
        },
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
    return ProjectConfig(rails=rails, ac_profiles=ac_profiles)


def load_config(filepath: str) -> List[PowerRail]:
    """Backward-compatible rail-only loader."""
    return load_project_config(filepath).rails


def load_ac_profiles(filepath: str):
    return load_project_config(filepath).ac_profiles

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
        "distribution_mode": load.distribution_mode
    }

def _dict_to_load(data: dict) -> UnifiedLoad:
    """Convert dictionary to UnifiedLoad."""
    return UnifiedLoad(
        component_ref=ComponentRef(ref_des=data["ref_des"]),
        total_current=data.get("total_current", 0.0),
        pad_names=data.get("pad_names", []),
        distribution_mode=data.get("distribution_mode", "UNIFORM")
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
        "efficiency": reg.efficiency
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
        efficiency=data.get("efficiency", 0.85)
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
