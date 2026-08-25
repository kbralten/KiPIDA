import json
from typing import List, Optional, Union
from pathlib import Path

try:
    from .models import PowerRail, UnifiedSource, UnifiedLoad, VoltageRegulator, ComponentRef
except (ImportError, ValueError):
    from models import PowerRail, UnifiedSource, UnifiedLoad, VoltageRegulator, ComponentRef

CONFIG_VERSION = "1.0"


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

def save_config(rails: List[PowerRail], filepath: str):
    """
    Save power network configuration to JSON file.
    
    Args:
        rails: List of PowerRail objects to serialize
        filepath: Path to save JSON file
    """
    config = {
        "version": CONFIG_VERSION,
        "rails": [_rail_to_dict(rail) for rail in rails]
    }
    
    with open(filepath, 'w') as f:
        json.dump(config, f, indent=2)

def load_config(filepath: str) -> List[PowerRail]:
    """
    Load power network configuration from JSON file.
    
    Args:
        filepath: Path to JSON config file
        
    Returns:
        List of PowerRail objects
        
    Raises:
        FileNotFoundError: If config file doesn't exist
        ValueError: If config file is invalid
    """
    if not Path(filepath).exists():
        raise FileNotFoundError(f"Config file not found: {filepath}")
    
    with open(filepath, 'r') as f:
        config = json.load(f)
    
    # Validate version
    version = config.get("version", "unknown")
    if version != CONFIG_VERSION:
        raise ValueError(f"Unsupported config version: {version}")
    
    # Deserialize rails
    rails = [_dict_to_rail(rail_dict) for rail_dict in config.get("rails", [])]
    return rails

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
