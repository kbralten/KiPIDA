from dataclasses import dataclass, field
from typing import Dict, List, Optional

@dataclass
class ComponentRef:
    """Reference to a specific footprint component on the board."""
    ref_des: str
    
    def __hash__(self):
        return hash(self.ref_des)
    
    def __eq__(self, other):
        if not isinstance(other, ComponentRef):
            return False
        return self.ref_des == other.ref_des

@dataclass
class UnifiedSource:
    """
    Represents a component acting as a voltage source.
    Voltage is defined by the parent PowerRail.
    """
    component_ref: ComponentRef
    pad_names: List[str] = field(default_factory=list)

@dataclass
class UnifiedLoad:
    """
    Represents a component acting as a current load.
    distribution_mode: 'UNIFORM' divides current equally among enabled pads.
    """
    component_ref: ComponentRef
    total_current: float = 0.0
    pad_names: List[str] = field(default_factory=list)
    distribution_mode: str = "UNIFORM" 

@dataclass
class VoltageRegulator:
    """
    Represents a voltage regulator connecting two PowerRails.
    Input/Output are defined by component RefDes and specific pads.
    """
    name: str  # Name of the regulator instance (e.g. "Buck 1")
    
    input_rail_name: str
    input_ref_des: str
    input_pad_names: List[str]
    
    output_rail_name: str
    output_ref_des: str
    output_pad_names: List[str]
    
    reg_type: str = "LINEAR"  # "LINEAR" or "SWITCHING"
    efficiency: float = 0.85  # Only used if SWITCHING. 0.0-1.0

@dataclass
class PowerRail:
    """
    High-level representation of a power domain.
    """
    net_name: str
    nominal_voltage: float = 0.0
    sources: List[UnifiedSource] = field(default_factory=list)
    loads: List[UnifiedLoad] = field(default_factory=list)
    # Regulators where this rail is the INPUT
    child_regulators: List[VoltageRegulator] = field(default_factory=list)
    
    def add_source(self, source: UnifiedSource):
        self.sources.append(source)
    
    def add_load(self, load: UnifiedLoad):
        self.loads.append(load)
    
    def add_child_regulator(self, reg: VoltageRegulator):
        self.child_regulators.append(reg)


@dataclass
class MeshBranch:
    """Physical branch retained for frequency-domain analysis."""
    node_a: int
    node_b: int
    resistance_ohm: float
    inductance_h: float = 0.0
    kind: str = "lateral"


@dataclass
class ACSourceModel:
    """Small-signal source impedance between a rail and its return net."""
    ref_des: str = ""
    rail_pad_names: List[str] = field(default_factory=list)
    ground_pad_names: List[str] = field(default_factory=list)
    resistance_ohm: float = 0.01
    inductance_h: float = 1.0e-9


@dataclass
class ACMeasurementPort:
    """Differential rail-to-ground observation port."""
    ref_des: str = ""
    rail_pad_names: List[str] = field(default_factory=list)
    ground_pad_names: List[str] = field(default_factory=list)


@dataclass
class CapacitorModel:
    """Lumped RLC model for a mounted decoupling capacitor."""
    ref_des: str
    rail_pad_names: List[str] = field(default_factory=list)
    ground_pad_names: List[str] = field(default_factory=list)
    capacitance_f: float = 0.0
    esr_ohm: float = 0.01
    esl_h: float = 0.8e-9
    enabled: bool = True
    candidate: bool = False
    model_source: str = "estimated"


@dataclass
class ACAnalysisSettings:
    """Persisted settings for one rail impedance profile."""
    rail_name: str = ""
    ground_net_name: str = "GND"
    frequency_start_hz: float = 1.0e3
    frequency_stop_hz: float = 1.0e8
    frequency_points: int = 121
    target_impedance_ohm: float = 0.05
    source: ACSourceModel = field(default_factory=ACSourceModel)
    measurement_port: ACMeasurementPort = field(default_factory=ACMeasurementPort)
    capacitors: List[CapacitorModel] = field(default_factory=list)
    optimizer_values_f: List[float] = field(
        default_factory=lambda: [10e-9, 47e-9, 100e-9, 470e-9, 1e-6, 4.7e-6, 10e-6]
    )
    optimizer_max_additions: int = 8


@dataclass
class ProjectConfig:
    """Complete persisted Ki-PIDA project configuration."""
    rails: List[PowerRail] = field(default_factory=list)
    ac_profiles: Dict[str, ACAnalysisSettings] = field(default_factory=dict)


@dataclass
class ImpedanceSweepResult:
    """Complex impedance sweep and derived summary values."""
    frequencies_hz: List[float]
    impedance_ohm: List[complex]
    target_impedance_ohm: float = 0.0
    worst_frequency_hz: float = 0.0
    worst_impedance_ohm: float = 0.0
    meets_target: bool = False


@dataclass
class OptimizationRecommendation:
    ref_des: str
    capacitance_f: float
    action: str = "populate"


@dataclass
class DecouplingOptimizationResult:
    baseline: ImpedanceSweepResult
    optimized: ImpedanceSweepResult
    recommendations: List[OptimizationRecommendation] = field(default_factory=list)
    reached_target: bool = False

def generate_regulator_name(input_ref_des: str, output_ref_des: str, output_rail_name: str = "") -> str:
    """
    Generate a regulator name based on input and output components and rails.
    If input == output, name is "input (output_rail)".
    Otherwise, name is "input -> output (output_rail)".
    """
    name = ""
    if not input_ref_des:
        name = output_ref_des if output_ref_des else ""
    elif not output_ref_des:
        name = input_ref_des
    elif input_ref_des == output_ref_des:
        name = input_ref_des
    else:
        name = f"{input_ref_des} -> {output_ref_des}"
        
    if output_rail_name:
        name += f" ({output_rail_name})"
        
    return name

