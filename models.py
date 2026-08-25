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
    thermal_profile: Optional["ThermalAnalysisSettings"] = None
    cfd_profile: Optional["EnclosureCFDSettings"] = None


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


@dataclass
class DCSolveResult:
    """Detailed DC result used by electro-thermal analysis."""
    voltages: Dict[int, float] = field(default_factory=dict)
    branch_currents_a: List[float] = field(default_factory=list)
    branch_losses_w: List[float] = field(default_factory=list)
    total_loss_w: float = 0.0


@dataclass
class AirflowSettings:
    """Convective boundary settings for the exposed PCB surfaces."""
    mode: str = "NATURAL"  # NATURAL, FORCED, CUSTOM
    velocity_m_s: float = 0.0
    direction_deg: float = 0.0
    custom_h_w_m2k: float = 10.0
    expose_top: bool = True
    expose_bottom: bool = True
    expose_edges: bool = True


@dataclass
class ThermalComponentModel:
    """Compact component-to-board thermal model and heat source."""
    ref_des: str
    power_w: float = 0.0
    width_mm: float = 3.0
    depth_mm: float = 3.0
    height_mm: float = 1.0
    theta_jb_c_per_w: float = 20.0
    max_junction_c: float = 125.0
    enabled: bool = True
    model_source: str = "estimated"


@dataclass
class ThermalAnalysisSettings:
    """Persisted settings for steady-state 3D thermal analysis."""
    ambient_c: float = 25.0
    grid_size_mm: float = 1.0
    airflow: AirflowSettings = field(default_factory=AirflowSettings)
    include_radiation: bool = True
    emissivity: float = 0.9
    include_dc_copper_losses: bool = True
    coupled_iterations: int = 6
    convergence_c: float = 0.1
    relaxation: float = 0.6
    copper_temp_coefficient_per_c: float = 0.00393
    components: List[ThermalComponentModel] = field(default_factory=list)


@dataclass
class ThermalHotspot:
    node_id: int
    x_mm: float
    y_mm: float
    z_mm: float
    temperature_c: float


@dataclass
class ComponentThermalResult:
    ref_des: str
    board_temperature_c: float
    junction_temperature_c: float
    power_w: float
    max_junction_c: float
    margin_c: float
    model_source: str = "estimated"


@dataclass
class ThermalResult:
    temperatures_c: Dict[int, float] = field(default_factory=dict)
    hotspot: Optional[ThermalHotspot] = None
    component_results: List[ComponentThermalResult] = field(default_factory=list)
    total_input_power_w: float = 0.0
    total_boundary_power_w: float = 0.0
    energy_balance_error_pct: float = 0.0
    convection_coefficient_w_m2k: float = 0.0
    iterations: int = 1
    converged: bool = True


@dataclass
class ElectroThermalResult:
    thermal: ThermalResult
    dc_results: Dict[str, DCSolveResult] = field(default_factory=dict)
    iterations: int = 1
    converged: bool = True


@dataclass
class FluidProperties:
    """Air properties used by the enclosure CFD solver (SI units)."""
    density_kg_m3: float = 1.184
    dynamic_viscosity_pa_s: float = 1.85e-5
    heat_capacity_j_kgk: float = 1007.0
    conductivity_w_mk: float = 0.0262
    thermal_expansion_per_k: float = 0.00335


@dataclass
class CFDBoundaryPatch:
    """Rectangular boundary patch on one of the six enclosure faces."""
    name: str
    kind: str = "VENT"  # WALL, INLET, OUTLET, VENT, FAN
    face: str = "XMIN"  # XMIN/XMAX/YMIN/YMAX/ZMIN/ZMAX
    center_u: float = 0.5
    center_v: float = 0.5
    size_u: float = 0.25
    size_v: float = 0.25
    velocity_m_s: float = 0.0
    temperature_c: float = 25.0
    pressure_pa: float = 0.0


@dataclass
class EnclosureGeometrySettings:
    """Axis-aligned enclosure and PCB placement settings."""
    width_mm: float = 120.0
    depth_mm: float = 100.0
    height_mm: float = 50.0
    board_orientation: str = "XY"  # XY, XZ, YZ
    board_offset_x_mm: float = 0.0
    board_offset_y_mm: float = 0.0
    board_offset_z_mm: float = 15.0
    wall_heat_transfer_w_m2k: float = 5.0


@dataclass
class CFDSolverSettings:
    """Numerical controls for steady incompressible enclosure flow."""
    cell_size_mm: float = 5.0
    max_iterations: int = 250
    tolerance: float = 1.0e-4
    relaxation: float = 0.45
    pseudo_time_step_s: float = 0.02
    pressure_iterations: int = 60
    include_buoyancy: bool = True
    gravity_x_m_s2: float = 0.0
    gravity_y_m_s2: float = 0.0
    gravity_z_m_s2: float = -9.81
    max_cells: int = 250000


@dataclass
class EnclosureCFDSettings:
    """Persisted Phase 4 enclosure, fluid, boundary, and solver settings."""
    ambient_c: float = 25.0
    geometry: EnclosureGeometrySettings = field(default_factory=EnclosureGeometrySettings)
    fluid: FluidProperties = field(default_factory=FluidProperties)
    solver: CFDSolverSettings = field(default_factory=CFDSolverSettings)
    patches: List[CFDBoundaryPatch] = field(default_factory=list)
    use_phase3_heat_sources: bool = True
    include_dc_copper_losses: bool = True


@dataclass
class CFDResidualHistory:
    continuity: List[float] = field(default_factory=list)
    momentum: List[float] = field(default_factory=list)
    energy: List[float] = field(default_factory=list)


@dataclass
class EnclosureCFDResult:
    """Volumetric CFD fields and conservation diagnostics."""
    pressure_pa: List[float] = field(default_factory=list)
    velocity_u_m_s: List[float] = field(default_factory=list)
    velocity_v_m_s: List[float] = field(default_factory=list)
    velocity_w_m_s: List[float] = field(default_factory=list)
    air_temperature_c: List[float] = field(default_factory=list)
    solid_temperature_c: List[float] = field(default_factory=list)
    residuals: CFDResidualHistory = field(default_factory=CFDResidualHistory)
    iterations: int = 0
    converged: bool = False
    mass_balance_error_pct: float = 0.0
    energy_balance_error_pct: float = 0.0
    maximum_velocity_m_s: float = 0.0
    maximum_air_temperature_c: float = 0.0
    maximum_solid_temperature_c: float = 0.0
    total_heat_w: float = 0.0

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

