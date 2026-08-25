import unittest
import json
import tempfile
import os
from pathlib import Path

import sys
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from config_manager import get_project_config_path, save_config, load_config, load_project_config
from models import (
    ACAnalysisSettings, ACMeasurementPort, ACSourceModel, AirflowSettings, CapacitorModel,
    CFDBoundaryPatch, CFDSolverSettings, EnclosureCFDSettings, EnclosureGeometrySettings,
    FluidProperties,
    PowerRail, UnifiedSource, UnifiedLoad, VoltageRegulator, ComponentRef,
    ThermalAnalysisSettings, ThermalComponentModel,
)

class TestConfigManager(unittest.TestCase):

    def test_kicad_10_project_file_config_path(self):
        project_file = Path(r"C:\projects\p02_alimentation\p02_alimentation.kicad_pro")

        self.assertEqual(
            get_project_config_path(project_file, "p02_alimentation"),
            Path(r"C:\projects\p02_alimentation\p02_alimentation.kipida.json"),
        )

    def test_project_directory_config_path(self):
        project_dir = Path(r"C:\projects\p02_alimentation")

        self.assertEqual(
            get_project_config_path(project_dir, "p02_alimentation"),
            project_dir / "p02_alimentation.kipida.json",
        )
    
    def setUp(self):
        """Create sample power network configuration."""
        # Create rails
        self.rail_12v = PowerRail(net_name="12V", nominal_voltage=12.0)
        self.rail_5v = PowerRail(net_name="5V", nominal_voltage=5.0)
        self.rail_3v3 = PowerRail(net_name="3V3", nominal_voltage=3.3)
        
        # Add source to 12V
        self.rail_12v.add_source(UnifiedSource(
            component_ref=ComponentRef(ref_des="J1"),
            pad_names=["1", "2"]
        ))
        
        # Add load to 5V
        self.rail_5v.add_load(UnifiedLoad(
            component_ref=ComponentRef(ref_des="U1"),
            total_current=1.5,
            pad_names=["VDD"],
            distribution_mode="UNIFORM",
            thermal_mode="LOCAL",
        ))
        
        # Add regulator from 12V to 5V
        self.rail_12v.add_child_regulator(VoltageRegulator(
            name="Buck1",
            input_rail_name="12V",
            input_ref_des="U2",
            input_pad_names=["VIN"],
            output_rail_name="5V",
            output_ref_des="U2",
            output_pad_names=["VOUT"],
            reg_type="SWITCHING",
            efficiency=0.90,
            thermal_ref_des="U2",
        ))
        
        # Add regulator from 5V to 3V3
        self.rail_5v.add_child_regulator(VoltageRegulator(
            name="LDO1",
            input_rail_name="5V",
            input_ref_des="U3",
            input_pad_names=["IN"],
            output_rail_name="3V3",
            output_ref_des="U3",
            output_pad_names=["OUT"],
            reg_type="LINEAR",
            efficiency=1.0
        ))
        
        self.rails = [self.rail_12v, self.rail_5v, self.rail_3v3]
    
    def test_serialization(self):
        """Test that PowerRail objects can be serialized to JSON."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
        
        try:
            save_config(self.rails, filepath)
            
            # Verify file exists
            self.assertTrue(Path(filepath).exists())
            
            # Load and verify JSON structure
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            self.assertEqual(data["version"], "1.4")
            self.assertEqual(len(data["rails"]), 3)
            
            # Verify 12V rail
            rail_12v = data["rails"][0]
            self.assertEqual(rail_12v["net_name"], "12V")
            self.assertEqual(rail_12v["nominal_voltage"], 12.0)
            self.assertEqual(len(rail_12v["sources"]), 1)
            self.assertEqual(len(rail_12v["child_regulators"]), 1)
            
            # Verify source
            src = rail_12v["sources"][0]
            self.assertEqual(src["ref_des"], "J1")
            self.assertEqual(src["pad_names"], ["1", "2"])
            
            # Verify regulator
            reg = rail_12v["child_regulators"][0]
            self.assertEqual(reg["name"], "Buck1")
            self.assertEqual(reg["reg_type"], "SWITCHING")
            self.assertEqual(reg["efficiency"], 0.90)
            self.assertEqual(reg["thermal_ref_des"], "U2")

            load = data["rails"][1]["loads"][0]
            self.assertEqual(load["thermal_mode"], "LOCAL")
            
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)
    
    def test_deserialization(self):
        """Test that JSON can be deserialized to PowerRail objects."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
        
        try:
            save_config(self.rails, filepath)
            loaded_rails = load_config(filepath)
            
            # Verify rail count
            self.assertEqual(len(loaded_rails), 3)
            
            # Verify 12V rail
            rail_12v = loaded_rails[0]
            self.assertEqual(rail_12v.net_name, "12V")
            self.assertEqual(rail_12v.nominal_voltage, 12.0)
            self.assertEqual(len(rail_12v.sources), 1)
            self.assertEqual(len(rail_12v.child_regulators), 1)
            
            # Verify source
            src = rail_12v.sources[0]
            self.assertEqual(src.component_ref.ref_des, "J1")
            self.assertEqual(src.pad_names, ["1", "2"])
            
            # Verify load on 5V
            rail_5v = loaded_rails[1]
            self.assertEqual(len(rail_5v.loads), 1)
            load = rail_5v.loads[0]
            self.assertEqual(load.component_ref.ref_des, "U1")
            self.assertEqual(load.total_current, 1.5)
            self.assertEqual(load.pad_names, ["VDD"])
            self.assertEqual(load.thermal_mode, "LOCAL")
            
            # Verify regulator
            reg = rail_12v.child_regulators[0]
            self.assertEqual(reg.name, "Buck1")
            self.assertEqual(reg.input_rail_name, "12V")
            self.assertEqual(reg.output_rail_name, "5V")
            self.assertEqual(reg.reg_type, "SWITCHING")
            self.assertEqual(reg.efficiency, 0.90)
            self.assertEqual(reg.thermal_ref_des, "U2")
            
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)
    
    def test_round_trip(self):
        """Test that save → load → save produces identical results."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath1 = f.name
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath2 = f.name
        
        try:
            # Save original
            save_config(self.rails, filepath1)
            
            # Load and save again
            loaded_rails = load_config(filepath1)
            save_config(loaded_rails, filepath2)
            
            # Compare JSON files
            with open(filepath1, 'r') as f:
                data1 = json.load(f)
            with open(filepath2, 'r') as f:
                data2 = json.load(f)
            
            self.assertEqual(data1, data2)
            
        finally:
            if Path(filepath1).exists():
                os.unlink(filepath1)
            if Path(filepath2).exists():
                os.unlink(filepath2)
    
    def test_empty_rails(self):
        """Test saving and loading empty rails list."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
        
        try:
            save_config([], filepath)
            loaded_rails = load_config(filepath)
            self.assertEqual(len(loaded_rails), 0)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_ac_profile_round_trip(self):
        """AC settings are stored beside the existing power-tree configuration."""
        profile = ACAnalysisSettings(
            rail_name="5V",
            ground_net_name="GND",
            frequency_start_hz=100.0,
            frequency_stop_hz=20e6,
            frequency_points=77,
            target_impedance_ohm=0.025,
            source=ACSourceModel(
                ref_des="U2", rail_pad_names=["VOUT"], ground_pad_names=["GND"],
                resistance_ohm=0.008, inductance_h=0.7e-9,
            ),
            measurement_port=ACMeasurementPort(
                ref_des="U1", rail_pad_names=["VDD"], ground_pad_names=["VSS"],
            ),
            capacitors=[CapacitorModel(
                ref_des="C12", rail_pad_names=["1"], ground_pad_names=["2"],
                capacitance_f=4.7e-6, esr_ohm=0.015, esl_h=0.6e-9,
            )],
            optimizer_values_f=[100e-9, 1e-6],
            optimizer_max_additions=3,
        )
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name

        try:
            save_config(self.rails, filepath, {"5V": profile})
            project = load_project_config(filepath)
            loaded = project.ac_profiles["5V"]

            self.assertEqual(len(project.rails), 3)
            self.assertEqual(loaded.source.ref_des, "U2")
            self.assertEqual(loaded.measurement_port.ground_pad_names, ["VSS"])
            self.assertAlmostEqual(loaded.capacitors[0].capacitance_f, 4.7e-6)
            self.assertEqual(loaded.optimizer_max_additions, 3)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_legacy_v1_config_migrates_without_ac_profiles(self):
        """Existing Phase 1 files remain loadable after the schema extension."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
            json.dump({"version": "1.0", "rails": []}, f)

        try:
            project = load_project_config(filepath)
            self.assertEqual(project.rails, [])
            self.assertEqual(project.ac_profiles, {})
            self.assertIsNone(project.thermal_profile)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_thermal_profile_round_trip(self):
        profile = ThermalAnalysisSettings(
            ambient_c=35.0,
            grid_size_mm=1.5,
            airflow=AirflowSettings(
                mode="FORCED", velocity_m_s=2.5, direction_deg=0.0,
                expose_top=True, expose_bottom=False, expose_edges=True,
            ),
            include_radiation=True,
            emissivity=0.82,
            components=[ThermalComponentModel(
                ref_des="U2", power_w=1.25, width_mm=4.0, depth_mm=4.0,
                theta_jb_c_per_w=8.0, max_junction_c=150.0,
                model_source="user",
            )],
        )
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name

        try:
            save_config(self.rails, filepath, thermal_profile=profile)
            loaded = load_project_config(filepath).thermal_profile

            self.assertIsNotNone(loaded)
            self.assertEqual(loaded.airflow.mode, "FORCED")
            self.assertAlmostEqual(loaded.airflow.velocity_m_s, 2.5)
            self.assertFalse(loaded.airflow.expose_bottom)
            self.assertTrue(loaded.include_radiation)
            self.assertEqual(loaded.components[0].ref_des, "U2")
            self.assertAlmostEqual(loaded.components[0].power_w, 1.25)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_legacy_v11_config_migrates_without_thermal_profile(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
            json.dump({"version": "1.1", "rails": [], "ac_profiles": {}}, f)

        try:
            project = load_project_config(filepath)
            self.assertIsNone(project.thermal_profile)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_cfd_profile_round_trip(self):
        profile = EnclosureCFDSettings(
            ambient_c=32.0,
            geometry=EnclosureGeometrySettings(
                width_mm=180.0, depth_mm=120.0, height_mm=65.0,
                board_orientation="XZ", board_offset_x_mm=3.0,
                wall_heat_transfer_w_m2k=7.5,
            ),
            fluid=FluidProperties(density_kg_m3=1.16, conductivity_w_mk=0.027),
            solver=CFDSolverSettings(
                cell_size_mm=4.0, max_iterations=350, tolerance=2e-5,
                relaxation=0.35, include_buoyancy=False, max_cells=123456,
            ),
            patches=[CFDBoundaryPatch(
                "Front fan", "FAN", "XMIN", 0.5, 0.6, 0.3, 0.4,
                velocity_m_s=1.4, temperature_c=29.0,
            )],
            use_phase3_heat_sources=True,
            include_dc_copper_losses=False,
        )
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name

        try:
            save_config(self.rails, filepath, cfd_profile=profile)
            loaded = load_project_config(filepath).cfd_profile

            self.assertIsNotNone(loaded)
            self.assertEqual(loaded.geometry.board_orientation, "XZ")
            self.assertAlmostEqual(loaded.geometry.width_mm, 180.0)
            self.assertAlmostEqual(loaded.fluid.density_kg_m3, 1.16)
            self.assertFalse(loaded.solver.include_buoyancy)
            self.assertEqual(loaded.solver.max_cells, 123456)
            self.assertEqual(loaded.patches[0].kind, "FAN")
            self.assertAlmostEqual(loaded.patches[0].velocity_m_s, 1.4)
            self.assertFalse(loaded.include_dc_copper_losses)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_legacy_v12_config_migrates_without_cfd_profile(self):
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
            json.dump({
                "version": "1.2", "rails": [], "ac_profiles": {},
                "thermal_profile": None,
            }, f)

        try:
            project = load_project_config(filepath)
            self.assertIsNone(project.cfd_profile)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

    def test_legacy_v13_power_tree_uses_safe_thermal_defaults(self):
        """Existing configs gain safe semantics without manual migration."""
        legacy = {
            "version": "1.3",
            "rails": [{
                "net_name": "5V",
                "nominal_voltage": 5.0,
                "sources": [],
                "loads": [{
                    "ref_des": "J6",
                    "total_current": 2.0,
                    "pad_names": ["1"],
                    "distribution_mode": "UNIFORM",
                }],
                "child_regulators": [{
                    "name": "Buck",
                    "input_rail_name": "12V",
                    "input_ref_des": "U4",
                    "input_pad_names": ["VIN"],
                    "output_rail_name": "5V",
                    "output_ref_des": "L1",
                    "output_pad_names": ["1"],
                    "reg_type": "SWITCHING",
                    "efficiency": 0.9,
                }],
            }],
        }
        with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as f:
            json.dump(legacy, f)
            filepath = f.name
        try:
            project = load_project_config(filepath)
            load = project.rails[0].loads[0]
            regulator = project.rails[0].child_regulators[0]
            self.assertEqual(load.thermal_mode, "AUTO")
            self.assertEqual(regulator.thermal_ref_des, "")
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)
    
    def test_file_not_found(self):
        """Test that FileNotFoundError is raised for missing file."""
        with self.assertRaises(FileNotFoundError):
            load_config("/nonexistent/path/config.json")
    
    def test_invalid_version(self):
        """Test that ValueError is raised for unsupported version."""
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json') as f:
            filepath = f.name
            json.dump({"version": "99.0", "rails": []}, f)
        
        try:
            with self.assertRaises(ValueError):
                load_config(filepath)
        finally:
            if Path(filepath).exists():
                os.unlink(filepath)

if __name__ == '__main__':
    unittest.main()
