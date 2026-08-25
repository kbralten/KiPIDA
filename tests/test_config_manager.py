import unittest
import json
import tempfile
import os
from pathlib import Path

import sys
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from config_manager import get_project_config_path, save_config, load_config
from models import PowerRail, UnifiedSource, UnifiedLoad, VoltageRegulator, ComponentRef

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
            distribution_mode="UNIFORM"
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
            efficiency=0.90
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
            
            self.assertEqual(data["version"], "1.0")
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
            
            # Verify regulator
            reg = rail_12v.child_regulators[0]
            self.assertEqual(reg.name, "Buck1")
            self.assertEqual(reg.input_rail_name, "12V")
            self.assertEqual(reg.output_rail_name, "5V")
            self.assertEqual(reg.reg_type, "SWITCHING")
            self.assertEqual(reg.efficiency, 0.90)
            
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
