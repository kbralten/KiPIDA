import os
import sys
import unittest

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from models import ComponentRef, PowerRail, UnifiedLoad, VoltageRegulator
from thermal_model import PowerLossEstimator


class TestPowerLossEstimator(unittest.TestCase):
    def _load(self, ref_des, current):
        return UnifiedLoad(
            component_ref=ComponentRef(ref_des=ref_des),
            total_current=current,
        )

    def test_switching_regulator_loss(self):
        rail_12v = PowerRail(net_name="12V", nominal_voltage=12.0)
        rail_5v = PowerRail(net_name="5V", nominal_voltage=5.0)
        rail_5v.add_load(self._load("U_LOAD", 1.0))
        rail_12v.add_child_regulator(VoltageRegulator(
            name="BUCK", input_rail_name="12V", input_ref_des="U1",
            input_pad_names=["VIN"],
            output_rail_name="5V", output_ref_des="U1",
            output_pad_names=["VOUT"],
            reg_type="SWITCHING", efficiency=0.8,
        ))

        models = {item.ref_des: item for item in PowerLossEstimator.estimate([rail_12v, rail_5v])}

        self.assertAlmostEqual(models["U_LOAD"].power_w, 5.0)
        self.assertAlmostEqual(models["U1"].power_w, 1.25)

    def test_linear_regulator_loss(self):
        rail_5v = PowerRail(net_name="5V", nominal_voltage=5.0)
        rail_3v3 = PowerRail(net_name="3V3", nominal_voltage=3.3)
        rail_3v3.add_load(self._load("U_LOAD", 0.5))
        rail_5v.add_child_regulator(VoltageRegulator(
            name="LDO", input_rail_name="5V", input_ref_des="U2",
            input_pad_names=["IN"],
            output_rail_name="3V3", output_ref_des="U2",
            output_pad_names=["OUT"],
            reg_type="LINEAR", efficiency=1.0,
        ))

        models = {item.ref_des: item for item in PowerLossEstimator.estimate([rail_5v, rail_3v3])}

        self.assertAlmostEqual(models["U2"].power_w, 0.85)


if __name__ == "__main__":
    unittest.main()
