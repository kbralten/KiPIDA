import os
import sys
import unittest
from types import SimpleNamespace

plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from ac_model import ACModelBuilder, format_capacitance, parse_capacitance


def pad(number, net_name):
    return SimpleNamespace(number=number, net=SimpleNamespace(name=net_name))


class TestCapacitorDiscovery(unittest.TestCase):
    def test_parses_common_kicad_values(self):
        self.assertAlmostEqual(parse_capacitance("100n"), 100e-9)
        self.assertAlmostEqual(parse_capacitance("4u7"), 4.7e-6)
        self.assertAlmostEqual(parse_capacitance("10uF 16V"), 10e-6)
        self.assertAlmostEqual(parse_capacitance("1.5pF"), 1.5e-12)
        self.assertIsNone(parse_capacitance("DNP"))

    def test_formats_engineering_units(self):
        self.assertEqual(format_capacitance(4.7e-6), "4.7uF")
        self.assertEqual(format_capacitance(100e-9), "100nF")

    def test_discovers_populated_and_dnp_rail_capacitors(self):
        populated = SimpleNamespace(
            reference="C1", value="100n", footprint_name="Capacitor_SMD:C_0402",
            pads=[pad("1", "3V3"), pad("2", "GND")], dnp=False,
        )
        candidate = SimpleNamespace(
            reference="C2", value="DNP", footprint_name="Capacitor_SMD:C_0805",
            pads=[pad("1", "3V3"), pad("2", "GND")], dnp=True,
        )
        unrelated = SimpleNamespace(
            reference="C3", value="1u", footprint_name="Capacitor_SMD:C_0603",
            pads=[pad("1", "5V"), pad("2", "GND")], dnp=False,
        )
        builder = ACModelBuilder(SimpleNamespace(footprints=[populated, candidate, unrelated]))

        capacitors = builder.discover_capacitors("3V3", "GND")

        self.assertEqual([cap.ref_des for cap in capacitors], ["C1", "C2"])
        self.assertTrue(capacitors[0].enabled)
        self.assertFalse(capacitors[0].candidate)
        self.assertFalse(capacitors[1].enabled)
        self.assertTrue(capacitors[1].candidate)
        self.assertLess(capacitors[0].esl_h, capacitors[1].esl_h)
        self.assertEqual(builder.discover_ground_nets(), ["GND"])


if __name__ == "__main__":
    unittest.main()
