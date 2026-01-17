
import unittest
import sys
import os

# Ensure the plugin root is in the path
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

class MockPcbnew:
    class ActionPlugin: pass

if 'pcbnew' not in sys.modules:
    sys.modules['pcbnew'] = MockPcbnew

class TestPluginStructure(unittest.TestCase):
    def test_imports(self):
        """test that the main modules can be imported"""
        try:
            import extractor
            import mesh
            import solver
            import plugin
        except ImportError as e:
            self.fail(f"Failed to import Core Modules: {e}")
        
    def test_ui_import(self):
        """test that the UI modules can be imported (requires wx)"""
        try:
            import wx
        except ImportError:
            self.skipTest("wxPython not available")

        try:
            from ui import main_dialog
        except ImportError as e:
            self.fail(f"Failed to import UI Modules: {e}")

if __name__ == '__main__':
    unittest.main()
