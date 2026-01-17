try:
    from .plugin import KiPIDA_Plugin
except (ImportError, ValueError):
    import os
    import sys
    # Add current dir to path if needed for flat loading
    plugin_dir = os.path.dirname(os.path.abspath(__file__))
    if plugin_dir not in sys.path:
        sys.path.insert(0, plugin_dir)
    from plugin import KiPIDA_Plugin

def register_plugin():
    try:
        KiPIDA_Plugin().register()
    except Exception as e:
        print(f"KiPIDA registration failed: {e}")

register_plugin()
