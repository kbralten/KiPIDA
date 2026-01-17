from .plugin import KiPIDA_Plugin

def register_plugin():
    KiPIDA_Plugin().register()

register_plugin()
