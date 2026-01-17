import pcbnew
import os
import wx
try:
    from .ui.main_dialog import KiPIDA_MainDialog
except (ImportError, ValueError):
    from ui.main_dialog import KiPIDA_MainDialog

class KiPIDA_Plugin(pcbnew.ActionPlugin):
    def defaults(self):
        self.name = "Ki-PIDA"
        self.category = "Simulation"
        self.description = "Power Integrity & Delivery Analyzer for KiCad"
        self.show_toolbar_button = True
        self.icon_file_name = os.path.join(os.path.dirname(__file__), 'resources', 'kipida.png')


    def Run(self):
        # The entry function of the plugin that is executed on user action
        board = pcbnew.GetBoard()
        dlg = KiPIDA_MainDialog(None)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()
