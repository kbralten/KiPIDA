import wx
import wx.dataview

try:
    from ..models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef
    from .component_selector import ComponentSelectorDialog
    from ..discovery import NetDiscoverer
except (ImportError, ValueError):
    from models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef
    from ui.component_selector import ComponentSelectorDialog
    from discovery import NetDiscoverer

class NetSelectionDialog(wx.Dialog):
    def __init__(self, parent, nets):
        super().__init__(parent, title="Select Net", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.nets = sorted(nets)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.lb = wx.ListBox(self, choices=self.nets, style=wx.LB_SINGLE | wx.LB_NEEDED_SB)
        sizer.Add(self.lb, 1, wx.EXPAND | wx.ALL, 10)
        
        btns = self.CreateButtonSizer(wx.OK | wx.CANCEL)
        sizer.Add(btns, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(sizer)
        self.SetSize((300, 400))

    def GetSelectedNet(self):
        sel = self.lb.GetSelection()
        if sel == -1: return None
        return self.lb.GetString(sel)

class PowerTreePanel(wx.Panel):
    def __init__(self, parent, board):
        super().__init__(parent)
        self.board = board
        self.discoverer = NetDiscoverer(board)
        
        # Data
        self.rails = [] # List of PowerRail objects
        self.active_rail = None
        
        self._init_ui()
        
    def _init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Toolbar / Top Buttons
        top_hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_scan = wx.Button(self, label="Scan Board")
        top_hbox.Add(self.btn_scan, 0, wx.ALL, 5)
        self.btn_refresh = wx.Button(self, label="Refresh Net")
        top_hbox.Add(self.btn_refresh, 0, wx.ALL, 5)
        main_sizer.Add(top_hbox, 0, wx.EXPAND)
        
        # Splitter: Rail List (Left) vs Rail Details (Right)
        # Actually, let's use a master-detail approach. 
        # ListBox of Rails -> Details Panel.
        
        splitter = wx.SplitterWindow(self)
        
        # LEFT: Rail List
        left_panel = wx.Panel(splitter)
        left_sizer = wx.BoxSizer(wx.VERTICAL)
        left_sizer.Add(wx.StaticText(left_panel, label="Detected Power Rails"), 0, wx.ALL, 2)
        
        self.rail_list = wx.ListBox(left_panel, style=wx.LB_SINGLE)
        left_sizer.Add(self.rail_list, 1, wx.EXPAND | wx.ALL, 2)
        
        self.btn_add_rail = wx.Button(left_panel, label="Add Manual Net")
        left_sizer.Add(self.btn_add_rail, 0, wx.EXPAND | wx.ALL, 2)
        
        left_panel.SetSizer(left_sizer)
        
        # RIGHT: Details
        self.detail_panel = wx.Panel(splitter)
        self.detail_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Rail Properties
        prop_sizer = wx.StaticBoxSizer(wx.VERTICAL, self.detail_panel, "Rail Properties")
        
        h_volt = wx.BoxSizer(wx.HORIZONTAL)
        h_volt.Add(wx.StaticText(prop_sizer.GetStaticBox(), label="Nominal Voltage (V):"), 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
        self.txt_voltage = wx.TextCtrl(prop_sizer.GetStaticBox())
        h_volt.Add(self.txt_voltage, 1, wx.EXPAND)
        
        prop_sizer.Add(h_volt, 0, wx.EXPAND | wx.ALL, 5)
        self.detail_sizer.Add(prop_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        # Components List (Sources/Loads)
        # We'll use a ListView here
        self.comp_list = wx.ListCtrl(self.detail_panel, style=wx.LC_REPORT)
        self.comp_list.InsertColumn(0, "Role", width=80)
        self.comp_list.InsertColumn(1, "Ref", width=60)
        self.comp_list.InsertColumn(2, "Value", width=80)
        self.comp_list.InsertColumn(3, "Pads", width=100)
        
        self.detail_sizer.Add(self.comp_list, 1, wx.EXPAND | wx.ALL, 5)
        
        # Component Actions
        act_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_add_src = wx.Button(self.detail_panel, label="+ Source")
        self.btn_add_load = wx.Button(self.detail_panel, label="+ Load")
        self.btn_del_comp = wx.Button(self.detail_panel, label="- Remove")
        
        act_sizer.Add(self.btn_add_src, 0, wx.RIGHT, 5)
        act_sizer.Add(self.btn_add_load, 0, wx.RIGHT, 5)
        act_sizer.Add(self.btn_del_comp, 0, wx.RIGHT, 5)
        
        self.detail_sizer.Add(act_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        self.detail_panel.SetSizer(self.detail_sizer)
        self.detail_panel.Disable() # Disable until rail selected
        
        splitter.SplitVertically(left_panel, self.detail_panel, 200)
        main_sizer.Add(splitter, 1, wx.EXPAND | wx.ALL, 5)
        
        self.SetSizer(main_sizer)
        
        self.btn_add_rail.Bind(wx.EVT_BUTTON, self.on_add_rail)
        self.btn_scan.Bind(wx.EVT_BUTTON, self.on_scan)
        self.rail_list.Bind(wx.EVT_LISTBOX, self.on_rail_select)
        self.txt_voltage.Bind(wx.EVT_TEXT, self.on_voltage_change)
        
        self.btn_add_src.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("SOURCE"))
        self.btn_add_load.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("LOAD"))
        self.btn_del_comp.Bind(wx.EVT_BUTTON, self.on_del_component)

    def on_scan(self, event):
        self.rails = self.discoverer.discover_power_nets()
        self.rail_list.Clear()
        for r in self.rails:
            lbl = f"{r.net_name}"
            if r.nominal_voltage:
                lbl += f" ({r.nominal_voltage}V)"
            self.rail_list.Append(lbl)
            
        if self.rails:
            self.rail_list.SetSelection(0)
            self.on_rail_select(None)

    def on_add_rail(self, event):
        # Find all nets not already in self.rails
        existing_nets = {r.net_name for r in self.rails}
        all_nets = self.board.GetNetsByName()
        candidate_nets = [str(n) for n in all_nets.keys() if str(n) and str(n) not in existing_nets]
        
        dlg = NetSelectionDialog(self, candidate_nets)
        if dlg.ShowModal() == wx.ID_OK:
            net_name = dlg.GetSelectedNet()
            if net_name:
                rail = PowerRail(net_name=net_name)
                # Try to estimate voltage
                v = self.discoverer._estimate_voltage(net_name)
                if v: rail.nominal_voltage = v
                
                self.rails.append(rail)
                lbl = f"{rail.net_name}"
                if rail.nominal_voltage:
                    lbl += f" ({rail.nominal_voltage}V)"
                self.rail_list.Append(lbl)
                self.rail_list.SetSelection(self.rail_list.GetCount() - 1)
                self.on_rail_select(None)
        dlg.Destroy()

    def on_rail_select(self, event):
        idx = self.rail_list.GetSelection()
        if idx == -1: 
            self.active_rail = None
            self.detail_panel.Disable()
            return
            
        self.active_rail = self.rails[idx]
        self.detail_panel.Enable()
        
        # Populate Details
        self.txt_voltage.ChangeValue(str(self.active_rail.nominal_voltage) if self.active_rail.nominal_voltage else "0.0")
        
        self.refresh_comp_list()
        
    def refresh_comp_list(self):
        self.comp_list.DeleteAllItems()
        if not self.active_rail: return
        
        # Sources
        for s in self.active_rail.sources:
            idx = self.comp_list.InsertItem(self.comp_list.GetItemCount(), "SOURCE")
            self.comp_list.SetItem(idx, 1, s.component_ref.ref_des)
            self.comp_list.SetItem(idx, 2, "---") # Voltage is rail level
            self.comp_list.SetItem(idx, 3, str(len(s.pad_names)))
            
        # Loads
        for l in self.active_rail.loads:
            idx = self.comp_list.InsertItem(self.comp_list.GetItemCount(), "LOAD")
            self.comp_list.SetItem(idx, 1, l.component_ref.ref_des)
            self.comp_list.SetItem(idx, 2, f"{l.total_current} A")
            self.comp_list.SetItem(idx, 3, str(len(l.pad_names)))

    def on_voltage_change(self, event):
        if not self.active_rail: return
        try:
            val = float(self.txt_voltage.GetValue())
            self.active_rail.nominal_voltage = val
        except: pass

    def on_add_component(self, mode):
        if not self.active_rail: return
        
        # Get components on this net
        comps = self.discoverer.get_components_on_net(self.active_rail.net_name)
        if not comps:
            wx.MessageBox("No components found on this net.")
            return

        dlg = ComponentSelectorDialog(self, "Select Component", self.active_rail.net_name, comps)
        dlg.set_mode(mode)
        
        if dlg.ShowModal() == wx.ID_OK:
            ref_des, val, pads = dlg.GetSelection()
            ref = ComponentRef(ref_des)
            
            if mode == "SOURCE":
                # Remove check if exists?
                s = UnifiedSource(ref, pads)
                self.active_rail.add_source(s)
            else:
                l = UnifiedLoad(ref, val, pads)
                self.active_rail.add_load(l)
                
            self.refresh_comp_list()
            
        dlg.Destroy()

    def on_del_component(self, event):
        sel = self.comp_list.GetFirstSelected()
        if sel == -1 or not self.active_rail: return
        
        # This is tricky because list is flattened sources + loads
        # Better strategy: keep a parallel list of objects or verify index
        # For MVP, simpler: Sources are first, then Loads
        
        n_src = len(self.active_rail.sources)
        if sel < n_src:
            self.active_rail.sources.pop(sel)
        else:
            self.active_rail.loads.pop(sel - n_src)
            
        self.refresh_comp_list()
