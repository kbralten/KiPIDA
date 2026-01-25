import wx
import wx.dataview

try:
    from models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef
    from ui.component_selector import ComponentSelectorDialog
    from discovery import NetDiscoverer
except (ImportError, ValueError):
    # Fallback or just re-raise if absolute fails (shouldn't happen with sys.path set)
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
    def __init__(self, parent, board, log_callback=None):
        super().__init__(parent)
        self.board = board
        self.log_callback = log_callback
        self.discoverer = NetDiscoverer(board, log_callback=log_callback)
        
        # Data
        self.rails = [] # List of PowerRail objects
        self.active_rail = None
        
        self._init_ui()

    def log(self, msg):
        if self.log_callback:
            self.log_callback(msg)
        else:
            print(msg)
        
    def _init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Auto-scan will happen after UI is built
        
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
        self.rail_list.Bind(wx.EVT_LISTBOX, self.on_rail_select)
        self.txt_voltage.Bind(wx.EVT_TEXT, self.on_voltage_change)
        
        self.btn_add_src.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("SOURCE"))
        self.btn_add_load.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("LOAD"))
        self.btn_del_comp.Bind(wx.EVT_BUTTON, self.on_del_component)
        self.comp_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_edit_component)

    def auto_scan(self):
        """Auto-scan board on startup"""
        self.rails = self.discoverer.discover_power_nets()
        self.log(f"Scan complete. Discovered {len(self.rails)} power rail candidates.")
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
        # Use discoverer to get nets via IPC
        all_nets = self.discoverer.get_all_net_names()
        self.log(f"Found {len(all_nets)} total nets on board.")
        candidate_nets = [n for n in all_nets if n and n not in existing_nets]
        
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
        all_comps = self.discoverer.get_components_on_net(self.active_rail.net_name)
        if not all_comps:
            wx.MessageBox("No components found on this net.")
            return
        
        # Filter out components where all pads are already assigned
        # Build set of all used pads (ref-pad format)
        used_pads = set()
        for src in self.active_rail.sources:
            for pad in src.pad_names:
                used_pads.add(f"{src.component_ref.ref_des}-{pad}")
        for load in self.active_rail.loads:
            for pad in load.pad_names:
                used_pads.add(f"{load.component_ref.ref_des}-{pad}")
        
        # Filter components
        available_comps = {}
        for ref, pads in all_comps.items():
            # Check if any pad is not used
            has_available_pad = False
            for pad in pads:
                pad_name = getattr(pad, 'number', getattr(pad, 'name', ''))
                if f"{ref}-{pad_name}" not in used_pads:
                    has_available_pad = True
                    break
            
            if has_available_pad:
                available_comps[ref] = pads
        
        if not available_comps:
            wx.MessageBox("All components on this net have been assigned.")
            return

        dlg = ComponentSelectorDialog(self, "Select Component", self.active_rail.net_name, available_comps)
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

    def on_edit_component(self, event):
        sel = self.comp_list.GetFirstSelected()
        if sel == -1 or not self.active_rail: return
        
        n_src = len(self.active_rail.sources)
        comp_obj = None
        mode = ""
        if sel < n_src:
            comp_obj = self.active_rail.sources[sel]
            mode = "SOURCE"
        else:
            comp_obj = self.active_rail.loads[sel - n_src]
            mode = "LOAD"
            
        # Get components on this net
        all_comps = self.discoverer.get_components_on_net(self.active_rail.net_name)
        if not all_comps: return
        
        # Filter out components where all pads are already assigned, 
        # but KEEP the pads of the component being edited.
        used_pads = set()
        for src in self.active_rail.sources:
            if src is not comp_obj:
                for pad in src.pad_names:
                    used_pads.add(f"{src.component_ref.ref_des}-{pad}")
        for load in self.active_rail.loads:
            if load is not comp_obj:
                for pad in load.pad_names:
                    used_pads.add(f"{load.component_ref.ref_des}-{pad}")
                    
        # Filter components
        available_comps = {}
        for ref, pads in all_comps.items():
            has_available_pad = False
            for pad in pads:
                pad_name = getattr(pad, 'number', getattr(pad, 'name', ''))
                if f"{ref}-{pad_name}" not in used_pads:
                    has_available_pad = True
                    break
            if has_available_pad:
                available_comps[ref] = pads

        dlg = ComponentSelectorDialog(self, "Edit Component", self.active_rail.net_name, available_comps)
        dlg.set_mode(mode)
        
        # Prepopulate
        val = 0.0
        if mode == "LOAD":
            val = comp_obj.total_current
        dlg.prepopulate(comp_obj.component_ref.ref_des, val, comp_obj.pad_names)
        
        if dlg.ShowModal() == wx.ID_OK:
            ref_des, val, pads = dlg.GetSelection()
            comp_obj.component_ref.ref_des = ref_des
            comp_obj.pad_names = pads
            if mode == "LOAD":
                comp_obj.total_current = val
                
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
