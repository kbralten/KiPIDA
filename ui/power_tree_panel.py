import wx
import wx.dataview

try:
    from models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef, VoltageRegulator
    from ui.component_selector import ComponentSelectorDialog
    from ui.regulator_dialog import RegulatorDialog
    from discovery import NetDiscoverer
    from config_manager import get_project_config_path, save_config, load_config
except (ImportError, ValueError):
    # Fallback or just re-raise if absolute fails (shouldn't happen with sys.path set)
    from models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef, VoltageRegulator
    from ui.component_selector import ComponentSelectorDialog
    from ui.regulator_dialog import RegulatorDialog
    from discovery import NetDiscoverer
    from config_manager import get_project_config_path, save_config, load_config

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
    def __init__(self, parent, board, project=None, log_callback=None):
        super().__init__(parent)
        self.board = board
        self.project = project
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
        
        # Config buttons
        config_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_save_config = wx.Button(left_panel, label="Save Config")
        self.btn_load_config = wx.Button(left_panel, label="Load Config")
        config_sizer.Add(self.btn_save_config, 1, wx.EXPAND | wx.RIGHT, 2)
        config_sizer.Add(self.btn_load_config, 1, wx.EXPAND)
        left_sizer.Add(config_sizer, 0, wx.EXPAND | wx.ALL, 2)
        
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
        self.btn_add_reg = wx.Button(self.detail_panel, label="+ Regulator")
        self.btn_del_comp = wx.Button(self.detail_panel, label="- Remove")
        
        act_sizer.Add(self.btn_add_src, 0, wx.RIGHT, 5)
        act_sizer.Add(self.btn_add_load, 0, wx.RIGHT, 5)
        act_sizer.Add(self.btn_add_reg, 0, wx.RIGHT, 5)
        act_sizer.Add(self.btn_del_comp, 0, wx.RIGHT, 5)
        
        self.detail_sizer.Add(act_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        self.detail_panel.SetSizer(self.detail_sizer)
        self.detail_panel.Disable() # Disable until rail selected
        
        splitter.SplitVertically(left_panel, self.detail_panel, 200)
        main_sizer.Add(splitter, 1, wx.EXPAND | wx.ALL, 5)
        
        self.SetSizer(main_sizer)
        
        self.btn_add_rail.Bind(wx.EVT_BUTTON, self.on_add_rail)
        self.btn_save_config.Bind(wx.EVT_BUTTON, self.on_save_config)
        self.btn_load_config.Bind(wx.EVT_BUTTON, self.on_load_config)
        self.rail_list.Bind(wx.EVT_LISTBOX, self.on_rail_select)
        self.txt_voltage.Bind(wx.EVT_TEXT, self.on_voltage_change)
        
        self.btn_add_src.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("SOURCE"))
        self.btn_add_load.Bind(wx.EVT_BUTTON, lambda e: self.on_add_component("LOAD"))
        self.btn_add_reg.Bind(wx.EVT_BUTTON, self.on_add_regulator)
        self.btn_del_comp.Bind(wx.EVT_BUTTON, self.on_del_component)
        self.comp_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_edit_component)

    def auto_scan(self):
        """Auto-scan board on startup or load config if it exists."""
        # Try to load config first
        config_path = self._get_config_path()
        if config_path and config_path.exists():
            try:
                self.rails = load_config(str(config_path))
                self.log(f"Loaded configuration from {config_path.name} ({len(self.rails)} rails)")
            except Exception as e:
                self.log(f"Failed to load config: {e}. Running auto-scan instead.")
                self.rails = self.discoverer.discover_power_nets()
                self.log(f"Scan complete. Discovered {len(self.rails)} power rail candidates.")
        else:
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

        # Regulators (Outputting to another rail)
        for r in self.active_rail.child_regulators:
            idx = self.comp_list.InsertItem(self.comp_list.GetItemCount(), f"REG -> {r.output_rail_name}")
            self.comp_list.SetItem(idx, 1, r.name)
            self.comp_list.SetItem(idx, 2, "---") # V determined by output rail
            self.comp_list.SetItem(idx, 3, r.reg_type)

        # Incoming Regulators (Sources from another rail)
        for r in self.rails:
            if r == self.active_rail: continue
            for reg in r.child_regulators:
                if reg.output_rail_name == self.active_rail.net_name:
                    idx = self.comp_list.InsertItem(self.comp_list.GetItemCount(), f"REG FROM {r.net_name}")
                    self.comp_list.SetItem(idx, 1, reg.name)
                    self.comp_list.SetItem(idx, 2, "SOURCE") 
                    self.comp_list.SetItem(idx, 3, reg.reg_type)

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
        """Handle double-click on component list item."""
        sel = self.comp_list.GetFirstSelected()
        if sel == -1 or not self.active_rail:
            return
        
        n_src = len(self.active_rail.sources)
        n_load = len(self.active_rail.loads)
        n_reg = len(self.active_rail.child_regulators)
        
        # Determine which type of component was clicked
        if sel < n_src:
            # Edit Source
            self._edit_source(sel)
        elif sel < n_src + n_load:
            # Edit Load
            self._edit_load(sel - n_src)
        elif sel < n_src + n_load + n_reg:
            # Edit outgoing Regulator
            self._edit_regulator(sel - n_src - n_load)
        else:
            # Incoming regulator - find it
            incoming_idx = sel - n_src - n_load - n_reg
            self._edit_incoming_regulator(incoming_idx)
    
    def _edit_source(self, idx):
        """Edit a source component."""
        comp_obj = self.active_rail.sources[idx]
        
        # Get components on this net
        all_comps = self.discoverer.get_components_on_net(self.active_rail.net_name)
        if not all_comps:
            wx.MessageBox("No components found on this net.")
            return
        
        # Build set of used pads (excluding current component)
        used_pads = set()
        for i, src in enumerate(self.active_rail.sources):
            if i != idx:  # Skip current component
                for pad in src.pad_names:
                    used_pads.add(f"{src.component_ref.ref_des}-{pad}")
        for load in self.active_rail.loads:
            for pad in load.pad_names:
                used_pads.add(f"{load.component_ref.ref_des}-{pad}")
        
        # Filter available components
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
        
        dlg = ComponentSelectorDialog(self, "Edit Source", self.active_rail.net_name, available_comps)
        dlg.set_mode("SOURCE")
        dlg.prepopulate(comp_obj.component_ref.ref_des, 0.0, comp_obj.pad_names)
        
        if dlg.ShowModal() == wx.ID_OK:
            ref_des, val, pads = dlg.GetSelection()
            comp_obj.component_ref.ref_des = ref_des
            comp_obj.pad_names = pads
            self.refresh_comp_list()
        
        dlg.Destroy()
    
    def _edit_load(self, idx):
        """Edit a load component."""
        comp_obj = self.active_rail.loads[idx]
        
        # Get components on this net
        all_comps = self.discoverer.get_components_on_net(self.active_rail.net_name)
        if not all_comps:
            wx.MessageBox("No components found on this net.")
            return
        
        # Build set of used pads (excluding current component)
        used_pads = set()
        for src in self.active_rail.sources:
            for pad in src.pad_names:
                used_pads.add(f"{src.component_ref.ref_des}-{pad}")
        for i, load in enumerate(self.active_rail.loads):
            if i != idx:  # Skip current component
                for pad in load.pad_names:
                    used_pads.add(f"{load.component_ref.ref_des}-{pad}")
        
        # Filter available components
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
        
        dlg = ComponentSelectorDialog(self, "Edit Load", self.active_rail.net_name, available_comps)
        dlg.set_mode("LOAD")
        dlg.prepopulate(comp_obj.component_ref.ref_des, comp_obj.total_current, comp_obj.pad_names)
        
        if dlg.ShowModal() == wx.ID_OK:
            ref_des, val, pads = dlg.GetSelection()
            comp_obj.component_ref.ref_des = ref_des
            comp_obj.pad_names = pads
            comp_obj.total_current = val
            self.refresh_comp_list()
        
        dlg.Destroy()
    
    def _edit_regulator(self, idx):
        """Edit an outgoing regulator."""
        reg = self.active_rail.child_regulators[idx]
        
        rail_names = [r.net_name for r in self.rails]
        dlg = RegulatorDialog(
            self, 
            "Edit Regulator", 
            rail_names, 
            self.discoverer,
            input_rail=reg.input_rail_name,
            output_rail=reg.output_rail_name
        )
        
        # Prepopulate the dialog with existing regulator data
        dlg.prepopulate(reg)
        
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetValue()
            
            # Update regulator object
            reg.name = data['name']
            reg.input_rail_name = data['input_rail']
            reg.input_ref_des = data['input_ref_des']
            reg.input_pad_names = data['input_pads']
            reg.output_rail_name = data['output_rail']
            reg.output_ref_des = data['output_ref_des']
            reg.output_pad_names = data['output_pads']
            reg.reg_type = data['type']
            reg.efficiency = data['efficiency']
            
            # If input rail changed, move regulator to new rail
            if data['input_rail'] != self.active_rail.net_name:
                self.active_rail.child_regulators.pop(idx)
                target_rail = next((r for r in self.rails if r.net_name == data['input_rail']), None)
                if target_rail:
                    target_rail.add_child_regulator(reg)
                    self.log(f"Regulator moved to {target_rail.net_name}")
            
            self.refresh_comp_list()
        
        dlg.Destroy()
    
    def _edit_incoming_regulator(self, incoming_idx):
        """Edit an incoming regulator (one that outputs to this rail)."""
        # Find the incoming regulator
        count = 0
        source_rail = None
        reg = None
        
        for r in self.rails:
            if r == self.active_rail:
                continue
            for candidate_reg in r.child_regulators:
                if candidate_reg.output_rail_name == self.active_rail.net_name:
                    if count == incoming_idx:
                        source_rail = r
                        reg = candidate_reg
                        break
                    count += 1
            if reg:
                break
        
        if not reg or not source_rail:
            wx.MessageBox("Could not locate regulator.")
            return
        
        # Edit it (same as _edit_regulator but from source_rail context)
        rail_names = [r.net_name for r in self.rails]
        dlg = RegulatorDialog(
            self,
            "Edit Regulator",
            rail_names,
            self.discoverer,
            input_rail=reg.input_rail_name,
            output_rail=reg.output_rail_name
        )
        
        dlg.prepopulate(reg)
        
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetValue()
            
            # Update regulator object
            reg.name = data['name']
            reg.input_rail_name = data['input_rail']
            reg.input_ref_des = data['input_ref_des']
            reg.input_pad_names = data['input_pads']
            reg.output_rail_name = data['output_rail']
            reg.output_ref_des = data['output_ref_des']
            reg.output_pad_names = data['output_pads']
            reg.reg_type = data['type']
            reg.efficiency = data['efficiency']
            
            # If input rail changed, move regulator
            if data['input_rail'] != source_rail.net_name:
                source_rail.child_regulators.remove(reg)
                target_rail = next((r for r in self.rails if r.net_name == data['input_rail']), None)
                if target_rail:
                    target_rail.add_child_regulator(reg)
                    self.log(f"Regulator moved to {target_rail.net_name}")
            
            self.refresh_comp_list()
        
        dlg.Destroy()

    def on_add_regulator(self, event):
        if not self.active_rail: return
        
        rail_names = [r.net_name for r in self.rails]
        # Pass discoverer to dialog
        dlg = RegulatorDialog(self, "Add Regulator", rail_names, self.discoverer, input_rail=self.active_rail.net_name)
        
        if dlg.ShowModal() == wx.ID_OK:
            data = dlg.GetValue()
            
            # Create Regulator Object
            reg = VoltageRegulator(
                name=data['name'],
                input_rail_name=data['input_rail'],
                input_ref_des=data['input_ref_des'],
                input_pad_names=data['input_pads'],
                output_rail_name=data['output_rail'],
                output_ref_des=data['output_ref_des'],
                output_pad_names=data['output_pads'],
                reg_type=data['type'],
                efficiency=data['efficiency']
            )
            
            # Add to the INPUT rail (which is self.active_rail)
            # Actually, the user might have changed input_rail in the dialog!
            # We need to find the correct input rail object.
            
            target_rail = None
            for r in self.rails:
                if r.net_name == data['input_rail']:
                    target_rail = r
                    break
            
            if target_rail:
                target_rail.add_child_regulator(reg)
                self.refresh_comp_list() # Only refreshes active rail list
                
                # If we added to a different rail than active, maybe switch? or just log?
                if target_rail != self.active_rail:
                    self.log(f"Regulator added to {target_rail.net_name}")
            else:
                wx.MessageBox("Error locating input rail.")
            
        dlg.Destroy()

    def on_del_component(self, event):
        sel = self.comp_list.GetFirstSelected()
        if sel == -1 or not self.active_rail: return
        
        n_src = len(self.active_rail.sources)
        n_load = len(self.active_rail.loads)
        n_reg = len(self.active_rail.child_regulators)
        
        if sel < n_src:
            self.active_rail.sources.pop(sel)
        elif sel < n_src + n_load:
            self.active_rail.loads.pop(sel - n_src)
        elif sel < n_src + n_load + n_reg:
            self.active_rail.child_regulators.pop(sel - n_src - n_load)
            
        self.refresh_comp_list()
    
    def _get_config_path(self):
        """Get the path to the config file based on project name."""
        try:
            from pathlib import Path
            
            # Use project object if available (IPC API)
            if self.project and hasattr(self.project, 'path'):
                return get_project_config_path(
                    self.project.path,
                    getattr(self.project, 'name', None),
                )
            
            # Fallback: try to get from board (legacy)
            if hasattr(self.board, 'filename'):
                project_file = self.board.filename
            elif hasattr(self.board, 'get_filename'):
                project_file = self.board.get_filename()
            else:
                self.log("ERROR: Cannot determine project path - no project object and board has no filename attribute.")
                return None
            
            if not project_file:
                return None
            
            # Board filenames are stored beside the generated configuration.
            project_path = Path(project_file)
            config_filename = f"{project_path.stem}.kipida.json"
            config_path = project_path.parent / config_filename
            return config_path
        except Exception as e:
            self.log(f"Error getting config path: {e}")
            return None
    
    def on_save_config(self, event):
        """Save current configuration to JSON file."""
        config_path = self._get_config_path()
        if not config_path:
            wx.MessageBox("Unable to determine project path. Cannot save config.", "Error", wx.OK | wx.ICON_ERROR)
            return
        
        try:
            save_config(self.rails, str(config_path))
            self.log(f"Configuration saved to {config_path.name}")
            wx.MessageBox(f"Configuration saved successfully to:\n{config_path}", "Success", wx.OK | wx.ICON_INFORMATION)
        except Exception as e:
            self.log(f"Failed to save config: {e}")
            wx.MessageBox(f"Failed to save configuration:\n{e}", "Error", wx.OK | wx.ICON_ERROR)
    
    def on_load_config(self, event):
        """Load configuration from JSON file."""
        config_path = self._get_config_path()
        if not config_path:
            wx.MessageBox("Unable to determine project path. Cannot load config.", "Error", wx.OK | wx.ICON_ERROR)
            return
        
        if not config_path.exists():
            wx.MessageBox(f"Config file not found:\n{config_path}", "Error", wx.OK | wx.ICON_ERROR)
            return
        
        try:
            self.rails = load_config(str(config_path))
            self.log(f"Loaded configuration from {config_path.name} ({len(self.rails)} rails)")
            
            # Refresh UI
            self.rail_list.Clear()
            for r in self.rails:
                lbl = f"{r.net_name}"
                if r.nominal_voltage:
                    lbl += f" ({r.nominal_voltage}V)"
                self.rail_list.Append(lbl)
            
            if self.rails:
                self.rail_list.SetSelection(0)
                self.on_rail_select(None)
            
            wx.MessageBox(f"Configuration loaded successfully:\n{len(self.rails)} rails", "Success", wx.OK | wx.ICON_INFORMATION)
        except Exception as e:
            self.log(f"Failed to load config: {e}")
            wx.MessageBox(f"Failed to load configuration:\n{e}", "Error", wx.OK | wx.ICON_ERROR)
