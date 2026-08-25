import wx
try:
    from models import generate_regulator_name
except ImportError:
    from ..models import generate_regulator_name

class RegulatorDialog(wx.Dialog):
    def __init__(self, parent, title, available_rails, discoverer, input_rail=None, output_rail=None):
        super(RegulatorDialog, self).__init__(parent, title=title, size=(500, 650))
        
        self.available_rails = sorted(available_rails)
        self.discoverer = discoverer
        self.input_rail = input_rail
        self.output_rail = output_rail
        
        self.input_comps = {} # ref -> [pads]
        self.output_comps = {} # ref -> [pads]
        self._last_input_ref = ""
        
        self._init_ui()
        self.Centre()
        
        # Trigger initial population
        if self.input_rail: self._on_input_rail_change(None)
        if self.output_rail: self._on_output_rail_change(None)
        
    def _init_ui(self):
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        # Name
        hbox_name = wx.BoxSizer(wx.HORIZONTAL)
        hbox_name.Add(wx.StaticText(self, label="Regulator Name:"), 0, wx.ALIGN_CENTER_VERTICAL|wx.RIGHT, 5)
        self.txt_name = wx.TextCtrl(self, style=wx.TE_READONLY)
        hbox_name.Add(self.txt_name, 1, wx.EXPAND)
        vbox.Add(hbox_name, 0, wx.EXPAND | wx.ALL, 10)
        
        # --- Input Section ---
        sb_in = wx.StaticBoxSizer(wx.VERTICAL, self, "Input Side")
        
        # Input Rail
        grid_in = wx.FlexGridSizer(3, 2, 5, 5)
        grid_in.Add(wx.StaticText(sb_in.GetStaticBox(), label="Input Rail:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_input_rail = wx.ComboBox(sb_in.GetStaticBox(), choices=self.available_rails, style=wx.CB_READONLY)
        if self.input_rail: self.cmb_input_rail.SetValue(self.input_rail)
        self.cmb_input_rail.Bind(wx.EVT_COMBOBOX, self._on_input_rail_change)
        grid_in.Add(self.cmb_input_rail, 1, wx.EXPAND)
        
        # Input Component
        grid_in.Add(wx.StaticText(sb_in.GetStaticBox(), label="Input Component:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_input_comp = wx.ComboBox(sb_in.GetStaticBox(), style=wx.CB_READONLY)
        self.cmb_input_comp.Bind(wx.EVT_COMBOBOX, self._on_input_comp_change)
        grid_in.Add(self.cmb_input_comp, 1, wx.EXPAND)
        
        # Input Pads
        grid_in.Add(wx.StaticText(sb_in.GetStaticBox(), label="Input Pads:"), 0, wx.ALIGN_TOP)
        self.lst_input_pads = wx.CheckListBox(sb_in.GetStaticBox(), size=(-1, 80))
        grid_in.Add(self.lst_input_pads, 1, wx.EXPAND)
        
        grid_in.AddGrowableCol(1, 1)
        sb_in.Add(grid_in, 1, wx.EXPAND | wx.ALL, 5)
        vbox.Add(sb_in, 0, wx.EXPAND | wx.ALL, 5)
        
        # --- Output Section ---
        sb_out = wx.StaticBoxSizer(wx.VERTICAL, self, "Output Side")
        
        # Output Rail
        grid_out = wx.FlexGridSizer(3, 2, 5, 5)
        grid_out.Add(wx.StaticText(sb_out.GetStaticBox(), label="Output Rail:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_output_rail = wx.ComboBox(sb_out.GetStaticBox(), choices=self.available_rails, style=wx.CB_READONLY)
        if self.output_rail: self.cmb_output_rail.SetValue(self.output_rail)
        self.cmb_output_rail.Bind(wx.EVT_COMBOBOX, self._on_output_rail_change)
        grid_out.Add(self.cmb_output_rail, 1, wx.EXPAND)
        
        # Output Component
        grid_out.Add(wx.StaticText(sb_out.GetStaticBox(), label="Output Component:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_output_comp = wx.ComboBox(sb_out.GetStaticBox(), style=wx.CB_READONLY)
        self.cmb_output_comp.Bind(wx.EVT_COMBOBOX, self._on_output_comp_change)
        grid_out.Add(self.cmb_output_comp, 1, wx.EXPAND)
        
        # Output Pads
        grid_out.Add(wx.StaticText(sb_out.GetStaticBox(), label="Output Pads:"), 0, wx.ALIGN_TOP)
        self.lst_output_pads = wx.CheckListBox(sb_out.GetStaticBox(), size=(-1, 80))
        grid_out.Add(self.lst_output_pads, 1, wx.EXPAND)
        
        grid_out.AddGrowableCol(1, 1)
        sb_out.Add(grid_out, 1, wx.EXPAND | wx.ALL, 5)
        vbox.Add(sb_out, 0, wx.EXPAND | wx.ALL, 5)
        
        # --- Params Section ---
        sb_param = wx.StaticBoxSizer(wx.VERTICAL, self, "Parameters")
        grid_p = wx.FlexGridSizer(2, 2, 5, 5)
        
        grid_p.Add(wx.StaticText(sb_param.GetStaticBox(), label="Type:"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_type = wx.ComboBox(sb_param.GetStaticBox(), choices=["LINEAR", "SWITCHING"], style=wx.CB_READONLY)
        self.cmb_type.SetValue("LINEAR")
        self.cmb_type.Bind(wx.EVT_COMBOBOX, self.on_type_change)
        grid_p.Add(self.cmb_type, 1, wx.EXPAND)
        
        self.lbl_eff = wx.StaticText(sb_param.GetStaticBox(), label="Efficiency (0.0-1.0):")
        grid_p.Add(self.lbl_eff, 0, wx.ALIGN_CENTER_VERTICAL)
        self.txt_eff = wx.TextCtrl(sb_param.GetStaticBox(), value="0.85")
        self.txt_eff.Disable()
        grid_p.Add(self.txt_eff, 1, wx.EXPAND)

        grid_p.Add(wx.StaticText(
            sb_param.GetStaticBox(), label="Thermal loss component:"
        ), 0, wx.ALIGN_CENTER_VERTICAL)
        self.cmb_thermal_comp = wx.ComboBox(
            sb_param.GetStaticBox(), choices=[], style=wx.CB_READONLY
        )
        grid_p.Add(self.cmb_thermal_comp, 1, wx.EXPAND)
        
        grid_p.AddGrowableCol(1, 1)
        sb_param.Add(grid_p, 1, wx.EXPAND | wx.ALL, 5)
        vbox.Add(sb_param, 0, wx.EXPAND | wx.ALL, 5)
        
        # Buttons
        btns = wx.StdDialogButtonSizer()
        self.btn_ok = wx.Button(self, wx.ID_OK)
        self.btn_ok.Bind(wx.EVT_BUTTON, self.on_ok)
        btns.AddButton(self.btn_ok)
        btns.AddButton(wx.Button(self, wx.ID_CANCEL))
        btns.Realize()
        vbox.Add(btns, 0, wx.EXPAND | wx.ALL, 10)
        
        self.SetSizer(vbox)
        
    def _on_input_rail_change(self, event):
        rail = self.cmb_input_rail.GetValue()
        self.input_comps = self.discoverer.get_components_on_net(rail)
        refs = sorted(self.input_comps.keys())
        self.cmb_input_comp.Set(refs)
        if refs:
            self.cmb_input_comp.SetSelection(0)
            self._on_input_comp_change(None)
        else:
            self.cmb_input_comp.Clear()
            self.lst_input_pads.Clear()

    def _on_input_comp_change(self, event):
        ref = self.cmb_input_comp.GetValue()
        pads = self.input_comps.get(ref, [])
        pad_names = [getattr(p, 'number', getattr(p, 'name', str(p))) for p in pads]
        self.lst_input_pads.Set(pad_names) 
        self._update_auto_name()
        thermal_ref = self.cmb_thermal_comp.GetValue()
        preferred = ref if not thermal_ref or thermal_ref == self._last_input_ref else None
        self._last_input_ref = ref
        self._update_thermal_components(preferred)

    def _on_output_rail_change(self, event):
        rail = self.cmb_output_rail.GetValue()
        self.output_comps = self.discoverer.get_components_on_net(rail)
        refs = sorted(self.output_comps.keys())
        self.cmb_output_comp.Set(refs)
        if refs:
            self.cmb_output_comp.SetSelection(0)
            self._on_output_comp_change(None)
        else:
            self.cmb_output_comp.Clear()
            self.lst_output_pads.Clear()
            self._update_auto_name()
            
    def _on_output_comp_change(self, event):
        ref = self.cmb_output_comp.GetValue()
        pads = self.output_comps.get(ref, [])
        pad_names = [getattr(p, 'number', getattr(p, 'name', str(p))) for p in pads]
        self.lst_output_pads.Set(pad_names)
        self._update_auto_name()
        self._update_thermal_components()

    def _update_thermal_components(self, preferred=None):
        """Keep thermal placement independent from connectivity endpoints."""
        current = preferred or self.cmb_thermal_comp.GetValue()
        refs = set(self.input_comps) | set(self.output_comps)
        refs.update(filter(None, (
            self.cmb_input_comp.GetValue(),
            self.cmb_output_comp.GetValue(),
            preferred,
        )))
        choices = sorted(refs)
        self.cmb_thermal_comp.Set(choices)
        if current in choices:
            self.cmb_thermal_comp.SetValue(current)
        elif self.cmb_input_comp.GetValue() in choices:
            self.cmb_thermal_comp.SetValue(self.cmb_input_comp.GetValue())
        elif choices:
            self.cmb_thermal_comp.SetSelection(0)

    def on_type_change(self, event):
        if self.cmb_type.GetValue() == "SWITCHING":
            self.txt_eff.Enable()
        else:
            self.txt_eff.Disable()


    def _update_auto_name(self):
        inp_ref = self.cmb_input_comp.GetValue()
        out_ref = self.cmb_output_comp.GetValue()
        out_rail = self.cmb_output_rail.GetValue()
        
        name = generate_regulator_name(inp_ref, out_ref, out_rail)
        self.txt_name.ChangeValue(name)
    
    def prepopulate(self, regulator):
        """Prepopulate dialog with existing regulator data."""
        self.txt_name.SetValue(regulator.name)
        
        # Set rails
        self.cmb_input_rail.SetValue(regulator.input_rail_name)
        self.cmb_output_rail.SetValue(regulator.output_rail_name)
        
        # Trigger component population
        self._on_input_rail_change(None)
        self._on_output_rail_change(None)
        
        # Set components
        self.cmb_input_comp.SetValue(regulator.input_ref_des)
        self.cmb_output_comp.SetValue(regulator.output_ref_des)
        
        # Trigger pad population
        self._on_input_comp_change(None)
        self._on_output_comp_change(None)
        
        # Check the appropriate pads
        for i in range(self.lst_input_pads.GetCount()):
            pad_name = self.lst_input_pads.GetString(i)
            if pad_name in regulator.input_pad_names:
                self.lst_input_pads.Check(i, True)
        
        for i in range(self.lst_output_pads.GetCount()):
            pad_name = self.lst_output_pads.GetString(i)
            if pad_name in regulator.output_pad_names:
                self.lst_output_pads.Check(i, True)
        
        # Set type and efficiency
        self.cmb_type.SetValue(regulator.reg_type)
        self.txt_eff.SetValue(str(regulator.efficiency))
        self._update_thermal_components(
            regulator.thermal_ref_des or regulator.input_ref_des
        )
        
        # Enable/disable efficiency based on type
        if regulator.reg_type == "SWITCHING":
            self.txt_eff.Enable()
        else:
            self.txt_eff.Disable()
            
    def on_ok(self, event):
        name = self.txt_name.GetValue().strip()
        if not name:
            wx.MessageBox("Name is required.", "Error", wx.OK | wx.ICON_ERROR)
            return

        inp = self.cmb_input_rail.GetValue()
        out = self.cmb_output_rail.GetValue()
        
        if not inp or not out:
             wx.MessageBox("Select both Input and Output rails.", "Error", wx.OK | wx.ICON_ERROR)
             return
        if inp == out:
            wx.MessageBox("Input and Output rails cannot be the same.", "Error", wx.OK | wx.ICON_ERROR)
            return
            
        if not self.cmb_input_comp.GetValue() or not self.lst_input_pads.GetCheckedItems():
            wx.MessageBox("Select Input Component and at least one Pad.", "Error", wx.OK | wx.ICON_ERROR)
            return
            
        if not self.cmb_output_comp.GetValue() or not self.lst_output_pads.GetCheckedItems():
            wx.MessageBox("Select Output Component and at least one Pad.", "Error", wx.OK | wx.ICON_ERROR)
            return

        if not self.cmb_thermal_comp.GetValue():
            wx.MessageBox("Select the component that dissipates the regulator loss.", "Error", wx.OK | wx.ICON_ERROR)
            return

        if self.cmb_type.GetValue() == "SWITCHING":
            try:
                eff = float(self.txt_eff.GetValue())
                if not (0.0 < eff <= 1.0): raise ValueError()
            except ValueError:
                wx.MessageBox("Efficiency must be between 0.0 and 1.0", "Error", wx.OK | wx.ICON_ERROR)
                return
                
        event.Skip()
        
    def GetValue(self):
        # Get selected pads names
        in_pads = [self.lst_input_pads.GetString(i) for i in self.lst_input_pads.GetCheckedItems()]
        out_pads = [self.lst_output_pads.GetString(i) for i in self.lst_output_pads.GetCheckedItems()]
        
        return {
            'name': self.txt_name.GetValue().strip(),
            'input_rail': self.cmb_input_rail.GetValue(),
            'input_ref_des': self.cmb_input_comp.GetValue(),
            'input_pads': in_pads,
            'output_rail': self.cmb_output_rail.GetValue(),
            'output_ref_des': self.cmb_output_comp.GetValue(),
            'output_pads': out_pads,
            'type': self.cmb_type.GetValue(),
            'efficiency': float(self.txt_eff.GetValue()) if self.cmb_type.GetValue() == "SWITCHING" else 1.0,
            'thermal_ref_des': self.cmb_thermal_comp.GetValue(),
        }
