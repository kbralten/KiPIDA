import wx
import pcbnew

class ComponentSelectorDialog(wx.Dialog):
    def __init__(self, parent, title, net_name, component_dict):
        super(ComponentSelectorDialog, self).__init__(parent, title=title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        self.net_name = net_name
        self.comp_dict = component_dict # { ref: [pads] }
        self.sorted_refs = sorted(list(self.comp_dict.keys()))
        
        self.selected_ref = None
        self.selected_pads = []
        self.target_value = 0.0 # Voltage or Current
        self.mode = "LOAD" # or SOURCE
        
        self._init_ui()
        self.Center()
        
    def _init_ui(self):
        sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 1. Component Dropdown
        h1 = wx.BoxSizer(wx.HORIZONTAL)
        h1.Add(wx.StaticText(self, label="Component:"), 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.cb_comp = wx.ComboBox(self, choices=self.sorted_refs, style=wx.CB_READONLY)
        h1.Add(self.cb_comp, 1, wx.EXPAND | wx.ALL, 5)
        sizer.Add(h1, 0, wx.EXPAND | wx.ALL, 5)
        
        # 2. Pads Selector (CheckListBox)
        sizer.Add(wx.StaticText(self, label="Enabled Pads:"), 0, wx.LEFT, 10)
        self.clb_pads = wx.CheckListBox(self, size=(-1, 100))
        sizer.Add(self.clb_pads, 1, wx.EXPAND | wx.ALL, 5)
        
        # 3. Value Input
        h2 = wx.BoxSizer(wx.HORIZONTAL)
        self.lbl_value = wx.StaticText(self, label="Value:")
        h2.Add(self.lbl_value, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.txt_value = wx.TextCtrl(self)
        h2.Add(self.txt_value, 1, wx.EXPAND | wx.ALL, 5)
        self.lbl_unit = wx.StaticText(self, label="A")
        h2.Add(self.lbl_unit, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        sizer.Add(h2, 0, wx.EXPAND | wx.ALL, 5)
        
        # Buttons
        btns = self.CreateButtonSizer(wx.OK | wx.CANCEL)
        sizer.Add(btns, 0, wx.EXPAND | wx.ALL, 5)
        
        self.SetSizer(sizer)
        self.Fit()
        
        # Bindings
        self.cb_comp.Bind(wx.EVT_COMBOBOX, self.on_comp_select)
        
        # Init first if exists
        if(self.sorted_refs):
            self.cb_comp.SetSelection(0)
            self.on_comp_select(None)
            
    def set_mode(self, mode):
        self.mode = mode
        if mode == "SOURCE":
            self.lbl_value.SetLabel("Voltage (Not used)")
            self.txt_value.Disable() # Voltage is rail-level now
            self.lbl_unit.SetLabel("V")
            self.SetTitle(f"Add Source to {self.net_name}")
        else:
            self.lbl_value.SetLabel("Total Current:")
            self.txt_value.Enable()
            self.lbl_unit.SetLabel("A")
            self.SetTitle(f"Add Load to {self.net_name}")

    def on_comp_select(self, event):
        ref = self.cb_comp.GetValue()
        pads = self.comp_dict.get(ref, [])
        
        self.clb_pads.Clear()
        
        # Store raw pad objects to retrieve data later
        # Use p.GetName() for sorting as it's the pad 'number' string
        self.current_pads = sorted(pads, key=lambda p: p.GetName())
        
        display_names = []
        for p in self.current_pads:
            num = p.GetName()
            # Try to get pin function (name in schematic)
            try:
                name = p.GetPinFunction()
            except:
                name = ""
                
            if name:
                display_names.append(f"{num} - {name}")
            else:
                display_names.append(num)
        
        self.clb_pads.InsertItems(display_names, 0)
        
        # Check all by default
        for i in range(len(display_names)):
            self.clb_pads.Check(i, True)
            
    def GetSelection(self):
        """Returns (ref, value, [enabled_pad_names])"""
        ref = self.cb_comp.GetValue()
        
        # Get checked pad names from the original objects
        selected_pad_names = []
        for i in range(self.clb_pads.GetCount()):
            if self.clb_pads.IsChecked(i):
                selected_pad_names.append(self.current_pads[i].GetName())
                
        val_str = self.txt_value.GetValue().strip()
        try:
            val = float(val_str) if val_str else 0.0
        except:
            val = 0.0
            
        return ref, val, selected_pad_names
