import wx

class ComponentSelectorDialog(wx.Dialog):
    def __init__(self, parent, title, net_name, components):
        super().__init__(parent, title=title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.net_name = net_name
        self.components = components
        self.mode = "LOAD"
        
        self.selected_ref = None
        self.selected_value = 0.0
        self.selected_pads = []
        
        self._init_ui()
        self.SetSize((400, 500))
        self.Center()
        
    def set_mode(self, mode):
        self.mode = mode
        if mode == "SOURCE":
            self.lbl_instruction.SetLabel("Select Component to use as Voltage Source:")
            # Hide value input for sources (voltage is set at rail level)
            self.lbl_val.Hide()
            self.txt_val.Hide()
        else:
            self.lbl_instruction.SetLabel("Select Component to use as Load:")
            self.lbl_val.SetLabel("Current (A):")
            self.txt_val.SetValue("1.0")
            self.lbl_val.Show()
            self.txt_val.Show()
            
    def _init_ui(self):
        panel = wx.Panel(self)
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        self.lbl_instruction = wx.StaticText(panel, label="Select Component:")
        vbox.Add(self.lbl_instruction, 0, wx.ALL, 5)
        
        self.lst_comps = wx.ListBox(panel, style=wx.LB_SINGLE | wx.LB_SORT)
        self.comp_map = {}
        
        # self.components is a dict { ref: [pad_objs] }
        if isinstance(self.components, dict):
             for ref, pads in self.components.items():
                 self.lst_comps.Append(ref)
                 self.comp_map[ref] = pads
        else:
             # Fallback for list of objects (legacy)
             for c in self.components:
                 ref = getattr(c, 'ref_des', str(c))
                 self.lst_comps.Append(ref)
                 self.comp_map[ref] = c
            
        vbox.Add(self.lst_comps, 1, wx.EXPAND | wx.ALL, 5)
        
        h_val = wx.BoxSizer(wx.HORIZONTAL)
        self.lbl_val = wx.StaticText(panel, label="Value:")
        self.txt_val = wx.TextCtrl(panel)
        h_val.Add(self.lbl_val, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        h_val.Add(self.txt_val, 1, wx.EXPAND)
        vbox.Add(h_val, 0, wx.EXPAND | wx.ALL, 5)
        
        # Pad selection
        vbox.Add(wx.StaticText(panel, label="Select Pads:"), 0, wx.ALL, 5)
        
        self.chk_all_pads = wx.CheckBox(panel, label="Use all connected pads")
        self.chk_all_pads.SetValue(True)
        vbox.Add(self.chk_all_pads, 0, wx.ALL, 5)
        
        # Pad list (shown when "use all" is unchecked)
        self.pad_list = wx.CheckListBox(panel)
        self.pad_list.Enable(False)  # Disabled until component selected and "use all" unchecked
        vbox.Add(self.pad_list, 0, wx.EXPAND | wx.ALL, 5)
        
        h_btns = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_ok = wx.Button(panel, wx.ID_OK, "OK")
        self.btn_cancel = wx.Button(panel, wx.ID_CANCEL, "Cancel")
        
        # Helper to make OK default
        self.btn_ok.SetDefault()
        
        h_btns.AddStretchSpacer()
        h_btns.Add(self.btn_cancel, 0, wx.RIGHT, 10)
        h_btns.Add(self.btn_ok, 0)
        
        vbox.Add(h_btns, 0, wx.EXPAND | wx.ALL, 5)
        
        panel.SetSizer(vbox)
        
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        main_sizer.Add(panel, 1, wx.EXPAND)
        self.SetSizer(main_sizer)
        
        self.lst_comps.Bind(wx.EVT_LISTBOX, self.on_select)
        self.chk_all_pads.Bind(wx.EVT_CHECKBOX, self.on_all_pads_toggle)
        
    def on_all_pads_toggle(self, event):
        """Enable/disable pad list based on checkbox"""
        use_all = self.chk_all_pads.GetValue()
        self.pad_list.Enable(not use_all)
        
    def on_select(self, event):
        """Update pad list when component is selected"""
        sel_idx = self.lst_comps.GetSelection()
        if sel_idx != wx.NOT_FOUND:
            ref = self.lst_comps.GetString(sel_idx)
            comp_data = self.comp_map.get(ref)
            
            # Populate pad list
            self.pad_list.Clear()
            if isinstance(comp_data, list):
                for p in comp_data:
                    p_name = getattr(p, 'number', getattr(p, 'name', ''))
                    if p_name:
                        self.pad_list.Append(str(p_name))
                        self.pad_list.Check(self.pad_list.GetCount() - 1, True)  # Check all by default
            elif hasattr(comp_data, 'pads'):
                for p in comp_data.pads:
                    p_name = getattr(p, 'number', getattr(p, 'name', ''))
                    if p_name:
                        self.pad_list.Append(str(p_name))
                        self.pad_list.Check(self.pad_list.GetCount() - 1, True)

    def GetSelection(self):
        sel_idx = self.lst_comps.GetSelection()
        if sel_idx != wx.NOT_FOUND:
            ref = self.lst_comps.GetString(sel_idx)
            try:
                val = float(self.txt_val.GetValue()) if self.txt_val.IsShown() else 0.0
            except:
                val = 0.0
                
            pads = []
            
            # Check if using all pads or selected pads
            if self.chk_all_pads.GetValue():
                # Use all pads
                comp_data = self.comp_map.get(ref)
                
                if isinstance(comp_data, list):
                    for p in comp_data:
                        p_name = getattr(p, 'number', getattr(p, 'name', ''))
                        if p_name:
                            pads.append(str(p_name))
                elif hasattr(comp_data, 'pads'):
                    for p in comp_data.pads:
                        p_name = getattr(p, 'number', getattr(p, 'name', ''))
                        if p_name:
                            pads.append(str(p_name))
            else:
                # Use only checked pads
                for i in range(self.pad_list.GetCount()):
                    if self.pad_list.IsChecked(i):
                        pads.append(self.pad_list.GetString(i))
            
            return ref, val, pads
        return None, 0.0, []

    def prepopulate(self, ref_des, value, pads):
        """Pre-set the dialog state for editing"""
        idx = self.lst_comps.FindString(ref_des)
        if idx != wx.NOT_FOUND:
            self.lst_comps.SetSelection(idx)
            self.on_select(None) # Load pads
            
            # Set value
            if self.mode == "LOAD":
                self.txt_val.SetValue(str(value))
            
            # Set pads
            all_comp_pads = [self.pad_list.GetString(i) for i in range(self.pad_list.GetCount())]
            if set(all_comp_pads) == set(pads):
                self.chk_all_pads.SetValue(True)
                self.pad_list.Enable(False)
            else:
                self.chk_all_pads.SetValue(False)
                self.pad_list.Enable(True)
                for i in range(self.pad_list.GetCount()):
                    p_name = self.pad_list.GetString(i)
                    self.pad_list.Check(i, p_name in pads)