import wx
# import pcbnew # Removed to avoid SWIG dependency crash
class SourceLoadDialog(wx.Dialog):
    def __init__(self, parent, title, net_name, board, used_pads=None, initial_value=0.0, initial_pads=None):
        super(SourceLoadDialog, self).__init__(parent, title=title, size=(400, 350))
        
        self.net_name = net_name
        self.board = board
        self.used_pads = used_pads if used_pads else []
        self.initial_pads = initial_pads if initial_pads else []
        self.selected_pads = []
        self.value = initial_value
        
        self._init_ui()
        if initial_value != 0.0:
            self.txt_val.SetValue(str(initial_value))
        self.Centre()
        
    def _init_ui(self):
        vbox = wx.BoxSizer(wx.VERTICAL)
        
        # 1. Inputs
        grid = wx.FlexGridSizer(2, 2, 10, 10)
        
        lbl_val = wx.StaticText(self, label="Value (V or A):")
        self.txt_val = wx.TextCtrl(self)
        
        lbl_pads = wx.StaticText(self, label="Select Pads (Disabled = Already Used):")
        self.lst_pads = wx.CheckListBox(self, choices=[])
        self._populate_pads()
        
        grid.Add(lbl_val, 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(self.txt_val, 1, wx.EXPAND)
        grid.Add(lbl_pads, 0, wx.ALIGN_TOP)
        grid.Add(self.lst_pads, 1, wx.EXPAND | wx.ALL, 5)
        
        grid.AddGrowableCol(1, 1)
        grid.AddGrowableRow(1, 1)
        
        vbox.Add(grid, 1, wx.EXPAND | wx.ALL, 10)
        
        # 2. Buttons
        self.btn_ok = wx.Button(self, wx.ID_OK)
        self.btn_cancel = wx.Button(self, wx.ID_CANCEL)
        
        btns = wx.StdDialogButtonSizer()
        btns.AddButton(self.btn_ok)
        btns.AddButton(self.btn_cancel)
        btns.Realize()
        vbox.Add(btns, 0, wx.EXPAND | wx.ALL, 10)
        
        self.SetSizer(vbox)
        
        # Bind OK button for validation
        self.btn_ok.Bind(wx.EVT_BUTTON, self.on_ok)

    def on_ok(self, event):
        val_str = self.txt_val.GetValue().strip()
        
        # 1. Value Validation
        try:
            val = float(val_str)
            if val == 0:
                wx.MessageBox("Value cannot be zero.", "Validation Error", wx.OK | wx.ICON_ERROR)
                return
        except ValueError:
            wx.MessageBox("Please enter a valid numeric value.", "Validation Error", wx.OK | wx.ICON_ERROR)
            return
            
        # 2. Pad Selection Validation
        if not self.GetSelectedPads():
            wx.MessageBox("Please select at least one pad.", "Validation Error", wx.OK | wx.ICON_ERROR)
            return
            
        event.Skip() # Allow default processing (closes dialog with ID_OK)
        
    def _populate_pads(self):
        # Find all pads on this net
        net = self.board.FindNet(self.net_name)
        if not net: return
        
        net_code = net.GetNetCode()
        pads = []
        
        # Helper to check if pad is used
        def get_used_ref(pad):
            # If it's in our INITIAL pads, it's NOT "already used" by someone else
            for ip in self.initial_pads:
                if ip.GetPosition() == pad.GetPosition() and ip.GetName() == pad.GetName():
                    return False
            for up in self.used_pads:
                if up.GetPosition() == pad.GetPosition() and up.GetName() == pad.GetName():
                    return True
            return False
            
        def is_initial(pad):
            for ip in self.initial_pads:
                if ip.GetPosition() == pad.GetPosition() and ip.GetName() == pad.GetName():
                    return True
            return False

        for fp in self.board.GetFootprints():
            for pad in fp.Pads():
                if pad.GetNetCode() == net_code:
                    name = f"{fp.GetReference()} - Pad {pad.GetName()}"
                    is_p_used = get_used_ref(pad)
                    if is_p_used:
                        name = "[ALREADY USED] " + name
                    pads.append((name, pad, is_p_used))
                    
        # Sort by name
        pads.sort(key=lambda x: x[0])
        
        self.pad_data = pads # list of (name, pad_obj, is_used)
        names = [p[0] for p in pads]
        self.lst_pads.Set(names)
        
        # Pre-check initial pads
        for i, (name, pad, used) in enumerate(pads):
            if is_initial(pad):
                self.lst_pads.Check(i, True)
        
    def GetSelectedPads(self):
        indices = self.lst_pads.GetCheckedItems()
        selected = []
        for i in indices:
            name, pad, used = self.pad_data[i]
            if not used:
                selected.append(pad)
            else:
                # Uncheck it if they clicked a used one (dialog stays open)
                self.lst_pads.Check(i, False)
        return selected
        
    def GetValue(self):
        try:
            return float(self.txt_val.GetValue())
        except ValueError:
            return 0.0
            
    def GetSelectedPads(self):
        indices = self.lst_pads.GetCheckedItems()
        return [self.pad_data[i][1] for i in indices]

class VoltageSourceDialog(SourceLoadDialog):
    def _init_ui(self):
        super()._init_ui()
        # Customize labels if needed? 
        # For MVP generic is fine, but let's be clearer
        # We can't easily change static text after Create but we can in init
        pass

