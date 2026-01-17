import wx
import wx.dataview
import pcbnew
import sys
import os

# Ensure plugin dir is in path to import modules
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

from extractor import GeometryExtractor
from mesh import Mesher

from solver import Solver
from ui.setup_dialogs import SourceLoadDialog

class KiPIDA_MainDialog(wx.Dialog):
    def __init__(self, parent):
        super(KiPIDA_MainDialog, self).__init__(parent, title="Ki-PIDA: Power Integrity Analyzer", 
                                                style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        self.SetSize((900, 700))
        self.SetMinSize((600, 500))
        
        self.board = pcbnew.GetBoard()
        self.sl_items = [] # list of dict {'pads': [pad_objs], 'value': float, 'type': 'V' or 'I'}
        
        self._init_ui()
        self.Center()
        
        self._load_nets()
    
    def _init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 1. Notebook for tabs
        self.notebook = wx.Notebook(self)
        
        # Tab 1: Configuration / Net Management
        self.tab_config = wx.Panel(self.notebook)
        self._init_config_tab(self.tab_config)
        self.notebook.AddPage(self.tab_config, "Configuration")
        
        # Tab 2: Results
        self.tab_results = wx.Panel(self.notebook)
        self._init_results_tab(self.tab_results)
        self.notebook.AddPage(self.tab_results, "Results")
        
        # Tab 3: Log/Debug
        self.tab_log = wx.Panel(self.notebook)
        self._init_log_tab(self.tab_log)
        self.notebook.AddPage(self.tab_log, "Log")
        
        main_sizer.Add(self.notebook, 1, wx.EXPAND | wx.ALL, 5)
        
        # 2. Action Buttons (Bottom)
        btn_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.btn_run = wx.Button(self, label="Run Simulation")
        self.btn_cancel = wx.Button(self, wx.ID_CANCEL, "Close")
        
        btn_sizer.AddStretchSpacer()
        btn_sizer.Add(self.btn_run, 0, wx.ALL, 5)
        btn_sizer.Add(self.btn_cancel, 0, wx.ALL, 5)
        
        main_sizer.Add(btn_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        self.SetSizer(main_sizer)
        
        # Bind events
        self.btn_run.Bind(wx.EVT_BUTTON, self.on_run)
        self.btn_cancel.Bind(wx.EVT_BUTTON, self.on_close)
    
    def _init_config_tab(self, parent):
        # Splitter: Left = Nets, Right = Sources/Loads
        splitter = wx.BoxSizer(wx.HORIZONTAL)
        
        # Left: Net Selection
        left_sizer = wx.BoxSizer(wx.VERTICAL)
        lbl_net = wx.StaticText(parent, label="1. Select Net:")
        self.net_list = wx.ListCtrl(parent, style=wx.LC_REPORT | wx.LC_SINGLE_SEL)
        self.net_list.InsertColumn(0, "Net Name", width=250)
        left_sizer.Add(lbl_net, 0, wx.ALL, 5)
        left_sizer.Add(self.net_list, 1, wx.EXPAND | wx.ALL, 5)
        
        splitter.Add(left_sizer, 1, wx.EXPAND | wx.ALL, 5)
        
        # Right: Sources and Loads
        right_sizer = wx.BoxSizer(wx.VERTICAL)
        lbl_sl = wx.StaticText(parent, label="2. Define Sources & Loads:")
        right_sizer.Add(lbl_sl, 0, wx.ALL, 5)
        
        self.sl_list = wx.ListCtrl(parent, style=wx.LC_REPORT)
        self.sl_list.InsertColumn(0, "Type", width=60)
        self.sl_list.InsertColumn(1, "Value", width=60)
        self.sl_list.InsertColumn(2, "Location", width=150)
        right_sizer.Add(self.sl_list, 1, wx.EXPAND | wx.ALL, 5)
        
        # Buttons
        btn_sl_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.btn_add_src = wx.Button(parent, label="+ Voltage Src")
        self.btn_add_load = wx.Button(parent, label="+ Current Load")
        self.btn_remove_sl = wx.Button(parent, label="- Remove Selected")
        self.btn_clear_sl = wx.Button(parent, label="Clear All")
        
        btn_sl_sizer.Add(self.btn_add_src, 0, wx.RIGHT, 5)
        btn_sl_sizer.Add(self.btn_add_load, 0, wx.RIGHT, 5)
        btn_sl_sizer.Add(self.btn_remove_sl, 0, wx.RIGHT, 5)
        btn_sl_sizer.Add(self.btn_clear_sl, 0, wx.RIGHT, 5)
        right_sizer.Add(btn_sl_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        # Grid Size Setting
        grid_sett_sizer = wx.BoxSizer(wx.HORIZONTAL)
        lbl_grid = wx.StaticText(parent, label="Mesh Grid Size (mm):")
        self.txt_grid_size = wx.TextCtrl(parent, value="0.1", size=(60, -1))
        grid_sett_sizer.Add(lbl_grid, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        grid_sett_sizer.Add(self.txt_grid_size, 0, wx.ALIGN_CENTER_VERTICAL)
        right_sizer.Add(grid_sett_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        # Debug Mode Setting
        debug_sett_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.chk_debug = wx.CheckBox(parent, label="Enable Debug Logging")
        self.chk_debug.SetValue(False)
        debug_sett_sizer.Add(self.chk_debug, 0, wx.ALIGN_CENTER_VERTICAL)
        right_sizer.Add(debug_sett_sizer, 0, wx.EXPAND | wx.ALL, 5)
        
        splitter.Add(right_sizer, 1, wx.EXPAND | wx.ALL, 5)
        
        parent.SetSizer(splitter)
        
        # Events
        self.btn_add_src.Bind(wx.EVT_BUTTON, self.on_add_source)
        self.btn_add_load.Bind(wx.EVT_BUTTON, self.on_add_load)
        self.btn_remove_sl.Bind(wx.EVT_BUTTON, self.on_remove_sl)
        self.btn_clear_sl.Bind(wx.EVT_BUTTON, self.on_clear_sl)
        
        self.sl_list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.on_edit_sl)

    def _init_results_tab(self, parent):
        sizer = wx.BoxSizer(wx.VERTICAL)
        lbl = wx.StaticText(parent, label="Simulation Results")
        sizer.Add(lbl, 0, wx.ALL, 5)
        
        # Splitter: Top=Text Stats, Bottom=Image
        self.result_splitter = wx.SplitterWindow(parent)
        
        self.pnl_text = wx.Panel(self.result_splitter)
        text_sizer = wx.BoxSizer(wx.VERTICAL)
        self.result_text = wx.TextCtrl(self.pnl_text, style=wx.TE_MULTILINE | wx.TE_READONLY)
        text_sizer.Add(self.result_text, 1, wx.EXPAND | wx.ALL, 5)
        self.pnl_text.SetSizer(text_sizer)
        
        self.pnl_image = wx.Panel(self.result_splitter)
        image_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # Scrolled Window for image
        self.scrolled_img = wx.ScrolledWindow(self.pnl_image, style=wx.HSCROLL | wx.VSCROLL)
        self.scrolled_img.SetScrollRate(10, 10)
        
        # result_image MUST be a child of scrolled_img
        self.result_image = wx.StaticBitmap(self.scrolled_img)
        
        scrolled_sizer = wx.BoxSizer(wx.VERTICAL)
        scrolled_sizer.Add(self.result_image, 1, wx.CENTER | wx.ALL, 5)
        self.scrolled_img.SetSizer(scrolled_sizer)
        
        image_sizer.Add(self.scrolled_img, 1, wx.EXPAND | wx.ALL, 0)
        self.pnl_image.SetSizer(image_sizer)
        
        self.result_splitter.SplitHorizontally(self.pnl_text, self.pnl_image, 100)
        self.result_splitter.SetMinimumPaneSize(50)
        
        sizer.Add(self.result_splitter, 1, wx.EXPAND | wx.ALL, 5)
        parent.SetSizer(sizer)

    def _init_log_tab(self, parent):
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.log_ctrl = wx.TextCtrl(parent, style=wx.TE_MULTILINE | wx.TE_READONLY)
        sizer.Add(self.log_ctrl, 1, wx.EXPAND | wx.ALL, 5)
        parent.SetSizer(sizer)

    def _load_nets(self):
        self.net_list.DeleteAllItems()
        nets = self.board.GetNetsByName()
        valid_nets = []
        for name, net in nets.items():
            s_name = str(name)
            if s_name == "": continue
            valid_nets.append(s_name)
        valid_nets.sort()
        for name in valid_nets:
            self.net_list.Append([name])

    def log(self, msg):
        self.log_ctrl.AppendText(msg + "\n")
        
    def get_selected_net(self):
        sel = self.net_list.GetFirstSelected()
        if sel == -1: return None
        return self.net_list.GetItemText(sel)

    def on_add_source(self, event):
        net = self.get_selected_net()
        if not net: 
            wx.MessageBox("Select a Net first.")
            return
            
        used = []
        for item in self.sl_items: used.extend(item['pads'])
            
        dlg = SourceLoadDialog(self, f"Add Voltage Source for {net}", net, self.board, used_pads=used)
        if dlg.ShowModal() == wx.ID_OK:
            val = dlg.GetValue()
            pads = dlg.GetSelectedPads()
            if not pads:
                wx.MessageBox("No pads selected.")
                return
            
            self.sl_items.append({'type': 'V', 'value': val, 'pads': pads})
            self._update_sl_list()
        dlg.Destroy()
        
    def on_add_load(self, event):
        net = self.get_selected_net()
        if not net: 
            wx.MessageBox("Select a Net first.")
            return

        used = []
        for item in self.sl_items: used.extend(item['pads'])

        dlg = SourceLoadDialog(self, f"Add Current Load for {net}", net, self.board, used_pads=used)
        if dlg.ShowModal() == wx.ID_OK:
            val = dlg.GetValue()
            pads = dlg.GetSelectedPads()
            if not pads: return
            
            self.sl_items.append({'type': 'I', 'value': val, 'pads': pads})
            self._update_sl_list()
        dlg.Destroy()

    def on_edit_sl(self, event):
        sel = self.sl_list.GetFirstSelected()
        if sel == -1: return
        
        net = self.get_selected_net()
        item = self.sl_items[sel]
        
        is_src = (item['type'] == 'V')
        title = f"Edit {'Voltage Source' if is_src else 'Current Load'} for {net}"
            
        used = []
        for it in self.sl_items: used.extend(it['pads'])
        
        dlg = SourceLoadDialog(self, title, net, self.board, 
                               used_pads=used, 
                               initial_value=item['value'], 
                               initial_pads=item['pads'])
        
        if dlg.ShowModal() == wx.ID_OK:
            val = dlg.GetValue()
            pads = dlg.GetSelectedPads()
            if not pads:
                wx.MessageBox("No pads selected.")
                return
            
            # Update memory - type stays the same!
            item['value'] = val
            item['pads'] = pads
            self._update_sl_list()
            
        dlg.Destroy()

    def on_remove_sl(self, event):
        sel = self.sl_list.GetFirstSelected()
        if sel == -1: return
        self.sl_items.pop(sel)
        self._update_sl_list()
        
    def on_clear_sl(self, event):
        self.sl_items = []
        self._update_sl_list()

    def _update_sl_list(self):
        self.sl_list.DeleteAllItems()
        for item in self.sl_items:
            is_src = (item['type'] == 'V')
            type_str = "Voltage" if is_src else "Load"
            unit = "V" if is_src else "A"
            loc = f"{len(item['pads'])} pads"
            self.sl_list.Append([type_str, f"{item['value']} {unit}", loc])

        # Lock net selection if items exist
        has_items = (len(self.sl_items) > 0)
        self.net_list.Enable(not has_items)

    def on_run(self, event):
        net_name = self.get_selected_net()
        if not net_name: return
        
        self.notebook.SetSelection(2) # Show Log
        self.log(f"--- Starting Simulation for {net_name} ---")
        
        if not self.sl_items:
            self.log("WARNING: No sources or loads defined. Mesh will be generated but results will be zero.")
        
        try:
            # 1. Extraction
            debug_mode = self.chk_debug.GetValue()
            if debug_mode:
                self.log("Debug mode enabled - detailed diagnostics will follow...")
            extractor = GeometryExtractor(self.board, debug=debug_mode, log_callback=self.log)
            stackup = extractor.get_board_stackup()
            geo = extractor.get_net_geometry(net_name)
            if not geo:
                self.log("No geometry found.")
                return

            self.log(f"Geometry extracted. Layers: {list(geo.keys())}")
            
            # Visualize extracted geometry if debug mode
            if debug_mode:
                self.log("Generating geometry visualization...")
                try:
                    import wx
                    geo_buf = extractor.plot_geometry(geo)
                    if geo_buf:
                        image = wx.Image(geo_buf, wx.BITMAP_TYPE_PNG)
                        geo_bitmap = wx.Bitmap(image)
                        # Show in results tab temporarily
                        self.result_image.SetBitmap(geo_bitmap)
                        self.scrolled_img.GetSizer().Layout()
                        self.scrolled_img.FitInside()
                        self.notebook.SetSelection(1)  # Show Results tab
                        self.log("Geometry plot displayed in Results tab.")
                except Exception as e:
                    self.log(f"Geometry plot failed: {e}")
            
            # 2. Meshing
            self.log("Generating Mesh...")
            try:
                gs = float(self.txt_grid_size.GetValue())
                if gs < 0.01: gs = 0.01 # Safety
            except:
                gs = 0.1
                
            mesher = Mesher(self.board, debug=debug_mode, log_callback=self.log)
            mesh = mesher.generate_mesh(net_name, geo, stackup, grid_size_mm=gs)
            self.log(f"Mesh generated: {len(mesh.nodes)} nodes, {len(mesh.edges)} edges.")
            
            # Show mesh visualization before solving if debug mode
            if debug_mode:
                self.log("Generating pre-solve mesh visualization...")
                try:
                    import wx
                    mesh_bitmap = mesher.debug_plot(mesh, stackup=stackup)
                    if mesh_bitmap:
                        self.result_image.SetBitmap(mesh_bitmap)
                        self.scrolled_img.GetSizer().Layout()
                        self.scrolled_img.FitInside()
                        self.notebook.SetSelection(1)  # Show Results tab
                        self.log("Pre-solve mesh plot displayed (no voltage data yet).")
                except Exception as e:
                    self.log(f"Pre-solve mesh plot failed: {e}")
            
            # 3. Solver Setup
            # Need to map Source/Load Pads to Mesh Nodes
            # Mesh has `node_coords` { nid: (x,y,layer) }
            # and `node_map` { (ix, iy, layer): nid }
            
            solver_sources = []
            solver_loads = []
            
            grid_step = mesh.grid_step
            origin = mesh.grid_origin
            
            # Helper to find node closest to pad
            def find_nodes_for_pads(pads):
                found_nodes = []
                for pad in pads:
                    pos = pad.GetPosition()
                    px = pcbnew.ToMM(pos.x)
                    py = pcbnew.ToMM(pos.y)
                    
                    # Target grid index
                    tx = int(round((px - origin[0]) / grid_step))
                    ty = int(round((py - origin[1]) / grid_step))
                    
                    # Search 3x3 neighborhood for ANY node on ANY layer
                    node_found = False
                    for dx in [-1, 0, 1]:
                        for dy in [-1, 0, 1]:
                            ix, iy = tx + dx, ty + dy
                            for layer in range(32): # All logical copper layers
                                nid = mesh.node_map.get((ix, iy, layer))
                                if nid is not None:
                                    found_nodes.append(nid)
                                    node_found = True
                    
                    if not node_found:
                        self.log(f"Warning: Pad at ({px:.2f}, {py:.2f}) [Grid {tx},{ty}] not on mesh nodes.")
                    else:
                        # Log success for debug
                        pass
                        
                return list(set(found_nodes)) # Unique nodes

            # Prepare Sources/Loads
            for item in self.sl_items:
                nodes = find_nodes_for_pads(item['pads'])
                if not nodes: continue
                
                if item['type'] == 'V':
                    # Voltage source: Apply V to all found nodes
                    for nid in nodes:
                        solver_sources.append({'node_id': nid, 'voltage': item['value']})
                else:
                    # Current source: Divide I by number of nodes found
                    i_per_node = item['value'] / len(nodes)
                    for nid in nodes:
                        solver_loads.append({'node_id': nid, 'current': i_per_node})

            # 4. Solve
            if not solver_sources:
                 self.log("Error: No Voltage Sources defined! Singular matrix.")
            else:
                self.log("Solving Linear System...")
                solver = Solver()
                results = solver.solve(mesh, solver_sources, solver_loads)
                self.log("Solution Complete.")
                
                # Show Results
                if results:
                    v_vals = list(results.values())
                    v_min = min(v_vals)
                    v_max = max(v_vals)
                    drop = v_max - v_min
                    r_str = f"Min V: {v_min:.4f} V\nMax V: {v_max:.4f} V\nDrop: {drop:.4f} V"
                    self.result_text.SetValue(r_str)
                    
                # 5. Visualize
                mesh.results = results
                self.log("Generating Plot...")
                bitmap = mesher.debug_plot(mesh, stackup=stackup)
                
                if bitmap:
                    self.result_image.SetBitmap(bitmap)
                    self.scrolled_img.GetSizer().Layout()
                    self.scrolled_img.FitInside()
                    self.log("Results plot generated in-memory.")
                else:
                    self.log("Failed to generate plot bitmap.")
                
                self.notebook.SetSelection(1) # Results tab
            
        except Exception as e:
            self.log(f"Error: {e}")
            import traceback
            traceback.print_exc()

    def on_close(self, event):
        self.EndModal(wx.ID_CANCEL)
