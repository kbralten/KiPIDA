import wx
import wx.dataview
import pcbnew
import sys
import os

# Ensure plugin dir is in path to import modules
plugin_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if plugin_dir not in sys.path:
    sys.path.insert(0, plugin_dir)

try:
    from ..extractor import GeometryExtractor
    from ..mesh import Mesher
    from ..solver import Solver
    from .power_tree_panel import PowerTreePanel
except (ImportError, ValueError):
    from extractor import GeometryExtractor
    from mesh import Mesher
    from solver import Solver
    from ui.power_tree_panel import PowerTreePanel

class KiPIDA_MainDialog(wx.Dialog):
    def __init__(self, parent):
        super(KiPIDA_MainDialog, self).__init__(parent, title="Ki-PIDA: Power Integrity Analyzer", 
                                                style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        
        self.SetSize((1000, 700))
        self.SetMinSize((800, 500))
        
        self.board = pcbnew.GetBoard()
        
        self._init_ui()
        self.Center()
        
    def _init_ui(self):
        main_sizer = wx.BoxSizer(wx.VERTICAL)
        
        # 1. Notebook for tabs
        self.notebook = wx.Notebook(self)
        
        # Tab 1: Configuration (New Power Tree Panel)
        self.tab_config = wx.Panel(self.notebook)
        self.power_tree = PowerTreePanel(self.tab_config, self.board)
        
        # Config Tab Layout
        config_sizer = wx.BoxSizer(wx.VERTICAL)
        config_sizer.Add(self.power_tree, 1, wx.EXPAND | wx.ALL, 5)
        
        # Global Settings (Grid Size, Debug)
        sett_sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        lbl_grid = wx.StaticText(self.tab_config, label="Mesh Resolution (mm):")
        self.txt_grid_size = wx.TextCtrl(self.tab_config, value="0.1", size=(60, -1))
        sett_sizer.Add(lbl_grid, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 5)
        sett_sizer.Add(self.txt_grid_size, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, 20)
        
        self.chk_debug = wx.CheckBox(self.tab_config, label="Enable Debug Log")
        sett_sizer.Add(self.chk_debug, 0, wx.ALIGN_CENTER_VERTICAL)
        
        config_sizer.Add(sett_sizer, 0, wx.EXPAND | wx.ALL, 5)
        self.tab_config.SetSizer(config_sizer)
        
        self.notebook.AddPage(self.tab_config, "Power Tree & Config")
        
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
    
    def _init_results_tab(self, parent):
        sizer = wx.BoxSizer(wx.VERTICAL)
        
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

    def log(self, msg):
        self.log_ctrl.AppendText(msg + "\n")
        
    def on_run(self, event):
        rail = self.power_tree.active_rail
        if not rail:
            wx.MessageBox("Please select a Power Rail to simulate.")
            return
            
        net_name = rail.net_name
        self.notebook.SetSelection(2) # Switch to Log
        self.log(f"--- Starting Simulation for {net_name} ---")
        self.log(f"Nominal Voltage: {rail.nominal_voltage} V")
        
        if not rail.sources:
             self.log("WARNING: No Voltage Sources defined. Simulation will likely fail or return all zeros.")

        try:
            debug_mode = self.chk_debug.GetValue()
            
            # --- 1. Extraction ---
            extractor = GeometryExtractor(self.board, debug=debug_mode, log_callback=self.log)
            stackup = extractor.get_board_stackup()
            geo = extractor.get_net_geometry(net_name)
            
            if not geo:
                self.log("No geometry found for net.")
                return

            # Plot Geometry (Debug)
            if debug_mode:
                 self._debug_plot_geo(extractor, geo)

            # --- 2. Meshing ---
            try:
                gs = float(self.txt_grid_size.GetValue())
                if gs < 0.01: gs = 0.01
            except: gs = 0.1
            
            mesher = Mesher(self.board, debug=debug_mode, log_callback=self.log)
            mesh = mesher.generate_mesh(net_name, geo, stackup, grid_size_mm=gs)
            self.log(f"Mesh: {len(mesh.nodes)} nodes, {len(mesh.edges)} edges.")
            
            if len(mesh.nodes) == 0:
                self.log("Error: Mesh has 0 nodes.")
                return

            if debug_mode:
                self._debug_plot_mesh(mesher, mesh, stackup)

            # --- 3. Map Sources/Loads to Nodes ---
            solver_sources = []
            solver_loads = []
            
            grid_step = mesh.grid_step
            origin = mesh.grid_origin
            
            def get_mesh_nodes_for_component_pads(ref_des, pad_names):
                """Helper to find mesh nodes corresponding to specific pads of a component."""
                nodes = []
                fp = self.board.FindFootprintByReference(ref_des)
                if not fp:
                    self.log(f"Warning: Footprint {ref_des} not found.")
                    return []
                
                for pad_name in pad_names:
                    pad = fp.FindPadByNumber(pad_name)
                    if not pad: continue
                    
                    pos = pad.GetPosition()
                    px, py = pcbnew.ToMM(pos.x), pcbnew.ToMM(pos.y)
                    
                    # Convert to grid
                    tx = int(round((px - origin[0]) / grid_step))
                    ty = int(round((py - origin[1]) / grid_step))
                    
                    # Search neighborhood
                    # TODO: Properly use pad bounding box? MVP: 3x3 search
                    found_for_pad = False
                    for dx in [-1, 0, 1]:
                        for dy in [-1, 0, 1]:
                            ix, iy = tx + dx, ty + dy
                            # Check all layers
                            for layer in range(32): # iterate reasonably
                                nid = mesh.node_map.get((ix, iy, layer))
                                if nid is not None:
                                    nodes.append(nid)
                                    found_for_pad = True
                    if not found_for_pad and debug_mode:
                        self.log(f"  Pad {ref_des}.{pad_name} at ({px:.2f},{py:.2f}) not on mesh.")
                        
                return list(set(nodes))

            # Assemble Sources
            for src in rail.sources:
                nodes = get_mesh_nodes_for_component_pads(src.component_ref.ref_des, src.pad_names)
                if not nodes:
                    self.log(f"Warning: Source {src.component_ref.ref_des} has no connected mesh nodes.")
                    continue
                
                # Voltage defined at Rail level
                v_set = rail.nominal_voltage
                for nid in nodes:
                    solver_sources.append({'node_id': nid, 'voltage': v_set})
                    
            # Assemble Loads
            for load in rail.loads:
                nodes = get_mesh_nodes_for_component_pads(load.component_ref.ref_des, load.pad_names)
                if not nodes:
                    self.log(f"Warning: Load {load.component_ref.ref_des} has no connected mesh nodes.")
                    continue
                
                # Distribute current
                # Default UNIFORM: Total / Node Count
                i_per_node = load.total_current / len(nodes)
                for nid in nodes:
                    solver_loads.append({'node_id': nid, 'current': i_per_node})

            # --- 4. Solve ---
            if not solver_sources:
                self.log("Error: No valid source nodes found on mesh. Cannot solve.")
            else:
                self.log("Solving...")
                solver = Solver()
                results = solver.solve(mesh, solver_sources, solver_loads)
                
                # Stats
                v_vals = list(results.values())
               
                v_min = min(v_vals)
                v_max = max(v_vals)
                drop_abs = v_max - v_min
                drop_pct = (drop_abs / rail.nominal_voltage * 100) if rail.nominal_voltage > 0 else 0
                
                res_str = f"Rail: {net_name}\n"
                res_str += f"Nominal: {rail.nominal_voltage} V\n"
                res_str += f"Min V: {v_min:.4f} V\n"
                res_str += f"Max V: {v_max:.4f} V\n"
                res_str += f"Max Drop: {drop_abs:.4f} V ({drop_pct:.2f}%)"
                
                self.result_text.SetValue(res_str)
                self.log("Simulation Success.")
                
                # --- 5. Visualize ---
                mesh.results = results
                bitmap = mesher.debug_plot(mesh, stackup)
                if bitmap:
                    self.result_image.SetBitmap(bitmap)
                    self.scrolled_img.FitInside()
                    self.notebook.SetSelection(1) # Results Tab

        except Exception as e:
            self.log(f"Exception: {e}")
            import traceback
            traceback.print_exc()

    def _debug_plot_geo(self, extractor, geo):
        try:
            buf = extractor.plot_geometry(geo)
            if buf:
                img = wx.Image(buf, wx.BITMAP_TYPE_PNG)
                self.result_image.SetBitmap(wx.Bitmap(img))
                self.scrolled_img.FitInside()
        except: pass

    def _debug_plot_mesh(self, mesher, mesh, stackup):
        try:
            bmp = mesher.debug_plot(mesh, stackup)
            if bmp:
                self.result_image.SetBitmap(bmp)
                self.scrolled_img.FitInside()
        except: pass

    def on_close(self, event):
        self.EndModal(wx.ID_CANCEL)
