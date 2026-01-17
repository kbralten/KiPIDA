
import pcbnew
import sys
import math

try:
    import numpy as np
    from shapely.geometry import Point, box
    from shapely.prepared import prep
    import matplotlib.path
except ImportError:
    np = None
    Point = box = prep = None
    matplotlib = None

class Mesh:
    def __init__(self):
        self.nodes = [] # List of node_ids (integers)
        self.node_coords = {} # { node_id: (x_mm, y_mm, layer_id) }
        self.edges = [] # List of (node_a, node_b, conductance_G) - Legacy support
        self.node_map = {} # { (x_idx, y_idx, layer_id): node_id }
        self.grid_origin = (0, 0)
        self.grid_step = 0
        
        # New sparse matrix components
        self.G_coo_data = [] # [g, g, -g, -g, ...]
        self.G_coo_row = []
        self.G_coo_col = []
        self.G_final_csr = None # To be filled by solver or mesher if configured
        
    def add_edge_direct(self, u, v, g):
        """Adds an edge directly to the sparse data arrays."""
        # G[u,u] += g
        self.G_coo_row.append(u)
        self.G_coo_col.append(u)
        self.G_coo_data.append(g)
        
        # G[v,v] += g
        self.G_coo_row.append(v)
        self.G_coo_col.append(v)
        self.G_coo_data.append(g)
        
        # G[u,v] -= g
        self.G_coo_row.append(u)
        self.G_coo_col.append(v)
        self.G_coo_data.append(-g)
        
        # G[v,u] -= g
        self.G_coo_row.append(v)
        self.G_coo_col.append(u)
        self.G_coo_data.append(-g)

class Mesher:
    def __init__(self, board, debug=False, log_callback=None):
        self.board = board
        self.debug = debug
        self.log_callback = log_callback
        if np is None or Point is None or matplotlib is None:
            raise ImportError("NumPy, Shapely, and Matplotlib are required for Meshing.")
    
    def _log(self, msg):
        """Helper to log debug messages."""
        if self.debug and self.log_callback:
            self.log_callback(f"[MESH] {msg}")

    def generate_mesh(self, net_name, geometry_by_layer, stackup, grid_size_mm=0.5):
        """
        Generates a resistive mesh from the geometry using vectorized operations.
        """
        mesh = Mesh()
        mesh.grid_step = grid_size_mm
        
        # 1. Calculate Bounding Box
        min_x, min_y, max_x, max_y = float('inf'), float('inf'), float('-inf'), float('-inf')
        
        has_geometry = False
        valid_layers = []
        for lid, poly in geometry_by_layer.items():
            if poly.is_empty: continue
            has_geometry = True
            valid_layers.append(lid)
            b = poly.bounds
            min_x = min(min_x, b[0])
            min_y = min(min_y, b[1])
            max_x = max(max_x, b[2])
            max_y = max(max_y, b[3])
            
        if not has_geometry:
            return mesh

        # Pad bounds slightly
        pad = grid_size_mm
        min_x -= pad
        min_y -= pad
        max_x += pad
        max_y += pad
        
        mesh.grid_origin = (min_x, min_y)
        
        # 2. Create Grid Coordinates
        width_mm = max_x - min_x
        height_mm = max_y - min_y
        
        nx = int(math.ceil(width_mm / grid_size_mm))
        ny = int(math.ceil(height_mm / grid_size_mm))
        
        x_coords = np.linspace(min_x, min_x + (nx * grid_size_mm), nx + 1)
        y_coords = np.linspace(min_y, min_y + (ny * grid_size_mm), ny + 1)
        
        xv, yv = np.meshgrid(x_coords, y_coords) # shape (ny+1, nx+1)
        # Flatten for vectorized checks logic
        # But we need structure for neighbor identifying.
        
        grid_points = np.column_stack((xv.ravel(), yv.ravel()))
        
        if self.debug:
            self._log(f"Grid setup: {nx+1}x{ny+1} points, bounds ({min_x:.1f},{min_y:.1f}) to ({max_x:.1f},{max_y:.1f})")

        # 3. Vectorized Rasterization
        # We will build a 3D boolean mask: presence[layer_idx, y_idx, x_idx]
        # But layer IDs are sparse (e.g. 0, 1, 31). So we map them.
        
        sorted_layers = sorted(valid_layers)
        layer_map = { lid: idx for idx, lid in enumerate(sorted_layers) }
        
        # node_indices: [layer, y, x] -> node_id (or -1 if empty)
        node_grid = np.full((len(sorted_layers), ny + 1, nx + 1), -1, dtype=int)
        
        node_counter = 0
        
        for lid in sorted_layers:
            poly = geometry_by_layer[lid]
            if poly.is_empty: continue
            
            # Simple buffering to ensure boundary inclusion - though check 'contains' logic
            # Using matplotlib path for speed
            # Matplotlib Path uses vertices. 
            # If poly is MultiPolygon, iterate parts.
            
            polys_to_check = [poly] if poly.geom_type == 'Polygon' else list(poly.geoms)
            
            # Create a combined boolean mask for this layer
            layer_mask = np.zeros(len(grid_points), dtype=bool)
            
            for p in polys_to_check:
                # Buffer slightly to include points on edges
                pb = p.buffer(1e-5) 
                
                # Extract exterior coords
                
                codes = []
                verts = []
                
                # Exterior
                ext_coords = list(pb.exterior.coords)
                verts.extend(ext_coords)
                codes.append(matplotlib.path.Path.MOVETO)
                codes.extend([matplotlib.path.Path.LINETO] * (len(ext_coords) - 2))
                codes.append(matplotlib.path.Path.CLOSEPOLY)
                
                # Interiors (Holes)
                for interior in pb.interiors:
                    int_coords = list(interior.coords)
                    verts.extend(int_coords)
                    codes.append(matplotlib.path.Path.MOVETO)
                    codes.extend([matplotlib.path.Path.LINETO] * (len(int_coords) - 2))
                    codes.append(matplotlib.path.Path.CLOSEPOLY)
                
                path = matplotlib.path.Path(verts, codes)
                
                # Check points
                # radius=0 means exact point check. Could use small radius for tolerance.
                mask = path.contains_points(grid_points, radius=1e-9)
                layer_mask |= mask
            
            # Reshape back to grid
            mask_2d = layer_mask.reshape((ny + 1, nx + 1))
            
            # Assign Node IDs
            count_on_layer = np.count_nonzero(mask_2d)
            if count_on_layer > 0:
                # Get indices where mask is true
                y_idxs, x_idxs = np.nonzero(mask_2d)
                
                # Generate new IDs
                new_ids = np.arange(node_counter, node_counter + count_on_layer)
                node_grid[layer_map[lid], y_idxs, x_idxs] = new_ids
                
                # Save to mesh.nodes and mesh.node_coords
                
                # For `mesh.nodes` (list of ints)
                mesh.nodes.extend(new_ids)

                
                for i in range(count_on_layer):
                    nid = new_ids[i]
                    xi = x_idxs[i]
                    yi = y_idxs[i]
                    mesh.node_map[(xi, yi, lid)] = nid
                    mesh.node_coords[nid] = (
                        min_x + xi * grid_size_mm,
                        min_y + yi * grid_size_mm,
                        lid
                    )
                
                node_counter += count_on_layer
                
                # 4. Generate Lateral Edges (Vectorized)
                # Physical props
                copper_info = stackup.get('copper', {}).get(lid, {})
                thick = copper_info.get('thickness_mm', 0.035)
                rho = stackup.get('resistivity', 1.7e-5)
                g_lat_val = thick / rho
                
                # Horizontal Neighbors (x, y) <-> (x+1, y)
                # Check where node and right-neighbor both exist
                
                # mask_2d is boolean. node_grid has IDs.
                # Valid H edges: mask[:, :-1] & mask[:, 1:]
                
                # Right neighbors
                right_mask = mask_2d[:, :-1] & mask_2d[:, 1:]
                if np.any(right_mask):
                    y_r, x_r = np.nonzero(right_mask)
                    # Nodes at (y,x)
                    u_ids = node_grid[layer_map[lid], y_r, x_r]
                    # Nodes at (y, x+1)
                    v_ids = node_grid[layer_map[lid], y_r, x_r + 1]
                    
                    self._bulk_add_edges(mesh, u_ids, v_ids, g_lat_val)
                
                # Vertical (Top) Neighbors (x, y) <-> (x, y+1)
                top_mask = mask_2d[:-1, :] & mask_2d[1:, :]
                if np.any(top_mask):
                    y_t, x_t = np.nonzero(top_mask)
                    # Nodes at (y, x)
                    u_ids = node_grid[layer_map[lid], y_t, x_t]
                    # Nodes at (y+1, x)
                    v_ids = node_grid[layer_map[lid], y_t + 1, x_t]
                    
                    self._bulk_add_edges(mesh, u_ids, v_ids, g_lat_val)

            if self.debug:
                self._log(f"  Layer {lid} vectorized mesh: {count_on_layer} nodes.")

        # 5. Vertical Connections (Vias & PTH) - Unchanged logic, but updated to use new map
        
        if self.log_callback:
            self.log_callback("Adding vertical interconnects...")
        net = self.board.FindNet(net_name)
        if net:
            net_code = net.GetNetCode()
            for track in self.board.GetTracks():
                if track.GetNetCode() == net_code and isinstance(track, pcbnew.PCB_VIA):
                    self._add_vertical_link(mesh, track, stackup)
                    
            for fp in self.board.GetFootprints():
                for pad in fp.Pads():
                    if pad.GetNetCode() == net_code and pad.GetAttribute() == pcbnew.PAD_ATTRIB_PTH:
                        self._add_vertical_stack(mesh, pad.GetPosition(), layers=None, diameter=pcbnew.ToMM(pad.GetDrillSize().x), stackup=stackup)

        return mesh

    def _bulk_add_edges(self, mesh, u_ids, v_ids, g):
        """Adds multiple edges at once to sparse arrays."""
        # This is where we gain massive speed in construction
        n = len(u_ids)
        # We need to replicate g for all edges
        gs = np.full(n, g)
        neg_gs = np.full(n, -g)
        
        # Prepare arrays
        rows = []
        cols = []
        data = []
        
        # u,u
        rows.append(u_ids); cols.append(u_ids); data.append(gs)
        # v,v
        rows.append(v_ids); cols.append(v_ids); data.append(gs)
        # u,v
        rows.append(u_ids); cols.append(v_ids); data.append(neg_gs)
        # v,u
        rows.append(v_ids); cols.append(u_ids); data.append(neg_gs)
        
        # Concatenate and extend
        mesh.G_coo_row.extend(np.concatenate(rows))
        mesh.G_coo_col.extend(np.concatenate(cols))
        mesh.G_coo_data.extend(np.concatenate(data))

    def _calculate_vertical_g(self, layer_a, layer_b, stackup, diameter_mm):
        # ... logic unchanged from previous implementation ...
        if layer_a == layer_b: return 1.0e9
        plating_thick = 0.025 
        area = math.pi * (diameter_mm * plating_thick - plating_thick**2)
        if area <= 0: return 1000.0 
        rho = stackup.get('resistivity', 1.68e-5)
        h = 0
        l_min, l_max = min(layer_a, layer_b), max(layer_a, layer_b)
        for sub in stackup.get('substrate', []):
            sb = sub['between']
            if sb[0] is not None and sb[1] is not None:
                if min(sb) >= l_min and max(sb) <= l_max:
                    h += sub['thickness_mm']
        if h <= 0: h = 0.5 
        return area / (rho * h)

    def _add_vertical_link(self, mesh, via, stackup):
        layers = sorted(list(via.GetLayerSet().Seq()))
        if not layers: return
        start_id, end_id = min(layers), max(layers)
        pos = via.GetPosition()
        x_mm, y_mm = pcbnew.ToMM(pos.x), pcbnew.ToMM(pos.y)
        ix = int(round((x_mm - mesh.grid_origin[0]) / mesh.grid_step))
        iy = int(round((y_mm - mesh.grid_origin[1]) / mesh.grid_step))
        
        all_cu_layers = sorted(stackup['copper'].keys())
        nodes_in_stack = []
        for lid in all_cu_layers:
            if start_id <= lid <= end_id:
                nid = mesh.node_map.get((ix, iy, lid))
                if nid is not None:
                    nodes_in_stack.append(nid)
        
        drill_dia = pcbnew.ToMM(via.GetDrillValue())
        for i in range(len(nodes_in_stack) - 1):
            nid_a = nodes_in_stack[i]
            nid_b = nodes_in_stack[i+1]
            la = mesh.node_coords[nid_a][2]
            lb = mesh.node_coords[nid_b][2]
            g_via = self._calculate_vertical_g(la, lb, stackup, drill_dia)
            mesh.add_edge_direct(nid_a, nid_b, g_via)

    def _add_vertical_stack(self, mesh, pos, layers, diameter, stackup):
        if layers is None:
            layers = sorted(stackup['copper'].keys())
        x_mm, y_mm = pcbnew.ToMM(pos.x), pcbnew.ToMM(pos.y)
        ix = int(round((x_mm - mesh.grid_origin[0]) / mesh.grid_step))
        iy = int(round((y_mm - mesh.grid_origin[1]) / mesh.grid_step))
        
        nodes_in_stack = []
        for layer in layers:
            nid = mesh.node_map.get((ix, iy, layer))
            if nid is not None:
                nodes_in_stack.append(nid)
                
        for i in range(len(nodes_in_stack) - 1):
            nid_a = nodes_in_stack[i]
            nid_b = nodes_in_stack[i+1]
            la = mesh.node_coords[nid_a][2]
            lb = mesh.node_coords[nid_b][2]
            g_via = self._calculate_vertical_g(la, lb, stackup, diameter)
            mesh.add_edge_direct(nid_a, nid_b, g_via)

    def debug_plot(self, mesh, stackup=None):
        # Unchanged legacy method, kept for UI visualization
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            import io
            import wx
            
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            
            xs, ys, zs, c = [], [], [], []
            has_results = hasattr(mesh, 'results') and mesh.results
            
            layer_to_z = {}
            if stackup and 'layer_order' in stackup:
                for idx, layer_id in enumerate(stackup['layer_order']):
                    layer_to_z[layer_id] = 10.0 - idx
            else:
                 copper_layers = sorted(stackup['copper'].keys()) if stackup else []
                 for idx, layer_id in enumerate(copper_layers):
                    layer_to_z[layer_id] = 10.0 - idx
            
            for nid, (x, y, layer) in mesh.node_coords.items():
                xs.append(x)
                ys.append(y)
                zs.append(layer_to_z.get(layer, 10 - layer * 0.5))
                c.append(mesh.results.get(nid, 0.0) if has_results else layer)
                
            sc = ax.scatter(xs, ys, zs, c=c, cmap='viridis')
            if has_results:
                plt.colorbar(sc, label='Voltage (V)')
            
            ax.set_xlabel('X (mm)'); ax.set_ylabel('Y (mm)'); ax.set_zlabel('L (pseudo)')
            ax.invert_yaxis()
            
            # Set equal aspect ratio for X and Y to show correct geometry proportions
            x_limits = ax.get_xlim3d()
            y_limits = ax.get_ylim3d()
            
            x_range = x_limits[1] - x_limits[0]
            y_range = y_limits[1] - y_limits[0]
            max_range = max(x_range, y_range)
            
            x_mid = (x_limits[0] + x_limits[1]) / 2.0
            y_mid = (y_limits[0] + y_limits[1]) / 2.0
            
            # Force limits to be the same size
            ax.set_xlim3d([x_mid - max_range/2, x_mid + max_range/2])
            ax.set_ylim3d([y_mid - max_range/2, y_mid + max_range/2])
            
            # Force consistent tick intervals on both axes
            import matplotlib.ticker as ticker
            locator = ticker.MaxNLocator(nbins=6)
            tick_values = locator.tick_values(x_mid - max_range/2, x_mid + max_range/2)
            if len(tick_values) > 1:
                step = tick_values[1] - tick_values[0]
                ax.xaxis.set_major_locator(ticker.MultipleLocator(step))
                ax.yaxis.set_major_locator(ticker.MultipleLocator(step))
            
            import os
            filename = 'debug_mesh_solved.png' if has_results else 'debug_mesh_presolved.png'
            output_path = os.path.join(os.path.dirname(__file__), filename)
            plt.savefig(output_path, format='png', dpi=150, bbox_inches='tight')
            
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=100)
            plt.close(fig)
            buf.seek(0)
            image = wx.Image(buf, wx.BITMAP_TYPE_PNG)
            return wx.Bitmap(image)
        except Exception:
            return None
