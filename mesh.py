
import pcbnew
import sys
import math

try:
    import numpy as np
    from shapely.geometry import Point, box
    from shapely.prepared import prep
except ImportError:
    np = None
    Point = box = prep = None

class Mesh:
    def __init__(self):
        self.nodes = [] # List of node_ids (integers)
        self.node_coords = {} # { node_id: (x_mm, y_mm, layer_id) }
        self.edges = [] # List of (node_a, node_b, conductance_G)
        self.node_map = {} # { (x_idx, y_idx, layer_id): node_id }
        self.grid_origin = (0, 0)
        self.grid_step = 0

class Mesher:
    def __init__(self, board):
        self.board = board
        if np is None or Point is None:
            raise ImportError("NumPy and Shapely are required for Meshing.")

    def generate_mesh(self, net_name, geometry_by_layer, stackup, grid_size_mm=0.5):
        """
        Generates a resistive mesh from the geometry.
        
        Args:
            net_name (str): The net being meshed.
            geometry_by_layer (dict): { layer_id: shapely.Polygon }
            stackup (dict): { layer_id: {thickness, resistivity} }
            grid_size_mm (float): Mesh resolution.
            
        Returns:
            Mesh object.
        """
        mesh = Mesh()
        mesh.grid_step = grid_size_mm
        
        # 1. Calculate Bounding Box
        min_x, min_y, max_x, max_y = float('inf'), float('inf'), float('-inf'), float('-inf')
        
        has_geometry = False
        for poly in geometry_by_layer.values():
            if poly.is_empty: continue
            has_geometry = True
            b = poly.bounds
            min_x = min(min_x, b[0])
            min_y = min(min_y, b[1])
            max_x = max(max_x, b[2])
            max_y = max(max_y, b[3])
            
        if not has_geometry:
            return mesh

        # Pad bounds slightly
        min_x -= grid_size_mm
        min_y -= grid_size_mm
        max_x += grid_size_mm
        max_y += grid_size_mm
        
        mesh.grid_origin = (min_x, min_y)
        
        # 2. Create Grid Arrays
        width_mm = max_x - min_x
        height_mm = max_y - min_y
        
        nx = int(math.ceil(width_mm / grid_size_mm))
        ny = int(math.ceil(height_mm / grid_size_mm))
        
        x_coords = np.linspace(min_x, min_x + (nx * grid_size_mm), nx + 1)
        y_coords = np.linspace(min_y, min_y + (ny * grid_size_mm), ny + 1)
        
        # 3. Rasterization (Lateral Mesh)
        node_counter = 0
        
        for layer_id, poly in geometry_by_layer.items():
            if poly.is_empty: continue
            
            buffered_poly = poly.buffer(1e-5)
            
            # Get physical props
            copper_info = stackup.get('copper', {}).get(layer_id, {})
            thick = copper_info.get('thickness_mm', 0.035)
            rho = stackup.get('resistivity', 1.7e-5)
            
            g_lat_val = thick / rho
            
            layer_nodes = [] # List of (x_idx, y_idx, node_id)
            
            # Let's try to be somewhat efficient: only check points in poly bounds
            pb = poly.bounds
            ix_min = int((pb[0] - min_x) / grid_size_mm)
            ix_max = int((pb[2] - min_x) / grid_size_mm) + 2
            iy_min = int((pb[1] - min_y) / grid_size_mm)
            iy_max = int((pb[3] - min_y) / grid_size_mm) + 2
            
            for ix in range(max(0, ix_min), min(nx+1, ix_max)):
                for iy in range(max(0, iy_min), min(ny+1, iy_max)):
                    px = min_x + ix * grid_size_mm
                    py = min_y + iy * grid_size_mm
                    
                    # Direct interaction with buffered polygon
                    if buffered_poly.intersects(Point(px, py)):
                        # Found a node
                        nid = node_counter
                        node_counter += 1
                        
                        mesh.nodes.append(nid)
                        mesh.node_coords[nid] = (px, py, layer_id)
                        mesh.node_map[(ix, iy, layer_id)] = nid
                        layer_nodes.append((ix, iy, nid))

            # Connect Lateral Neighbors
            for ix, iy, nid in layer_nodes:
                # Right Neighbor
                nid_right = mesh.node_map.get((ix + 1, iy, layer_id))
                if nid_right is not None:
                    mesh.edges.append((nid, nid_right, g_lat_val))
                    
                # Top Neighbor
                nid_top = mesh.node_map.get((ix, iy + 1, layer_id))
                if nid_top is not None:
                    mesh.edges.append((nid, nid_top, g_lat_val))

        # 4. Vertical Connections (Vias & PTH)
        net = self.board.FindNet(net_name)
        if not net:
            return mesh
        net_code = net.GetNetCode()
        
        # 4a. Vias
        for track in self.board.GetTracks():
            if track.GetNetCode() == net_code and isinstance(track, pcbnew.PCB_VIA):
                self._add_vertical_link(mesh, track, stackup)

        # 4b. Plated Through Holes (Pads)
        for fp in self.board.GetFootprints():
            for pad in fp.Pads():
                # Check Net
                if pad.GetNetCode() != net_code:
                    continue
                
                # Check Attribute: PTH
                # pcbnew.PAD_ATTRIB_PTH is the enum, usually 0. 
                # PAD_ATTRIB_SMD = 1, CONN = 2, NPTH = 3.
                # Use getattr safely for cross-version compat if needed, but standard is PAD_ATTRIB_PTH
                if pad.GetAttribute() == pcbnew.PAD_ATTRIB_PTH:
                    # PTH connects F_Cu to B_Cu (and all in between usually)
                    # We need to know specific layer span? 
                    # Usually PTH goes through all layers.
                    # Or check pad.GetLayerSet().
                    
                    # Assume all enabled copper layers for PTH
                    # We need to link adjacent stacked layers.
                    self._add_vertical_stack(mesh, pad.GetPosition(), layers=None, diameter=pcbnew.ToMM(pad.GetDrillSize().x), stackup=stackup)

        return mesh
    def _calculate_vertical_g(self, layer_a, layer_b, stackup, diameter_mm):
        """
        Calculates conductance between two copper layers based on via/PTH physical properties.
        $G = A / (rho * L)$ where $A = pi * d_avg * t_plating$
        """
        if layer_a == layer_b: return 1.0e9 # short
        
        # 1. Constant plating thickness (typically 20-25um)
        plating_thick = 0.025 
        
        # 2. Cross sectional area of the tube
        # A = pi * ( (d/2)^2 - (d/2 - w)^2 )
        # A = pi * ( d_outer*w - w^2 ) 
        area = math.pi * (diameter_mm * plating_thick - plating_thick**2)
        if area <= 0: return 1000.0 # fallback
        
        # 3. Resistivity
        rho = stackup.get('resistivity', 1.68e-5)
        
        # 4. Total height h (sum of dielectric segments between these layers)
        h = 0
        l_min, l_max = min(layer_a, layer_b), max(layer_a, layer_b)
        
        for sub in stackup.get('substrate', []):
            sb = sub['between']
            if sb[0] is not None and sb[1] is not None:
                if min(sb) >= l_min and max(sb) <= l_max:
                    h += sub['thickness_mm']
        
        if h <= 0: h = 0.5 # fallback dielectric 0.5mm
        
        # G = A / (rho * h)
        return area / (rho * h)

    def _add_vertical_link(self, mesh, via, stackup):
        """Links nodes for a via between its start/end layers."""
        layers = sorted(list(via.GetLayerSet().Seq()))
        if not layers: return
        
        start_id = min(layers)
        end_id = max(layers)
        
        pos = via.GetPosition()
        x_mm = pcbnew.ToMM(pos.x)
        y_mm = pcbnew.ToMM(pos.y)
        
        ix = int(round((x_mm - mesh.grid_origin[0]) / mesh.grid_step))
        iy = int(round((y_mm - mesh.grid_origin[1]) / mesh.grid_step))
        
        all_cu_layers = sorted(stackup['copper'].keys())
        
        # Identify nodes at this (ix, iy) on these layers
        nodes_in_stack = []
        for lid in all_cu_layers:
            if start_id <= lid <= end_id:
                nid = mesh.node_map.get((ix, iy, lid))
                if nid is not None:
                    nodes_in_stack.append(nid)
        
        # Link them sequentially
        drill_dia = pcbnew.ToMM(via.GetDrillValue())
        
        for i in range(len(nodes_in_stack) - 1):
            nid_a = nodes_in_stack[i]
            nid_b = nodes_in_stack[i+1]
            la = mesh.node_coords[nid_a][2]
            lb = mesh.node_coords[nid_b][2]
            g_via = self._calculate_vertical_g(la, lb, stackup, drill_dia)
            mesh.edges.append((nid_a, nid_b, g_via))

    def _add_vertical_stack(self, mesh, pos, layers, diameter, stackup):
        # Similar to _add_vertical_link but for PTH (all layers)
        if layers is None:
            layers = sorted(stackup['copper'].keys())
            
        x_mm = pcbnew.ToMM(pos.x)
        y_mm = pcbnew.ToMM(pos.y)
        
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
            mesh.edges.append((nid_a, nid_b, g_via))

    def debug_plot(self, mesh):
        """
        Generates a 3D plot of the mesh and returns it as a wx.Bitmap.
        """
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            import io
            import wx
            
            # Plot nodes
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
            
            xs = []
            ys = []
            zs = []
            c = []
            
            has_results = hasattr(mesh, 'results') and mesh.results
            
            for nid, (x, y, layer) in mesh.node_coords.items():
                xs.append(x)
                ys.append(y)
                # Map layer ID to Z height pseudo
                if layer == 0: z = 10
                elif layer == 31: z = 0
                else: z = 5
                zs.append(z)
                
                if has_results:
                    val = mesh.results.get(nid, 0.0)
                    c.append(val)
                else:
                    c.append(layer)
                
            sc = ax.scatter(xs, ys, zs, c=c, cmap='viridis')
            
            if has_results:
                plt.colorbar(sc, label='Voltage (V)')
                ax.set_title(f"Voltage Distribution ({min(c):.3f}V - {max(c):.3f}V)")
            else:
                ax.set_title("Mesh Nodes (Color=Layer)")
                
            ax.set_xlabel('X (mm)')
            ax.set_ylabel('Y (mm)')
            ax.set_zlabel('Layer (pseudo)')
            ax.invert_yaxis()
            
            # Save to memory buffer
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=100)
            plt.close(fig) # Close figure to free memory
            buf.seek(0)
            
            # Convert to wx.Image then wx.Bitmap
            image = wx.Image(buf, wx.BITMAP_TYPE_PNG)
            return wx.Bitmap(image)
            
        except Exception as e:
            print(f"Plotting failed: {e}")
            return None
