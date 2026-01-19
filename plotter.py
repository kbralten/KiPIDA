
import matplotlib
# Use Agg backend to avoid GUI requirement for matplotlib, since we just want images
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import io
import wx
import numpy as np

class Plotter:
    def __init__(self, debug=False):
        self.debug = debug

    def plot_3d_mesh(self, mesh, stackup=None):
        """
        Generates a 3D scatter plot of the mesh nodes.
        Returns a wx.Bitmap.
        """
        try:
            fig = plt.figure(figsize=(7, 5), constrained_layout=True)
            ax = fig.add_subplot(111, projection='3d')
            
            xs, ys, zs, c = [], [], [], []
            has_results = hasattr(mesh, 'results') and mesh.results
            
            layer_to_z = {}
            if stackup and 'layer_order' in stackup:
                for idx, layer_id in enumerate(stackup['layer_order']):
                    layer_to_z[layer_id] = 10.0 - idx
            else:
                copper_layers = sorted(stackup['copper'].keys()) if stackup and 'copper' in stackup else []
                for idx, layer_id in enumerate(copper_layers):
                    layer_to_z[layer_id] = 10.0 - idx
            
            # If no node_coords populated yet (edge case), return None
            if not mesh.node_coords:
                plt.close(fig)
                return None

            for nid, (x, y, layer) in mesh.node_coords.items():
                xs.append(x)
                ys.append(-y)  # Invert Y to match KiCad
                zs.append(layer_to_z.get(layer, 10 - layer * 0.5))
                c.append(mesh.results.get(nid, 0.0) if has_results else layer)
                
            sc = ax.scatter(xs, ys, zs, c=c, cmap='viridis')
            if has_results:
                plt.colorbar(sc, label='Voltage (V)', shrink=0.8)
            
            ax.set_xlabel('X (mm)'); ax.set_ylabel('Y (mm)'); ax.set_zlabel('L (pseudo)')
            
            # Equal aspect ratio
            x_limits = ax.get_xlim3d()
            y_limits = ax.get_ylim3d()
            x_range = x_limits[1] - x_limits[0]
            y_range = y_limits[1] - y_limits[0]
            max_range = max(x_range, y_range)
            x_mid = (x_limits[0] + x_limits[1]) / 2.0
            y_mid = (y_limits[0] + y_limits[1]) / 2.0
            ax.set_xlim3d([x_mid - max_range/2, x_mid + max_range/2])
            ax.set_ylim3d([y_mid - max_range/2, y_mid + max_range/2])
            
            return self._fig_to_bitmap(fig)
        except Exception as e:
            if self.debug: print(f"Plotter 3D Error: {e}")
            return None

    def plot_layer_2d(self, mesh, layer_id, stackup=None, vmin=None, vmax=None):
        """
        Generates a 2D plot (heatmap) for a specific layer.
        Returns a wx.Bitmap.
        """
        try:
            fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
            
            xs, ys, vs = [], [], []
            has_results = hasattr(mesh, 'results') and mesh.results
            
            # Filter nodes for this layer
            nodes_on_layer = [nid for nid in mesh.nodes if mesh.node_coords[nid][2] == layer_id]
            
            if not nodes_on_layer:
                plt.close(fig)
                return None

            for nid in nodes_on_layer:
                coords = mesh.node_coords[nid]
                xs.append(coords[0])
                ys.append(-coords[1]) # Invert Y
                val = mesh.results.get(nid, 0.0) if has_results else 0.0
                vs.append(val)
                
            if not xs:
                plt.close(fig)
                return None

            # Scatter plot for now - tripcolor or imshow is better if we have regular grid, 
            # but scatter is robust for sparse nodes.
            # Using a fixed marker size might be tricky, let's try a reasonable default.
            # Ideally s should relate to grid_size, but scatter size is in points^2.
            # Let's just use a standard size for visibility.
            sc = ax.scatter(xs, ys, c=vs, cmap='viridis', vmin=vmin, vmax=vmax, s=20)
            
            if has_results:
                plt.colorbar(sc, label='Voltage (V)')
            
            layer_name = str(layer_id)
            if stackup and 'copper' in stackup and layer_id in stackup['copper']:
                 # Try to get layer name? currently stackup dict structure in test is simple
                 pass

            ax.set_title(f"Layer {layer_name}")
            ax.set_xlabel('X (mm)')
            ax.set_ylabel('Y (mm)')
            ax.set_aspect('equal', 'box')
            
            return self._fig_to_bitmap(fig)

        except Exception as e:
            if self.debug: print(f"Plotter 2D Error: {e}")
            return None

    def _fig_to_bitmap(self, fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=100)
        plt.close(fig)
        buf.seek(0)
        image = wx.Image(buf, wx.BITMAP_TYPE_PNG)
        return wx.Bitmap(image)
