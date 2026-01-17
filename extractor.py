import pcbnew
import logging
try:
    from shapely.geometry import LineString, Polygon, MultiPolygon, Point, box
    from shapely.ops import unary_union
except ImportError:
    LineString = Polygon = MultiPolygon = Point = box = unary_union = None

class GeometryExtractor:
    def __init__(self, board, debug=False, log_callback=None):
        self.board = board
        self.debug = debug
        self.log_callback = log_callback
        self.logger = logging.getLogger('KiPIDA.GeometryExtractor')
        if debug:
            self.logger.setLevel(logging.DEBUG)
            # Clear any existing handlers to avoid duplicates
            self.logger.handlers.clear()
            if log_callback:
                # Use custom handler that calls the UI log callback
                class UILogHandler(logging.Handler):
                    def __init__(self, callback):
                        super().__init__()
                        self.callback = callback
                    def emit(self, record):
                        msg = self.format(record)
                        self.callback(msg)
                
                handler = UILogHandler(log_callback)
                handler.setFormatter(logging.Formatter('[%(levelname)s] %(message)s'))
                self.logger.addHandler(handler)
            else:
                # Fallback to stderr
                handler = logging.StreamHandler()
                handler.setFormatter(logging.Formatter('[%(name)s] %(levelname)s: %(message)s'))
                self.logger.addHandler(handler)
        if LineString is None:
            raise ImportError("The 'shapely' library is required for Ki-PIDA but was not found.")
        
        # Detect KiCad version for API compatibility
        self._detect_kicad_version()
    
    def _detect_kicad_version(self):
        """Detect KiCad version for API compatibility."""
        try:
            version_str = pcbnew.GetBuildVersion()
            self.kicad_version = version_str
            if self.debug:
                self.logger.debug(f"KiCad version detected: {version_str}")
        except:
            self.kicad_version = "unknown"
            if self.debug:
                self.logger.debug("Could not detect KiCad version")

    def get_board_stackup(self):
        """
        Extracts physical stackup properties.
        Returns a dict with copper data and substrate segments.
        """
        rho_copper = 1.68e-5 # Ohm-mm
        copper_data = {}
        layer_order = []  # Preserve physical layer order top-to-bottom
        substrates = []
        
        try:
            settings = self.board.GetDesignSettings()
            if self.debug:
                self.logger.debug(f"DesignSettings type: {type(settings)}")
                if hasattr(settings, 'm_HasStackup'):
                    self.logger.debug(f"m_HasStackup = {settings.m_HasStackup}")
            
            stackup = settings.GetStackupDescriptor()
            if self.debug:
                self.logger.debug(f"StackupDescriptor type: {type(stackup)}")
                # Show ALL methods (not just Stackup-related)
                all_methods = [m for m in dir(stackup) if not m.startswith('_')]
                self.logger.debug(f"All stackup methods: {all_methods}")
                
                # Try various size methods
                for method_name in ['size', 'Size', 'count', 'Count', 'length', 'GetCount']:
                    if hasattr(stackup, method_name):
                        try:
                            result = getattr(stackup, method_name)()
                            self.logger.debug(f"stackup.{method_name}() = {result}")
                        except Exception as e:
                            self.logger.debug(f"stackup.{method_name}() failed: {e}")
                
                self.logger.debug("Extracting board stackup from GetStackupDescriptor()...")
            
            last_copper = None
            # Use GetList() to get stackup items
            try:
                try:
                    items_list = stackup.GetList()
                    if self.debug:
                        self.logger.debug(f"stackup.GetList() returned: {type(items_list)}, len={len(items_list)}")
                except AttributeError as e:
                    if self.debug:
                        self.logger.warning(f"GetList() not available: {e}")
                    raise e
                    
                for item in items_list:
                    lid = item.GetLayerId()
                    thickness = pcbnew.ToMM(item.GetThickness())
                    
                    if pcbnew.IsCopperLayer(lid):
                        layer_name = self.board.GetLayerName(lid)
                        copper_data[lid] = {
                            'name': layer_name,
                            'thickness_mm': thickness if thickness > 0 else 0.035
                        }
                        layer_order.append(lid)
                        if self.debug:
                            self.logger.debug(f"  Copper layer {lid} ({layer_name}): thickness={copper_data[lid]['thickness_mm']:.4f}mm")
                        last_copper = lid
                    elif lid == -1:  # substrate/dielectric
                        if last_copper is not None:
                            if len(substrates) == 0 or substrates[-1]['between'][1] is not None:
                                substrates.append({'thickness_mm': thickness, 'between': [last_copper, None]})
                                if self.debug:
                                    self.logger.debug(f"  Substrate after layer {last_copper}: thickness={thickness:.4f}mm")
                            else:
                                substrates[-1]['thickness_mm'] += thickness
                                if self.debug:
                                    self.logger.debug(f"  Added {thickness:.4f}mm to previous substrate, total={substrates[-1]['thickness_mm']:.4f}mm")
                
                # Link substrates
                for sub in substrates:
                    start_lid = sub['between'][0]
                    found_start = False
                    for item in items_list:
                        lid = item.GetLayerId()
                        if lid == start_lid:
                            found_start = True
                        elif found_start and pcbnew.IsCopperLayer(lid):
                            sub['between'][1] = lid
                            if self.debug:
                                self.logger.debug(f"  Substrate links layer {start_lid} to {lid}")
                            break
            except Exception as list_error:
                if self.debug:
                    self.logger.warning(f"GetList() extraction failed: {list_error}, using fallback")
                raise list_error
                
        except Exception as e:
            # Fallback: Just get copper layers and assume default dielectric
            print(f"Detailed stackup extraction failed: {e}. Using fallbacks.")
            seq = self.board.GetEnabledLayers()
            enabled_layers = sorted(list(seq.Seq()))
            
            last_lid = None
            for lid in enabled_layers:
                if pcbnew.IsCopperLayer(lid):
                    # Try to get individual thickness if board.GetLayerThickness exists (KiCad 7+)
                    try:
                        thickness = pcbnew.ToMM(self.board.GetLayerThickness(lid))
                    except:
                        thickness = 0.035
                        
                    copper_data[lid] = {
                        'name': self.board.GetLayerName(lid),
                        'thickness_mm': thickness if thickness > 0 else 0.035
                    }
                    layer_order.append(lid)  # Add to layer_order in fallback path too
                    if last_lid is not None:
                        # Add a default 1.0mm dielectric between layers
                        substrates.append({'thickness_mm': 1.0, 'between': [last_lid, lid]})
                    last_lid = lid
        
        # Ensure back copper is at the end of layer_order
        # Layer 31 is standard B.Cu, but some boards use other IDs like layer 2
        if layer_order:
            # Find back copper layer (highest physical layer)
            back_copper = None
            for lid in layer_order:
                layer_name = copper_data.get(lid, {}).get('name', '')
                if 'B.Cu' in layer_name or 'Back' in layer_name:
                    back_copper = lid
                    break
            
            # If found and not already last, move to end
            if back_copper and layer_order[-1] != back_copper:
                layer_order.remove(back_copper)
                layer_order.append(back_copper)

        return {
            'copper': copper_data,
            'layer_order': layer_order,
            'substrate': substrates,
            'resistivity': rho_copper
        }

    def get_net_geometry(self, net_name):
        """
        Extracts and merges geometry for a specific net.
        Returns a dictionary: { layer_id: shapely.geometry.Polygon }
        """
        net = self.board.FindNet(net_name)
        if net is None:
            return {}
        
        net_code = net.GetNetCode()
        layer_shapes = {} # { layer_id: [shapely_objects] }

        def add_shape(layer, shape):
            if layer not in layer_shapes:
                layer_shapes[layer] = []
            layer_shapes[layer].append(shape)

        # 1. Process Tracks
        track_count = 0
        for track in self.board.GetTracks():
            if track.GetNetCode() != net_code:
                continue
            
            # Tracks and Vias are in the same list. 
            # Vias usually have layer = Undefined or encompass multple.
            # We only care about conductive traces on specific layers here.
            # Vias are handled separately in mesh generation (usually). 
            # NOTE: But for geometry, a via pad is effectively a circle on all connected layers.
            
            if isinstance(track, pcbnew.PCB_VIA):
                # Add via annular rings to all connected layers
                via = track
                top_layer = via.TopLayer()
                bottom_layer = via.BottomLayer()
                # Get all copper layers between top and bottom
                # We can iterate our stackup_data keys
                
                # Check layers
                radius_mm = pcbnew.ToMM(via.GetWidth()) / 2.0
                pos = via.GetPosition()
                x_mm = pcbnew.ToMM(pos.x)
                y_mm = pcbnew.ToMM(pos.y)
                via_shape = Point(x_mm, y_mm).buffer(radius_mm)
                
                # Add via copper to all layers it touches
                # Get layer set and iterate through it
                via_count = 0
                layer_set = via.GetLayerSet()
                for layer_id in layer_set.Seq():
                    if pcbnew.IsCopperLayer(layer_id):
                        add_shape(layer_id, via_shape)
                        via_count += 1
                
                if self.debug and via_count > 0:
                    self.logger.debug(f"Via at ({x_mm:.2f},{y_mm:.2f}), radius={radius_mm:.2f}mm, added to {via_count} layers")
                
            elif track.GetClass() == "PCB_TRACK" or track.GetClass() == "PCB_ARC":
                # Standard Track
                track_count += 1
                layer = track.GetLayer()
                width_mm = pcbnew.ToMM(track.GetWidth())
                
                start = track.GetStart()
                end = track.GetEnd()
                
                p0 = (pcbnew.ToMM(start.x), pcbnew.ToMM(start.y))
                p1 = (pcbnew.ToMM(end.x), pcbnew.ToMM(end.y))
                
                if self.debug and track_count <= 5:  # Log first 5 tracks
                    length = ((p1[0]-p0[0])**2 + (p1[1]-p0[1])**2)**0.5
                    self.logger.debug(f"Track #{track_count}: layer={layer}, width={width_mm:.2f}mm, length={length:.2f}mm, from ({p0[0]:.1f},{p0[1]:.1f}) to ({p1[0]:.1f},{p1[1]:.1f})")
                
                # Create LineString and buffer it to get a Polygon (track width)
                # Cap style 1 is ROUND
                line = LineString([p0, p1])
                poly = line.buffer(width_mm / 2, cap_style=1)
                
                add_shape(layer, poly)
        
        if self.debug and track_count > 0:
            self.logger.debug(f"Processed {track_count} track(s) for net '{net_name}'")

        # 2. Process Pads
        # Pads are part of footprints
        for footprint in self.board.GetFootprints():
            for pad in footprint.Pads():
                if pad.GetNetCode() != net_code:
                    continue
                
                # Pads can be on F_Cu, B_Cu, or all (Through Hole)
                layers_to_process = []
                # API check for pad layers
                # Setup a list of target layers for this pad
                pad_layers = pad.GetLayerSet()
                
                # Iterate all copper layers to see if pad is on them
                stackup = self.get_board_stackup()
                for layer_id in stackup['copper'].keys():
                    if pad_layers.Contains(layer_id):
                        layers_to_process.append(layer_id)
                
                if not layers_to_process:
                    continue
                    
                # Create Shape
                pos = pad.GetPosition()
                x_mm = pcbnew.ToMM(pos.x)
                y_mm = pcbnew.ToMM(pos.y)
                size = pad.GetSize()
                w_mm = pcbnew.ToMM(size.x)
                h_mm = pcbnew.ToMM(size.y)
                shape_type = pad.GetShape()
                rot_deg = pad.GetOrientation().AsDegrees()
                
                pad_poly = None
                
                if shape_type == pcbnew.PAD_SHAPE_CIRCLE:
                    pad_poly = Point(x_mm, y_mm).buffer(w_mm / 2)
                    
                elif shape_type == pcbnew.PAD_SHAPE_RECT:
                    # Create box centered at 0,0 then translate/rotate
                    minx, miny = -w_mm/2, -h_mm/2
                    maxx, maxy = w_mm/2, h_mm/2
                    base_box = box(minx, miny, maxx, maxy)
                    # Rotate and translate provided by Shapely affinity is one way, 
                    # but check if we need to do it manually or if pad.GetPosition() is center.
                    # Yes, pad position is center.
                    from shapely import affinity
                    rotated = affinity.rotate(base_box, -rot_deg, origin=(0,0)) # KiCad rotation is CCW? checking..
                    pad_poly = affinity.translate(rotated, x_mm, y_mm)
                    
                elif shape_type == pcbnew.PAD_SHAPE_ROUNDRECT:
                    # Approximate with box buffer or similar
                    # For now treat as RECT
                    minx, miny = -w_mm/2, -h_mm/2
                    maxx, maxy = w_mm/2, h_mm/2
                    base_box = box(minx, miny, maxx, maxy)
                    from shapely import affinity
                    rotated = affinity.rotate(base_box, -rot_deg, origin=(0,0))
                    pad_poly = affinity.translate(rotated, x_mm, y_mm)

                elif shape_type == pcbnew.PAD_SHAPE_OVAL:
                    # Segment with rounded ends
                    # If w > h, horizontal oval. length = w - h. radius = h/2.
                    # Simplified: Treat as point buffer if close to circle, or thick line.
                    # Implementing as a thick line segment
                    # For now, fallback to Rect approximation for simplicity or Skip if complex
                    # Just use the bounding box logic for MVP
                    minx, miny = -w_mm/2, -h_mm/2
                    maxx, maxy = w_mm/2, h_mm/2
                    base_box = box(minx, miny, maxx, maxy)
                    from shapely import affinity
                    rotated = affinity.rotate(base_box, -rot_deg, origin=(0,0))
                    pad_poly = affinity.translate(rotated, x_mm, y_mm)
                
                if pad_poly:
                    for layer in layers_to_process:
                        add_shape(layer, pad_poly)

        # 3. Process Zones
        zone_count = 0
        for zone in self.board.Zones():
            if zone.GetNetCode() != net_code:
                continue
            if not zone.GetIsRuleArea(): # Ensure it's a copper zone
                zone_count += 1
                
                # Check if zone is filled
                is_filled = zone.IsFilled()
                if self.debug:
                    self.logger.debug(f"Processing zone #{zone_count} for net '{net_name}', Filled: {is_filled}")
                
                if not is_filled:
                    if self.debug:
                        self.logger.warning(f"  Zone #{zone_count} is NOT filled - may result in missing geometry!")
                    # Continue anyway to see what we can extract
                
                # Multi-layer handling
                layers_to_process = []
                layer = zone.GetLayer()
                
                if layer < 0:
                    # Multi-layer zone
                    try:
                        lset = zone.GetLayerSet()
                        for id in lset.Seq():
                            layers_to_process.append(id)
                        if self.debug:
                            self.logger.debug(f"  Multi-layer zone, layers: {layers_to_process}")
                    except Exception as e:
                        if self.debug:
                            self.logger.warning(f"  Failed to get layer set: {e}")
                else:
                    layers_to_process.append(layer)
                    if self.debug:
                        self.logger.debug(f"  Single-layer zone on layer {layer}")
                    
                for target_layer in layers_to_process:
                    poly_set = None
                    api_method_used = None
                    
                    # Try multiple API methods for getting filled polys
                    try:
                        # Method 1: GetFilledPolysList with layer parameter (KiCad 9.0+)
                        try:
                            poly_set = zone.GetFilledPolysList(target_layer)
                            api_method_used = "GetFilledPolysList(layer)"
                        except (TypeError, AttributeError) as e1:
                            # Method 2: GetFilledPolysList without parameter (fallback)
                            try:
                                poly_set = zone.GetFilledPolysList()
                                api_method_used = "GetFilledPolysList()"
                            except Exception as e2:
                                # Method 3: Try FilledPolysList property (alternative API)
                                try:
                                    poly_set = zone.FilledPolysList(target_layer)
                                    api_method_used = "FilledPolysList(layer)"
                                except:
                                    try:
                                        poly_set = zone.FilledPolysList()
                                        api_method_used = "FilledPolysList()"
                                    except Exception as e3:
                                        if self.debug:
                                            self.logger.error(f"  All API methods failed for layer {target_layer}")
                                            self.logger.error(f"    Method 1 error: {e1}")
                                            self.logger.error(f"    Method 2 error: {e2}")
                                            self.logger.error(f"    Method 3 error: {e3}")
                                        continue
                        
                        if poly_set is None:
                            if self.debug:
                                self.logger.warning(f"  No polygon set retrieved for layer {target_layer}")
                            continue
                        
                        if self.debug:
                            self.logger.debug(f"  Successfully got polygons using {api_method_used}")
                        
                        count = poly_set.OutlineCount()
                        if self.debug:
                            self.logger.debug(f"  Polygon set has {count} outline(s)")
                        
                        for i in range(count):
                            try:
                                outline = poly_set.Outline(i)
                                
                                # SHAPE_LINE_CHAIN extraction
                                points = []
                                pc = outline.PointCount()
                                
                                if self.debug:
                                    self.logger.debug(f"    Outline {i} has {pc} points")
                                
                                if pc < 3:
                                    if self.debug:
                                        self.logger.warning(f"    Outline {i} has less than 3 points, skipping")
                                    continue
                                
                                point_method_used = None
                                for j in range(pc):
                                    pt = None
                                    # Try multiple point access methods
                                    try:
                                        # Method 1: CPoint (KiCad 9.0)
                                        pt = outline.CPoint(j)
                                        if point_method_used is None:
                                            point_method_used = "CPoint"
                                    except AttributeError:
                                        try:
                                            # Method 2: Point (older versions)
                                            pt = outline.Point(j)
                                            if point_method_used is None:
                                                point_method_used = "Point"
                                        except AttributeError:
                                            try:
                                                # Method 3: GetPoint (potential alternative)
                                                pt = outline.GetPoint(j)
                                                if point_method_used is None:
                                                    point_method_used = "GetPoint"
                                            except Exception as e_pt:
                                                if self.debug and j == 0:
                                                    self.logger.error(f"    Failed to access point {j}: {e_pt}")
                                                continue
                                    
                                    if pt is None:
                                        continue
                                    
                                    x = getattr(pt, 'x', None)
                                    y = getattr(pt, 'y', None)
                                    if x is not None and y is not None:
                                        points.append((pcbnew.ToMM(x), pcbnew.ToMM(y)))
                                    else:
                                        if self.debug and j == 0:
                                            self.logger.warning(f"    Point {j} missing x or y attribute")
                                
                                if self.debug and point_method_used:
                                    self.logger.debug(f"    Successfully extracted points using {point_method_used}")
                                        
                                if len(points) >= 3:
                                    p = Polygon(points)
                                    if not p.is_valid:
                                        if self.debug:
                                            self.logger.debug(f"    Polygon invalid, buffering to fix")
                                        p = p.buffer(0)
                                    if p.is_valid:
                                        # Log bounds of individual outline BEFORE adding
                                        if self.debug:
                                            bounds = p.bounds
                                            width = bounds[2] - bounds[0]
                                            height = bounds[3] - bounds[1]
                                            self.logger.debug(f"    Outline {i} bounds: {width:.1f}x{height:.1f}mm at ({bounds[0]:.1f},{bounds[1]:.1f}) to ({bounds[2]:.1f},{bounds[3]:.1f})")
                                        
                                        add_shape(target_layer, p)
                                        if self.debug:
                                            self.logger.debug(f"    Added polygon with {len(points)} points to layer {target_layer}")
                                    else:
                                        if self.debug:
                                            self.logger.warning(f"    Polygon still invalid after buffer, skipping")
                                else:
                                    if self.debug:
                                        self.logger.warning(f"    Only extracted {len(points)} points (need 3+), skipping")
                            except Exception as e:
                                if self.debug:
                                    self.logger.error(f"    Error parsing outline {i}: {e}")
                                else:
                                    print(f"Error parsing outline {i}: {e}")
                    except Exception as e_zone:
                        if self.debug:
                            self.logger.error(f"  Error processing zone on layer {target_layer}: {e_zone}")
                        else:
                            print(f"Error processing zone layer {target_layer}: {e_zone}")
        
        if self.debug and zone_count > 0:
            self.logger.debug(f"Processed {zone_count} zone(s) for net '{net_name}'") 
                
        # 4. Merge
        merged_geometry = {}
        for layer, shapes in layer_shapes.items():
            if not shapes:
                continue
            merged = unary_union(shapes)
            merged_geometry[layer] = merged
            
            # Debug: Log geometry details
            if self.debug:
                bounds = merged.bounds  # (minx, miny, maxx, maxy)
                area = merged.area
                width = bounds[2] - bounds[0]
                height = bounds[3] - bounds[1]
                self.logger.debug(f"Layer {layer} geometry: area={area:.2f} mm², size={width:.1f}x{height:.1f} mm, bounds=({bounds[0]:.1f}, {bounds[1]:.1f}) to ({bounds[2]:.1f}, {bounds[3]:.1f})")
            
        return merged_geometry

    def plot_geometry(self, geometry_by_layer):
        """Visualize the extracted geometry as a 2D plot."""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches
            from matplotlib.collections import PatchCollection
            import io
            
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Color map for layers
            layer_colors = {0: 'red', 2: 'blue', 31: 'green'}
            
            for layer_id, geom in geometry_by_layer.items():
                if geom.is_empty:
                    continue
                
                color = layer_colors.get(layer_id, 'gray')
                
                # Handle MultiPolygon or Polygon
                if geom.geom_type == 'Polygon':
                    polys = [geom]
                elif geom.geom_type == 'MultiPolygon':
                    polys = list(geom.geoms)
                else:
                    continue
                
                for poly in polys:
                    x, y = poly.exterior.xy
                    ax.fill(x, y, alpha=0.5, fc=color, ec='black', linewidth=0.5, label=f'Layer {layer_id}')
                    
                    # Plot holes
                    for interior in poly.interiors:
                        x, y = interior.xy
                        ax.fill(x, y, alpha=1.0, fc='white', ec='black', linewidth=0.5)
            
            ax.set_xlabel('X (mm)')
            ax.set_ylabel('Y (mm)')
            ax.set_title('Extracted Geometry (Pre-Mesh)')
            ax.set_aspect('equal')
            ax.grid(True, alpha=0.3)
            ax.invert_yaxis()  # Match KiCad coordinate system
            
            # Remove duplicate labels
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            ax.legend(by_label.values(), by_label.keys())
            
            plt.tight_layout()
            
            # Save to PNG file
            import os
            output_path = os.path.join(os.path.dirname(__file__), 'debug_geometry.png')
            plt.savefig(output_path, format='png', dpi=150, bbox_inches='tight')
            if self.debug:
                self.logger.debug(f"Saved geometry plot to {output_path}")
            
            # Also return as buffer for UI display
            buf = io.BytesIO()
            plt.savefig(buf, format='png', dpi=100)
            plt.close(fig)
            buf.seek(0)
            return buf
            
        except Exception as e:
            if self.debug:
                self.logger.error(f"Failed to plot geometry: {e}")
            return None
