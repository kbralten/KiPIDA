import pcbnew
try:
    from shapely.geometry import LineString, Polygon, MultiPolygon, Point, box
    from shapely.ops import unary_union
except ImportError:
    LineString = Polygon = MultiPolygon = Point = box = unary_union = None

class GeometryExtractor:
    def __init__(self, board):
        self.board = board
        if LineString is None:
            raise ImportError("The 'shapely' library is required for Ki-PIDA but was not found.")

    def get_board_stackup(self):
        """
        Extracts physical stackup properties.
        Returns a dict with copper data and substrate segments.
        """
        rho_copper = 1.68e-5 # Ohm-mm
        copper_data = {}
        substrates = []
        
        try:
            settings = self.board.GetDesignSettings()
            stackup = settings.GetStackupDescriptor()
            
            last_copper = None
            # Check if stackup is valid and has GetCount
            if hasattr(stackup, 'GetCount'):
                for i in range(stackup.GetCount()):
                    item = stackup.GetItem(i)
                    lid = item.GetLayerId()
                    thickness = pcbnew.ToMM(item.GetThickness())
                    
                    if pcbnew.IsCopperLayer(lid):
                        copper_data[lid] = {
                            'name': self.board.GetLayerName(lid),
                            'thickness_mm': thickness if thickness > 0 else 0.035
                        }
                        last_copper = lid
                    elif lid == -1: # substrate/dielectric
                        if last_copper is not None:
                            if len(substrates) == 0 or substrates[-1]['between'][1] is not None:
                                substrates.append({'thickness_mm': thickness, 'between': [last_copper, None]})
                            else:
                                substrates[-1]['thickness_mm'] += thickness
                
                # Link substrates
                for sub in substrates:
                    start_lid = sub['between'][0]
                    found_start = False
                    for j in range(stackup.GetCount()):
                        lid = stackup.GetItem(j).GetLayerId()
                        if lid == start_lid:
                            found_start = True
                        elif found_start and pcbnew.IsCopperLayer(lid):
                            sub['between'][1] = lid
                            break
            else:
                raise AttributeError("StackupDescriptor has no GetCount")
                
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
                    if last_lid is not None:
                        # Add a default 1.0mm dielectric between layers
                        substrates.append({'thickness_mm': 1.0, 'between': [last_lid, lid]})
                    last_lid = lid

        return {
            'copper': copper_data,
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
                
                # TODO: Accurately determine which layers the via hits
                # For now, let's just add to top and bottom to prove concept, 
                # or iterate all copper layers and check if between top/bottom.
                pass 
                
            elif track.GetClass() == "PCB_TRACK" or track.GetClass() == "PCB_ARC":
                # Standard Track
                layer = track.GetLayer()
                width_mm = pcbnew.ToMM(track.GetWidth())
                
                start = track.GetStart()
                end = track.GetEnd()
                
                p0 = (pcbnew.ToMM(start.x), pcbnew.ToMM(start.y))
                p1 = (pcbnew.ToMM(end.x), pcbnew.ToMM(end.y))
                
                # Create LineString and buffer it to get a Polygon (track width)
                # Cap style 1 is ROUND
                line = LineString([p0, p1])
                poly = line.buffer(width_mm / 2, cap_style=1)
                
                add_shape(layer, poly)

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
        for zone in self.board.Zones():
            if zone.GetNetCode() != net_code:
                continue
            if not zone.GetIsRuleArea(): # Ensure it's a copper zone
                # Multi-layer handling
                layers_to_process = []
                layer = zone.GetLayer()
                
                if layer < 0:
                    # Multi-layer zone
                    try:
                        lset = zone.GetLayerSet()
                        for id in lset.Seq():
                            layers_to_process.append(id)
                    except:
                        pass
                else:
                    layers_to_process.append(layer)
                    
                for target_layer in layers_to_process:
                    try:
                        try:
                            poly_set = zone.GetFilledPolysList(target_layer)
                        except TypeError:
                            poly_set = zone.GetFilledPolysList()
                        
                        count = poly_set.OutlineCount()
                        
                        for i in range(count):
                            try:
                                outline = poly_set.Outline(i)
                                
                                # SHAPE_LINE_CHAIN extraction
                                points = []
                                pc = outline.PointCount()
                                if pc < 3: continue
                                
                                for j in range(pc):
                                    try:
                                        # KiCad 9.0 uses CPoint
                                        pt = outline.CPoint(j)
                                    except AttributeError:
                                        # Fallback
                                        pt = outline.Point(j)
                                        
                                    x = getattr(pt, 'x', None)
                                    y = getattr(pt, 'y', None)
                                    if x is not None and y is not None:
                                        points.append((pcbnew.ToMM(x), pcbnew.ToMM(y)))
                                        
                                if len(points) >= 3:
                                    p = Polygon(points)
                                    if not p.is_valid:
                                        p = p.buffer(0)
                                    add_shape(target_layer, p)
                            except Exception as e:
                                print(f"Error parsing outline {i}: {e}")
                    except Exception as e_zone:
                        print(f"Error processing zone layer {target_layer}: {e_zone}")
                        pass 
                
        # 4. Merge
        merged_geometry = {}
        for layer, shapes in layer_shapes.items():
            if not shapes:
                continue
            merged = unary_union(shapes)
            merged_geometry[layer] = merged
            
        return merged_geometry
