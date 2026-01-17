import re
import pcbnew

try:
    from .models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef
except (ImportError, ValueError):
    from models import PowerRail, UnifiedSource, UnifiedLoad, ComponentRef

class NetDiscoverer:
    """
    Scans the board to identify likely power rails and their properties.
    """
    def __init__(self, board):
        self.board = board
        
        # Regex for voltage parsing
        # Matches: +3V3, 3V3, +5V, 12V, VCC_3V3
        self.voltage_patterns = [
            r'(\d+)V(\d+)',  # 3V3 -> 3.3
            r'[+]?(\d+\.?\d*)V', # +5V, 12V -> 5.0, 12.0
            r'(\d+)V', # 5V -> 5.0
        ]
        
        # Regex for likely power names to prioritize
        self.power_name_patterns = [
            r'(?i)VCC', r'(?i)VDD', r'(?i)PWR', r'(?i)\+.*V', r'(?i)GND'
        ]

    def discover_power_nets(self):
        """
        Returns a list of PowerRail objects found on the board.
        Does NOT fully populate sources/loads automatically yet, just sets up candidates.
        """
        candidates = []
        nets = self.board.GetNetsByName()
        
        for netname, net in nets.items():
            s_name = str(netname)
            if not s_name: continue
            
            # 1. Check if name looks like power
            is_power_name = any(re.search(p, s_name) for p in self.power_name_patterns)
            
            # 2. Check Pin Types (TODO: This requires iterating all pads, can be slow on huge boards)
            # For now, rely heavily on name + voltage detection for the "Quick Scan"
            
            est_voltage = self._estimate_voltage(s_name)
            
            # If it has a voltage in name, or is explicitly VCC/GND-like, add it.
            if est_voltage is not None or is_power_name:
                rail = PowerRail(net_name=s_name)
                if est_voltage:
                    rail.nominal_voltage = est_voltage
                candidates.append(rail)
                
        # Sort candidates: Voltage found first, then name
        candidates.sort(key=lambda r: (r.nominal_voltage if r.nominal_voltage else -1, r.net_name), reverse=True)
        return candidates

    def _estimate_voltage(self, name):
        """
        Tries to extract a float voltage from the net name.
        Returns float or None.
        """
        # Exclude common current signals like "ISENSE_5V" if needed, 
        # but usually "5V" means 5V domain.
        
        # Pattern 1: 3V3 -> 3.3
        m = re.search(r'(\d+)V(\d+)', name, re.IGNORECASE)
        if m:
            try:
                major = int(m.group(1))
                minor = int(m.group(2))
                # infer decimal place. 'V' usually separates integer and decimal
                return float(f"{major}.{minor}")
            except: pass
            
        # Pattern 2: +12V or 5.0V
        m = re.search(r'plus(\d+\.?\d*)V', name, re.IGNORECASE) or \
            re.search(r'\+(\d+\.?\d*)V', name, re.IGNORECASE) 
        if m:
             try: return float(m.group(1))
             except: pass
             
        # Pattern 3: Simple "5V" suffix or prefix
        # Be careful not to match "V1" (version 1)
        m = re.search(r'\b(\d+\.?\d*)V\b', name, re.IGNORECASE)
        if m:
            try: return float(m.group(1))
            except: pass
            
        return None

    def get_components_on_net(self, net_name):
        """
        Returns a dict of { ref_des: [pad_obj, ...] } for all components on this net.
        """
        components = {}
        net = self.board.FindNet(net_name)
        if not net: return {}
        
        # Iterate all footprints, then pads
        # This is safer than GetPads() for grouping by component
        for fp in self.board.GetFootprints():
            ref = fp.GetReference()
            relevant_pads = []
            for pad in fp.Pads():
                if pad.GetNet().GetNetname() == net_name:
                    relevant_pads.append(pad)
            
            if relevant_pads:
                components[ref] = relevant_pads
                
        return components
