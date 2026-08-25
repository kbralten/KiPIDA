
import sys

try:
    from .models import DCSolveResult
except (ImportError, ValueError):
    from models import DCSolveResult

# Try imports, handle if missing (though they should be in KiCad 9)
try:
    import numpy as np
    import scipy
    import scipy.sparse
    import scipy.sparse.linalg
except ImportError:
    np = None
    scipy = None

try:
    import pypardiso
except ImportError:
    pypardiso = None

class Solver:
    def __init__(self, debug=False, log_callback=None):
        self.debug = debug
        self.log_callback = log_callback
        if np is None or scipy is None:
            raise ImportError("NumPy and SciPy are required for Solver backend.")

    def _log(self, msg):
        if self.log_callback:
            self.log_callback(f"[SOLVER] {msg}")

    def solve(self, mesh, sources, loads, branch_resistance_scales=None):
        """
        Solves the DC circuit Mesh: [G][V] = [I]
        
        Args:
            mesh: Mesh object with .nodes (list) and .edges (list of u,v,g)
            sources: list of dicts { 'node_id': int, 'voltage': float } (Dirichlet)
            loads: list of dicts { 'node_id': int, 'current': float } (Neumann/Source)
            
        Returns:
            dict: { node_id: voltage_float }
        """
        if not mesh.nodes:
            self._log("Mesh has no nodes. Returning empty result.")
            return {}
            
        # 1. Map Node IDs to Matrix Indices (0..N-1)
        nodes = list(mesh.nodes)
        N = len(nodes)
        id_to_idx = { nid: i for i, nid in enumerate(nodes) }
        idx_to_id = { i: nid for i, nid in enumerate(nodes) }
        
        # 2. Build Matrix G
        if branch_resistance_scales is not None and getattr(mesh, 'branches', None):
            G = scipy.sparse.lil_matrix((N, N))
            scales = list(branch_resistance_scales)
            if len(scales) != len(mesh.branches):
                raise ValueError("One resistance scale is required for every mesh branch.")
            for branch, scale in zip(mesh.branches, scales):
                resistance = max(branch.resistance_ohm * max(float(scale), 1e-9), 1e-15)
                conductance = 1.0 / resistance
                u = id_to_idx.get(branch.node_a)
                v = id_to_idx.get(branch.node_b)
                if u is None or v is None:
                    continue
                G[u, u] += conductance
                G[v, v] += conductance
                G[u, v] -= conductance
                G[v, u] -= conductance
        elif hasattr(mesh, 'G_coo_data') and len(mesh.G_coo_data) > 0:
            if self.debug:
                self._log(f"Using pre-computed sparse matrix ({len(mesh.G_coo_data)} entries).")
            # We have raw node IDs in G_coo_row/col, need to map to indices 0..N-1
            # Optimally, Mesher should produce 0..N indices, but it deals with arbitrary node IDs.
            # If node IDs are 0..N-1 sequential (which they are in Mesher implementation), we can skip mapping?
            # Let's verify: Mesher.generate_mesh starts node_counter=0 and increments.
            # So if mesh.nodes is sorted 0..N-1, id_to_idx is identity.
            
            # Check if mapping is needed (heuristic: if max node id < N, probably okay, but safer to map)
            # Vectorized mapping using numpy is fast.
            
            # Convert to numpy arrays if not already
            row_ids = np.array(mesh.G_coo_row)
            col_ids = np.array(mesh.G_coo_col)
            data = np.array(mesh.G_coo_data)
            
            # Check if simple identity mapping checks out
            # Mesher logic guarantees 0...N-1 sequentially.
            # But let's be robust:
            # We can use a fast lookup array if max(nodes) isn't huge.
            G = scipy.sparse.coo_matrix((data, (row_ids, col_ids)), shape=(N, N))
            G = G.tolil()
            
        else:
            if self.debug:
                self._log("Using legacy edge iteration for matrix build.")
            # Legacy Path (Slow)
            G = scipy.sparse.lil_matrix((N, N))
            for u_id, v_id, g in mesh.edges:
                if u_id not in id_to_idx or v_id not in id_to_idx:
                    continue
                
                u = id_to_idx[u_id]
                v = id_to_idx[v_id]
                
                G[u, u] += g
                G[v, v] += g
                G[u, v] -= g
                G[v, u] -= g
            
        # Determine electrical connectivity before applying Dirichlet rows.
        # Replacing a source row with an identity row changes the matrix graph
        # and can create artificial one-node islands in diagnostics.  It also
        # leaves genuinely floating copper islands in a singular system.
        valid_node_mask = np.ones(N, dtype=bool)
        floating_representatives = []
        try:
            from scipy.sparse.csgraph import connected_components

            connectivity_matrix = G.tocsr()
            n_components, component_labels = connected_components(
                csgraph=connectivity_matrix,
                directed=False,
                return_labels=True,
            )
            source_components = set()
            for source in sources:
                idx = id_to_idx.get(source.get('node_id'))
                if idx is not None:
                    source_components.add(int(component_labels[idx]))

            valid_node_mask = np.array(
                [int(label) in source_components for label in component_labels],
                dtype=bool,
            )
            if n_components > 1:
                self._log(f"Detected {n_components} copper islands before source constraints.")

            for component in range(n_components):
                if component in source_components:
                    continue
                component_indices = np.flatnonzero(component_labels == component)
                if len(component_indices) == 0:
                    continue
                floating_representatives.append(int(component_indices[0]))
                load_nodes = {
                    load.get('node_id') for load in loads
                    if load.get('node_id') in id_to_idx
                    and int(component_labels[id_to_idx[load.get('node_id')]]) == component
                    and abs(float(load.get('current', 0.0))) > 0.0
                }
                if load_nodes:
                    self._log(
                        f"ERROR: Island #{component} ({len(component_indices)} nodes) has "
                        f"{len(load_nodes)} load node(s) but no voltage source; "
                        "its loads and voltages are excluded."
                    )
                else:
                    self._log(
                        f"Ignoring floating island #{component} ({len(component_indices)} nodes): "
                        "no voltage source or load."
                    )
        except Exception as e:
            self._log(f"Connectivity diagnostic failed: {e}")

        # Initialize Vector I
        I = np.zeros(N)

        # 4. Apply Loads (Current Sources)
        for load in loads:
            nid = load.get('node_id')
            current = load.get('current', 0.0)
            if nid in id_to_idx and valid_node_mask[id_to_idx[nid]]:
                idx = id_to_idx[nid]
                I[idx] -= current
                
        # 5. Apply Voltage Sources (Dirichlet BCs)
        for source in sources:
            nid = source.get('node_id')
            voltage = source.get('voltage', 0.0)
            if nid in id_to_idx:
                idx = id_to_idx[nid]
                
                # Zero out the row efficiently
                G.rows[idx] = [idx]
                G.data[idx] = [1.0]
                
                I[idx] = voltage

        # Anchor one node in each excluded component so the full sparse matrix
        # remains nonsingular.  These reference values are never returned and
        # therefore cannot contaminate voltage-drop statistics or plots.
        for idx in floating_representatives:
            G.rows[idx] = [idx]
            G.data[idx] = [1.0]
            I[idx] = 0.0
                
        # 6. Solve System
        # Convert to CSR for solving efficiency
        G_csr = G.tocsr()
        
        try:
            if pypardiso is not None:
                if self.debug: self._log("Using high-performance PyPardiso solver.")
                V_solution = pypardiso.spsolve(G_csr, I)
            else:
                if self.debug: self._log("Using standard SciPy spsolve (PyPardiso not found).")
                V_solution = scipy.sparse.linalg.spsolve(G_csr, I)
            
            if np.any(np.isnan(V_solution)):
                self._log("Warning: Solution contains NaN values.")
        except Exception as e:
            self._log(f"Solver Exception: {e}")
            return {}
            
        # 7. Map results back
        results = {}
        for i, v_val in enumerate(V_solution):
            if not valid_node_mask[i]:
                continue
            nid = idx_to_id[i]
            results[nid] = float(v_val)
            
        return results

    def solve_detailed(self, mesh, sources, loads, branch_resistance_scales=None):
        """Solve DC and retain branch currents/losses for thermal coupling."""
        voltages = self.solve(
            mesh,
            sources,
            loads,
            branch_resistance_scales=branch_resistance_scales,
        )
        scales = (
            list(branch_resistance_scales)
            if branch_resistance_scales is not None
            else [1.0] * len(getattr(mesh, 'branches', []))
        )
        if len(scales) != len(getattr(mesh, 'branches', [])):
            raise ValueError("One resistance scale is required for every mesh branch.")
        currents = []
        losses = []
        for branch, scale in zip(getattr(mesh, 'branches', []), scales):
            resistance = max(branch.resistance_ohm * max(float(scale), 1e-9), 1e-15)
            voltage_delta = voltages.get(branch.node_a, 0.0) - voltages.get(branch.node_b, 0.0)
            current = voltage_delta / resistance
            currents.append(float(current))
            losses.append(float(current * current * resistance))
        return DCSolveResult(
            voltages=voltages,
            branch_currents_a=currents,
            branch_losses_w=losses,
            total_loss_w=float(sum(losses)),
        )

