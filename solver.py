
import sys

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
    def __init__(self, log_callback=None):
        self.log_callback = log_callback
        if np is None or scipy is None:
            raise ImportError("NumPy and SciPy are required for Solver backend.")

    def _log(self, msg):
        if self.log_callback:
            self.log_callback(f"[SOLVER] {msg}")

    def solve(self, mesh, sources, loads):
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
        if hasattr(mesh, 'G_coo_data') and len(mesh.G_coo_data) > 0:
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
            
        # Initialize Vector I
        I = np.zeros(N)

        # 4. Apply Loads (Current Sources)
        for load in loads:
            nid = load.get('node_id')
            current = load.get('current', 0.0)
            if nid in id_to_idx:
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
                
        # 6. Solve System
        # Convert to CSR for solving efficiency
        G_csr = G.tocsr()
        
        # Diagnostics: Island Detection
        try:
            from scipy.sparse.csgraph import connected_components
            n_components, labels = connected_components(csgraph=G_csr, directed=False, return_labels=True)
            
            if n_components > 1:
                self._log(f"Detected {n_components} isolated copper islands.")
                # Check if each island has a source
                island_has_source = [False] * n_components
                for source in sources:
                    idx = id_to_idx.get(source['node_id'])
                    if idx is not None:
                        island_has_source[labels[idx]] = True
                
                for i, has_src in enumerate(island_has_source):
                    if not has_src:
                        n_nodes = np.count_nonzero(labels == i)
                        self._log(f"  Warning: Island #{i} ({n_nodes} nodes) has no voltage source. Results may be undefined.")
        except Exception as e:
            self._log(f"Connectivity diagnostic failed: {e}")
            
        try:
            if pypardiso is not None:
                self._log("Using high-performance PyPardiso solver.")
                V_solution = pypardiso.spsolve(G_csr, I)
            else:
                self._log("Using standard SciPy spsolve (PyPardiso not found).")
                V_solution = scipy.sparse.linalg.spsolve(G_csr, I)
            
            if np.any(np.isnan(V_solution)):
                self._log("Warning: Solution contains NaN values.")
        except Exception as e:
            self._log(f"Solver Exception: {e}")
            return {}
            
        # 7. Map results back
        results = {}
        for i, v_val in enumerate(V_solution):
            nid = idx_to_id[i]
            results[nid] = float(v_val)
            
        return results

