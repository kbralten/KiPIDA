
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

class Solver:
    def __init__(self):
        if np is None or scipy is None:
            raise ImportError("NumPy and SciPy are required for Solver backend.")

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
            return {}
            
        # 1. Map Node IDs to Matrix Indices (0..N-1)
        # We need a stable mapping.
        # mesh.nodes might be arbitrary integers.
        nodes = list(mesh.nodes)
        N = len(nodes)
        id_to_idx = { nid: i for i, nid in enumerate(nodes) }
        idx_to_id = { i: nid for i, nid in enumerate(nodes) }
        
        # 2. Initialize Matrix G (LiL for speed of filling) and Vector I
        G = scipy.sparse.lil_matrix((N, N))
        I = np.zeros(N)
        
        # 3. Fill G Matrix (Stamps Method)
        # G_ii = Sum(g_connected)
        # G_ij = -g_connected
        for u_id, v_id, g in mesh.edges:
            if u_id not in id_to_idx or v_id not in id_to_idx:
                continue
            
            u = id_to_idx[u_id]
            v = id_to_idx[v_id]
            
            G[u, u] += g
            G[v, v] += g
            G[u, v] -= g
            G[v, u] -= g
            
        # 4. Apply Loads (Current Sources)
        # KCL: Sum(G*V) = Sum(I_in) - Sum(I_out)
        # Standard form: G*V = I_vector
        # I_vector[n] represents current injecting INTO node n.
        # A Load draws current OUT of node n. So I_inj = -I_load.
        for load in loads:
            nid = load.get('node_id')
            current = load.get('current', 0.0)
            if nid in id_to_idx:
                idx = id_to_idx[nid]
                I[idx] -= current
                
        # 5. Apply Voltage Sources (Dirichlet BCs)
        # We enforce V_n = V_set.
        # Method: Row Replacement.
        # Zero out row 'n' in G.
        # Set G[n, n] = 1.
        # Set I[n] = V_set.
        
        # Note: If multiple sources are connected to same node, last one wins (short circuit logic not handled).
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
        
        # Diagnostics: Check for rows that are all-zero (isolated nodes)
        # In G, a node should have at least one non-zero on the diagonal or in the row.
        # But if it's strictly isolated, G[idx, idx] = 0.
        row_sums = np.abs(G_csr).sum(axis=1).A1
        zero_rows = np.where(row_sums == 0)[0]
        if len(zero_rows) > 0:
            pass
            # We can't solve for these. 
            # Strategy: Set them to 0V or something to avoid NaN? 
            # Or just warn.
            
        # Diagnostics: Island Detection
        try:
            from scipy.sparse.csgraph import connected_components
            n_components, labels = connected_components(csgraph=G_csr, directed=False, return_labels=True)
            
            if n_components > 1:
                pass
                # Check if each island has a source
                island_has_source = [False] * n_components
                for source in sources:
                    idx = id_to_idx.get(source['node_id'])
                    if idx is not None:
                        island_has_source[labels[idx]] = True
                
                for i, has_src in enumerate(island_has_source):
                    if not has_src:
                        # Count nodes in this floating island
                        n_nodes = np.count_nonzero(labels == i)
                        pass
        except Exception as e:
            print(f"Connectivity diagnostic failed: {e}")
            
        try:
            V_solution = scipy.sparse.linalg.spsolve(G_csr, I)
            
            if np.any(np.isnan(V_solution)):
                pass
        except Exception as e:
            print(f"Solver Exception: {e}")
            return {}
            
        # 7. Map results back
        results = {}
        for i, v_val in enumerate(V_solution):
            nid = idx_to_id[i]
            results[nid] = float(v_val)
            
        return results
