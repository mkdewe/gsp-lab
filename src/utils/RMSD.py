#!/usr/bin/env python3
"""
single_rmsd.py - Oblicza RMSD między dwoma plikami PDB
Użycie: python single_rmsd.py reference.pdb model.pdb
"""

import sys
import numpy as np

def parse_pdb_coords_fast(pdb_file, atom_name="C3'"):
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                current_atom = line[12:16].strip()
                if current_atom == atom_name:
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except ValueError:
                        continue
    return np.array(coords)

def kabsch_rmsd(P, Q):
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    
    C = np.dot(P_centered.T, Q_centered)
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    
    U = np.dot(V, W)
    P_rotated = np.dot(P_centered, U)
    
    diff = P_rotated - Q_centered
    rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
    
    return rmsd

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python single_rmsd.py reference.pdb model.pdb")
        sys.exit(1)
    
    ref_file = sys.argv[1]
    model_file = sys.argv[2]
    
    ref_coords = parse_pdb_coords_fast(ref_file)
    model_coords = parse_pdb_coords_fast(model_file)
    
    if len(ref_coords) == 0 or len(model_coords) == 0:
        print("Error: No C3' atoms found in one or both files")
        sys.exit(1)
    
    if len(ref_coords) != len(model_coords):
        print(f"Warning: Different number of C3' atoms (ref={len(ref_coords)}, model={len(model_coords)})")
    
    rmsd = kabsch_rmsd(model_coords, ref_coords)
    print(f"RMSD: {rmsd:.3f} Å")