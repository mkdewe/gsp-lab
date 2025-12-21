#!/usr/bin/env python3
import sys
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial.distance import cdist

def main():
    if len(sys.argv) != 3:
        print("Użycie: python3 map_residues.py target.pdb model.pdb")
        return
    
    target_file, model_file = sys.argv[1], sys.argv[2]
    
    parser = PDBParser()
    target = parser.get_structure("target", target_file)
    model = parser.get_structure("model", model_file)
    
    # Pobierz sekwencje
    def get_sequence(structure):
        return [(res.id[1], res.resname.strip()) for res in structure.get_residues() 
                if res.id[0] == ' ' and res.resname.strip() in ['A','C','G','U']]
    
    target_seq = get_sequence(target)
    model_seq = get_sequence(model)
    
    print("=== PORÓWNANIE SEKWENCJI ===")
    print(f"Target: {[f'{r[1]}{r[0]}' for r in target_seq]}")
    print(f"Model:  {[f'{r[1]}{r[0]}' for r in model_seq]}")
    
    # Automatyczne mapowanie przez współrzędne P
    print("\n=== AUTOMATYCZNE MAPOWANIE ===")
    # ... (wstaw kod z KROKU 3 powyżej)

if __name__ == "__main__":
    main()