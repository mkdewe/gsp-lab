#!/usr/bin/env python3
"""
Skrypt do naprawy numeracji atomów w plikach PDB z obsługą rekordów TER
"""

import sys
from Bio import PDB
from Bio.PDB import PDBIO

def fix_pdb_atom_numbers(input_file, output_file):
    """Naprawia numerację atomów w pliku PDB z obsługą TER"""
    # Statystyki
    atom_count = 0
    ter_count = 0
    
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                atom_count += 1
                # Zachowujemy resztę linii, zmieniamy tylko numer atomu
                new_line = f"ATOM  {atom_count:>5}{line[11:]}"
                f_out.write(new_line)
                
            elif line.startswith('TER'):
                # Dla TER: numerujemy tylko jeśli oryginał miał numer
                if line[6:11].strip().isdigit():
                    ter_count += 1
                    new_line = f"TER   {ter_count:>5}{line[11:]}"
                else:
                    new_line = "TER\n"  # Bez numeru
                f_out.write(new_line)
                
            else:
                f_out.write(line)
    
    print(f"Naprawiono numerację:")
    print(f"- Atomów: {atom_count}")
    print(f"- Rekordów TER: {ter_count}")
    print(f"Zapisano do: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Użycie: python fix_pdb_ter.py input.pdb output.pdb")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    
    fix_pdb_atom_numbers(input_pdb, output_pdb)