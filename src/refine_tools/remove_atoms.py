#!/usr/bin/env python3
"""
Skrypt do usuwania nadmiarowych atomów z modelu PDB na podstawie target
"""

import sys
import os
from Bio import PDB

def normalize_atom_name(atom_name):
    """Normalizuje nazwę atomu do postaci kanonicznej"""
    name = atom_name.strip().replace(" ", "").upper()
    if "*" in name:
        name = name.replace("*", "'")
    return name

def build_reference_set(reference_pdb):
    """Buduje zestaw atomów referencyjnych z unikalnymi identyfikatorami"""
    atom_set = set()
    
    with open(reference_pdb, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
                
            # Ekstrakcja kluczowych identyfikatorów
            chain_id = line[21]
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                continue
            
            residue_name = line[17:20].strip()
            atom_name = normalize_atom_name(line[12:16].strip())
            
            # Utwórz unikalny klucz atomu
            atom_key = (chain_id, res_seq, residue_name, atom_name)
            atom_set.add(atom_key)
            
    return atom_set

def remove_extra_atoms(reference_pdb, model_pdb, output_pdb):
    """Usuwa z modelu atomy, których nie ma w referencji"""
    # Krok 1: Budowanie zestawu atomów referencyjnych
    reference_atoms = build_reference_set(reference_pdb)
    
    # Krok 2: Filtrowanie atomów w modelu
    atom_counter = 1
    output_lines = []
    
    with open(model_pdb, 'r') as fin:
        model_lines = fin.readlines()
    
    removed_count = 0
    for line in model_lines:
        # Zachowaj linie nie-ATOM
        if not line.startswith('ATOM'):
            output_lines.append(line)
            continue
            
        # Ekstrakcja klucza identyfikującego atom
        chain_id = line[21]
        try:
            res_seq = int(line[22:26].strip())
        except ValueError:
            # Zachowaj linię jeśli nie możemy sparsować
            output_lines.append(line)
            continue
        
        residue_name = line[17:20].strip()
        atom_name = normalize_atom_name(line[12:16].strip())
        
        # Utwórz unikalny klucz atomu
        atom_key = (chain_id, res_seq, residue_name, atom_name)
        
        # Sprawdź czy atom istnieje w referencji
        if atom_key in reference_atoms:
            # Aktualizuj numer atomu
            new_line = line[:6] + f"{atom_counter:>5}" + line[11:]
            output_lines.append(new_line)
            atom_counter += 1
        else:
            removed_count += 1
            print(f"Usunięto atom: {atom_key}")
    
    # Dodaj końcowy rekord TER
    output_lines.append("TER\n")
    
    # Zapisz do pliku
    with open(output_pdb, 'w') as fout:
        fout.writelines(output_lines)
    
    print(f"Usunięto {removed_count} atomów")
    return True

def main():
    if len(sys.argv) != 4:
        print("Użycie: python remove_atoms.py reference.pdb model.pdb output.pdb")
        sys.exit(1)
    
    reference_pdb = sys.argv[1]
    model_pdb = sys.argv[2]
    output_pdb = sys.argv[3]
    
    if not os.path.exists(reference_pdb):
        print(f"❌ Plik referencyjny nie istnieje: {reference_pdb}")
        sys.exit(1)
    
    if not os.path.exists(model_pdb):
        print(f"❌ Plik model nie istnieje: {model_pdb}")
        sys.exit(1)
    
    try:
        success = remove_extra_atoms(reference_pdb, model_pdb, output_pdb)
        if success:
            print(f"✅ Pomyślnie usunięto nadmiarowe atomy")
        else:
            print("❌ Nie udało się usunąć atomów")
            sys.exit(1)
    except Exception as e:
        print(f"❌ Błąd podczas usuwania atomów: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()