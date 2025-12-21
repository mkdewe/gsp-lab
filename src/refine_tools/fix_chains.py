#!/usr/bin/env python3
"""
Uniwersalny skrypt do naprawy łańcuchów w plikach PDB
Dzieli łańcuch modelu na taką samą liczbę łańcuchów jak w target
"""

import sys
import os
from Bio import PDB
from collections import defaultdict

def analyze_chain_structure(pdb_file):
    """
    Analizuje strukturę łańcuchów w pliku PDB
    Zwraca: listę łańcuchów i liczbę reszt w każdym łańcuchu
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    chain_info = []
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            residues = list(chain.get_residues())
            residue_count = len(residues)
            
            chain_info.append({
                'id': chain_id,
                'residue_count': residue_count,
                'residues': residues
            })
    
    return chain_info

def get_atom_lines_from_pdb(pdb_file):
    """
    Zwraca wszystkie linie ATOM/HETATM z pliku PDB
    """
    atom_lines = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_lines.append(line)
    return atom_lines

def split_model_by_target_structure(target_chain_info, model_atom_lines, output_pdb):
    """
    Dzieli atomy modelu na łańcuchy zgodnie ze strukturą target
    """
    # Oblicz całkowitą liczbę reszt w target
    total_target_residues = sum(chain['residue_count'] for chain in target_chain_info)
    total_model_atoms = len(model_atom_lines)
    
    print(f"Target: {len(target_chain_info)} łańcuchów, {total_target_residues} reszt")
    print(f"Model: {total_model_atoms} atomów")
    
    # Jeśli model ma mniej atomów niż target ma reszt, to może być problem
    if total_model_atoms < total_target_residues:
        print(f"⚠️  Ostrzeżenie: Model ma mniej atomów ({total_model_atoms}) niż target ma reszt ({total_target_residues})")
    
    # Przygotowanie danych wyjściowych
    output_lines = []
    atom_counter = 1
    
    # Używamy standardowych identyfikatorów łańcuchów
    standard_chain_ids = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    
    # Dla każdego łańcucha w target
    current_atom_index = 0
    for i, target_chain in enumerate(target_chain_info):
        if i >= len(standard_chain_ids):
            print(f"❌ Zbyt wiele łańcuchów w target (max 26)")
            return False
            
        new_chain_id = standard_chain_ids[i]
        residues_needed = target_chain['residue_count']
        
        print(f"Tworzenie łańcucha {new_chain_id} z {residues_needed} reszt...")
        
        # Zbierz atomy dla tego łańcucha
        chain_atoms = []
        residues_found = 0
        current_residue_num = None
        
        while residues_found < residues_needed and current_atom_index < total_model_atoms:
            line = model_atom_lines[current_atom_index]
            
            # Sprawdź numer reszty
            try:
                res_seq = int(line[22:26].strip())
            except ValueError:
                current_atom_index += 1
                continue
            
            # Jeśli to nowa reszta, zwiększ licznik
            if current_residue_num != res_seq:
                current_residue_num = res_seq
                residues_found += 1
            
            # Zmodyfikuj linię: nowy łańcuch i numer atomu
            new_line = line[:21] + new_chain_id + line[22:]  # Zmiana łańcucha
            new_line = new_line[:6] + f"{atom_counter:>5}" + new_line[11:]  # Nowy numer atomu
            
            chain_atoms.append(new_line)
            atom_counter += 1
            current_atom_index += 1
        
        # Dodaj atomy tego łańcucha do wyniku
        output_lines.extend(chain_atoms)
        
        # Dodaj rekord TER na końcu łańcucha
        if chain_atoms:
            output_lines.append(f"TER   {atom_counter:>5}      {new_chain_id}\n")
            atom_counter += 1
    
    # Dodaj pozostałe atomy (jeśli model miał więcej atomów)
    remaining_atoms = total_model_atoms - current_atom_index
    if remaining_atoms > 0:
        print(f"⚠️  Pominięto {remaining_atoms} nadmiarowych atomów z modelu")
    
    # Dodaj rekord END na końcu
    output_lines.append("END\n")
    
    # Zapisz wynik
    with open(output_pdb, 'w') as f_out:
        f_out.writelines(output_lines)
    
    return True

def main():
    if len(sys.argv) != 4:
        print("Użycie: python fix_chains.py target.pdb model.pdb output.pdb")
        print("")
        print("Opis:")
        print("  target.pdb - plik referencyjny z pożądaną strukturą łańcuchów")
        print("  model.pdb  - plik do naprawy (zwykle z jednym łańcuchem)")
        print("  output.pdb - wynikowy plik z naprawionymi łańcuchami")
        sys.exit(1)
    
    target_pdb = sys.argv[1]
    model_pdb = sys.argv[2]
    output_pdb = sys.argv[3]
    
    if not os.path.exists(target_pdb):
        print(f"❌ Plik target nie istnieje: {target_pdb}")
        sys.exit(1)
    
    if not os.path.exists(model_pdb):
        print(f"❌ Plik model nie istnieje: {model_pdb}")
        sys.exit(1)
    
    try:
        print("🔍 Analiza struktury target...")
        target_chain_info = analyze_chain_structure(target_pdb)
        
        print("Struktura target:")
        for chain in target_chain_info:
            print(f"  Łańcuch {chain['id']}: {chain['residue_count']} reszt")
        
        print("\n🔍 Wczytywanie atomów z modelu...")
        model_atom_lines = get_atom_lines_from_pdb(model_pdb)
        
        print("\n🛠️  Rozpoczynanie podziału łańcuchów...")
        success = split_model_by_target_structure(target_chain_info, model_atom_lines, output_pdb)
        
        if success:
            print(f"\n✅ Pomyślnie utworzono plik: {output_pdb}")
            
            # Weryfikacja
            print("\n🔍 Weryfikacja wyniku:")
            result_chain_info = analyze_chain_structure(output_pdb)
            print("Struktura wynikowa:")
            for chain in result_chain_info:
                print(f"  Łańcuch {chain['id']}: {chain['residue_count']} reszt")
                
        else:
            print("❌ Nie udało się naprawić łańcuchów")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Błąd podczas naprawy łańcuchów: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()