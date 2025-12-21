#!/usr/bin/env python3
"""
Rozszerzony skrypt do porównywania plików PDB - wykrywa wszystkie różnice strukturalne
"""

import sys
from Bio import PDB
from collections import defaultdict

def normalize_atom_name(atom_name):
    """Normalizuje nazwę atomu do postaci kanonicznej"""
    name = atom_name.strip().replace(" ", "").upper()
    if "*" in name:
        name = name.replace("*", "'")
    return name

def build_structure_map(pdb_file):
    """Buduje szczegółową mapę struktury PDB bez koordynatów"""
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    structure_map = {
        'chains': {},
        'residues': defaultdict(dict),
        'atoms': defaultdict(lambda: defaultdict(dict)),
        'atom_count': 0,
        'residue_count': 0
    }
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            structure_map['chains'][chain_id] = {
                'residues': [],
                'atom_count': 0
            }
            
            for residue in chain:
                res_id = residue.id
                res_name = residue.resname.strip()
                res_num = res_id[1]
                res_icode = res_id[2] if len(res_id) > 2 else " "
                
                residue_key = (chain_id, res_num, res_icode)
                structure_map['residues'][chain_id][residue_key] = {
                    'name': res_name,
                    'atoms': {},
                    'full_id': res_id
                }
                
                for atom in residue:
                    atom_name = normalize_atom_name(atom.name)
                    structure_map['atoms'][chain_id][residue_key][atom_name] = {
                        'element': atom.element,
                        'serial_number': atom.serial_number,
                    }
                    structure_map['atom_count'] += 1
                    structure_map['chains'][chain_id]['atom_count'] += 1
                
                structure_map['residue_count'] += 1
                structure_map['chains'][chain_id]['residues'].append(residue_key)
    
    return structure_map

def compare_structures(target_pdb, model_pdb):
    """Porównuje dwie struktury PDB i znajduje wszystkie różnice strukturalne"""
    
    print("=" * 80)
    print("KOMPLEKSOWE PORÓWNANIE STRUKTUR PDB")
    print("=" * 80)
    print(f"Target: {target_pdb}")
    print(f"Model:  {model_pdb}")
    print("=" * 80)
    
    # Budowanie map struktur
    print("Wczytywanie struktur...")
    target_map = build_structure_map(target_pdb)
    model_map = build_structure_map(model_pdb)
    
    # Statystyki podstawowe
    print("\nSTATYSTYKI PODSTAWOWE:")
    print(f"Target: {target_map['atom_count']} atomów, {target_map['residue_count']} reszt, {len(target_map['chains'])} łańcuchów")
    print(f"Model:  {model_map['atom_count']} atomów, {model_map['residue_count']} reszt, {len(model_map['chains'])} łańcuchów")
    
    differences = {
        'chain_differences': [],
        'residue_differences': [],
        'atom_differences': [],
        'mutations': [],
        'numbering_differences': []
    }
    
    # 1. PORÓWNANIE ŁAŃCUCHÓW
    print("\n" + "=" * 50)
    print("PORÓWNANIE ŁAŃCUCHÓW")
    print("=" * 50)
    
    target_chains = set(target_map['chains'].keys())
    model_chains = set(model_map['chains'].keys())
    
    missing_chains = target_chains - model_chains
    extra_chains = model_chains - target_chains
    common_chains = target_chains & model_chains
    
    if missing_chains:
        print(f"Łańcuchy BRAKUJĄCE w modelu: {sorted(missing_chains)}")
        differences['chain_differences'].extend([f"Brakujący łańcuch: {chain}" for chain in missing_chains])
    
    if extra_chains:
        print(f"Łańcuchy DODATKOWE w modelu: {sorted(extra_chains)}")
        differences['chain_differences'].extend([f"Dodatkowy łańcuch: {chain}" for chain in extra_chains])
    
    if not missing_chains and not extra_chains:
        print("Łańcuchy: IDENTYCZNE")
    
    # 2. PORÓWNANIE RESZT
    print("\n" + "=" * 50)
    print("PORÓWNANIE RESZT")
    print("=" * 50)
    
    residue_differences_found = False
    
    for chain in common_chains:
        target_residues = set(target_map['residues'][chain].keys())
        model_residues = set(model_map['residues'][chain].keys())
        
        missing_residues = target_residues - model_residues
        extra_residues = model_residues - target_residues
        
        if missing_residues:
            residue_differences_found = True
            print(f"Łańcuch {chain}: Brakujące reszty w modelu:")
            for res_key in sorted(missing_residues):
                res_info = target_map['residues'][chain][res_key]
                print(f"  - {chain}{res_key[1]} ({res_info['name']})")
                differences['residue_differences'].append(f"Brakująca reszta: {chain}{res_key[1]} ({res_info['name']})")
        
        if extra_residues:
            residue_differences_found = True
            print(f"Łańcuch {chain}: Dodatkowe reszty w modelu:")
            for res_key in sorted(extra_residues):
                res_info = model_map['residues'][chain][res_key]
                print(f"  - {chain}{res_key[1]} ({res_info['name']})")
                differences['residue_differences'].append(f"Dodatkowa reszta: {chain}{res_key[1]} ({res_info['name']})")
        
        # Sprawdzanie mutacji (ta sama pozycja, inna reszta)
        common_residues = target_residues & model_residues
        for res_key in common_residues:
            target_res_name = target_map['residues'][chain][res_key]['name']
            model_res_name = model_map['residues'][chain][res_key]['name']
            
            if target_res_name != model_res_name:
                residue_differences_found = True
                print(f"MUTACJA: {chain}{res_key[1]} {target_res_name} -> {model_res_name}")
                differences['mutations'].append(f"{chain}{res_key[1]}: {target_res_name} -> {model_res_name}")
    
    if not residue_differences_found:
        print("Reszty: IDENTYCZNE")
    
    # 3. PORÓWNANIE ATOMÓW
    print("\n" + "=" * 50)
    print("PORÓWNANIE ATOMÓW")
    print("=" * 50)
    
    atom_differences_found = False
    
    for chain in common_chains:
        for res_key in set(target_map['residues'][chain].keys()) & set(model_map['residues'][chain].keys()):
            target_atoms = set(target_map['atoms'][chain][res_key].keys())
            model_atoms = set(model_map['atoms'][chain][res_key].keys())
            
            missing_atoms = target_atoms - model_atoms
            extra_atoms = model_atoms - target_atoms
            
            if missing_atoms:
                atom_differences_found = True
                res_name = target_map['residues'][chain][res_key]['name']
                print(f"Łańcuch {chain}, Reszta {res_key[1]} ({res_name}): Brakujące atomy w modelu:")
                for atom in sorted(missing_atoms):
                    print(f"  - {atom}")
                    differences['atom_differences'].append(f"Brakujący atom: {chain}{res_key[1]}.{atom}")
            
            if extra_atoms:
                atom_differences_found = True
                res_name = model_map['residues'][chain][res_key]['name']
                print(f"Łańcuch {chain}, Reszta {res_key[1]} ({res_name}): Dodatkowe atomy w modelu:")
                for atom in sorted(extra_atoms):
                    print(f"  - {atom}")
                    differences['atom_differences'].append(f"Dodatkowy atom: {chain}{res_key[1]}.{atom}")
    
    if not atom_differences_found:
        print("Atomy: IDENTYCZNE")
    
    # 4. PORÓWNANIE NUMERACJI ATOMÓW
    print("\n" + "=" * 50)
    print("PORÓWNANIE NUMERACJI ATOMÓW")
    print("=" * 50)
    
    all_target_atoms = []
    all_model_atoms = []
    
    for chain in target_map['atoms']:
        for res_key in target_map['atoms'][chain]:
            for atom_name, atom_data in target_map['atoms'][chain][res_key].items():
                all_target_atoms.append(atom_data['serial_number'])
    
    for chain in model_map['atoms']:
        for res_key in model_map['atoms'][chain]:
            for atom_name, atom_data in model_map['atoms'][chain][res_key].items():
                all_model_atoms.append(atom_data['serial_number'])
    
    if sorted(all_target_atoms) != sorted(all_model_atoms):
        print("UWAGA: Różnice w numeracji atomów!")
        print(f"Target: {len(all_target_atoms)} atomów (numery: {min(all_target_atoms)}-{max(all_target_atoms)})")
        print(f"Model:  {len(all_model_atoms)} atomów (numery: {min(all_model_atoms)}-{max(all_model_atoms)})")
        
        # Sprawdź czy to tylko problem z numeracją czy z liczbą atomów
        if len(all_target_atoms) != len(all_model_atoms):
            print("RÓŻNA LICZBA ATOMÓW - to nie tylko problem z numeracją!")
        else:
            print("Ta sama liczba atomów, ale różna numeracja")
            
        differences['numbering_differences'].append("Różna numeracja atomów")
    else:
        print("Numeracja atomów: SPÓJNA")
    
    # 5. PODSUMOWANIE
    print("\n" + "=" * 80)
    print("PODSUMOWANIE RÓŻNIC")
    print("=" * 80)
    
    total_differences = (len(differences['chain_differences']) + 
                        len(differences['residue_differences']) + 
                        len(differences['atom_differences']) + 
                        len(differences['mutations']) + 
                        len(differences['numbering_differences']))
    
    print(f"Łączna liczba różnic: {total_differences}")
    print(f"- Różnice w łańcuchach: {len(differences['chain_differences'])}")
    print(f"- Różnice w resztach: {len(differences['residue_differences'])}")
    print(f"- Różnice w atomach: {len(differences['atom_differences'])}")
    print(f"- Mutacje: {len(differences['mutations'])}")
    print(f"- Różnice w numeracji: {len(differences['numbering_differences'])}")
    
    if total_differences == 0:
        print("\n✓ Struktury są IDENTYCZNE pod względem strukturalnym!")
    else:
        print("\n✗ Znaleziono różnice między strukturami!")
    
    return differences

def analyze_first_residue(target_pdb, model_pdb):
    """Specjalna funkcja do analizy pierwszej reszty"""
    print("\n" + "=" * 80)
    print("ANALIZA PIERWSZEJ RESZTY")
    print("=" * 80)
    
    parser = PDB.PDBParser(QUIET=True)
    
    for name, pdb_file in [("Target", target_pdb), ("Model", model_pdb)]:
        print(f"\n{name}: {pdb_file}")
        structure = parser.get_structure(name, pdb_file)
        
        for model in structure:
            for chain in model:
                print(f"\nŁańcuch {chain.id}:")
                for residue in chain:
                    res_num = residue.id[1]
                    res_name = residue.resname
                    atoms = [atom.name for atom in residue]
                    
                    print(f"  Reszta {res_num} ({res_name}): {', '.join(sorted(atoms))}")
                    
                    # Szczegółowa analiza pierwszej reszty
                    if res_num == 1:
                        print(f"  >>> SZCZEGÓŁY RESZTY 1:")
                        for atom in residue:
                            print(f"      Atom: {atom.name} (numer: {atom.serial_number})")
                    
                    # Przerwij po kilku pierwszych resztach dla czytelności
                    if res_num >= 5:
                        print(f"  ... i dalsze reszty")
                        break
                break  # Tylko pierwszy łańcuch
            break  # Tylko pierwszy model

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Użycie: python compare_pdb.py target.pdb model.pdb")
        print("\nPrzykład: python compare_pdb.py structure1.pdb structure2.pdb")
        sys.exit(1)
    
    target_pdb = sys.argv[1]
    model_pdb = sys.argv[2]
    
    try:
        # Główne porównanie struktur
        differences = compare_structures(target_pdb, model_pdb)
        
        # Specjalna analiza pierwszej reszty
        analyze_first_residue(target_pdb, model_pdb)
        
    except Exception as e:
        print(f"\n❌ Błąd podczas analizy: {e}")
        import traceback
        traceback.print_exc()