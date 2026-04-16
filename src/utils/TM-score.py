#!/usr/bin/env python3
"""
Skrypt do obliczania TM-score dla modeli RNA w zestawach PZ9 i PZ20.
Wykorzystuje bibliotekę 'tmtools' i 'BioPython'.
Instalacja: pip install tmtools biopython numpy
"""

import os
import csv
import re
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser

# === KONFIGURACJA ===
PUZZLE_SETS = {
    "PZ9": Path(r"D:\Studia\Projekty\gsp-lab\data\finished\PZ9\pdb"),
    "PZ20": Path(r"D:\Studia\Projekty\gsp-lab\data\finished\PZ20\pdb")
}
RESULTS_BASE = Path(r"D:\Studia\Projekty\gsp-lab\results\other-metrics\rnapuzzles.github.io")

def extract_lab_and_num(filename):
    """
    Wyodrębnia nazwę laboratorium (Lab) i numer modelu (Num) z nazwy pliku.
    Np. '20_Bujnicki_1.pdb' -> ('Bujnicki', 1)
        '9_Bujnicki_3_rpr.pdb' -> ('Bujnicki', 3)
        '20_RNAComposerHuman_1.pdb' -> ('RNAComposerHuman', 1)
    """
    name = Path(filename).stem  # Usuwa '.pdb'

    # Usuń prefix puzzla (np. '20_', '9_')
    name = re.sub(r'^\d+_', '', name)
    
    # Usuń suffix '_rpr', jeśli istnieje
    if name.endswith('_rpr'):
        name = name[:-4]
    
    # GŁÓWNA LOGIKA: Znajdź OSTATNI ciąg znaków i numer
    match = re.search(r'([^0-9]+)_(\d+)$', name)
    
    if match:
        lab = match.group(1).strip('_')  # Usuń ewentualne podkreślniki z krańców
        num = int(match.group(2))
        return lab, num
    else:
        # Fallback dla plików niepasujących do wzorca
        return name, 0

def extract_coords_and_seq_from_pdb(pdb_file, atom_name="C3'"):
    """
    Ekstrahuje współrzędne atomów i sekwencję jednoliterową z pliku PDB.
    Dla RNA standardowo używamy atomów C3' jako reprezentatywnych dla reszty.
    Jeśli atom C3' nie istnieje, funkcja spróbuje użyć atomu P.
    
    Args:
        pdb_file: ścieżka do pliku PDB
        atom_name: nazwa atomu do ekstrakcji (domyślnie "C3'" dla RNA)
    
    Returns:
        krotka (coords, seq) gdzie:
        - coords: numpy array o kształcie (N, 3)
        - seq: string jednoliterowych kodów reszt (np. 'AAAU' dla RNA)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('model', str(pdb_file))
    
    coords_list = []
    seq_list = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Pobierz hetero-flag i nazwę reszty
                hetero_flag, resseq, icode = residue.get_id()
                resname = residue.get_resname().strip()
                
                # Pomijaj reszty, które nie są standardowymi nukleotydami/wodą
                if hetero_flag.startswith('H_') or resname in ['HOH', 'WAT']:
                    continue
                
                # Mapowanie nazw reszt RNA na jednoliterowe kody
                # Możesz rozszerzyć to mapowanie według potrzeb
                resname_to_code = {
                    'A': 'A', 'ADE': 'A',
                    'U': 'U', 'URA': 'U',
                    'G': 'G', 'GUA': 'G',
                    'C': 'C', 'CYT': 'C'
                }
                
                res_code = resname_to_code.get(resname, 'X')  # 'X' dla nieznanych
                
                # Spróbuj wybrać żądany atom
                if atom_name in residue:
                    atom = residue[atom_name]
                    coords_list.append(atom.get_coord())
                    seq_list.append(res_code)
                # Fallback: jeśli nie ma C3', spróbuj użyć atomu P
                elif atom_name == "C3'" and 'P' in residue:
                    atom = residue['P']
                    coords_list.append(atom.get_coord())
                    seq_list.append(res_code)
    
    if not coords_list:
        raise ValueError(f"Nie znaleziono atomów '{atom_name}' w pliku {pdb_file}")
    
    return np.array(coords_list), ''.join(seq_list)

def calculate_tmscore_tmtools(model_file, solution_file):
    """
    Oblicza TM-score używając biblioteki tmtools.
    
    Args:
        model_file: ścieżka do pliku modelu PDB
        solution_file: ścieżka do pliku rozwiązania PDB
    
    Returns:
        Wartość TM-score (float) lub None w przypadku błędu
    """
    try:
        # Import tmtools - główna funkcja
        from tmtools import tm_align
        
        # Ekstrahuj współrzędne i sekwencje
        coords_model, seq_model = extract_coords_and_seq_from_pdb(model_file)
        coords_solution, seq_solution = extract_coords_and_seq_from_pdb(solution_file)
        
        # Oblicz TM-score
        result = tm_align(coords_model, coords_solution, seq_model, seq_solution)
        
        # Użyj TM-score znormalizowanego względem pierwszego łańcucha (modelu)
        return result.tm_norm_chain1
        
    except ImportError as e:
        print(f"Błąd importu tmtools: {e}")
        print("Upewnij się, że biblioteka jest zainstalowana: pip install tmtools")
        return None
    except Exception as e:
        print(f"Błąd obliczania TM-score: {str(e)[:100]}")
        return None

def find_solution_file(puzzle_dir):
    """Znajduje plik rozwiązania (target) w katalogu puzzla."""
    # Szukaj plików zawierających 'solution' (case-insensitive)
    for pattern in ['*solution*.pdb', '*target*.pdb']:
        solution_files = list(puzzle_dir.glob(pattern))
        if solution_files:
            # Zwróć pierwsze znalezione rozwiązanie
            return solution_files[0]
    raise FileNotFoundError(f"Nie znaleziono pliku rozwiązania w {puzzle_dir}")

def calculate_tmscore_for_set(puzzle_name, puzzle_path):
    """Oblicza TM-score dla wszystkich modeli w zestawie."""
    print(f"\n--- Przetwarzanie {puzzle_name} ---")
    
    # Ścieżka do wynikowego CSV
    output_dir = RESULTS_BASE / puzzle_name
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "TM-score.csv"
    
    # Znajdź plik rozwiązania
    try:
        solution_file = find_solution_file(puzzle_path)
        print(f"Plik rozwiązania: {solution_file.name}")
    except FileNotFoundError as e:
        print(f"BŁĄD: {e}")
        return []
    
    # Znajdź wszystkie pliki modeli (pomijając plik rozwiązania)
    model_files = []
    for pdb_file in puzzle_path.glob('*.pdb'):
        if pdb_file != solution_file and 'solution' not in pdb_file.name.lower():
            model_files.append(pdb_file)
    
    print(f"Znaleziono {len(model_files)} modeli do porównania.")
    
    # Wyniki
    results = []
    
    for i, model_file in enumerate(model_files, 1):
        lab, num = extract_lab_and_num(model_file.name)
        
        print(f"  [{i}/{len(model_files)}] {model_file.name}...", end=" ")
        
        try:
            # OBLICZ TM-SCORE używając tmtools
            tm_score = calculate_tmscore_tmtools(model_file, solution_file)
            
            if tm_score is not None:
                tm_score_rounded = round(tm_score, 4)
                print(f"TM-score = {tm_score_rounded}")
                results.append([lab, num, tm_score_rounded])
            else:
                print("BŁĄD")
                results.append([lab, num, None])
                
        except Exception as e:
            print(f"BŁĄD: {str(e)[:100]}")
            results.append([lab, num, None])
    
    # Sortuj wyniki według TM-score (malejąco)
    results.sort(key=lambda x: x[2] if x[2] is not None else 0, reverse=True)
    
    # Zapisz do CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Lab', 'Num', 'TM-score'])
        writer.writerows(results)
    
    print(f"Wyniki zapisano do: {output_file}")
    
    # Podsumowanie
    valid_scores = [r[2] for r in results if r[2] is not None]
    if valid_scores:
        print(f"Obliczono {len(valid_scores)} TM-score.")
        if len(valid_scores) > 1:
            print(f"Zakres: {min(valid_scores):.4f} - {max(valid_scores):.4f}")
        print(f"Średnia: {sum(valid_scores)/len(valid_scores):.4f}")
    
    return results

def main():
    """Główna funkcja skryptu."""
    print("="*60)
    print("SKRYPT DO OBLICZANIA TM-SCORE DLA RNA (tmtools)")
    print("="*60)
    
    # Sprawdź dostępność wymaganych bibliotek
    try:
        # Test importu tmtools
        from tmtools import tm_align
        print("✓ Biblioteka 'tmtools' załadowana pomyślnie.")
    except ImportError as e:
        print(f"\n✗ BŁĄD: Nie znaleziono biblioteki 'tmtools'.")
        print(f"  Szczegóły: {e}")
        print("  Zainstaluj: pip install tmtools biopython numpy")
        return
    
    try:
        # Test importu BioPython
        from Bio.PDB import PDBParser
        print("✓ Biblioteka 'BioPython' załadowana pomyślnie.")
    except ImportError as e:
        print(f"\n✗ BŁĄD: Nie znaleziono biblioteki 'BioPython'.")
        print(f"  Zainstaluj: pip install biopython")
        return
    
    # Przetwórz każdy zestaw
    all_results = {}
    for puzzle_name, puzzle_path in PUZZLE_SETS.items():
        if not puzzle_path.exists():
            print(f"\n✗ Katalog nie istnieje: {puzzle_path}")
            continue
        
        results = calculate_tmscore_for_set(puzzle_name, puzzle_path)
        all_results[puzzle_name] = results
    
    # Podsumowanie globalne
    print("\n" + "="*60)
    print("PODSUMOWANIE")
    print("="*60)
    for puzzle_name, results in all_results.items():
        valid = sum(1 for r in results if r[2] is not None)
        total = len(results)
        print(f"{puzzle_name}: {valid}/{total} pomyślnie obliczonych TM-score")

if __name__ == "__main__":
    main()