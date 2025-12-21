import os
import sys
import re
from Bio.PDB import PDBParser

def normalize_atom_name(atom_name):
    """Normalizuje nazwę atomu do postaci kanonicznej"""
    # Usuń spacje i konwersja na wielkie litery
    name = atom_name.replace(" ", "").upper()
    
    # Standaryzacja nazw z primami
    if "*" in name:
        name = name.replace("*", "'")
    
    # Usuń oznaczenia konformacji (jak 'A', 'B')
    if len(name) > 1 and name[-1] in "ABCDEFGHIJKLMNOPQRSTUVWXYZ":
        if name[-2] in "0123456789":
            name = name[:-1]
    
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

def exact_match_clean_model(reference_pdb, model_pdb, output_pdb):
    """Czyści model, zachowując tylko atomy obecne w referencji"""
    # Krok 1: Budowanie zestawu atomów referencyjnych
    reference_atoms = build_reference_set(reference_pdb)
    
    # Krok 2: Filtrowanie atomów w modelu
    atom_counter = 1
    output_lines = []
    
    # Zbierz wszystkie linie, aby zachować nagłówki
    with open(model_pdb, 'r') as fin:
        model_lines = fin.readlines()
    
    for line in model_lines:
        # Zachowaj linie nie-ATOM (HEADER, TITLE, itp.)
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
            # Aktualizuj numer atomu (zachowując ciągłość)
            new_line = line[:6] + f"{atom_counter:>5}" + line[11:]
            output_lines.append(new_line)
            atom_counter += 1
        else:
            # Logowanie usuwanych atomów
            print(f"Usunięto atom: {atom_key} - {line[12:26]}")
    
    # Dodaj końcowy rekord TER
    output_lines.append("TER\n")
    
    # Zapisz do pliku
    with open(output_pdb, 'w') as fout:
        fout.writelines(output_lines)

def select_files_from_list(file_list):
    """Wyświetla listę plików i umożliwia ich wybór"""
    print("\nNumer | Plik")
    print("------|------")
    for i, file in enumerate(file_list, 1):
        print(f"{i:5} | {file}")
    
    while True:
        selection = input("\nWybierz pliki (np. 1, 2-4, 1,3): ").strip()
        
        if not selection:
            return []
        
        try:
            selected_indices = set()
            parts = re.split(r'[,\s]+', selection)
            
            for part in parts:
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    if start < 1 or end > len(file_list) or start > end:
                        raise ValueError("Nieprawidłowy zakres")
                    selected_indices.update(range(start, end + 1))
                else:
                    index = int(part)
                    if index < 1 or index > len(file_list):
                        raise ValueError("Numer poza zakresem")
                    selected_indices.add(index)
            
            return [file_list[i-1] for i in sorted(selected_indices)]
        
        except ValueError as e:
            print(f"Błąd: {e}. Spróbuj ponownie.")
        except Exception:
            print("Nieprawidłowy format. Spróbuj ponownie.")

def process_directory(directory, reference_pdb):
    """Przetwarza wszystkie wybrane pliki PDB w katalogu"""
    pdb_files = [f for f in os.listdir(directory) 
                 if f.lower().endswith('.pdb') and os.path.isfile(os.path.join(directory, f))]
    
    if not pdb_files:
        print("Brak plików PDB w katalogu.")
        return
    
    print("\nZnaleziono pliki PDB:")
    selected_files = select_files_from_list(pdb_files)
    
    if not selected_files:
        print("Nie wybrano żadnych plików.")
        return
    
    for filename in selected_files:
        input_path = os.path.join(directory, filename)
        base_name = os.path.splitext(filename)[0]
        output_path = os.path.join(directory, f"{base_name}_cleaned.pdb")
        
        print(f"\nPrzetwarzanie: {filename} -> {base_name}_cleaned.pdb")
        exact_match_clean_model(reference_pdb, input_path, output_path)
    
    print("\nPrzetwarzanie zakończone!")

def main():
    if len(sys.argv) < 2:
        print("Sposób użycia:")
        print("  Pojedynczy plik: python clean_model.py <reference.pdb> <input.pdb> <output.pdb>")
        print("  Katalog:        python clean_model.py -d <sciezka/do/katalogu> <reference.pdb>")
        sys.exit(1)
    
    # Tryb katalogu
    if sys.argv[1] in ['-d', '--directory']:
        if len(sys.argv) < 3:
            print("Proszę podać ścieżkę do katalogu")
            sys.exit(1)
        
        directory = sys.argv[2]
        
        if len(sys.argv) < 4:
            print("Proszę podać plik referencyjny")
            sys.exit(1)
        
        reference_pdb = sys.argv[3]
        
        if not os.path.isdir(directory):
            print(f"Błąd: {directory} nie jest katalogiem")
            sys.exit(1)
            
        if not os.path.isfile(reference_pdb):
            print(f"Błąd: Plik referencyjny {reference_pdb} nie istnieje")
            sys.exit(1)
            
        print(f"\nPrzetwarzanie katalogu: {directory}")
        print(f"Plik referencyjny: {reference_pdb}")
        process_directory(directory, reference_pdb)
    
    # Tryb pojedynczego pliku
    else:
        if len(sys.argv) < 4:
            print("Proszę podać plik referencyjny, wejściowy i wyjściowy")
            sys.exit(1)
            
        reference_pdb = sys.argv[1]
        input_pdb = sys.argv[2]
        output_pdb = sys.argv[3]
        
        if not os.path.isfile(reference_pdb):
            print(f"Błąd: Plik referencyjny {reference_pdb} nie istnieje")
            sys.exit(1)
            
        if not os.path.isfile(input_pdb):
            print(f"Błąd: Plik wejściowy {input_pdb} nie istnieje")
            sys.exit(1)
            
        print(f"\nPrzetwarzanie pojedynczego pliku")
        print(f"Plik referencyjny: {reference_pdb}")
        print(f"Plik wejściowy:    {input_pdb}")
        print(f"Plik wyjściowy:    {output_pdb}")
        exact_match_clean_model(reference_pdb, input_pdb, output_pdb)
        print("Przetwarzanie zakończone!")

if __name__ == "__main__":
    main()