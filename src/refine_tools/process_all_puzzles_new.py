import os
import subprocess
import glob
from complete_rna_mapper_new import SemanticRNAMapper

def parse_selection(selection, max_index):
    """Parsuje wybór użytkownika (np. '1,3,5-7')"""
    selection = selection.replace(" ", "")
    chosen = set()
    for part in selection.split(","):
        if "-" in part:
            try:
                start, end = part.split("-")
                chosen.update(range(int(start), int(end) + 1))
            except ValueError:
                print(f"Nieprawidłowy zakres: {part}")
                continue
        else:
            try:
                chosen.add(int(part))
            except ValueError:
                print(f"Nieprawidłowy numer: {part}")
                continue
    return [i - 1 for i in chosen if 0 < i <= max_index]

def select_puzzle_sets(external_base):
    """Wyświetla listę zestawów i pozwala użytkownikowi wybrać"""
    # Znajdź wszystkie zestawy
    puzzle_sets = [d for d in os.listdir(external_base) 
                  if os.path.isdir(os.path.join(external_base, d))]
    puzzle_sets.sort()
    
    if not puzzle_sets:
        print("Nie znaleziono żadnych zestawów!")
        return []
    
    print(f"\nZnaleziono {len(puzzle_sets)} zestawów:")
    for i, puzzle_set in enumerate(puzzle_sets, 1):
        print(f"{i:3}. {puzzle_set}")
    
    while True:
        print("\n" + "="*50)
        print("Wybierz zestawy do przetworzenia:")
        print("  - Podaj numery (np. 1,3,5-7)")
        print("  - 'all' dla wszystkich zestawów")
        print("  - 'q' aby wyjść")
        
        choice = input("\nTwój wybór: ").strip().lower()
        
        if choice == 'q':
            return []
        elif choice == 'all':
            print(f"\nWybrano WSZYSTKIE zestawy ({len(puzzle_sets)}). Rozpoczynam przetwarzanie...")
            return puzzle_sets
        else:
            try:
                indices = parse_selection(choice, len(puzzle_sets))
                if indices:
                    selected_sets = [puzzle_sets[i] for i in indices]
                    print(f"\nWybrano {len(selected_sets)} zestawów. Rozpoczynam przetwarzanie...")
                    return selected_sets
                else:
                    print("Błędny wybór. Spróbuj ponownie.")
            except Exception as e:
                print(f"Błąd: {e}. Spróbuj ponownie.")

def process_selected_puzzles(puzzle_sets, remove_extra_atoms=False):
    """Przetwarzaj wybrane zestawy - BEZ subprocess."""
    # Ścieżki bazowe
    external_base = r"D:\Studia\Projekty\gsp-lab\data\external\rnapuzzles.github.io\data"
    processed_base = r"D:\Studia\Projekty\gsp-lab\data\processed\rnapuzzles.github.io\data"
    
    total_success = 0
    total_errors = 0
    total_skipped = 0
    
    # Utwórz mapper
    mapper = SemanticRNAMapper()
    
    for idx, puzzle_set in enumerate(puzzle_sets, 1):
        print(f"\n{'='*60}")
        print(f"Zestaw {idx}/{len(puzzle_sets)}: {puzzle_set}")
        print(f"{'='*60}")
        
        set_external_path = os.path.join(external_base, puzzle_set, "pdb")
        set_processed_path = os.path.join(processed_base, puzzle_set, "pdb")
        
        if not os.path.exists(set_external_path):
            print("  Brak folderu pdb - pomijam")
            total_skipped += 1
            continue
        
        os.makedirs(set_processed_path, exist_ok=True)
        
        # Znajdź pliki PDB
        pdb_files = glob.glob(os.path.join(set_external_path, "*.pdb"))
        if not pdb_files:
            print("  Brak plików PDB")
            total_skipped += 1
            continue
        
        print(f"  Pliki PDB: {len(pdb_files)}")
        
        # Znajdź solution
        solution_files = [f for f in pdb_files if 'solution' in os.path.basename(f).lower()]
        if not solution_files:
            solution_files = [f for f in pdb_files if 'target' in os.path.basename(f).lower()]
        if not solution_files:
            solution_files = [f for f in pdb_files if 'native' in os.path.basename(f).lower()]
        
        if solution_files:
            target_file = solution_files[0]
            print(f"  Target: {os.path.basename(target_file)}")
        else:
            print("  BRAK SOLUTION - pomijam zestaw")
            total_skipped += 1
            continue
        
        # Przetwarzaj modele
        model_files = [f for f in pdb_files if f != target_file]
        success = 0
        errors = 0
        
        for model_idx, model_file in enumerate(model_files, 1):
            model_name = os.path.splitext(os.path.basename(model_file))[0]
            output_file = os.path.join(set_processed_path, f"{model_name}_refined.pdb")
            
            print(f"    [{model_idx}/{len(model_files)}] {model_name}...", end=" ", flush=True)
            
            try:
                # Użyj bezpośrednio klasy
                result = mapper.process_all(target_file, model_file, output_file, remove_extra_atoms=remove_extra_atoms)
                
                (residue_mapping, missing_atoms_report, mean_dist, changed_count, 
                 removed_atoms_count, removed_residues_count, atom_comparison_log, 
                 total_extra_atoms, total_missing_atoms, solution_output) = result
                
                if residue_mapping:
                    print(f"OK (zmapowano {len(residue_mapping)} reszt, średnia odległość: {mean_dist:.3f}Å)")
                    success += 1
                    
                    if solution_output and os.path.exists(solution_output):
                        print(f"      Utworzono solution: {os.path.basename(solution_output)}")
                        
                    # Logowanie szczegółów jeśli są błędy
                    if missing_atoms_report:
                        print(f"      Brakujące atomy: {total_missing_atoms}")
                    if removed_atoms_count > 0:
                        print(f"      Usunięto {removed_atoms_count} nadmiarowych atomów")
                        
                else:
                    print(f"Brak mapowania")
                    errors += 1
                    
            except Exception as e:
                print(f"BŁĄD: {e}")
                import traceback
                traceback.print_exc()
                errors += 1
        
        print(f"  Wynik: {success} OK, {errors} błędów")
        total_success += success
        total_errors += errors
    
    print(f"\n{'='*60}")
    print("PODSUMOWANIE:")
    print(f"{'='*60}")
    print(f"Łącznie przetworzonych zestawów: {len(puzzle_sets)}")
    print(f"Łącznie modeli OK: {total_success}")
    print(f"Łącznie błędów: {total_errors}")
    print(f"Łącznie pominiętych zestawów: {total_skipped}")
    print(f"{'='*60}")

def main():
    """Główna funkcja programu"""
    # Sprawdź czy complete_rna_mapper_semantic.py istnieje
        # Sprawdź czy complete_rna_mapper_new.py istnieje
    if not os.path.exists("complete_rna_mapper_new.py"):  
        print("BŁĄD: complete_rna_mapper_new.py nie znaleziony!")  
        print("Upewnij się, że skrypt znajduje się w bieżącym katalogu.")
        exit(1)
    
    # Ścieżka bazowa
    external_base = r"D:\Studia\Projekty\gsp-lab\data\external\rnapuzzles.github.io\data"
    
    if not os.path.exists(external_base):
        print(f"BŁĄD: Ścieżka nie istnieje: {external_base}")
        exit(1)
    
    # Zapytaj o opcję remove_extra_atoms
    print("\n" + "="*50)
    print("OPCJE PRZETWARZANIA:")
    print("  Czy usunąć nadmiarowe atomy z modeli?")
    print("  - 'y' (tak) - usunie atomy nieobecne w target")
    print("  - 'n' (nie) - zachowa wszystkie atomy (domyślnie)")
    remove_choice = input("\nUsunąć nadmiarowe atomy? (y/N): ").strip().lower()
    remove_extra_atoms = remove_choice == 'y'
    
    if remove_extra_atoms:
        print("Tryb: usuwanie nadmiarowych atomów")
    else:
        print("Tryb: zachowywanie wszystkich atomów")
    
    # Wybierz zestawy
    selected_sets = select_puzzle_sets(external_base)
    
    if not selected_sets:
        print("\nNie wybrano żadnych zestawów. Koniec programu.")
        return
    
    # Przetwarzaj wybrane zestawy
    process_selected_puzzles(selected_sets, remove_extra_atoms=remove_extra_atoms)

if __name__ == "__main__":
    main()