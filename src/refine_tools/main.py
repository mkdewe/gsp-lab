#!/usr/bin/env python3
"""
Główny skrypt koordynujący diagnostykę i naprawę plików PDB
"""

import sys
import os
import subprocess
import shutil

# Import skryptu diagnostycznego
from compare_pdb import compare_structures

def run_diagnosis(target_pdb, model_pdb):
    """
    Uruchamia diagnostykę i zwraca szczegółowe informacje o różnicach
    """
    print("🔍 Uruchamianie diagnostyki...")
    differences = compare_structures(target_pdb, model_pdb)
    return differences

def check_problems(differences):
    """
    Sprawdza jakie problemy zostały wykryte i zwraca listę akcji do wykonania
    """
    actions = []
    
    # NOWA KOLEJNOŚĆ: najpierw łańcuchy, potem atomy, potem numeracja
    problem_priority = [
        ('chain_differences', 'fix_chains'),
        ('atom_differences', 'fix_atoms'), 
        ('numbering_differences', 'fix_numbering'),
        ('residue_differences', 'fix_residues'),
        ('mutations', 'fix_mutations')
    ]
    
    for problem_type, action_name in problem_priority:
        if differences.get(problem_type):
            print(f"❌ Wykryto problemy z {problem_type.replace('_', ' ')}")
            actions.append(action_name)
    
    return actions

def fix_chains(target_pdb, model_pdb, output_pdb):
    """
    Naprawia problemy z łańcuchami/niciami
    """
    print("🛠️  Naprawianie problemów z łańcuchami/niciami...")
    
    try:
        # Uruchom skrypt i przechwyć output
        result = subprocess.run([
            'python', 'fix_chains.py',
            target_pdb, model_pdb, output_pdb
        ], check=True, capture_output=True, text=True, timeout=30)
        
        # Wyświetl output ze skryptu
        if result.stdout:
            print("📝 Output z fix_chains.py:")
            print(result.stdout)
        if result.stderr:
            print("⚠️  Błędy z fix_chains.py:")
            print(result.stderr)
            
        # Sprawdź czy plik wyjściowy został utworzony
        if os.path.exists(output_pdb):
            # Sprawdź czy zmiana się udała - przeczytaj pierwszy atom
            with open(output_pdb, 'r') as f:
                for line in f:
                    if line.startswith('ATOM'):
                        chain = line[21]
                        print(f"🔍 Przykładowy atom po naprawie: łańcuch = '{chain}'")
                        if chain == 'A':
                            print("✅ Łańcuch został poprawnie zmieniony na 'A'")
                        else:
                            print(f"❌ Łańcuch nadal to '{chain}' zamiast 'A'")
                        break
            return True
        else:
            print("❌ Plik wyjściowy nie został utworzony")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Błąd podczas naprawy łańcuchów: {e}")
        if e.stdout:
            print("STDOUT:", e.stdout)
        if e.stderr:
            print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
        print("❌ Skrypt fix_chains.py nie został znaleziony")
        return False
    except Exception as e:
        print(f"❌ Nieoczekiwany błąd: {e}")
        return False

def fix_residues(target_pdb, model_pdb, output_pdb):
    """
    Naprawia problemy z resztami - usuwa nadmiarowe reszty
    """
    print("🛠️  Naprawianie problemów z resztami...")
    print("📝 Uwaga: Automatyczna naprawa reszt nie jest jeszcze zaimplementowana")
    return False

def fix_atoms(target_pdb, model_pdb, output_pdb):
    """
    Naprawia problemy z atomami (brakujące/zbędne atomy)
    """
    print("🛠️  Naprawianie problemów z atomami...")
    
    try:
        result = subprocess.run([
            'python', 'remove_atoms.py',
            target_pdb, model_pdb, output_pdb
        ], check=True, capture_output=True, text=True, timeout=30)
        
        if result.stdout:
            print("📝 Output z remove_atoms.py:")
            print(result.stdout)
            
        return os.path.exists(output_pdb)
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Błąd podczas usuwania atomów: {e}")
        return False
    except FileNotFoundError:
        print("❌ Skrypt remove_atoms.py nie został znaleziony")
        return False

def fix_mutations(target_pdb, model_pdb, output_pdb):
    """
    Naprawia mutacje (zmiany typów reszt)
    """
    print("🛠️  Naprawianie mutacji...")
    print("📝 Uwaga: Automatyczna naprawa mutacji nie jest jeszcze zaimplementowana")
    return False

def fix_numbering(target_pdb, model_pdb, output_pdb):
    """
    Naprawia problemy z numeracją atomów
    """
    print("🛠️  Naprawianie problemów z numeracją...")
    
    try:
        result = subprocess.run([
            'python', 'fix_pdb_numbers.py',
            model_pdb, output_pdb
        ], check=True, capture_output=True, text=True, timeout=30)
        
        if result.stdout:
            print("📝 Output z fix_pdb_numbers.py:")
            print(result.stdout)
            
        return os.path.exists(output_pdb)
            
    except subprocess.CalledProcessError as e:
        print(f"❌ Błąd podczas naprawy numeracji: {e}")
        return False
    except FileNotFoundError:
        print("❌ Skrypt fix_pdb_numbers.py nie został znaleziony")
        return False

def create_backup(file_path):
    """Tworzy kopię zapasową pliku"""
    backup_path = file_path + '.backup'
    shutil.copy2(file_path, backup_path)
    return backup_path

def verify_chain_fix(file_path):
    """Sprawdza czy łańcuchy zostały naprawione"""
    print(f"🔍 Weryfikacja naprawy łańcuchów w: {file_path}")
    
    chains_found = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21]
                chains_found.add(chain)
                # Sprawdzamy tylko pierwsze kilka atomów
                if len(chains_found) > 0:
                    break
    
    print(f"Znalezione łańcuchy: {chains_found}")
    if 'A' in chains_found and len(chains_found) == 1:
        print("✅ Weryfikacja łańcuchów: POPRAWNA")
        return True
    else:
        print(f"❌ Weryfikacja łańcuchów: NIEPOPRAWNA - znaleziono: {chains_found}")
        return False

def main_repair_cycle(target_pdb, model_pdb, max_cycles=5):
    """
    Główna pętla naprawcza - wykonuje diagnostykę i naprawy w cyklach
    """
    current_model = model_pdb
    cycle = 1
    
    print("=" * 80)
    print("🚀 URUCHAMIANIE AUTOMATYCZNEGO SYSTEMU NAPRAWCZEGO PDB")
    print("=" * 80)
    
    # Tworzenie kopii zapasowej oryginalnego modelu
    original_backup = create_backup(model_pdb)
    print(f"📦 Utworzono kopię zapasową: {original_backup}")
    
    while cycle <= max_cycles:
        print(f"\n🔄 CYKL NAPRAWCZY {cycle}/{max_cycles}")
        print("-" * 50)
        
        # Krok 1: Diagnostyka
        differences = run_diagnosis(target_pdb, current_model)
        
        # Krok 2: Sprawdzenie problemów
        actions = check_problems(differences)
        
        # Krok 3: Jeśli nie ma problemów - zakończ
        if not actions:
            print("✅ Nie wykryto problemów - proces zakończony pomyślnie!")
            return True
        
        # Krok 4: Wykonaj naprawy w określonej kolejności
        fixed_any = False
        temp_output = f"temp_cycle_{cycle}.pdb"
        
        # NOWA KOLEJNOŚĆ NAPRAW - łańcuchy pierwsze!
        repair_order = ['fix_chains', 'fix_atoms', 'fix_numbering', 'fix_residues', 'fix_mutations']
        
        for action in repair_order:
            if action in actions:
                print(f"\n🎯 Wykonywanie akcji: {action}")
                
                # Wywołaj odpowiednią funkcję naprawczą
                repair_function = globals()[action]
                success = repair_function(target_pdb, current_model, temp_output)
                
                if success:
                    # Jeśli naprawa się udała, użyj naprawionego pliku w kolejnym cyklu
                    if os.path.exists(temp_output):
                        # Zweryfikuj naprawę łańcuchów jeśli to była akcja fix_chains
                        if action == 'fix_chains':
                            chain_success = verify_chain_fix(temp_output)
                            if not chain_success:
                                print("❌ Naprawa łańcuchów nie powiodła się!")
                                continue
                        
                        if current_model != model_pdb and current_model.startswith('temp_cycle_'):
                            os.remove(current_model)  # Usuń stary plik tymczasowy
                        current_model = temp_output
                    fixed_any = True
                    print(f"✅ Sukces: {action}")
                    
                    # Przerwij obecny cykl i rozpocznij nową diagnostykę
                    break
                else:
                    print(f"❌ Niepowodzenie: {action}")
        
        # Jeśli w tym cyklu nic nie naprawiono, zakończ
        if not fixed_any:
            print("\n💥 Nie udało się naprawić żadnego problemu w tym cyklu")
            print("   Wymagana ręczna interwencja!")
            return False
        
        cycle += 1
    
    print(f"\n⏰ Osiągnięto maksymalną liczbę cykli ({max_cycles})")
    print("   Proces naprawczy został zatrzymany")
    return False

def interactive_mode():
    """Tryb interaktywny"""
    print("\n🔍 TRYB INTERAKTYWNY")
    print("=" * 50)
    
    target_pdb = input("Podaj ścieżkę do pliku TARGET (referencyjnego) PDB: ").strip()
    model_pdb = input("Podaj ścieżkę do pliku MODEL (do naprawy) PDB: ").strip()
    
    if not os.path.exists(target_pdb):
        print(f"❌ Plik target nie istnieje: {target_pdb}")
        return
    
    if not os.path.exists(model_pdb):
        print(f"❌ Plik model nie istnieje: {model_pdb}")
        return
    
    success = main_repair_cycle(target_pdb, model_pdb)
    
    if success:
        print("\n🎉 PROCES NAPRAWCZY ZAKOŃCZONY SUKCESEM!")
    else:
        print("\n💔 PROCES NAPRAWCZY NIE POWIÓDŁ SIĘ!")
        print("   Wymagana ręczna interwencja.")

def direct_mode():
    """Tryb bezpośredni"""
    if len(sys.argv) != 3:
        print("Użycie: python main_repair.py target.pdb model.pdb")
        print("Lub:    python main_repair.py --interactive")
        sys.exit(1)
    
    target_pdb = sys.argv[1]
    model_pdb = sys.argv[2]
    
    if not os.path.exists(target_pdb):
        print(f"❌ Plik target nie istnieje: {target_pdb}")
        sys.exit(1)
    
    if not os.path.exists(model_pdb):
        print(f"❌ Plik model nie istnieje: {model_pdb}")
        sys.exit(1)
    
    success = main_repair_cycle(target_pdb, model_pdb)
    
    if success:
        print("\n🎉 PROCES NAPRAWCZY ZAKOŃCZONY SUKCESEM!")
        sys.exit(0)
    else:
        print("\n💔 PROCES NAPRAWCZY NIE POWIÓDŁ SIĘ!")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in ['-i', '--interactive']:
        interactive_mode()
    else:
        direct_mode()