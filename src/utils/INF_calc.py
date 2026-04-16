#!/usr/bin/env python3
"""
Ulepszony skrypt wsadowy do obliczania INF-all dla RNA Puzzles.
Poprawione parsowanie nazw, wyszukiwanie plików i debugowanie.
"""

import os
import sys
import json
import csv
import warnings
import re
import argparse
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime
from collections import defaultdict

# Wyłącz wszystkie warningi
warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')

# Dodaj utils do ścieżki
sys.path.insert(0, '/work/src/utils')

# Import INF modkułu
try:
    from INF import compute_INF
    INF_AVAILABLE = True
except ImportError as e:
    print(f"⚠️  Błąd importu INF: {e}")
    print("Upewnij się, że INF.py jest w /work/src/utils/")
    INF_AVAILABLE = False

class INFAnalyzer:
    """Klasa do analizy INF dla zestawów RNA Puzzles."""
    
    def __init__(self, base_dir: str = "/work/data/finished", 
                 output_dir: str = "/work/results",
                 debug: bool = False):
        """
        Inicjalizacja analizatora.
        
        Args:
            base_dir: Główny katalog z danymi
            output_dir: Katalog na wyniki
            debug: Tryb debugowania
        """
        self.base_dir = Path(base_dir)
        self.output_dir = Path(output_dir)
        self.debug = debug
        
        # Katalogi dla wyników
        self.inf_all_dir = self.output_dir / "INF_all"
        self.inf_all_dir.mkdir(parents=True, exist_ok=True)
        
        self.log_file = self.output_dir / "inf_analysis.log"
        self.setup_logging()
        
        # Mapa zestawów
        self.puzzle_sets = {}
        self.results = {}
        self.summary = {}
        
        # Statystyki
        self.stats = {
            'start_time': datetime.now(),
            'total_puzzles': 0,
            'total_models': 0,
            'successful': 0,
            'failed': 0,
            'errors': []
        }
        
    def setup_logging(self):
        """Konfiguruje logowanie."""
        import logging
        logging.basicConfig(
            level=logging.DEBUG if self.debug else logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def log(self, level: str, message: str):
        """Loguje wiadomość na określonym poziomie."""
        if self.debug or level in ['ERROR', 'WARNING', 'INFO']:
            print(f"[{level}] {message}")
    
    def find_pdb_files(self, directory: Path) -> List[Path]:
        """
        Rekurencyjnie znajduje wszystkie pliki PDB w katalogu.
        Rozszerzona wersja: szuka również plików z rozszerzeniami .ent, .pdb1, .cif
        oraz plików bez rozszerzenia, które mogą zawierać struktury PDB.
        """
        pdb_files = []
        
        # Wzorce do wyszukiwania
        patterns = [
            "*.pdb", "*.PDB", 
            "*.pdb1", "*.PDB1",
            "*.ent", "*.ENT",
            "*.cif", "*.CIF",
            "*_rpr", "*_refined", "*_solution"
        ]
        
        for pattern in patterns:
            try:
                found_files = list(directory.rglob(pattern))
                pdb_files.extend(found_files)
            except Exception as e:
                if self.debug:
                    self.log("DEBUG", f"    Błąd wyszukiwania {pattern} w {directory}: {e}")
        
        # Dodatkowo: sprawdź pliki bez rozszerzenia o odpowiedniej strukturze nazw
        try:
            all_files = list(directory.rglob("*"))
            for file in all_files:
                if file.is_file():
                    # Sprawdź czy plik może być PDB po nazwie
                    name_lower = file.name.lower()
                    if (any(keyword in name_lower for keyword in [
                        '_rpr', '_refined', '_solution', '_model', 
                        'cr', 'pz', 'ts', 'adamiak', 'bujnicki', 'chen', 'das'
                    ]) and file.suffix in ['.pdb', '.PDB', '.pdb1', '.ent', '.cif', '']):
                        # Sprawdź czy to plik tekstowy (PDB)
                        try:
                            with open(file, 'r', encoding='latin-1') as f:
                                first_line = f.readline()
                                if 'HEADER' in first_line or 'ATOM' in first_line or 'REMARK' in first_line:
                                    pdb_files.append(file)
                        except:
                            pass
        except Exception as e:
            if self.debug:
                self.log("DEBUG", f"    Błąd przeszukiwania wszystkich plików: {e}")
        
        # Usuń duplikaty
        pdb_files = list(set(pdb_files))
        
        # Filtruj tylko pliki o rozsądnym rozmiarze (więcej niż 100 bajtów)
        pdb_files = [f for f in pdb_files if f.stat().st_size > 100]
        
        return pdb_files
    
    def parse_filename(self, filename: str) -> Tuple[str, str, int]:
        """
        Parsuje nazwę pliku RNA Puzzle i wyodrębnia: lab (laboratorium), num (numer modelu).
        Poprawiona dla wszystkich formatów z Twoich przykładów.
        
        Args:
            filename: Nazwa pliku
            
        Returns:
            Tuple: (base, lab, num)
        """
        stem = Path(filename).stem
        
        if self.debug:
            self.log("DEBUG", f"    Parsowanie nazwy: {stem}")
        
        # Inicjalizacja
        base = "Unknown"
        lab = ""
        num = 0
        
        # 1. Dla formatów CR: CR1116TS029_GWxraylab_x_rpr
        if stem.startswith('CR'):
            # Wyciągnij bazę (CR1116)
            cr_match = re.match(r'(CR\d+)', stem)
            if cr_match:
                base = cr_match.group(1)
            
            # Wyciągnij TSxxx (TS029, TS232)
            ts_match = re.search(r'(TS\d+)', stem)
            if ts_match:
                lab = ts_match.group(1)
            
            # Wyciągnij nazwę laboratorium (część po TSxxx_)
            # Format: CR1116TS029_GWxraylab_x_rpr
            pattern1 = r'CR\d+(TS\d+)_([A-Za-z0-9\-]+)_(\d+)_rpr'
            match1 = re.match(pattern1, stem)
            if match1:
                ts, lab_name, num_str = match1.groups()
                if not lab:  # Jeśli jeszcze nie mamy lab
                    lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
            
            # Format: CR1116TS232_AIchemy_RNA2_x_rpr
            pattern2 = r'CR\d+(TS\d+)_([A-Za-z0-9\-]+)_RNA2_(\d+)_rpr'
            match2 = re.match(pattern2, stem)
            if match2:
                ts, lab_name, num_str = match2.groups()
                if not lab:  # Jeśli jeszcze nie mamy lab
                    lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
        
        # 2. Dla formatów PZ z numerem na początku: 10_Bujnicki_x_rpr
        elif re.match(r'^\d+_', stem):
            # Wyciągnij numer (10, 12, 14, itp.)
            num_match = re.match(r'(\d+)_', stem)
            if num_match:
                base = f"PZ{num_match.group(1)}"
            
            # Wyciągnij lab i num
            # Format: 10_Bujnicki_x_rpr
            pattern = r'^\d+_([A-Za-z0-9\-]+)_(\d+)_'
            match = re.match(pattern, stem)
            if match:
                lab_name, num_str = match.groups()
                lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
            
            # Format: 5_adamiak_x_rpr_refined
            pattern_refined = r'^\d+_([a-z]+)_(\d+)_rpr_refined'
            match_refined = re.match(pattern_refined, stem)
            if match_refined:
                lab_name, num_str = match_refined.groups()
                lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
        
        # 3. Dla formatów PZ z prefiksem: PZ16b_3dRNAAS2_x_refined
        elif stem.startswith('PZ'):
            # Wyciągnij bazę (PZ16b, PZ22, PZ34, itp.)
            pz_match = re.match(r'(PZ\d+[a-z]*)', stem)
            if pz_match:
                base = pz_match.group(1)
            
            # Format: PZ16b_3dRNAAS2_x_refined
            pattern1 = r'PZ\d+[a-z]*_([A-Za-z0-9]+)_(\d+)_'
            match1 = re.match(pattern1, stem)
            if match1:
                lab_name, num_str = match1.groups()
                lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
            
            # Format: PZ22Dimer_Adamiak_x_refined
            pattern2 = r'(PZ\d+[a-z]*)_([A-Za-z0-9]+)_(\d+)_'
            match2 = re.match(pattern2, stem)
            if match2:
                _, lab_name, num_str = match2.groups()
                lab = lab_name
                try:
                    num = int(num_str)
                except:
                    pass
            
            # Format: PZ16b_solution_x
            pattern_sol = r'PZ\d+[a-z]*_solution_(\d+)'
            match_sol = re.match(pattern_sol, stem)
            if match_sol:
                lab = "solution"
                try:
                    num = int(match_sol.group(1))
                except:
                    pass
        
        # 4. Dla formatów z solution w środku: solution_x_AdamiakPostExp_x_rpr
        elif 'solution' in stem.lower():
            # To jest plik rozwiązania
            lab = "solution"
            
            # Spróbuj wyciągnąć numer z solution
            sol_match = re.search(r'solution_(\d+)', stem)
            if sol_match:
                try:
                    num = int(sol_match.group(1))
                except:
                    pass
            
            # Spróbuj wyciągnąć bazę z początku
            if stem.startswith('CR'):
                base_match = re.match(r'(CR\d+)', stem)
                if base_match:
                    base = base_match.group(1)
            elif stem[0].isdigit():
                num_match = re.match(r'(\d+)', stem)
                if num_match:
                    base = f"PZ{num_match.group(1)}"
        
        # 5. Fallback: spróbuj wyciągnąć numer jako ostatnią liczbę
        if num == 0:
            numbers = re.findall(r'\d+', stem)
            if numbers:
                try:
                    num = int(numbers[-1])
                except:
                    pass
        
        # 6. Fallback dla lab: jeśli puste, spróbuj wyciągnąć z nazwy
        if not lab:
            # Sprawdź typowe nazwy laboratoriów
            lab_patterns = [
                'Adamiak', 'Bujnicki', 'Chen', 'Das', 'Dokholyan', 'Ding', 'Major',
                'Xiao', 'Flores', 'GWxraylab', 'Manifold', 'UltraFold', 'Rookie',
                'RNApolis', 'UNRES', 'Graphen', 'Kiharalab', 'SHT', 'nucE2E',
                'GinobiFold', 'Yang', 'SoutheRNA', 'rDP', 'FoldEver', 'PerezLab',
                'LCBio', 'AIchemy', 'Coqualia', 'CoMMiT', '3dRNA', 'RNAComposer',
                'SimRNA', 'RW3D', 'Lee', 'Nikolay', 'YagoubAli', 'santalucia',
                'wildauer', 'mikolajczak', 'Kevin', 'Ares', 'Dfold', 'Nithin',
                'Perez', 'ARES', 'SWM', 'DiMaio', 'Zhou', 'Yang', 'Vfold'
            ]
            
            for pattern in lab_patterns:
                if pattern.lower() in stem.lower():
                    lab = pattern
                    break
            
            # Jeśli nadal puste, użyj drugiej części po podziale
            if not lab:
                parts = stem.split('_')
                if len(parts) > 1:
                    # Unikaj części, które są liczbami lub słowami kluczowymi
                    candidate = parts[1]
                    if (not candidate.isdigit() and 
                        candidate.lower() not in ['rpr', 'refined', 'solution', 'model'] and
                        len(candidate) > 1):
                        lab = candidate
        
        # 7. Oczyszczenie lab
        if lab.lower() in ['rpr', 'refined']:
            lab = ""
        
        if self.debug:
            self.log("DEBUG", f"    Wynik: base={base}, lab={lab}, num={num}")
        
        return base, lab, num

    def find_puzzle_sets(self, selected_puzzles: Optional[List[str]] = None) -> Dict:
        """
        Znajduje wszystkie zestawy puzzli z obsługą wyboru.
        Ulepszona wersja: szuka również w różnych podstrukturach katalogów.
        """
        self.log("INFO", f"🔍 Szukanie zestawów puzzli w: {self.base_dir}")
        
        if not self.base_dir.exists():
            self.log("ERROR", f"Katalog {self.base_dir} nie istnieje!")
            return {}
        
        puzzle_sets = {}
        
        # Znajdź wszystkie katalogi w base_dir
        puzzle_dirs = [d for d in self.base_dir.iterdir() if d.is_dir()]
        
        if selected_puzzles:
            # Rozszerz o możliwe warianty nazw
            selected_set = set(selected_puzzles)
            puzzle_dirs = [d for d in puzzle_dirs if d.name in selected_set]
            self.log("INFO", f"Wybrano zestawy: {', '.join(selected_puzzles)}")
        
        for puzzle_dir in sorted(puzzle_dirs):
            puzzle_name = puzzle_dir.name
            
            self.log("DEBUG", f"Analizowanie zestawu: {puzzle_name}")
            
            # Znajdź wszystkie pliki PDB w całej strukturze katalogu
            pdb_files = self.find_pdb_files(puzzle_dir)
            
            if not pdb_files:
                self.log("WARNING", f"  Brak plików PDB w {puzzle_name}")
                
                # Sprawdź czy są pliki w podkatalogach
                subdirs = [d for d in puzzle_dir.rglob("*") if d.is_dir()]
                for subdir in subdirs[:3]:  # Sprawdź tylko pierwsze 3 podkatalogi
                    sub_files = self.find_pdb_files(subdir)
                    if sub_files:
                        self.log("INFO", f"    Znaleziono {len(sub_files)} plików w podkatalogu: {subdir.relative_to(puzzle_dir)}")
                        pdb_files.extend(sub_files)
                
                if not pdb_files:
                    puzzle_sets[puzzle_name] = {
                        'path': puzzle_dir,
                        'files': [],
                        'models': [],
                        'solutions': [],
                        'status': 'NO_PDB_FILES'
                    }
                    continue
            
            # Rozdziel pliki na modele i rozwiązania
            models = []
            solutions = []
            
            for file in pdb_files:
                name_lower = file.name.lower()
                
                # SPRAWDZENIE ROZWIĄZANIA PIERWSZE - PRIORYTET
                # Rozwiązania to pliki z '_solution', 'solution_' lub 'solution_0'
                is_solution = any(keyword in name_lower for keyword in [
                    '_solution', 'solution_', '.solution'
                ]) or 'solution_0' in name_lower
                
                # Modele to pliki z '_refined', '_rpr', '_model' lub nazwy laboratoriów BEZ 'solution'
                is_model = (any(keyword in name_lower for keyword in [
                    '_refined', '_rpr', '_model', '_prediction'
                ]) and not is_solution) or \
                    (any(lab_keyword in name_lower for lab_keyword in [
                        'adamiak', 'bujnicki', 'chen', 'das', 'dokholyan', 
                        'major', 'xiao', 'flores', 'ding', 'gw'
                    ]) and not is_solution)
                
                # Logika priorytetu: jeśli to rozwiązanie -> do solutions
                if is_solution:
                    solutions.append(file)
                # Jeśli to model -> do models
                elif is_model:
                    models.append(file)
                else:
                    # Domyślnie traktuj jako model (ostrożnie)
                    models.append(file)
            
            self.log("INFO", f"  {puzzle_name}: {len(pdb_files)} PDB ({len(models)} modele, {len(solutions)} rozwiązania)")
            
            if self.debug:
                if models:
                    self.log("DEBUG", "  Przykładowe modele:")
                    for model in models[:3]:
                        self.log("DEBUG", f"    - {model.name}")
                    if len(models) > 3:
                        self.log("DEBUG", f"    ... i {len(models) - 3} więcej")
                
                if solutions:
                    self.log("DEBUG", "  Rozwiązania:")
                    for solution in solutions[:3]:
                        self.log("DEBUG", f"    - {solution.name}")
                    if len(solutions) > 3:
                        self.log("DEBUG", f"    ... i {len(solutions) - 3} więcej")
            
            puzzle_sets[puzzle_name] = {
                'path': puzzle_dir,
                'files': pdb_files,
                'models': models,
                'solutions': solutions,
                'status': 'READY' if models else 'NO_MODELS'
            }
        
        return puzzle_sets
    
    def debug_directory_structure(self, puzzle_name: str):
        """
        Szczegółowe debugowanie struktury katalogów dla zestawu.
        """
        puzzle_dir = self.base_dir / puzzle_name
        
        if not puzzle_dir.exists():
            self.log("ERROR", f"Katalog {puzzle_dir} nie istnieje!")
            return
        
        self.log("INFO", f"\n🔍 DEBUG STRUKTURY KATALOGU: {puzzle_name}")
        self.log("INFO", f"Ścieżka: {puzzle_dir}")
        
        # Zbadaj strukturę katalogu
        import os
        for root, dirs, files in os.walk(puzzle_dir):
            level = root.replace(str(puzzle_dir), '').count(os.sep)
            indent = ' ' * 2 * level
            self.log("INFO", f"{indent}{os.path.basename(root)}/")
            subindent = ' ' * 2 * (level + 1)
            
            # Pokaż pliki PDB
            pdb_files = [f for f in files if f.lower().endswith(('.pdb', '.pdb1', '.ent', '.cif')) or 
                        any(keyword in f.lower() for keyword in ['_rpr', '_refined', '_solution'])]
            
            for file in pdb_files[:10]:  # Pokaż tylko pierwsze 10
                self.log("INFO", f"{subindent}{file}")
            
            if len(pdb_files) > 10:
                self.log("INFO", f"{subindent}... i {len(pdb_files) - 10} więcej")
        
        # Znajdź wszystkie pliki
        all_files = list(puzzle_dir.rglob("*"))
        pdb_files = self.find_pdb_files(puzzle_dir)
        
        self.log("INFO", f"\n📊 STATYSTYKI:")
        self.log("INFO", f"  Wszystkich plików: {len(all_files)}")
        self.log("INFO", f"  Plików PDB: {len(pdb_files)}")
        
        if pdb_files:
            self.log("INFO", f"\n📁 PRZYKŁADOWE PLIKI PDB:")
            for i, file in enumerate(pdb_files[:20]):
                base, lab, num = self.parse_filename(file.name)
                self.log("INFO", f"  {i+1:2d}. {file.name}")
                self.log("INFO", f"      → base={base}, lab={lab}, num={num}")
            
            if len(pdb_files) > 20:
                self.log("INFO", f"      ... i {len(pdb_files) - 20} więcej")
    
    def analyze_file_structure(self, puzzle_name: str):
        """
        Analizuje strukturę plików w zestawie, aby zrozumieć format nazw.
        
        Args:
            puzzle_name: Nazwa zestawu
        """
        if puzzle_name not in self.puzzle_sets:
            self.log("ERROR", f"Zestaw {puzzle_name} nie został znaleziony")
            return
        
        puzzle_info = self.puzzle_sets[puzzle_name]
        
        self.log("INFO", f"\n🔍 ANALIZA STRUKTURY PLIKÓW: {puzzle_name}")
        self.log("INFO", f"Liczba plików: {len(puzzle_info['files'])}")
        
        # Analizuj nazwy plików
        self.log("INFO", "\nPrzykładowe nazwy plików:")
        for i, file in enumerate(puzzle_info['files'][:10], 1):
            base, lab, num = self.parse_filename(file.name)
            self.log("INFO", f"  {i}. {file.name}")
            self.log("INFO", f"     → base={base}, lab={lab}, num={num}")
        
        if len(puzzle_info['files']) > 10:
            self.log("INFO", f"  ... i {len(puzzle_info['files']) - 10} więcej")
        
        # Statystyki lab
        self.log("INFO", "\nStatystyki laboratoriów:")
        labs = defaultdict(int)
        for file in puzzle_info['files']:
            _, lab, _ = self.parse_filename(file.name)
            labs[lab] += 1
        
        for lab, count in sorted(labs.items()):
            self.log("INFO", f"  {lab}: {count} plików")
        
        # Rozdzielenie na modele i rozwiązania
        models = puzzle_info['models']
        solutions = puzzle_info['solutions']
        
        self.log("INFO", f"\nPodział:")
        self.log("INFO", f"  Modele: {len(models)}")
        self.log("INFO", f"  Rozwiązania: {len(solutions)}")
        
        if solutions:
            self.log("INFO", "\nPrzykładowe rozwiązania:")
            for i, sol in enumerate(solutions[:5], 1):
                self.log("INFO", f"  {i}. {sol.name}")
        
        if models:
            self.log("INFO", "\nPrzykładowe modele:")
            for i, model in enumerate(models[:5], 1):
                base, lab, num = self.parse_filename(model.name)
                self.log("INFO", f"  {i}. {model.name}")
                self.log("INFO", f"     → lab={lab}, num={num}")
    
    def find_matching_solution(self, model_file: Path, solutions: List[Path]) -> Optional[Path]:
        """
        Znajduje odpowiednie rozwiązanie dla modelu.
        
        Args:
            model_file: Plik modelu
            solutions: Lista plików rozwiązania
            
        Returns:
            Pasujące rozwiązanie lub None
        """
        model_base, model_lab, model_num = self.parse_filename(model_file.name)
        
        if self.debug:
            self.log("DEBUG", f"    Szukanie rozwiązania dla: {model_file.name}")
            self.log("DEBUG", f"      base={model_base}, lab={model_lab}, num={model_num}")
        
        # Priorytet 1: dokładne dopasowanie
        for solution in solutions:
            sol_base, sol_lab, sol_num = self.parse_filename(solution.name)
            if sol_base == model_base and sol_lab == model_lab and sol_num == model_num:
                if self.debug:
                    self.log("DEBUG", f"      Znaleziono dokładne dopasowanie: {solution.name}")
                return solution
        
        # Priorytet 2: dopasowanie base i lab
        for solution in solutions:
            sol_base, sol_lab, sol_num = self.parse_filename(solution.name)
            if sol_base == model_base and sol_lab == model_lab:
                if self.debug:
                    self.log("DEBUG", f"      Znaleziono dopasowanie base+lab: {solution.name}")
                return solution
        
        # Priorytet 3: solution_0 dla tego samego base
        for solution in solutions:
            sol_base, sol_lab, sol_num = self.parse_filename(solution.name)
            if sol_base == model_base and 'solution' in sol_lab.lower() and sol_num == 0:
                if self.debug:
                    self.log("DEBUG", f"      Znaleziono solution_0: {solution.name}")
                return solution
        
        # Priorytet 4: jakiekolwiek solution dla tego samego base
        for solution in solutions:
            sol_base, sol_lab, sol_num = self.parse_filename(solution.name)
            if sol_base == model_base and 'solution' in solution.name.lower():
                if self.debug:
                    self.log("DEBUG", f"      Znaleziono solution: {solution.name}")
                return solution
        
        # Priorytet 5: pierwsze rozwiązanie
        if solutions:
            if self.debug:
                self.log("DEBUG", f"      Użycie pierwszego rozwiązania: {solutions[0].name}")
            return solutions[0]
        
        return None
    
    def compute_INF_safe(self, solution_path: Path, model_path: Path) -> Dict[str, Any]:
        """
        Bezpieczne obliczanie INF z obsługą błędów.
        
        Args:
            solution_path: Ścieżka do rozwiązania
            model_path: Ścieżka do modelu
            
        Returns:
            Słownik z wynikami
        """
        try:
            if not INF_AVAILABLE:
                raise ImportError("Moduł INF nie jest dostępny")
            
            if not solution_path.exists():
                raise FileNotFoundError(f"Rozwiązanie nie istnieje: {solution_path}")
            
            if not model_path.exists():
                raise FileNotFoundError(f"Model nie istnieje: {model_path}")
            
            result = compute_INF(str(solution_path), str(model_path))
            
            # Dodaj metadane
            result['solution_exists'] = True
            result['model_exists'] = True
            result['computation_success'] = True
            
            return result
            
        except Exception as e:
            if self.debug:
                self.log("ERROR", f"    Błąd compute_INF: {e}")
            
            return {
                'solution_exists': solution_path.exists(),
                'model_exists': model_path.exists(),
                'computation_success': False,
                'error': str(e),
                'inf_all': None,
                'inf_wc': None,
                'inf_nwc': None,
                'inf_stack': None,
                'rmsd': None,
                'aligned_bases': 0,
                'total_bases': 0
            }
    
    def process_single_pair(self, model_path: Path, solution_path: Path) -> Dict[str, Any]:
        """
        Przetwarza pojedynczą parę model-rozwiązanie.
        
        Args:
            model_path: Ścieżka do modelu
            solution_path: Ścieżka do rozwiązania
            
        Returns:
            Słownik z wynikami
        """
        if self.debug:
            self.log("DEBUG", f"    Porównanie: {model_path.name} vs {solution_path.name}")
        
        try:
            # Oblicz INF
            result = self.compute_INF_safe(solution_path, model_path)
            
            # Dodaj metadane
            result['model_file'] = model_path.name
            result['solution_file'] = solution_path.name
            
            # Parsuj nazwy dla lepszego formatowania
            base, lab, num = self.parse_filename(model_path.name)
            result['lab'] = lab
            result['num'] = num
            result['base'] = base
            
            # Dodaj ścieżki
            result['model_path'] = str(model_path)
            result['solution_path'] = str(solution_path)
            
            # Sprawdź czy INF został obliczony
            if result.get('inf_all') is None:
                result['computation_success'] = False
                result['error'] = result.get('error', 'INF computation returned None')
            
            if self.debug and result.get('computation_success'):
                self.log("DEBUG", f"      Wynik: INF-all={result.get('inf_all'):.4f}, RMSD={result.get('rmsd'):.4f}")
            
            return result
            
        except Exception as e:
            if self.debug:
                self.log("ERROR", f"    Błąd przetwarzania pary: {e}")
            
            base, lab, num = self.parse_filename(model_path.name)
            return {
                'model_file': model_path.name,
                'solution_file': solution_path.name if solution_path else None,
                'lab': lab,
                'num': num,
                'base': base,
                'error': str(e),
                'computation_success': False,
                'rmsd': None,
                'inf_all': None,
                'inf_wc': None,
                'inf_nwc': None,
                'inf_stack': None,
                'model_path': str(model_path),
                'solution_path': str(solution_path) if solution_path else None
            }
    
    def process_puzzle_set(self, puzzle_name: str, puzzle_info: Dict) -> Tuple[List[Dict], int, int]:
        """
        Przetwarza pojedynczy zestaw puzzli.
        
        Args:
            puzzle_name: Nazwa zestawu
            puzzle_info: Informacje o zestawie
            
        Returns:
            Tuple: (wyniki, liczba modeli, liczba sukcesów)
        """
        self.log("INFO", f"\n{'='*60}")
        self.log("INFO", f"Zestaw: {puzzle_name}")
        self.log("INFO", f"Ścieżka: {puzzle_info['path']}")
        self.log("INFO", f"Status: {puzzle_info['status']}")
        self.log("INFO", f"{'='*60}")
        
        if puzzle_info['status'] != 'READY':
            self.log("WARNING", f"  Pomijanie zestawu w stanie: {puzzle_info['status']}")
            return [], 0, 0
        
        models = puzzle_info['models']
        solutions = puzzle_info['solutions']
        
        if not models:
            self.log("WARNING", "  Brak modeli do przetworzenia")
            return [], 0, 0
        
        if not solutions:
            self.log("WARNING", "  Brak plików rozwiązania")
            return [], 0, 0
        
        self.log("INFO", f"  Modele: {len(models)}, Rozwiązania: {len(solutions)}")
        
        if self.debug:
            self.log("DEBUG", "  Modele:")
            for model in models[:5]:  # Pokaż tylko pierwsze 5
                self.log("DEBUG", f"    - {model.name}")
            if len(models) > 5:
                self.log("DEBUG", f"    ... i {len(models) - 5} więcej")
            
            self.log("DEBUG", "  Rozwiązania:")
            for solution in solutions[:5]:
                self.log("DEBUG", f"    - {solution.name}")
            if len(solutions) > 5:
                self.log("DEBUG", f"    ... i {len(solutions) - 5} więcej")
        
        results = []
        success_count = 0
        fail_count = 0
        
        for i, model in enumerate(models, 1):
            self.log("INFO", f"\n  [{i}/{len(models)}] Przetwarzanie: {model.name}")
            
            # Znajdź odpowiednie rozwiązanie
            solution = self.find_matching_solution(model, solutions)
            
            if not solution:
                self.log("WARNING", f"    ✗ Nie znaleziono odpowiedniego rozwiązania")
                fail_count += 1
                
                # Dodaj wynik z błędem
                base, lab, num = self.parse_filename(model.name)
                results.append({
                    'model_file': model.name,
                    'solution_file': None,
                    'lab': lab,
                    'num': num,
                    'base': base,
                    'error': 'No matching solution found',
                    'computation_success': False,
                    'puzzle_name': puzzle_name,
                    'model_path': str(model)
                })
                continue
            
            # Przetwórz parę
            result = self.process_single_pair(model, solution)
            
            # Dodaj informacje o zestawie
            result['puzzle_name'] = puzzle_name
            
            if result.get('computation_success') and result.get('inf_all') is not None:
                success_count += 1
                self.log("INFO", f"    ✓ INF-all: {result.get('inf_all', 0):.4f}, RMSD: {result.get('rmsd', 0):.4f}")
            else:
                fail_count += 1
                error_msg = result.get('error', 'Unknown error')
                self.log("WARNING", f"    ✗ Błąd: {error_msg}")
            
            results.append(result)
        
        # Podsumowanie zestawu
        self.log("INFO", f"\n  Podsumowanie {puzzle_name}:")
        self.log("INFO", f"    • Przetworzone modele: {len(models)}")
        self.log("INFO", f"    • Udane obliczenia: {success_count}")
        self.log("INFO", f"    • Nieudane obliczenia: {fail_count}")
        
        if success_count > 0:
            success_rate = (success_count / len(models)) * 100
            self.log("INFO", f"    • Wskaźnik sukcesu: {success_rate:.1f}%")
        
        return results, len(models), success_count
    
    def save_puzzle_results(self, puzzle_name: str, results: List[Dict]) -> Tuple[Optional[Path], int]:
        """
        Zapisuje wyniki dla pojedynczego zestawu do CSV.
        Poprawione formatowanie: lab nie może być "refined" ani "rpr".
        
        Args:
            puzzle_name: Nazwa zestawu
            results: Lista wyników
            
        Returns:
            Tuple: (ścieżka do pliku, liczba zapisanych wyników)
        """
        # Filtruj tylko wyniki z poprawnym INF-all
        valid_results = []
        for result in results:
            if result.get('computation_success') and result.get('inf_all') is not None:
                # Popraw lab jeśli to "refined" lub "rpr"
                lab = result.get('lab', '')
                if lab.lower() in ['refined', 'rpr']:
                    # Spróbuj wyciągnąć prawdziwe lab z nazwy pliku
                    model_name = result.get('model_file', '')
                    lab_patterns = [
                        r'TS\d+', r'TSR\d+', r'GW[a-z]+', r'Das', 
                        r'Adamiak', r'Manifold', r'AIchemy', r'SimRNA',
                        r'Rosetta', r'MC', r'Vienna', r'RNAComposer'
                    ]
                    for pattern in lab_patterns:
                        match = re.search(pattern, model_name, re.IGNORECASE)
                        if match:
                            lab = match.group(0)
                            result['lab'] = lab
                            break
                    else:
                        lab = ''  # Puste jeśli nie znaleziono
                
                # Upewnij się, że lab nie jest pusty
                if not lab or lab == 'Unknown':
                    lab = ''
                
                result['lab'] = lab
                valid_results.append(result)
        
        if not valid_results:
            self.log("WARNING", f"  Brak poprawnych wyników do zapisania dla {puzzle_name}")
            return None, 0
        
        # Sortuj wyniki
        valid_results.sort(key=lambda x: (
            x.get('lab', ''),
            x.get('num', 0)
        ))
        
        # Nazwa pliku CSV
        csv_file = self.inf_all_dir / f"{puzzle_name}.csv"
        
        try:
            with open(csv_file, 'w', newline='', encoding='utf-8') as fh:
                writer = csv.writer(fh)
                
                # Nagłówek
                writer.writerow(['Lab', 'Num', 'INF-all', 'INF-wc', 'INF-nwc', 'INF-stack', 'RMSD', 'Solution'])
                
                for result in valid_results:
                    lab = result.get('lab', '')
                    # Upewnij się, że lab nie jest "refined" ani "rpr"
                    if lab.lower() in ['refined', 'rpr']:
                        lab = ''
                    
                    writer.writerow([
                        lab,
                        result.get('num', ''),
                        f"{result.get('inf_all', 0):.4f}" if result.get('inf_all') else '',
                        f"{result.get('inf_wc', 0):.4f}" if result.get('inf_wc') and result.get('inf_wc') != -1 else '',
                        f"{result.get('inf_nwc', 0):.4f}" if result.get('inf_nwc') and result.get('inf_nwc') != -1 else '',
                        f"{result.get('inf_stack', 0):.4f}" if result.get('inf_stack') and result.get('inf_stack') != -1 else '',
                        f"{result.get('rmsd', 0):.4f}" if result.get('rmsd') else '',
                        result.get('solution_file', '')
                    ])
            
            self.log("INFO", f"  ✓ Zapisano {len(valid_results)} wyników do: {csv_file}")
            return csv_file, len(valid_results)
            
        except Exception as e:
            self.log("ERROR", f"  ✗ Błąd zapisu CSV dla {puzzle_name}: {e}")
            return None, 0
    
    def generate_summary(self) -> Dict:
        """Generuje podsumowanie wszystkich wyników."""
        summary = {
            'total_puzzles': len(self.results),
            'total_models': 0,
            'successful_calculations': 0,
            'failed_calculations': 0,
            'puzzle_details': {},
            'timestamp': datetime.now().isoformat(),
            'execution_time': str(datetime.now() - self.stats['start_time'])
        }
        
        for puzzle_name, puzzle_data in self.results.items():
            total = puzzle_data['total_models']
            successful = puzzle_data['successful']
            
            summary['total_models'] += total
            summary['successful_calculations'] += successful
            summary['failed_calculations'] += (total - successful)
            
            # Oblicz średnie dla poprawnych wyników
            if successful > 0:
                valid_results = [r for r in puzzle_data['results'] 
                               if r.get('computation_success') and r.get('inf_all') is not None]
                
                if valid_results:
                    avg_inf_all = sum(r.get('inf_all', 0) for r in valid_results) / len(valid_results)
                    # Filtruj tylko wyniki z RMSD
                    rmsd_results = [r for r in valid_results if r.get('rmsd') is not None]
                    avg_rmsd = sum(r.get('rmsd', 0) for r in rmsd_results) / len(rmsd_results) if rmsd_results else 0
                else:
                    avg_inf_all = 0
                    avg_rmsd = 0
            else:
                avg_inf_all = 0
                avg_rmsd = 0
            
            summary['puzzle_details'][puzzle_name] = {
                'total_models': total,
                'successful': successful,
                'failed': total - successful,
                'success_rate': f"{(successful/total*100):.1f}%" if total > 0 else "0%",
                'avg_inf_all': round(avg_inf_all, 4),
                'avg_rmsd': round(avg_rmsd, 4),
                'has_results': successful > 0
            }
        
        if summary['total_models'] > 0:
            summary['overall_success_rate'] = f"{(summary['successful_calculations']/summary['total_models']*100):.1f}%"
        else:
            summary['overall_success_rate'] = "0%"
        
        return summary
    
    def save_global_results(self) -> Tuple[Path, Path]:
        """
        Zapisuje wszystkie wyniki do globalnych plików.
        
        Returns:
            Tuple: (ścieżka do JSON, ścieżka do CSV)
        """
        # Zbierz wszystkie poprawne wyniki
        all_valid_results = []
        for puzzle_name, puzzle_data in self.results.items():
            for result in puzzle_data['results']:
                if result.get('computation_success') and result.get('inf_all') is not None:
                    all_valid_results.append(result)
        
        # Zapisz jako CSV
        csv_file = self.output_dir / "INF_all_results.csv"
        if all_valid_results:
            # Sortuj wyniki
            all_valid_results.sort(key=lambda x: (
                x.get('puzzle_name', ''),
                x.get('lab', ''),
                x.get('num', 0)
            ))
            
            try:
                with open(csv_file, 'w', newline='', encoding='utf-8') as fh:
                    writer = csv.writer(fh)
                    writer.writerow(['Puzzle', 'Lab', 'Num', 'INF-all', 'INF-wc', 'INF-nwc', 'INF-stack', 'RMSD', 'Solution'])
                    
                    for result in all_valid_results:
                        writer.writerow([
                            result.get('puzzle_name', ''),
                            result.get('lab', ''),
                            result.get('num', ''),
                            f"{result.get('inf_all', 0):.4f}" if result.get('inf_all') else '',
                            f"{result.get('inf_wc', 0):.4f}" if result.get('inf_wc') and result.get('inf_wc') != -1 else '',
                            f"{result.get('inf_nwc', 0):.4f}" if result.get('inf_nwc') and result.get('inf_nwc') != -1 else '',
                            f"{result.get('inf_stack', 0):.4f}" if result.get('inf_stack') and result.get('inf_stack') != -1 else '',
                            f"{result.get('rmsd', 0):.4f}" if result.get('rmsd') else '',
                            result.get('solution_file', '')
                        ])
                
                self.log("INFO", f"✓ Globalny CSV zapisany: {csv_file} ({len(all_valid_results)} wierszy)")
                
            except Exception as e:
                self.log("ERROR", f"✗ Błąd zapisu globalnego CSV: {e}")
                csv_file = None
        else:
            self.log("WARNING", "✗ Brak poprawnych wyników do globalnego CSV")
            csv_file = None
        
        # Zapisz jako JSON
        json_file = self.output_dir / "INF_all_results.json"
        try:
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(self.results, f, indent=2, ensure_ascii=False, default=str)
            self.log("INFO", f"✓ Globalny JSON zapisany: {json_file}")
        except Exception as e:
            self.log("ERROR", f"✗ Błąd zapisu globalnego JSON: {e}")
            json_file = None
        
        return json_file, csv_file
    
    def print_summary_table(self, summary: Dict):
        """Drukuje podsumowanie w formie tabeli."""
        self.log("INFO", "\n" + "="*80)
        self.log("INFO", "PODSUMOWANIE")
        self.log("INFO", "="*80)
        
        self.log("INFO", f"\nPrzetworzono zestawów: {summary['total_puzzles']}")
        self.log("INFO", f"Łączna liczba modeli: {summary['total_models']}")
        self.log("INFO", f"Udanych obliczeń: {summary['successful_calculations']}")
        self.log("INFO", f"Nieudanych obliczeń: {summary['failed_calculations']}")
        self.log("INFO", f"Ogólny wskaźnik sukcesu: {summary['overall_success_rate']}")
        
        # Tabela szczegółowa
        self.log("INFO", "\nSzczegółowe wyniki per zestaw:")
        self.log("INFO", "-" * 90)
        self.log("INFO", f"{'Zestaw':<15} {'Modele':<6} {'Udane':<6} {'Sukces':<10} {'Śr. INF':<10} {'Śr. RMSD':<10} {'Status':<10}")
        self.log("INFO", "-" * 90)
        
        for puzzle_name, details in sorted(summary['puzzle_details'].items()):
            status = "OK" if details['has_results'] else "BRAK"
            self.log("INFO", f"{puzzle_name:<15} {details['total_models']:<6} {details['successful']:<6} "
                  f"{details['success_rate']:<10} {details['avg_inf_all']:<10.4f} "
                  f"{details['avg_rmsd']:<10.4f} {status:<10}")
    
    def run(self, selected_puzzles: Optional[List[str]] = None, 
            skip_empty: bool = False, force: bool = False):
        """
        Główna metoda uruchamiająca analizę.
        
        Args:
            selected_puzzles: Lista nazw zestawów do przetworzenia
            skip_empty: Pomijaj zestawy bez wyników
            force: Wymuś przetwarzanie nawet jeśli INF nie jest dostępny
        """
        # Sprawdź dostępność INF
        if not INF_AVAILABLE and not force:
            self.log("ERROR", "Moduł INF nie jest dostępny! Użyj --force aby kontynuować.")
            return
        
        self.log("INFO", "="*80)
        self.log("INFO", "SKRYPT DO OBLICZANIA INF-all DLA RNA PUZZLES")
        self.log("INFO", f"Czas rozpoczęcia: {self.stats['start_time']}")
        self.log("INFO", f"Tryb debug: {self.debug}")
        self.log("INFO", "="*80)
        
        # Znajdź zestawy
        self.puzzle_sets = self.find_puzzle_sets(selected_puzzles)
        
        if not self.puzzle_sets:
            self.log("ERROR", "Nie znaleziono żadnych zestawów do przetworzenia!")
            return
        
        # Przetwarzaj zestawy
        self.log("INFO", "\n" + "="*80)
        self.log("INFO", "ROZPOCZĘCIE PRZETWARZANIA")
        self.log("INFO", "="*80)
        
        processed_count = 0
        
        for puzzle_name in sorted(self.puzzle_sets.keys()):
            puzzle_info = self.puzzle_sets[puzzle_name]
            
            # Pomijaj puste zestawy jeśli requested
            if skip_empty and (not puzzle_info['models'] or not puzzle_info['solutions']):
                self.log("INFO", f"\n⏭️  Pomijanie pustego zestawu: {puzzle_name}")
                continue
            
            try:
                results, total_models, successful = self.process_puzzle_set(puzzle_name, puzzle_info)
                
                # Zapisz wyniki per zestaw
                csv_file, saved_count = self.save_puzzle_results(puzzle_name, results)
                
                self.results[puzzle_name] = {
                    'results': results,
                    'total_models': total_models,
                    'successful': successful,
                    'csv_file': str(csv_file) if csv_file else None,
                    'saved_count': saved_count
                }
                
                processed_count += 1
                
                # Krótka przerwa
                time.sleep(0.1)
                
            except Exception as e:
                self.log("ERROR", f"\n❌ Błąd przetwarzania zestawu {puzzle_name}: {e}")
                import traceback
                if self.debug:
                    self.log("DEBUG", traceback.format_exc())
                
                self.results[puzzle_name] = {
                    'results': [],
                    'total_models': 0,
                    'successful': 0,
                    'error': str(e),
                    'csv_file': None,
                    'saved_count': 0
                }
        
        # Generuj i wyświetl podsumowanie
        summary = self.generate_summary()
        self.print_summary_table(summary)
        
        # Zapisz globalne wyniki
        self.log("INFO", "\n" + "="*80)
        self.log("INFO", "ZAPISYWANIE WYNIKÓW")
        self.log("INFO", "="*80)
        
        json_file, csv_file = self.save_global_results()
        
        # Zapisz podsumowanie
        summary_file = self.output_dir / "INF_summary.json"
        try:
            with open(summary_file, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
            self.log("INFO", f"✓ Podsumowanie zapisane: {summary_file}")
        except Exception as e:
            self.log("ERROR", f"✗ Błąd zapisu podsumowania: {e}")
        
        # Koniec
        end_time = datetime.now()
        execution_time = end_time - self.stats['start_time']
        
        self.log("INFO", "\n" + "="*80)
        self.log("INFO", "PRZETWARZANIE ZAKOŃCZONE")
        self.log("INFO", "="*80)
        self.log("INFO", f"Czas zakończenia: {end_time}")
        self.log("INFO", f"Czas wykonania: {execution_time}")
        
        self.log("INFO", "\n📊 FINALNE PODSUMOWANIE:")
        self.log("INFO", f"  • Przetworzone zestawy: {summary['total_puzzles']}")
        self.log("INFO", f"  • Łącznie modeli: {summary['total_models']}")
        self.log("INFO", f"  • Udane obliczenia: {summary['successful_calculations']}")
        self.log("INFO", f"  • Sukces: {summary['overall_success_rate']}")
        
        self.log("INFO", "\n📁 LOKALIZACJA WYNIKÓW:")
        self.log("INFO", f"  • Pliki CSV per zestaw: {self.inf_all_dir}/")
        if csv_file:
            self.log("INFO", f"  • Globalny CSV: {csv_file}")
        if json_file:
            self.log("INFO", f"  • Globalny JSON: {json_file}")
        self.log("INFO", f"  • Podsumowanie: {summary_file}")
        self.log("INFO", f"  • Logi: {self.log_file}")

def main():
    """Główna funkcja z obsługą argumentów wiersza poleceń."""
    parser = argparse.ArgumentParser(
        description='Obliczanie INF-all dla zestawów RNA Puzzles',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Przykłady użycia:
  %(prog)s                          # Przetwórz wszystkie zestawy
  %(prog)s --puzzles CR1116 CR1108  # Przetwórz tylko wybrane zestawy
  %(prog)s --debug                  # Włącz tryb debugowania
  %(prog)s --skip-empty             # Pomijaj puste zestawy
  %(prog)s --list                   # Wyświetl dostępne zestawy
  %(prog)s --debug-structure CR1116 # Debuguj strukturę katalogów
  %(prog)s --analyze-files PZ10     # Analizuj strukturę plików
        """
    )
    
    parser.add_argument('--puzzles', '-p', nargs='+', 
                       help='Lista zestawów do przetworzenia (np. CR1107 CR1108)')
    parser.add_argument('--debug', '-d', action='store_true',
                       help='Włącz tryb debugowania')
    parser.add_argument('--skip-empty', '-s', action='store_true',
                       help='Pomijaj zestawy bez modeli lub rozwiązań')
    parser.add_argument('--force', '-f', action='store_true',
                       help='Wymuś uruchomienie nawet bez modułu INF')
    parser.add_argument('--list', '-l', action='store_true',
                       help='Wyświetl dostępne zestawy i wyjście')
    parser.add_argument('--debug-structure', metavar='PUZZLE',
                       help='Szczegółowe debugowanie struktury katalogu zestawu')
    parser.add_argument('--analyze-files', metavar='PUZZLE',
                       help='Analizuj strukturę plików w zestawie')
    parser.add_argument('--base-dir', default="/work/data/finished",
                       help='Główny katalog z danymi (domyślnie: /work/data/finished)')
    parser.add_argument('--output-dir', default="/work/results",
                       help='Katalog na wyniki (domyślnie: /work/results)')
    
    args = parser.parse_args()
    
    # Sprawdź czy katalog istnieje
    if not Path(args.base_dir).exists():
        print(f"❌ Błąd: Katalog {args.base_dir} nie istnieje!")
        sys.exit(1)
    
    # Opcja listy
    if args.list:
        analyzer = INFAnalyzer(
            base_dir=args.base_dir,
            output_dir=args.output_dir,
            debug=True
        )
        analyzer.puzzle_sets = analyzer.find_puzzle_sets()
        print("\n📁 DOSTĘPNE ZESTAWY:")
        for name, info in sorted(analyzer.puzzle_sets.items()):
            status = info.get('status', 'UNKNOWN')
            models = len(info.get('models', []))
            solutions = len(info.get('solutions', []))
            print(f"  • {name:<15} modele={models:<3} rozwiązania={solutions:<3} status={status}")
        sys.exit(0)
    
    # Opcja debugowania struktury
    if args.debug_structure:
        analyzer = INFAnalyzer(
            base_dir=args.base_dir,
            output_dir=args.output_dir,
            debug=True
        )
        analyzer.debug_directory_structure(args.debug_structure)
        sys.exit(0)
    
    # Opcja analizy plików
    if args.analyze_files:
        analyzer = INFAnalyzer(
            base_dir=args.base_dir,
            output_dir=args.output_dir,
            debug=True
        )
        analyzer.puzzle_sets = analyzer.find_puzzle_sets([args.analyze_files])
        if args.analyze_files in analyzer.puzzle_sets:
            analyzer.analyze_file_structure(args.analyze_files)
        else:
            print(f"Zestaw {args.analyze_files} nie został znaleziony")
        sys.exit(0)
    
    # Uruchom analizę
    try:
        analyzer = INFAnalyzer(
            base_dir=args.base_dir,
            output_dir=args.output_dir,
            debug=args.debug
        )
        analyzer.run(
            selected_puzzles=args.puzzles,
            skip_empty=args.skip_empty,
            force=args.force
        )
    except KeyboardInterrupt:
        print("\n\n⏹️  Przerwano przez użytkownika")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Krytyczny błąd: {e}")
        if args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()