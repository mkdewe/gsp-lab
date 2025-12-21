import unittest
import os
import subprocess
import pandas as pd
import shutil
from pathlib import Path

class TestRNAMapperIntegration(unittest.TestCase):
    
    def setUp(self):
        """Przygotowanie ścieżek do testów."""
        self.current_dir = r"D:\Studia\Projekty\gsp-lab\src\refine_tools\pdb\Incorrects_examples"
        self.test_dir = os.path.join(self.current_dir, "test_output")
        
        # USUŃ STARY KATALOG TESTOWY I UTWÓRZ NOWY
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        os.makedirs(self.test_dir, exist_ok=True)
        
        print(f"📁 Katalog testowy: {self.test_dir}")
        
        # Ścieżki do skryptów
        self.rna_mapper_script = r"D:\Studia\Projekty\gsp-lab\src\refine_tools\rna_mapper.py"
        self.compute_script = r"D:\Studia\Projekty\gsp-lab\src\gsp\compute.py"
        
        # Pliki testowe
        self.target_pdb = os.path.join(self.current_dir, "2_solution_0.pdb")
        self.model_pdb = os.path.join(self.current_dir, "PZ2_Bujnicki_1.pdb")
        
        print(f"🔍 Target: {os.path.exists(self.target_pdb)}")
        print(f"🔍 Model: {os.path.exists(self.model_pdb)}")
        print(f"🔍 RNA Mapper: {os.path.exists(self.rna_mapper_script)}")
        print(f"🔍 Compute: {os.path.exists(self.compute_script)}")
        
        if not all([os.path.exists(self.target_pdb), os.path.exists(self.model_pdb),
                   os.path.exists(self.rna_mapper_script), os.path.exists(self.compute_script)]):
            self.fail("Brak wymaganych plików!")
    
    def run_command_and_copy_files(self, cmd, timeout=300):
        """Uruchamia komendę i KOPIUJE pliki do katalogu testowego."""
        try:
            # Uruchom w katalogu Incorrects_examples (tam gdzie skrypt tworzy pliki)
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                cwd=self.current_dir,
                timeout=timeout,
                encoding='utf-8',
                errors='replace'
            )
            
            # PO URUCHOMIENIU: Skopiuj wszystkie nowo utworzone pliki do test_dir
            self._copy_generated_files()
            
            return result
            
        except Exception as e:
            self.fail(f"Błąd podczas uruchamiania komendy: {e}")
    
    def _copy_generated_files(self):
        """Kopiuje wygenerowane pliki do katalogu testowego."""
        current_path = Path(self.current_dir)
        
        # Lista plików do skopiowania
        files_to_copy = [
            "PZ2_Bujnicki_1_refined.pdb",
            "complete_mapping.txt", 
            "missing_atoms.txt",
            "atom_comparison.log",
            "*gGSP.csv"
        ]
        
        # Dodaj pliki solution jeśli istnieją
        files_to_copy.extend(current_path.glob("solution_*.pdb"))
        
        for pattern in files_to_copy:
            if isinstance(pattern, str):
                # Szukaj plików według wzorca
                for file_path in current_path.glob(pattern):
                    if file_path.is_file():
                        shutil.copy2(file_path, self.test_dir)
                        print(f"📋 Skopiowano: {file_path.name} -> test_output/")
            else:
                # To już jest Path object
                if pattern.is_file():
                    shutil.copy2(pattern, self.test_dir)
                    print(f"📋 Skopiowano: {pattern.name} -> test_output/")
    
    def test_01_rna_mapper_creates_refined_no_solution(self):
        """Test 1: Sprawdza czy rna_mapper tworzy plik refined, ale NIE tworzy solution."""
        print("\n" + "="*50)
        print("TEST 1: RNA MAPPER")
        print("="*50)
        
        cmd = ["python", self.rna_mapper_script, "2_solution_0.pdb", "PZ2_Bujnicki_1.pdb"]
        print(f"Uruchamianie: {' '.join(cmd)}")
        print(f"Katalog roboczy: {self.current_dir}")
        
        result = self.run_command_and_copy_files(cmd)
        
        print("WYNIK RNA MAPPER:")
        print(result.stdout)
        if result.stderr:
            print("BŁĘDY:")
            print(result.stderr)
        
        # SPRAWDŹ CZY PLIKI RZECZYWIŚCIE SĄ W TEST_DIR
        refined_path = os.path.join(self.test_dir, "PZ2_Bujnicki_1_refined.pdb")
        mapping_path = os.path.join(self.test_dir, "complete_mapping.txt")
        atom_log_path = os.path.join(self.test_dir, "atom_comparison.log")
        
        print(f"\n🔍 SPRAWDZANIE PLIKÓW W: {self.test_dir}")
        print(f"   Refined: {os.path.exists(refined_path)}")
        print(f"   Mapping: {os.path.exists(mapping_path)}") 
        print(f"   Atom log: {os.path.exists(atom_log_path)}")
        
        # Wypisz zawartość katalogu testowego
        test_files = os.listdir(self.test_dir)
        print(f"📁 Zawartość test_output: {test_files}")
        
        self.assertTrue(os.path.exists(refined_path), "Plik refined nie został utworzony!")
        print("✅ Plik refined utworzony")
        
        solution_files = [f for f in test_files if f.startswith("solution_") and f.endswith(".pdb")]
        self.assertEqual(len(solution_files), 0, "Plik solution został utworzony, ale nie powinien!")
        print("✅ Plik solution nie utworzony (poprawnie)")
        
        # Zapisz ścieżkę dla następnych testów
        self.refined_pdb_path = refined_path
        
        return refined_path
    
    def test_02_compute_script_creates_csv_with_correct_gGSP(self):
        """Test 2: Sprawdza czy compute.py tworzy plik CSV z oczekiwaną wartością gGSP."""
        print("\n" + "="*50)
        print("TEST 2: COMPUTE SCRIPT")
        print("="*50)
        
        # Najpierw uruchom test_01, aby mieć plik refined
        refined_path = self.test_01_rna_mapper_creates_refined_no_solution()
        
        # Użyj ścieżki do pliku refined w test_dir
        refined_relative_path = os.path.join("test_output", "PZ2_Bujnicki_1_refined.pdb")
        
        cmd = ["python", self.compute_script, "-t", "2_solution_0.pdb", "-m", refined_relative_path]
        print(f"Uruchamianie: {' '.join(cmd)}")
        print(f"Katalog roboczy: {self.current_dir}")
        
        result = self.run_command_and_copy_files(cmd)
        
        print("WYNIK COMPUTE SCRIPT:")
        print(result.stdout)
        if result.stderr:
            print("BŁĘDY:")
            print(result.stderr)
        
        # Znajdź plik CSV w test_dir
        csv_files = [f for f in os.listdir(self.test_dir) if f.endswith("_C1'-gGSP.csv")]
        print(f"📊 Znalezione pliki CSV: {csv_files}")
        
        self.assertGreater(len(csv_files), 0, "Nie znaleziono pliku CSV!")
        
        csv_path = os.path.join(self.test_dir, csv_files[0])
        print(f"✅ Znaleziono plik CSV: {csv_files[0]}")
        
        # Sprawdź zawartość CSV
        df = pd.read_csv(csv_path)
        gGSP_value = df['gGSP'].iloc[0]
        expected_gGSP = 99.963
        
        print(f"Wartość gGSP: {gGSP_value}")
        print(f"Oczekiwana wartość: {expected_gGSP}")
        
        self.assertAlmostEqual(gGSP_value, expected_gGSP, delta=0.1, 
                              msg=f"gGSP {gGSP_value} nie jest bliskie oczekiwanej {expected_gGSP}")
        print("✅ Wartość gGSP poprawna")
        
        return csv_path


if __name__ == "__main__":
    print("🎯 TEST INTEGRACYJNY RNA MAPPER + COMPUTE SCRIPT")
    print("Testuje:")
    print("1. RNA Mapper → tworzy refined.pdb (bez solution)")
    print("2. Compute Script → tworzy CSV z gGSP ≈ 99.963") 
    print("3. RZECZYWISTE przenoszenie plików do test_output")
    print("="*60)
    
    # Utwórz test suite i uruchom tylko test_02 (który zawiera pełny przepływ)
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromName('test_02_compute_script_creates_csv_with_correct_gGSP', TestRNAMapperIntegration)
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Dodatkowy raport
    print("\n" + "="*60)
    print("PODSUMOWANIE TESTÓW")
    print("="*60)
    
    if result.wasSuccessful():
        print("✅ WSZYSTKIE TESTY ZALICZONE!")
        current_dir = r"D:\Studia\Projekty\gsp-lab\src\refine_tools\pdb\Incorrects_examples"
        test_output_dir = os.path.join(current_dir, "test_output")
        print(f"📁 Wszystkie pliki zostały zapisane w: {test_output_dir}")
        
        # Pokaż zawartość katalogu test_output
        if os.path.exists(test_output_dir):
            print("\n📋 Zawartość katalogu test_output:")
            for item in sorted(os.listdir(test_output_dir)):
                item_path = os.path.join(test_output_dir, item)
                if os.path.isfile(item_path):
                    size = os.path.getsize(item_path)
                    print(f"   📄 {item} ({size} bytes)")
    
    exit(0 if result.wasSuccessful() else 1)