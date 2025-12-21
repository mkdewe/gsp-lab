#!/usr/bin/env python3
"""Główny skrypt RNA Mapper"""

import argparse
import os
import sys
from pathlib import Path

# Dodaj ścieżkę do modułów
current_dir = Path(__file__).parent
sys.path.insert(0, str(current_dir))

# Importy bezwzględne
from core.mapper import RNAMapper
from io_utils.pdb_handler import PDBHandler, AtomChecker
from io_utils.file_writer import FileWriter
from config import DEFAULT_OUTPUT_DIR, SOLUTION_PREFIX, REFINED_PREFIX
from utils.logger import Logger
from exceptions import RNAMapperError

def create_output_names(model_path, output_dir):
    model_stem = Path(model_path).stem
    if model_stem.startswith('PZ2_'): 
        model_suffix = model_stem[4:]
    else: 
        model_suffix = model_stem
    
    return {
        'refined': output_dir / f"{REFINED_PREFIX}_{model_suffix}.pdb",
        'solution': output_dir / f"{SOLUTION_PREFIX}_{model_suffix}.pdb",
        'mapping': output_dir / f"{model_suffix}_mapping.txt",
        'missing_atoms': output_dir / f"{model_suffix}_missing_atoms.txt"
    }

def main():
    parser = argparse.ArgumentParser(description='RNA Mapper')
    parser.add_argument('-t', '--target', required=True, help='Plik target PDB')
    parser.add_argument('-m', '--model', required=True, help='Plik model PDB')
    parser.add_argument('-o', '--output', default=DEFAULT_OUTPUT_DIR, help='Folder wynikowy')
    parser.add_argument('-l', '--log', action='store_true', help='Włącz logowanie')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.target):
        print(f"ERROR: Target '{args.target}' nie istnieje!")
        sys.exit(1)
    if not os.path.exists(args.model):
        print(f"ERROR: Model '{args.model}' nie istnieje!")
        sys.exit(1)
    
    Logger.set_enabled(args.log)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"=== RNA MAPPER ===")
    print(f"Target: {args.target}")
    print(f"Model: {args.model}") 
    print(f"Output: {output_dir}")
    
    try:
        mapper = RNAMapper()
        results = mapper.process(args.target, args.model)
        
        output_files = create_output_names(args.model, output_dir)
        pdb_handler = PDBHandler()
        
        # Sprawdź brakujące atomy
        missing_atoms = []
        for (m_chain, m_res), (t_chain, t_res) in results['residue_mapping'].items():
            m_res_obj = next((r for r in results['model_structure'][0][m_chain] if r.id[1] == m_res), None)
            t_res_obj = next((r for r in results['target_structure'][0][t_chain] if r.id[1] == t_res), None)
            if m_res_obj and t_res_obj:
                missing = AtomChecker.check_missing_atoms(t_res_obj, m_res_obj)
                if missing: 
                    missing_atoms.append({
                        'model_chain': m_chain, 'model_residue': m_res,
                        'target_chain': t_chain, 'target_residue': t_res,
                        'missing_atoms': missing
                    })
        
        # Utwórz pliki
        chain_order = mapper.structure_parser.get_chain_order(results['target_structure'])
        refined_count = pdb_handler.create_refined_pdb(
            args.model, output_files['refined'], results['residue_mapping'], chain_order)
        
        # POPRAWIONY WARUNEK: Solution tylko gdy są brakujące atomy LUB nie wszystkie reszty są zmapowane
        target_rna_count = mapper.structure_parser.count_rna_residues(results['target_structure'])
        mapped_residues_count = len(results['residue_mapping'])
        
        # Solution tworzymy tylko gdy:
        # 1. Są brakujące atomy w modelu LUB
        # 2. Nie udało się zmapować wszystkich reszt RNA z targetu
        need_solution = missing_atoms or (mapped_residues_count < target_rna_count)
        
        if need_solution:
            solution_count = pdb_handler.create_solution_pdb(
                args.target, output_files['solution'], results['residue_mapping'])
            print(f"Solution PDB: {output_files['solution']} (utworzony z powodu: {'brakujące atomy' if missing_atoms else 'niepełne mapowanie reszt'})")
        else:
            print("Solution PDB: nie utworzono (brak brakujących atomów i pełne mapowanie reszt)")
        
        # Zapisz wyniki
        metadata = {'target': args.target, 'model': args.model}
        FileWriter.save_mapping(output_files['mapping'], results['residue_mapping'], metadata)
        FileWriter.save_missing_atoms(output_files['missing_atoms'], missing_atoms)
        
        # Podsumowanie
        print(f"\n=== PODSUMOWANIE ===")
        print(f"Zmapowane reszty: {mapped_residues_count}/{target_rna_count}")
        print(f"Naprawione atomy: {refined_count}")
        print(f"Brakujące atomy: {len(missing_atoms)}")
        print(f"Plik wyjściowy: {output_files['refined']}")
        
    except RNAMapperError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Nieoczekiwany błąd: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()