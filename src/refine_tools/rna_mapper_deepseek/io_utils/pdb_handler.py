"""Operacje na plikach PDB"""

from pathlib import Path
from collections import defaultdict
from config import ESSENTIAL_ATOMS, SOLUTION_PREFIX, REFINED_PREFIX
from utils.logger import Logger

class PDBHandler:
    @staticmethod
    def create_refined_pdb(model_path, output_path, residue_mapping, target_chain_order):
        """Tworzy naprawiony PDB, zachowując oryginalny format BEZ DODAWANIA NAGŁÓWKÓW"""
        Logger.info(f"Tworzenie: {output_path}")
        with open(model_path, 'r') as f:
            lines = f.readlines()
        
        # Przetwórz atomy - zachowujemy tylko zmapowane reszty
        atom_lines = []
        current_res, has_mapping = None, False
        
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain, resid = PDBHandler._extract_residue_info(line)
                
                if current_res != (chain, resid):
                    current_res = (chain, resid)
                    has_mapping = (chain, resid) in residue_mapping
                
                if has_mapping:
                    new_chain, new_resid = residue_mapping[(chain, resid)]
                    new_line = line[:21] + new_chain + line[22:]
                    new_resid_str = str(new_resid).rjust(4)
                    new_line = new_line[:22] + new_resid_str + new_line[26:]
                    atom_lines.append(new_line)
        
        # Posortuj i ponumeruj atomy
        sorted_atoms = PDBHandler._sort_atoms(atom_lines, target_chain_order)
        
        # Nadaj nowe numery atomów TYLKO dla linii ATOM/HETATM
        atom_count = 1
        final_lines = []
        for line in sorted_atoms:
            if line.startswith(('ATOM', 'HETATM')):
                atom_num_str = str(atom_count).rjust(5)
                new_line = line[:6] + atom_num_str + line[11:]
                final_lines.append(new_line)
                atom_count += 1
            else:
                final_lines.append(line)
        
        # Zapisz BEZ DODAWANIA ŻADNYCH NAGŁÓWKÓW
        with open(output_path, 'w') as f:
            f.writelines(final_lines)
        
        return atom_count - 1
    
    @staticmethod
    def create_solution_pdb(target_path, output_path, residue_mapping):
        """Tworzy solution PDB tylko z zmapowanymi resztami, BEZ DODAWANIA NAGŁÓWKÓW"""
        Logger.info(f"Tworzenie solution: {output_path}")
        with open(target_path, 'r') as f:
            lines = f.readlines()
        
        # Odwróć mapowanie
        reverse_mapping = {v: k for k, v in residue_mapping.items()}
        target_residues_to_keep = set(reverse_mapping.keys())
        
        atom_lines = []
        current_residue = None
        current_residue_should_keep = False
        
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                try:
                    resid = int(line[22:26].strip())
                except ValueError:
                    continue
                
                # Sprawdź czy zmieniła się reszta
                if current_residue != (chain_id, resid):
                    current_residue = (chain_id, resid)
                    current_residue_should_keep = (chain_id, resid) in target_residues_to_keep
                
                if current_residue_should_keep:
                    atom_lines.append(line)
        
        # Nadaj nowe numery atomów
        atom_count = 1
        for i, line in enumerate(atom_lines):
            atom_num_str = str(atom_count).rjust(5)
            atom_lines[i] = line[:6] + atom_num_str + line[11:]
            atom_count += 1
        
        # Zapisz BEZ DODAWANIA ŻADNYCH NAGŁÓWKÓW
        with open(output_path, 'w') as f:
            f.writelines(atom_lines)
        
        return atom_count - 1
    
    @staticmethod
    def _extract_residue_info(line):
        """Wyodrębnia informacje o reszcie z linii PDB"""
        chain_id = line[21]
        try:
            resid = int(line[22:26].strip())
        except ValueError:
            resid = -1
        return chain_id, resid
    
    @staticmethod
    def _sort_atoms(atom_lines, chain_order):
        """Sortuje atomy według łańcuchów i numerów reszt"""
        # Tylko linie ATOM/HETATM - pomijamy wszystkie inne
        atom_dict = defaultdict(list)
        
        for line in atom_lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain, resid = PDBHandler._extract_residue_info(line)
                atom_dict[(chain, resid)].append(line)
        
        sorted_atoms = []
        
        # Atomy w kolejności łańcuchów targetu
        for chain_id in chain_order:
            residues = [(c, r) for (c, r) in atom_dict.keys() if c == chain_id]
            residues.sort(key=lambda x: x[1])
            for residue_key in residues:
                sorted_atoms.extend(atom_dict[residue_key])
        
        # Pozostałe łańcuchy
        remaining_chains = set([c for (c, r) in atom_dict.keys()]) - set(chain_order)
        for chain_id in sorted(remaining_chains):
            residues = [(c, r) for (c, r) in atom_dict.keys() if c == chain_id]
            residues.sort(key=lambda x: x[1])
            for residue_key in residues:
                sorted_atoms.extend(atom_dict[residue_key])
        
        return sorted_atoms

class AtomChecker:
    @staticmethod
    def check_missing_atoms(target_res, model_res):
        """Sprawdza brakujące atomy w modelu względem targetu"""
        missing = []
        base = target_res.resname.strip()
        
        # Sprawdź atomy szkieletu
        for atom in ESSENTIAL_ATOMS:
            try:
                target_res[atom]
                try: 
                    model_res[atom]
                except: 
                    missing.append(atom)
            except: 
                pass
        
        # Sprawdź atomy zasad specyficzne dla bazy
        base_atoms = {
            'A': ['N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4'],
            'G': ['N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4'],
            'C': ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'],
            'U': ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
        }
        
        if base in base_atoms:
            for atom in base_atoms[base]:
                try:
                    target_res[atom]
                    try: 
                        model_res[atom]
                    except: 
                        missing.append(atom)
                except: 
                    pass
        
        return missing