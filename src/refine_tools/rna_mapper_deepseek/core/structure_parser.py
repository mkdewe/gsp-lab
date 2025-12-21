"""Parsowanie struktur PDB"""

from Bio.PDB import PDBParser
import warnings
from exceptions import StructureParseError
from utils.logger import Logger

class StructureParser:
    def __init__(self):
        warnings.filterwarnings("ignore")
        self.parser = PDBParser()
    
    def parse(self, file_path, name):
        try:
            structure = self.parser.get_structure(name, file_path)
            sequences, residues = self.extract_rna_data(structure)
            rna_count = self.count_rna_residues(structure)
            Logger.info(f"{name}: {len(sequences)} łańcuchów, {rna_count} reszt RNA")
            return structure, sequences, residues
        except Exception as e:
            raise StructureParseError(f"Błąd parsowania {file_path}: {e}")
    
    def extract_rna_data(self, structure):
        sequences, residues = {}, {}
        for model in structure:
            for chain in model:
                seq, res_list = "", []
                for residue in chain:
                    if residue.id[0] == ' ':
                        resname = residue.resname.strip()
                        if resname in ['A','C','G','U']:
                            seq += resname
                            res_list.append(residue)
                if seq: 
                    sequences[chain.id], residues[chain.id] = seq, res_list
        return sequences, residues
    
    def get_chain_order(self, structure):
        return [chain.id for model in structure for chain in model]
    
    def count_rna_residues(self, structure):
        """Liczy całkowitą liczbę reszt RNA w strukturze"""
        count = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':
                        resname = residue.resname.strip()
                        if resname in ['A','C','G','U']:
                            count += 1
        return count