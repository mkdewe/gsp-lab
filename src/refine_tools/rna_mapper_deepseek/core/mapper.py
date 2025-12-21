"""Główna logika mapowania RNA"""

from exceptions import RNAMapperError
from utils.logger import Logger
from core.structure_parser import StructureParser
from core.sequence_aligner import SequenceAligner

class RNAMapper:
    def __init__(self):
        self.structure_parser = StructureParser()
        self.sequence_aligner = SequenceAligner()
    
    def process(self, target_path, model_path):
        Logger.section("RNA MAPPER")
        
        # Wczytaj struktury
        target_struct, target_seqs, target_res = self.structure_parser.parse(target_path, "target")
        model_struct, model_seqs, model_res = self.structure_parser.parse(model_path, "model")
        
        # Mapuj łańcuchy
        chain_mapping = self.sequence_aligner.find_chain_mapping(target_seqs, model_seqs)
        if not chain_mapping:
            raise RNAMapperError("Nie znaleziono pasujących sekwencji")
        
        # Mapuj reszty z alignowaniem
        residue_mapping = self._create_residue_mapping_with_alignment(chain_mapping, target_seqs, target_res, model_seqs, model_res)
        if not residue_mapping:
            raise RNAMapperError("Nie udało się zmapować reszt")
        
        Logger.info(f"Zmapowano {len(residue_mapping)} reszt")
        return {
            'target_structure': target_struct, 
            'model_structure': model_struct,
            'residue_mapping': residue_mapping, 
            'target_seqs': target_seqs
        }
    
    def _create_residue_mapping_with_alignment(self, chain_mapping, target_seqs, target_res, model_seqs, model_res):
        """Tworzy mapowanie reszt z użyciem alignowania sekwencji"""
        mapping = {}
        
        for target_chain, model_chain in chain_mapping.items():
            if target_chain not in target_res or model_chain not in model_res:
                continue
            
            target_seq = target_seqs[target_chain]
            model_seq = model_seqs[model_chain]
            
            Logger.info(f"Alignowanie łańcucha {target_chain} (target) i {model_chain} (model)")
            Logger.info(f"Target seq: {target_seq}")
            Logger.info(f"Model seq:  {model_seq}")
            
            # Wyrównaj sekwencje
            alignments = self.sequence_aligner.aligner.align(target_seq, model_seq)
            if not alignments:
                Logger.warning(f"Brak alignowania dla {target_chain} i {model_chain}")
                continue
                
            best_alignment = alignments[0]
            aligned_target, aligned_model = best_alignment[0], best_alignment[1]
            
            Logger.info(f"Aligned target: {aligned_target}")
            Logger.info(f"Aligned model:  {aligned_model}")
            
            # Mapowanie na podstawie alignowania
            target_idx, model_idx = 0, 0
            mapped_in_chain = 0
            
            for t_char, m_char in zip(aligned_target, aligned_model):
                if t_char != '-' and m_char != '-':
                    # Oba znaki są nukleotydami - mapuj reszty
                    if (target_idx < len(target_res[target_chain]) and 
                        model_idx < len(model_res[model_chain])):
                        
                        target_residue = target_res[target_chain][target_idx]
                        model_residue = model_res[model_chain][model_idx]
                        
                        # Sprawdź czy nukleotydy się zgadzają
                        if (target_residue.resname.strip() == t_char and 
                            model_residue.resname.strip() == m_char):
                            
                            mapping[(model_chain, model_residue.id[1])] = (target_chain, target_residue.id[1])
                            Logger.debug(f"  {model_chain}{model_residue.id[1]} -> {target_chain}{target_residue.id[1]}")
                            mapped_in_chain += 1
                        else:
                            Logger.warning(f"  Nukleotydy nie zgadzają się: {model_chain}{model_residue.id[1]}({model_residue.resname}) != {target_chain}{target_residue.id[1]}({target_residue.resname})")
                    
                    target_idx += 1
                    model_idx += 1
                    
                elif t_char == '-':
                    # Gap w targetcie - pomiń resztę w modelu
                    if model_idx < len(model_res[model_chain]):
                        Logger.debug(f"  Gap w targetcie, pomijam {model_chain}{model_res[model_chain][model_idx].id[1]}")
                    model_idx += 1
                    
                elif m_char == '-':
                    # Gap w modelu - pomiń resztę w targetcie
                    if target_idx < len(target_res[target_chain]):
                        Logger.debug(f"  Gap w modelu, pomijam {target_chain}{target_res[target_chain][target_idx].id[1]}")
                    target_idx += 1
            
            Logger.info(f"Zmapowano {mapped_in_chain} reszt dla {model_chain}->{target_chain}")
        
        return mapping