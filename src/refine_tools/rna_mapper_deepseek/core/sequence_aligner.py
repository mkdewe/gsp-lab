"""Dopasowanie sekwencji RNA"""

from Bio.Align import PairwiseAligner
from config import ALIGNER_PARAMS, SEQUENCE_THRESHOLD
from utils.logger import Logger

class SequenceAligner:
    def __init__(self):
        self.aligner = PairwiseAligner()
        for param, value in ALIGNER_PARAMS.items():
            setattr(self.aligner, param, value)
    
    def find_chain_mapping(self, target_seqs, model_seqs):
        mapping, used = {}, set()
        for target_chain, target_seq in target_seqs.items():
            best_score, best_model = 0, None
            for model_chain, model_seq in model_seqs.items():
                if model_chain in used: 
                    continue
                score = self.calculate_similarity(target_seq, model_seq)
                Logger.info(f"{target_chain} vs {model_chain}: {score:.3f}")
                if score > best_score and score > SEQUENCE_THRESHOLD:
                    best_score, best_model = score, model_chain
            if best_model:
                mapping[target_chain] = best_model
                used.add(best_model)
                Logger.info(f"Mapowanie: {target_chain} -> {best_model}")
        return mapping
    
    def calculate_similarity(self, seq1, seq2):
        if seq1 == seq2: 
            return 1.0
        try:
            alignments = self.aligner.align(seq1, seq2)
            if alignments:
                aligned = alignments[0]
                matches = sum(1 for a, b in zip(aligned[0], aligned[1]) if a == b)
                return matches / len(aligned[0])
        except: 
            pass
        return 0.0