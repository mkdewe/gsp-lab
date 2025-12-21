"""Konfiguracja"""

from pathlib import Path

DEFAULT_OUTPUT_DIR = Path("output")
DEFAULT_LOGGING = False

ALIGNER_PARAMS = {
    'mode': 'global', 
    'match_score': 2, 
    'mismatch_score': -1,
    'open_gap_score': -0.5, 
    'extend_gap_score': -0.1
}

SEQUENCE_THRESHOLD = 0.3
ESSENTIAL_ATOMS = ['P', 'C4\'', 'C1\'']
SOLUTION_PREFIX = "solution"
REFINED_PREFIX = "refined"