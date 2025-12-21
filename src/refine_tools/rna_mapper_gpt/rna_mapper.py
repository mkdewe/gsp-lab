from Bio.Seq import Seq
from Bio import pairwise2
import MDAnalysis as mda
from MDAnalysis.core.universe import Merge
import numpy as np
from prody import calcTransformation

# -------------------------------
# Helper: Needleman-Wunsch for mapping sequences
# -------------------------------

def best_seq_align(seq1, seq2):
    """Global alignment, returns aligned strings"""
    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -2, -0.5)
    best = alignments[0]
    return best.seqA, best.seqB


# -------------------------------
# 1. Load example PDBs
# -------------------------------

u_target = mda.Universe("target.pdb")
u_model = mda.Universe("model.pdb")

# Extract target chains A,B,C,D
target_chains = {seg.segid: seg.atoms for seg in u_target.segments}

print("Target chains:", target_chains.keys())

# Sequence of target chains (hardcoded example):
seqs_target = {
    'A': "AAAA",
    'B': "CCCC",
    'C': "UUUU",
    'D': "GGGG",
}

# The model has one chain X: AAAAUUUUGGGGCCCC
seq_model = "AAAAUUUUGGGGCCCC"


# -------------------------------
# 2. Split model into 4 chains by best alignment
# -------------------------------

ordered_target = ['A', 'B', 'C', 'D']
ordered_tseq   = ["AAAA", "CCCC", "UUUU", "GGGG"]

# Ideal segmentation lengths:
ideal_lengths = [len(s) for s in ordered_tseq]

model_segments = []
offset = 0
for li in ideal_lengths:
    model_segments.append(seq_model[offset:offset+li])
    offset += li

print("Model segments:", model_segments)

# -------------------------------
# 3. Build new model with 4 chains
# -------------------------------

atoms = u_model.atoms

# We assume residues in model are in correct order corresponding to seq_model
residues = u_model.residues
if len(residues) != len(seq_model):
    print("WARNING: model residues != seq length")

# Assign new chain IDs
new_chain_ids = []
idx = 0
for chain_id, seg_seq in zip(ordered_target, model_segments):
    for _ in seg_seq:
        new_chain_ids.append(chain_id)
        idx += 1

# Assign chain IDs to residues
for res, cid in zip(residues, new_chain_ids):
    res.segid = cid

# New 4-chain universe
model_split = u_model

# -------------------------------
# 4. Align each chain separately
# -------------------------------

def superpose_atoms(target_atoms, model_atoms):
    """RMSD alignment ignoring missing atoms."""
    # Match atoms by name
    T = []
    M = []
    t_dict = {a.name: a for a in target_atoms}
    for a in model_atoms:
        if a.name in t_dict:
            T.append(t_dict[a.name].position)
            M.append(a.position)
    T = np.array(T)
    M = np.array(M)
    if len(T) < 3:
        print("Too few atoms for alignment")
        return model_atoms

    # ProDy transformation
    transf = calcTransformation(M, T)
    new_positions = transf.apply(M)
    idx = 0
    for a in model_atoms:
        if a.name in t_dict:
            a.position = new_positions[idx]
            idx += 1


for cid in ordered_target:
    t = target_chains[cid]
    m = model_split.select_atoms(f"segid {cid}")
    if len(m) == 0:
        continue
    print(f"Aligning chain {cid}")
    superpose_atoms(t, m)

model_split.atoms.write("model_fixed.pdb")

print("DONE → Saved model_fixed.pdb")
