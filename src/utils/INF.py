# INF.py
# Utility wrapper to compute INF_all / INF_wc / INF_nwc / INF_stack and RMSD
# using the legacy RNA_assessment (RNA_normalizer) code.
#
# Place this file in: D:\Studia\Projekty\gsp-lab\src\utils\INF.py
# and call from your project:
#   from utils.INF import compute_INF
#   res = compute_INF(native_pdb, model_pdb)
#
# CLI:
#   python INF.py native.pdb model.pdb [--csv out.csv]
#
import os
import sys
import tempfile
import json
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings('ignore', category=PDBConstructionWarning)

# ------------------ EDIT THIS IF RNA_assessment IS ELSEWHERE ------------------
RNA_ASSESSMENT_ROOT = os.environ.get("RNA_ASSESSMENT_ROOT")
if not RNA_ASSESSMENT_ROOT:
    raise RuntimeError(
        "RNA_ASSESSMENT_ROOT is not set. "
        "In Docker it must point to /rna_assessment"
    )

# ------------------------------------------------------------------------------

# ensure RNA_assessment is importable
if RNA_ASSESSMENT_ROOT not in sys.path:
    sys.path.insert(0, RNA_ASSESSMENT_ROOT)

try:
    import RNA_normalizer
except Exception as e:
    raise ImportError(
        "Cannot import RNA_normalizer from RNA_assessment. "
        f"Check RNA_ASSESSMENT_ROOT (currently {RNA_ASSESSMENT_ROOT}) and that "
        "RNA_normalizer.py exists there."
    ) from e

# data files used by PDBNormalizer in RNA_assessment
DATA_DIR = os.path.join(RNA_ASSESSMENT_ROOT, "data")
RESIDUES_LIST = os.path.join(DATA_DIR, "residues.list")
ATOMS_LIST = os.path.join(DATA_DIR, "atoms.list")
MCANNOTATE_BIN = "MC-Annotate"



def _mktempdir(prefix="rna_inf_"):
    return tempfile.mkdtemp(prefix=prefix)

def _ensure_normalized_and_index(pdb_path):
    """
    Ensure structure is normalized (via PDBNormalizer.parse) and return
    tuple (pdb_to_load, index_path_or_None, tmpdir_used_or_None).
    If an existing .index file next to pdb_path is found, it will be used.
    The function may create temporary files (normalized pdb) and returns the
    tmpdir path so caller can remove it later.
    """
    pdb_path = os.path.abspath(pdb_path)
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    idx_candidate = pdb_path + ".index"
    if os.path.exists(idx_candidate):
        return pdb_path, idx_candidate, None

    # Create normalized file in temp dir
    tmpdir = _mktempdir()
    normalized_name = os.path.join(tmpdir, os.path.basename(pdb_path).replace('.pdb', '_normalized.pdb'))

    normalizer = RNA_normalizer.PDBNormalizer(RESIDUES_LIST, ATOMS_LIST)
    ok = normalizer.parse(pdb_path, normalized_name)
    if not ok:
        # still attempt to return original file (caller will see load errors)
        return pdb_path, None, None

    # check for index files produced/available
    idx_norm = normalized_name + ".index"
    if os.path.exists(idx_norm):
        return normalized_name, idx_norm, tmpdir

    # maybe index exists next to original file
    if os.path.exists(idx_candidate):
        return normalized_name, idx_candidate, tmpdir

    # no index found; return normalized pdb and None index (PDBStruct.load may accept None)
    return normalized_name, None, tmpdir

def compute_INF(native_pdb, model_pdb, write_csv=None):
    """
    Compute INF and RMSD using RNA_assessment's RNA_normalizer.
    Returns dict with keys: rmsd, inf_all, inf_wc, inf_nwc, inf_stack.
    If write_csv path is provided, saves results there.
    """
    tmpdirs = []
    try:
        native_p, native_idx, tmpn = _ensure_normalized_and_index(native_pdb)
        model_p, model_idx, tmpm = _ensure_normalized_and_index(model_pdb)
        if tmpn:
            tmpdirs.append(tmpn)
        if tmpm:
            tmpdirs.append(tmpm)

        # load structures
        res_struct = RNA_normalizer.PDBStruct()
        res_struct.load(native_p, native_idx)
        sol_struct = RNA_normalizer.PDBStruct()
        sol_struct.load(model_p, model_idx)

        comparer = RNA_normalizer.PDBComparer()

        try:
            rmsd = comparer.rmsd(sol_struct, res_struct)
        except Exception:
            rmsd = None

        def safe_inf(t):
            try:
                return comparer.INF(sol_struct, res_struct, type=t)
            except Exception:
                return None

        INF_ALL = safe_inf("ALL")
        INF_WC = safe_inf("PAIR_2D")
        INF_NWC = safe_inf("PAIR_3D")
        INF_STACK = safe_inf("STACK")

        result = {
            "native_pdb": os.path.abspath(native_p),
            "native_index": os.path.abspath(native_idx) if native_idx else None,
            "model_pdb": os.path.abspath(model_p),
            "model_index": os.path.abspath(model_idx) if model_idx else None,
            "rmsd": rmsd,
            "inf_all": INF_ALL,
            "inf_wc": INF_WC,
            "inf_nwc": INF_NWC,
            "inf_stack": INF_STACK,
        }

        if write_csv:
            # simple CSV with header
            import csv
            header = ["native_pdb","model_pdb","rmsd","inf_all","inf_wc","inf_nwc","inf_stack"]
            with open(write_csv, "w", newline='', encoding='utf-8') as fh:
                writer = csv.writer(fh)
                writer.writerow(header)
                writer.writerow([result.get(h) for h in header])

        return result

    finally:
        # cleanup temp dirs produced during normalization
        for d in tmpdirs:
            try:
                import shutil
                shutil.rmtree(d)
            except Exception:
                pass

# ---------------------- CLI support ----------------------
if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Compute INF/RMSD using RNA_assessment (RNA_normalizer).")
    p.add_argument("native", help="native/reference pdb file")
    p.add_argument("model", help="model pdb file")
    p.add_argument("--csv", help="optional output CSV file", default=None)
    args = p.parse_args()

    try:
        out = compute_INF(args.native, args.model, write_csv=args.csv)
    except Exception as e:
        print("ERROR computing INF:", e, file=sys.stderr)
        raise

    # pretty print JSON for easy parsing
    print(json.dumps(out, indent=2, ensure_ascii=False))
