# semantic_rna_mapper.py
"""
SemanticRNAMapper
- Weighted central-atom semantic mapping (profile B)
- Atom-name normalization (C4* -> C4', etc.)
- Atom-level fallback with adaptive local threshold
- Does NOT change coordinates (no transforms)
- Writes:
    - refined PDB (model coords unchanged, semantics from target, optionally pruned atoms)
    - solution PDB (original target trimmed to atoms present in refined)
- Entrypoint: process_all(target_file, model_file, output_pdb, remove_extra_atoms=False)
"""
import os
import math
import numpy as np
from collections import defaultdict
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings

warnings.filterwarnings("ignore", category=PDBConstructionWarning)
warnings.filterwarnings("ignore")

try:
    from scipy.optimize import linear_sum_assignment
except Exception:
    linear_sum_assignment = None


class SemanticRNAMapper:
    def __init__(self, profile='B'):
        # profile B weights
        # P:5, C4':4, C1':4, N1/N9:4, others:1
        self.weights = self._weights_for_profile(profile)
        # canonical central atoms we care about (with normalized names)
        self.central_atoms = {'P', "C4'", "C1'", "N1", "N9"}
        # atom name alias mapping (normalize during extraction)
        # many PDBs use * instead of '
        self.atom_aliases = {
            "C4*": "C4'",
            "C1*": "C1'",
            "O2*": "O2'",
            "C2*": "C2'",
            "C3*": "C3'",
            "O3*": "O3'",
            "O4*": "O4'",
            "OP1": "OP1", "OP2": "OP2",
            "O1P": "OP1", "O2P": "OP2"
        }

    def _weights_for_profile(self, profile):
        if profile == 'A':
            return {'P': 10, "C4'": 8, "C1'": 8, 'N1': 8, 'N9': 8}
        if profile == 'C':
            return {'P': 3, "C4'": 2, "C1'": 2, 'N1': 2, 'N9': 2}
        # default B
        return {'P': 5, "C4'": 4, "C1'": 4, 'N1': 4, 'N9': 4}

    # ---------------- atom name normalization ----------------
    def normalize_atom_name(self, name):
        """Normalize common aliases: * -> ', O1P/O2P -> OP1/OP2, uppercase."""
        if not name:
            return name
        n = name.strip()
        # replace * with '
        n = n.replace('*', "'")
        # direct alias map (case-insensitive)
        key = n.upper()
        if key in self.atom_aliases:
            return self.atom_aliases[key]
        # also ensure standard formatting: C4' etc stay as-is
        return n

    # ---------------- load structure ----------------
    def load_structure(self, path, label=None):
        parser = PDBParser(QUIET=True)
        label = label or os.path.basename(path)
        return parser.get_structure(label, path)

    # ---------------- extract residue representative points & atom sets ----------------
    def extract_residue_points(self, structure):
        """
        Returns:
            keys: list of (chain_id, resid)
            info: dict key -> {
                'residue': residue_obj,
                'coord': np.array([x,y,z]) for chosen representative atom (pref central),
                'resname': 'A'/'C'/'G'/'U',
                'atom_names': set([...])  # normalized atom names present,
                'atom_coords': dict(name -> np.array(coord))
            }
        """
        info = {}
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    resname = residue.resname.strip()
                    if resname not in ('A', 'C', 'G', 'U'):
                        continue
                    atom_coords = {}
                    atom_names = set()
                    # collect atoms (normalize their names)
                    for atom in residue:
                        aname = self.normalize_atom_name(atom.get_name())
                        atom_names.add(aname)
                        atom_coords[aname] = np.array(atom.get_coord())
                    # choose representative atom preferentially: P, C4', C1', N1/N9, fallback to any heavy atom
                    rep = None
                    for prefer in ['P', "C4'", "C1'", "N1", "N9"]:
                        if prefer in atom_coords:
                            rep = atom_coords[prefer]
                            break
                    if rep is None:
                        # fallback any non-hydrogen atom in atom_coords
                        if atom_coords:
                            # pick first
                            rep = next(iter(atom_coords.values()))
                        else:
                            continue
                    key = (chain.id, residue.id[1])
                    info[key] = {
                        'residue': residue,
                        'coord': rep,
                        'resname': resname,
                        'atom_names': atom_names,
                        'atom_coords': atom_coords
                    }
        keys = list(info.keys())
        return keys, info

    # ---------------- weighted residue distance ----------------
    def weighted_residue_distance(self, t_entry, m_entry):
        """
        Compute weighted distance between two residues using matching atom names.
        Use central-atom weights from self.weights; missing atoms penalized by a large distance.
        """
        # consider union of atom names present in either residue (but focus on central atoms)
        t_atoms = t_entry['atom_coords']
        m_atoms = m_entry['atom_coords']
        # union
        atom_names = set(t_atoms.keys()).union(set(m_atoms.keys()))
        total_w = 0.0
        wdist = 0.0
        # if both have an atom, include weighted distance; if missing on one side, add penalty
        for an in atom_names:
            # normalized name used as-is
            w = self.weights.get(an, 1.0) if an in self.weights else (self.weights.get(an, 1.0) if an in self.weights else (5.0 if an in self.central_atoms else 1.0))
            # if central atom by name, ensure weight
            if an in self.central_atoms and an not in self.weights:
                w = 4.0
            total_w += w
            if an in t_atoms and an in m_atoms:
                d = np.linalg.norm(t_atoms[an] - m_atoms[an])
                wdist += w * d
            else:
                # penalty for missing atom - use a moderate penalty (e.g., 6 Å)
                # but scale penalty with weight
                penalty = 6.0
                wdist += w * penalty
        if total_w == 0:
            return float('inf')
        return wdist / total_w

    # ---------------- build semantic mapping (weighted Hungarian) ----------------
    def build_semantic_residue_mapping(self, target_structure, model_structure, max_cost=None):
        """
        Primary mapping: build cost matrix using weighted_residue_distance,
        then solve Hungarian assignment (min sum distances).
        Returns mapping dict: model_key -> target_key and stats.
        """
        if linear_sum_assignment is None:
            raise ImportError("scipy.optimize.linear_sum_assignment is required")

        tgt_keys, tgt_info = self.extract_residue_points(target_structure)
        mdl_keys, mdl_info = self.extract_residue_points(model_structure)

        if not tgt_keys or not mdl_keys:
            return {}, {}, {'mapped_pairs': 0, 'n_target': len(tgt_keys), 'n_model': len(mdl_keys)}

        # compute cost matrix (n_t x n_m)
        nt = len(tgt_keys)
        nm = len(mdl_keys)
        cost = np.zeros((nt, nm), dtype=float)
        for i, tk in enumerate(tgt_keys):
            for j, mk in enumerate(mdl_keys):
                cost[i, j] = self.weighted_residue_distance(tgt_info[tk], mdl_info[mk])

        # Optionally cap costs if max_cost provided
        if max_cost is not None:
            cost = np.minimum(cost, max_cost)

        row_ind, col_ind = linear_sum_assignment(cost)

        mapping = {}
        # build mapping but optionally skip absurdly large distances (we'll allow and later fallback)
        for r, c in zip(row_ind, col_ind):
            mapping[mdl_keys[c]] = tgt_keys[r]

        stats = {'mapped_pairs': len(mapping), 'n_target': nt, 'n_model': nm, 'cost_matrix_mean': float(np.mean(cost))}
        return mapping, (tgt_info, mdl_info), stats

    # ---------------- adaptive local threshold ----------------
    def adaptive_threshold_for_model_residue(self, m_entry, all_target_coords, k=5, scale=2.0):
        """
        Compute adaptive distance threshold for a given model residue based on local density:
        - compute distances from m_entry to k nearest target representative points
        - threshold = median_of_k * scale (fallback to fixed value if not enough)
        """
        if len(all_target_coords) == 0:
            return 6.0
        dists = np.linalg.norm(all_target_coords - m_entry['coord'], axis=1)
        if len(dists) >= k:
            nn = np.partition(dists, k - 1)[:k]
            med = float(np.median(nn))
            return max(2.0, med * scale)
        else:
            med = float(np.median(dists)) if len(dists) > 0 else 3.0
            return max(2.0, med * scale)

    # ---------------- atom-level fallback ----------------
    def atom_level_fallback(self, residue_mapping, tgt_info, mdl_info):
        """
        For unmapped model residues, attempt to assign by atom-level agreement:
        - For each unmapped model residue, compute candidate target residues within adaptive threshold.
        - Score candidate by weighted count of atom-name matches where distance < atom_distance_cutoff (1.8 A).
        - If best score exceeds threshold (e.g., >= 0.4 of possible weighted sum), accept mapping.
        Returns augmented mapping and list of added mappings.
        """
        # build reverse mapping target->model for fast check
        mapped_targets = set(residue_mapping.values())
        added = {}
        tgt_keys = list(tgt_info.keys())
        mdl_keys = list(mdl_info.keys())
        all_tgt_coords = np.vstack([tgt_info[k]['coord'] for k in tgt_keys])

        # build quick lookup for candidate scoring
        for mk in mdl_keys:
            if mk in residue_mapping:
                continue  # already mapped
            m_entry = mdl_info[mk]
            # adaptive threshold
            thr = self.adaptive_threshold_for_model_residue(m_entry, all_tgt_coords)
            # candidates: target residues whose rep-point distance <= thr
            dists = np.linalg.norm(all_tgt_coords - m_entry['coord'], axis=1)
            candidate_idxs = np.where(dists <= thr)[0]
            if candidate_idxs.size == 0:
                # try nearest few anyway
                candidate_idxs = np.argsort(dists)[:5]
            best_score = 0.0
            best_tk = None
            for idx in candidate_idxs:
                tk = tgt_keys[int(idx)]
                if tk in mapped_targets:
                    # allow mapping to already mapped target? Generally prefer unique mapping,
                    # so skip already mapped targets to avoid collisions
                    continue
                t_entry = tgt_info[tk]
                # score by weighted matching atoms within atomic cutoff
                score, max_possible = self._atom_match_score(t_entry, m_entry, atom_distance_cutoff=1.8)
                # normalized score
                if max_possible > 0:
                    norm_score = score / max_possible
                else:
                    norm_score = 0.0
                if norm_score > best_score:
                    best_score = norm_score
                    best_tk = tk
            # acceptance threshold: e.g., >= 0.45 (tunable)
            if best_score >= 0.45 and best_tk is not None:
                residue_mapping[mk] = best_tk
                mapped_targets.add(best_tk)
                added[mk] = best_tk
        return residue_mapping, added

    def _atom_match_score(self, t_entry, m_entry, atom_distance_cutoff=1.8):
        """
        Compute weighted match score between t_entry and m_entry by counting matching atom names
        whose euclidean distance < atom_distance_cutoff. Weight by central atom weights.
        Returns (score, max_possible_score)
        """
        t_atoms = t_entry['atom_coords']
        m_atoms = m_entry['atom_coords']
        # union of atom names that appear in either
        atom_names = set(t_atoms.keys()).union(set(m_atoms.keys()))
        score = 0.0
        max_possible = 0.0
        for an in atom_names:
            w = self.weights.get(an, 1.0) if an in self.weights else (4.0 if an in self.central_atoms else 1.0)
            max_possible += w
            if an in t_atoms and an in m_atoms:
                d = np.linalg.norm(t_atoms[an] - m_atoms[an])
                if d <= atom_distance_cutoff:
                    score += w
        return score, max_possible

    # ---------------- write refined PDB (model coords unchanged) ----------------
    def write_refined_pdb(self, model_structure, residue_mapping, target_info, output_path):
        """
        Write refined model:
        - Coordinates from model_structure
        - labels (chain, resnum, resname) from mapped target_info
        - Only residues present in residue_mapping are written
        """
        # index model residues
        model_index = {}
        for model in model_structure:
            for ch in model:
                for r in ch:
                    if r.id[0] != ' ':
                        continue
                    model_index[(ch.id, r.id[1])] = r

        # invert mapping: target_key -> [model_keys]
        inv = defaultdict(list)
        for m_key, t_key in residue_mapping.items():
            inv[t_key].append(m_key)

        lines = []
        atom_serial = 1
        lines.append("REMARK   Refined model created by SemanticRNAMapper\n")
        for t_chain in sorted({k[0] for k in inv.keys()}):
            # collect target residue keys for this chain sorted by resnum
            t_reskeys = sorted([k for k in inv.keys() if k[0] == t_chain], key=lambda x: x[1])
            for t_key in t_reskeys:
                t_resnum = t_key[1]
                t_resname = target_info[t_key]['resname'] if t_key in target_info else 'UNK'
                # there might be multiple model residues mapped to same target (rare) -> write them sequentially
                for m_key in inv[t_key]:
                    if m_key not in model_index:
                        continue
                    res = model_index[m_key]
                    for atom in res:
                        atom_name = self.normalize_atom_name(atom.get_name())
                        try:
                            bfactor = atom.get_bfactor()
                        except Exception:
                            bfactor = 0.00
                        try:
                            occ = atom.get_occupancy()
                        except Exception:
                            occ = 1.00
                        x, y, z = atom.get_coord().tolist()
                        an_fmt = atom_name.rjust(4) if len(atom_name) < 4 else atom_name[:4]
                        line = ("ATOM  {serial:5d} {name:4s} {resname:>3s} {chain:1s}{resnum:4d}"
                                "    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}\n").format(
                            serial=atom_serial, name=an_fmt, resname=t_resname,
                            chain=t_chain, resnum=t_resnum, x=x, y=y, z=z, occ=occ, bf=bfactor
                        )
                        lines.append(line)
                        atom_serial += 1
        lines.append("END\n")
        with open(output_path, 'w') as fh:
            fh.writelines(lines)
        return atom_serial - 1

    # ---------------- write solution by trimming original target to atoms present in refined model ----------------
    def write_solution_from_refined(self, target_structure, model_structure, residue_mapping, output_path):
        """
        Build allowed set: for every mapped model residue m_key -> target t_key,
        include triple (t_chain, t_resnum, atom_name) for atoms that are present in the
        refined model residue (i.e., atoms kept after pruning).
        Then iterate original target and write only atoms that are in allowed set.
        """
        # build model index
        model_index = {}
        for model in model_structure:
            for ch in model:
                for r in ch:
                    if r.id[0] != ' ':
                        continue
                    model_index[(ch.id, r.id[1])] = r

        allowed = set()
        for m_key, t_key in residue_mapping.items():
            if m_key not in model_index:
                continue
            res = model_index[m_key]
            t_chain, t_resnum = t_key
            for atom in res:
                aname = self.normalize_atom_name(atom.get_name())
                allowed.add((t_chain, t_resnum, aname))

        # iterate original target and write only allowed atoms
        lines = []
        atom_serial = 1
        lines.append("REMARK   Solution built by trimming target to refined-model atoms\n")
        for model in target_structure:
            for ch in model:
                for r in ch:
                    if r.id[0] != ' ':
                        continue
                    t_chain = ch.id
                    t_resnum = r.id[1]
                    for atom in r:
                        aname = self.normalize_atom_name(atom.get_name())
                        if (t_chain, t_resnum, aname) not in allowed:
                            continue
                        try:
                            bfactor = atom.get_bfactor()
                        except Exception:
                            bfactor = 0.00
                        try:
                            occ = atom.get_occupancy()
                        except Exception:
                            occ = 1.00
                        x, y, z = atom.get_coord().tolist()
                        an_fmt = aname.rjust(4) if len(aname) < 4 else aname[:4]
                        line = ("ATOM  {serial:5d} {name:4s} {resname:>3s} {chain:1s}{resnum:4d}"
                                "    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}\n").format(
                            serial=atom_serial, name=an_fmt, resname=r.resname,
                            chain=t_chain, resnum=t_resnum, x=x, y=y, z=z, occ=occ, bf=bfactor
                        )
                        lines.append(line)
                        atom_serial += 1
        lines.append("END\n")
        with open(output_path, 'w') as fh:
            fh.writelines(lines)
        return atom_serial - 1

    # ---------------- compute mean representative distance ----------------
    def compute_mean_rep_distance(self, tgt_info, mdl_info, residue_mapping):
        dists = []
        for m_key, t_key in residue_mapping.items():
            m_e = mdl_info.get(m_key)
            t_e = tgt_info.get(t_key)
            if m_e is None or t_e is None:
                continue
            d = np.linalg.norm(m_e['coord'] - t_e['coord'])
            dists.append(d)
        if not dists:
            return float('inf'), []
        return float(np.mean(dists)), dists

    # ---------------- main pipeline ----------------
    def process_all(self, target_file, model_file, output_pdb, remove_extra_atoms=False):
        """
        Main entrypoint - compatible signature.
        Returns similar tuple as earlier implementations.
        """
        print("=== START: SemanticRNAMapper (profile B) ===")
        print(f" target: {target_file}")
        print(f" model:  {model_file}")
        print(f" output: {output_pdb}")
        print(f" remove_extra_atoms: {remove_extra_atoms}")

        try:
            if linear_sum_assignment is None:
                raise ImportError("scipy.optimize.linear_sum_assignment is required")

            target_structure = self.load_structure(target_file, 'target')
            model_structure = self.load_structure(model_file, 'model')

            # 1) build primary mapping by weighted residue matching
            residue_mapping, (tgt_info, mdl_info), stats = self.build_semantic_residue_mapping(
                target_structure, model_structure
            )
            print(f"Primary mapped pairs: {stats.get('mapped_pairs')} (target: {stats.get('n_target')}, model: {stats.get('n_model')})")

            # 2) atom-level fallback to try to recover unmapped residues
            before_unmapped = len(mdl_info) - len(residue_mapping)
            residue_mapping, added = self.atom_level_fallback(residue_mapping, tgt_info, mdl_info)
            after_unmapped = len(mdl_info) - len(residue_mapping)
            print(f"Atom-level fallback added: {len(added)} mappings (unmapped before: {before_unmapped}, after: {after_unmapped})")

            # 3) Optionally remove extra atoms in model residues not present in mapped target residue
            # apply_semantic_mapping previously mutated model; here we actually remove atoms in model_structure
            missing_atoms_report = []
            removed_atoms_count = 0
            changed_count = 0
            # index model residues
            model_index = {}
            for model in model_structure:
                for ch in model:
                    for r in list(ch):
                        if r.id[0] != ' ':
                            continue
                        m_key = (ch.id, r.id[1])
                        model_index[m_key] = r

            for m_key, t_key in list(residue_mapping.items()):
                if m_key not in model_index:
                    continue
                m_res = model_index[m_key]
                t_entry = tgt_info.get(t_key)
                if t_entry is None:
                    continue
                t_atom_names = t_entry['atom_names']
                model_atom_names = set(self.normalize_atom_name(a.get_name()) for a in m_res)
                missing = sorted(list(t_atom_names - model_atom_names))
                if missing:
                    missing_atoms_report.append({
                        'model_chain': m_key[0], 'model_residue': m_key[1],
                        'target_chain': t_key[0], 'target_residue': t_key[1],
                        'missing_atoms': missing
                    })
                if remove_extra_atoms:
                    # remove atoms present in model but not in target atom list
                    for atom in list(m_res):
                        an = self.normalize_atom_name(atom.get_name())
                        if an not in t_atom_names:
                            try:
                                m_res.detach_child(atom.id)
                                removed_atoms_count += 1
                            except Exception:
                                try:
                                    m_res.child_list.remove(atom)
                                    removed_atoms_count += 1
                                except Exception:
                                    pass
                changed_count += 1

            # 4) write refined (model coords unchanged, semantics from target)
            refined_output = output_pdb
            refined_written = self.write_refined_pdb(model_structure, residue_mapping, tgt_info, refined_output)
            print(f"Wrote refined: {refined_output} atoms={refined_written}")

            # 5) write solution by trimming original target to atoms present in refined model
            model_basename = os.path.basename(model_file)
            model_name = os.path.splitext(model_basename)[0]
            solution_filename = f"solution_{model_name}.pdb"
            solution_output = os.path.join(os.path.dirname(output_pdb) or '.', solution_filename)
            sol_written = self.write_solution_from_refined(target_structure, model_structure, residue_mapping, solution_output)
            print(f"Wrote solution: {solution_output} atoms={sol_written}")

            # 6) metrics & return tuple
            mean_dist, dlist = self.compute_mean_rep_distance(tgt_info, mdl_info, residue_mapping)
            atom_comparison_log = [f"- MISSING in model {r['target_chain']}{r['target_residue']}: {r['missing_atoms']} mapped from {r['model_chain']}{r['model_residue']}" for r in missing_atoms_report]
            total_extra_atoms = removed_atoms_count
            total_missing_atoms = len(missing_atoms_report)
            removed_residues_count = len(mdl_info) - len(residue_mapping)

            print("=== DONE: process_all ===")
            return residue_mapping, missing_atoms_report, mean_dist, changed_count, removed_atoms_count, removed_residues_count, atom_comparison_log, total_extra_atoms, total_missing_atoms, solution_output

        except Exception as e:
            print("=== ERROR in process_all ===")
            import traceback
            traceback.print_exc()
            return {}, [], float('inf'), 0, 0, 0, [], 0, 0, None


# If executed as script
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="SemanticRNAMapper - semantic mapping without transforming coords")
    parser.add_argument("target_pdb")
    parser.add_argument("model_pdb")
    parser.add_argument("--output", "-o", required=True, help="Refined output PDB path")
    parser.add_argument("--remove-extra-atoms", action='store_true', help="Remove model atoms not present in mapped target residues")
    args = parser.parse_args()

    mapper = SemanticRNAMapper(profile='B')
    res = mapper.process_all(args.target_pdb, args.model_pdb, args.output, remove_extra_atoms=args.remove_extra_atoms)
    if res and res[9]:
        print("Solution file:", res[9])
