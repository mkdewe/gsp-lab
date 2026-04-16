#!/usr/bin/env python3
"""
RNA Puzzle - Interactive PDB Editor
Compares and edits solution/model files, records operations, and applies them to other models.
Output files contain REMARK lines describing the operations performed.
"""

import os
import sys
import pathlib
import re
import datetime
from collections import defaultdict, OrderedDict
import numpy as np

# -------------------------------------------------------------------
# PATH CONFIGURATION
# -------------------------------------------------------------------
BASE_DIR = pathlib.Path(__file__).parent.parent.parent.resolve()
BASE_PROCESSED_DIR = BASE_DIR / "data" / "processed"
BASE_EXTERNAL_DIR = BASE_DIR / "data" / "external"

# -------------------------------------------------------------------
# DICTIONARIES FOR CORRECT NUCLEOTIDE ATOM ORDER
# -------------------------------------------------------------------
NUCLEOTIDE_ATOM_ORDER = {
    'G': ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4"],
    'C': ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"],
    'A': ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"],
    'U': ["P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "N1", "C2", "O2", "N3", "C4", "O4", "C5", "C6"]
}

RESIDUE_NAME_MAP = {
    'G': 'G', 'GUA': 'G', 'GTP': 'G',
    'C': 'C', 'CYT': 'C', 'CTP': 'C',
    'A': 'A', 'ADE': 'A', 'ATP': 'A',
    'U': 'U', 'URA': 'U', 'UTP': 'U'
}

ATOMS_TO_REMOVE = { "H5'", "H5''", "H2'", "H3'", "H4'", "HO2'", "H1'", "H8", "H6", "H2", "H61", "H62", "H41", "H42" }

# -------------------------------------------------------------------
# ATOM NAME NORMALIZATION
# -------------------------------------------------------------------
def normalize_atom_name(name):
    """
    Normalize atom name:
        - Replace asterisk (*) with prime (')
        - Map O1P -> OP1, O2P -> OP2
        - Strip whitespace
    """
    name = name.strip()
    # Fix prime symbol
    name = name.replace('*', "'")
    # Fix phosphate naming
    if name in ('O1P', 'O2P'):
        name = 'OP' + name[-1]   # becomes OP1 or OP2
    return name

# -------------------------------------------------------------------
# PDB READ/WRITE FUNCTIONS
# -------------------------------------------------------------------
def read_pdb_lines(pdb_file):
    with open(pdb_file, 'r') as f:
        return f.readlines()

def write_pdb_lines(pdb_file, lines):
    with open(pdb_file, 'w') as f:
        f.writelines(lines)

def get_atom_lines(lines):
    """Extract ATOM/HETATM lines and parse them into a list of (line_num, line, info)."""
    atom_lines = []
    for idx, line in enumerate(lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                atom_num = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                atom_name = normalize_atom_name(atom_name)   # <-- normalisation
                alt_loc = line[16:17].strip()
                res_name = line[17:20].strip()
                chain_id = line[21:22].strip()
                res_seq = int(line[22:26].strip())
                icode = line[26:27].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occupancy = line[54:60].strip()
                temp_factor = line[60:66].strip()
                element = line[76:78].strip()

                atom_info = {
                    'line': line.rstrip('\n'),
                    'atom_num': atom_num,
                    'atom_name': atom_name,
                    'alt_loc': alt_loc,
                    'res_name': res_name,
                    'chain_id': chain_id,
                    'res_seq': res_seq,
                    'icode': icode,
                    'x': x, 'y': y, 'z': z,
                    'occupancy': occupancy,
                    'temp_factor': temp_factor,
                    'element': element,
                    'line_num': idx,
                    'original_line': line
                }
                atom_lines.append((idx, line, atom_info))
            except Exception as e:
                print(f"Error parsing line {idx}: {line.strip()}")
    return atom_lines

def update_full_pdb(original_lines, atom_lines):
    """Replace ATOM/HETATM section with sorted and formatted atom lines."""
    first_atom_idx = None
    last_atom_idx = None
    for i, line in enumerate(original_lines):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if first_atom_idx is None:
                first_atom_idx = i
            last_atom_idx = i
    if first_atom_idx is None:
        return original_lines
    before_atoms = original_lines[:first_atom_idx]
    after_atoms = original_lines[last_atom_idx + 1:]
    atom_lines_sorted = sorted(atom_lines, key=lambda x: x[2]['atom_num'])
    new_atom_lines = []
    for _, line, info in atom_lines_sorted:
        line = line.rstrip('\n')
        if len(line) < 80:
            line = line.ljust(80)
        elif len(line) > 80:
            line = line[:80]
        new_atom_lines.append(line + '\n')
    return before_atoms + new_atom_lines + after_atoms

def format_pdb_atom_line(serial,
                         atom_name,
                         alt_loc,
                         res_name,
                         chain_id,
                         res_seq,
                         i_code,
                         x, y, z,
                         occupancy='1.00',
                         temp_factor='0.00',
                         element=' '):
    """Create a formatted ATOM line."""
    try:
        serial_i = int(serial)
    except:
        serial_i = 0
    try:
        res_seq_i = int(res_seq)
    except:
        res_seq_i = 0
    try:
        x_f = float(x)
        y_f = float(y)
        z_f = float(z)
    except:
        x_f = y_f = z_f = 0.0
    try:
        occ_f = float(occupancy)
    except:
        occ_f = 1.0
    try:
        temp_f = float(temp_factor)
    except:
        temp_f = 0.0

    atom_name = (atom_name or "").strip()
    atom_name = normalize_atom_name(atom_name)   # <-- normalisation
    atom_name_field = atom_name.ljust(4)[:4]
    alt_loc_field = (alt_loc or ' ')[:1]
    res_name_field = str(res_name or '').strip().rjust(3)[:3]
    chain_id_field = (chain_id or ' ')[:1]
    i_code_field = (i_code or ' ')[:1]
    element_field = (str(element or '').strip()).rjust(2)[:2]

    buf = [' '] * 80
    label = "ATOM"
    for i, ch in enumerate(label):
        buf[i] = ch
    serial_str = f"{serial_i:>5d}"
    for i, ch in enumerate(serial_str):
        buf[6 + i] = ch
    buf[11] = ' '
    for i, ch in enumerate(atom_name_field):
        buf[12 + i] = ch
    buf[16] = alt_loc_field
    for i, ch in enumerate(res_name_field):
        buf[17 + i] = ch
    buf[20] = ' '
    buf[21] = chain_id_field
    res_seq_str = f"{res_seq_i:>4d}"
    for i, ch in enumerate(res_seq_str):
        buf[22 + i] = ch
    buf[26] = i_code_field
    x_str = f"{x_f:8.3f}"
    for i, ch in enumerate(x_str):
        buf[30 + i] = ch
    y_str = f"{y_f:8.3f}"
    for i, ch in enumerate(y_str):
        buf[38 + i] = ch
    z_str = f"{z_f:8.3f}"
    for i, ch in enumerate(z_str):
        buf[46 + i] = ch
    occ_str = f"{occ_f:6.2f}"
    for i, ch in enumerate(occ_str):
        buf[54 + i] = ch
    temp_str = f"{temp_f:6.2f}"
    for i, ch in enumerate(temp_str):
        buf[60 + i] = ch
    for i, ch in enumerate(element_field):
        buf[76 + i] = ch

    line = "".join(buf)
    if len(line) < 80:
        line = line.ljust(80)
    elif len(line) > 80:
        line = line[:80]
    return line + "\n"

# -------------------------------------------------------------------
# HELPER FUNCTIONS
# -------------------------------------------------------------------
def write_report_csv(results, output_dir):
    """Write CSV report with incremental update (merge with existing)."""
    import csv
    report_path = output_dir / "report.csv"

    existing = {}
    if report_path.exists():
        try:
            with open(report_path, 'r', newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    key = (row['model'], row['solution'])
                    existing[key] = row
        except Exception as e:
            print(f"Warning: Could not read existing report: {e}. Will overwrite.")

    for res in results:
        model = res['model']
        if res['overall_success'] is False and not res['solutions']:
            key = (model, '')
            row = {
                'model': model,
                'solution': '',
                'status': 'failed',
                'rmsd': '',
                'reason': res.get('error', 'No solution processed')
            }
            existing[key] = row
        else:
            for sol_status in res['solutions']:
                key = (model, sol_status['solution'])
                row = {
                    'model': model,
                    'solution': sol_status['solution'],
                    'status': sol_status['status'],
                    'rmsd': f"{sol_status['rmsd']:.4f}" if sol_status['rmsd'] is not None else '',
                    'reason': sol_status.get('reason', '')
                }
                existing[key] = row

    with open(report_path, 'w', newline='') as f:
        fieldnames = ['model', 'solution', 'status', 'rmsd', 'reason']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(existing.keys()):
            writer.writerow(existing[key])

    print(f"\nReport saved to: {report_path}")

def get_successful_solutions(model, report_path):
    """Return a set of solution names that have status 'success' for the given model."""
    successful = set()
    if not report_path.exists():
        return successful
    import csv
    try:
        with open(report_path, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['model'] == model and row['status'] == 'success':
                    successful.add(row['solution'])
    except Exception as e:
        print(f"Warning: Could not read report for successful solutions: {e}")
    return successful

def display_comparison(solution_atoms, model_atoms, solution_filename="SOLUTION", model_filename="MODEL", highlight_diffs=True):
    """Display side-by-side atom lists to highlight differences."""
    print("\n" + "="*120)
    print(f"{solution_filename:^60} | {model_filename:^60}")
    print("="*120)
    max_len = max(len(solution_atoms), len(model_atoms))
    for i in range(max_len):
        sol_display = ""
        mod_display = ""
        if i < len(solution_atoms):
            _, _, sol_info = solution_atoms[i]
            sol_res_name = sol_info['res_name']
            if len(sol_res_name) != 3:
                sol_res_name = sol_res_name.ljust(3)
            sol_display = f"{sol_info['atom_num']:4d} {sol_res_name} {sol_info['chain_id']:1s} {sol_info['res_seq']:4d} {sol_info['atom_name']:>4s}"
        if i < len(model_atoms):
            _, _, mod_info = model_atoms[i]
            mod_res_name = mod_info['res_name']
            if len(mod_res_name) != 3:
                mod_res_name = mod_res_name.ljust(3)
            mod_display = f"{mod_info['atom_num']:4d} {mod_res_name} {mod_info['chain_id']:1s} {mod_info['res_seq']:4d} {mod_info['atom_name']:>4s}"
        is_different = False
        if i < len(solution_atoms) and i < len(model_atoms):
            sol_info = solution_atoms[i][2]
            mod_info = model_atoms[i][2]
            sol_res_name_clean = sol_info['res_name'].strip()
            mod_res_name_clean = mod_info['res_name'].strip()
            is_different = (sol_info['res_seq'] != mod_info['res_seq'] or
                           sol_info['atom_name'] != mod_info['atom_name'] or
                           sol_res_name_clean != mod_res_name_clean or
                           sol_info['chain_id'] != mod_info['chain_id'])
        marker = " |" if is_different and highlight_diffs else "  "
        print(f"{sol_display:60s}{marker}{mod_display:60s}")
    print("="*120)
    print(f"Total atoms: {len(solution_atoms)} in {solution_filename}, {len(model_atoms)} in {model_filename}")
    print("Legend: '|' indicates a difference between corresponding atoms (including chain ID)")

def parse_selection(selection_str, max_value):
    """Parse user input like '1,3,5-10,END' into a sorted list of integers."""
    selected = set()
    parts = selection_str.replace(' ', '').split(',')
    for part in parts:
        if '-' in part:
            if 'END' in part.upper():
                try:
                    start_str = part.split('-')[0]
                    start = int(start_str.strip())
                    if start < 1:
                        start = 1
                    if start > max_value:
                        start = max_value
                    selected.update(range(start, max_value + 1))
                except:
                    print(f"Invalid range with END: {part}")
                    continue
            else:
                try:
                    start, end = map(int, part.split('-'))
                    if start < 1:
                        start = 1
                    if end > max_value:
                        end = max_value
                    selected.update(range(start, end + 1))
                except:
                    print(f"Invalid range: {part}")
                    continue
        else:
            if part.upper() == 'END':
                selected.add(max_value)
            else:
                try:
                    num = int(part)
                    if 1 <= num <= max_value:
                        selected.add(num)
                    else:
                        print(f"Number {num} out of range (1‑{max_value})")
                except:
                    print(f"Invalid number: {part}")
                    continue
    return sorted(selected)

def delete_atoms(atom_lines, selection):
    """Delete atoms by their index (1‑based)."""
    if not selection:
        return atom_lines
    indices_to_delete = sorted([idx - 1 for idx in selection], reverse=True)
    for idx in indices_to_delete:
        if 0 <= idx < len(atom_lines):
            atom_lines.pop(idx)
    return atom_lines

def renumber_atoms(atom_lines):
    """Renumber atoms sequentially starting from 1."""
    result = []
    for new_serial, (line_num, line, info) in enumerate(atom_lines, start=1):
        new_line = format_pdb_atom_line(
            new_serial,
            info.get('atom_name', ''),
            info.get('alt_loc', ' '),
            info.get('res_name', '').strip(),
            info.get('chain_id', ' '),
            info.get('res_seq', 0),
            info.get('icode', ' '),
            info.get('x', 0.0),
            info.get('y', 0.0),
            info.get('z', 0.0),
            info.get('occupancy', '1.00'),
            info.get('temp_factor', '0.00'),
            info.get('element', '')
        )
        info['atom_num'] = new_serial
        info['line'] = new_line.rstrip('\n')
        info['original_line'] = new_line
        result.append((line_num, new_line, info))
    return result

def change_residue_number(atom_lines, old_num, new_num):
    """Change residue number from old_num to new_num."""
    updated = []
    for line_num, line, info in atom_lines:
        if info.get('res_seq') == old_num:
            info['res_seq'] = new_num
            new_line = format_pdb_atom_line(
                info.get('atom_num', 0),
                info.get('atom_name', ''),
                info.get('alt_loc', ' '),
                info.get('res_name', '').strip(),
                info.get('chain_id', ' '),
                info.get('res_seq', new_num),
                info.get('icode', ' '),
                info.get('x', 0.0),
                info.get('y', 0.0),
                info.get('z', 0.0),
                info.get('occupancy', '1.00'),
                info.get('temp_factor', '0.00'),
                info.get('element', '')
            )
            info['line'] = new_line.rstrip('\n')
            info['original_line'] = new_line
            updated.append((line_num, new_line, info))
        else:
            updated.append((line_num, line, info))
    return updated

def change_chain_id(atom_lines, old_chain, new_chain, residue_range=None):
    """Change chain ID for atoms matching old_chain, optionally restricting to a residue range."""
    updated = []
    for line_num, line, info in atom_lines:
        if info.get('chain_id') == old_chain:
            if residue_range:
                start, end = residue_range
                if not (start <= info.get('res_seq', 0) <= end):
                    updated.append((line_num, line, info))
                    continue
            info['chain_id'] = new_chain
            new_line = format_pdb_atom_line(
                info.get('atom_num', 0),
                info.get('atom_name', ''),
                info.get('alt_loc', ' '),
                info.get('res_name', '').strip(),
                info.get('chain_id', ' '),
                info.get('res_seq', 0),
                info.get('icode', ' '),
                info.get('x', 0.0),
                info.get('y', 0.0),
                info.get('z', 0.0),
                info.get('occupancy', '1.00'),
                info.get('temp_factor', '0.00'),
                info.get('element', '')
            )
            info['line'] = new_line.rstrip('\n')
            info['original_line'] = new_line
            updated.append((line_num, new_line, info))
        else:
            updated.append((line_num, line, info))
    return updated

def auto_renumber_residues(atom_lines):
    """Renumber residues sequentially across all chains."""
    if not atom_lines:
        return atom_lines
    groups = OrderedDict()
    for _, _, info in atom_lines:
        key = (info['chain_id'], info['res_name'], info['res_seq'])
        if key not in groups:
            groups[key] = []
        groups[key].append(info)
    new_res_seq = 1
    updated_atoms = []
    for (chain_id, res_name, old_res_seq), atoms in groups.items():
        for info in atoms:
            info['res_seq'] = new_res_seq
            res_seq_str = f"{new_res_seq:>4}"
            line = info['original_line']
            new_line = line[:22] + res_seq_str + line[26:]
            info['line'] = new_line.rstrip('\n')
            info['original_line'] = new_line
            updated_atoms.append(info)
        new_res_seq += 1
    result = []
    for i, info in enumerate(updated_atoms, 1):
        info['atom_num'] = i
        line = info['original_line']
        atom_num_str = f"{i:>5}"
        new_line = line[:6] + atom_num_str + line[11:]
        info['line'] = new_line.rstrip('\n')
        info['original_line'] = new_line
        result.append((0, new_line, info))
    return result

def get_residue_info(atom_lines):
    """Return a dictionary of residue information."""
    residues = {}
    for _, _, atom_info in atom_lines:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        if key not in residues:
            residues[key] = {
                'chain_id': atom_info['chain_id'],
                'res_seq': atom_info['res_seq'],
                'res_name': atom_info['res_name'],
                'atom_count': 0,
                'atoms': []
            }
        residues[key]['atom_count'] += 1
        residues[key]['atoms'].append(atom_info)
    return residues

def display_residues(atom_lines, target_name="SOLUTION"):
    """Display residue table."""
    residues = get_residue_info(atom_lines)
    print(f"\n{'='*80}")
    print(f"RESIDUES IN {target_name}: {len(residues)} residue(s)")
    print(f"{'='*80}")
    print(f"{'No':>4} {'Chain':>5} {'ResSeq':>6} {'ResName':>8} {'Atoms':>6} {'NucType':>10}")
    print(f"{'-'*80}")
    for i, ((chain_id, res_seq), residue_info) in enumerate(sorted(residues.items(), key=lambda x: (x[0][0], x[0][1])), 1):
        res_name = residue_info['res_name'].strip()
        nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
        print(f"{i:4d} {chain_id:>5} {res_seq:6d} {res_name:>8} {residue_info['atom_count']:6d} {nuc_type:>10}")

def remove_residue(atom_lines, chain_id, res_seq_input):
    """Remove entire residue(s) by chain and residue number(s)."""
    if isinstance(res_seq_input, int):
        res_seq_list = [res_seq_input]
    elif isinstance(res_seq_input, str):
        max_res_seq = 0
        for _, _, atom_info in atom_lines:
            if atom_info['chain_id'] == chain_id:
                max_res_seq = max(max_res_seq, atom_info['res_seq'])
        if max_res_seq == 0:
            print(f"No residues found in chain {chain_id}")
            return atom_lines
        res_seq_list = parse_selection(res_seq_input, max_res_seq)
    elif isinstance(res_seq_input, list):
        res_seq_list = res_seq_input
    else:
        print(f"Invalid input format: {res_seq_input}")
        return atom_lines
    if not res_seq_list:
        return atom_lines
    atoms_to_keep = []
    removed_count = 0
    removed_residues = set()
    for line_num, line, atom_info in atom_lines:
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] in res_seq_list:
            removed_count += 1
            removed_residues.add(atom_info['res_seq'])
            continue
        atoms_to_keep.append((line_num, line, atom_info))
    if removed_count > 0:
        removed_list = sorted(removed_residues)
        print(f"Removed residue(s) {chain_id}:{removed_list} ({removed_count} atoms)")
    return atoms_to_keep

def simple_swap_residues(atom_lines, chain_id, res_seq1, res_seq2):
    """Simple swap of residue numbers (used when residues overlap)."""
    updated_atoms = []
    for line_num, line, atom_info in atom_lines:
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq1:
            atom_info['res_seq'] = res_seq2
            res_seq_str = f"{res_seq2:>4}"
            new_line = line[:22] + res_seq_str + line[26:]
            atom_info['line'] = new_line.rstrip('\n')
            atom_info['original_line'] = new_line
            updated_atoms.append((line_num, new_line, atom_info))
        elif atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq2:
            atom_info['res_seq'] = res_seq1
            res_seq_str = f"{res_seq1:>4}"
            new_line = line[:22] + res_seq_str + line[26:]
            atom_info['line'] = new_line.rstrip('\n')
            atom_info['original_line'] = new_line
            updated_atoms.append((line_num, new_line, atom_info))
        else:
            updated_atoms.append((line_num, line, atom_info))
    return updated_atoms

def swap_residues(atom_lines, chain_id, res_seq1, res_seq2):
    """Swap the positions of two residues."""
    atoms_res1 = []
    atoms_res2 = []
    other_atoms = []
    for idx, (line_num, line, atom_info) in enumerate(atom_lines):
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq1:
            atoms_res1.append((idx, line_num, line, atom_info.copy()))
        elif atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq2:
            atoms_res2.append((idx, line_num, line, atom_info.copy()))
        else:
            other_atoms.append((idx, line_num, line, atom_info.copy()))
    if not atoms_res1 or not atoms_res2:
        print("Error: One of the residues does not exist!")
        return atom_lines
    print(f"Found: residue {chain_id}:{res_seq1} - {len(atoms_res1)} atoms, residue {chain_id}:{res_seq2} - {len(atoms_res2)} atoms")
    res1_indices = [idx for idx, _, _, _ in atoms_res1]
    res2_indices = [idx for idx, _, _, _ in atoms_res2]
    min1, max1 = min(res1_indices), max(res1_indices)
    min2, max2 = min(res2_indices), max(res2_indices)
    if (min1 <= max2 and max1 >= min2):
        print("WARNING: Residues overlap! Using simple swap.")
        return simple_swap_residues(atom_lines, chain_id, res_seq1, res_seq2)
    all_items = []
    current_idx = 0
    while current_idx < len(atom_lines):
        if current_idx == min1:
            for _, line_num, line, atom_info in atoms_res2:
                new_info = atom_info.copy()
                new_info['res_seq'] = res_seq1
                res_seq_str = f"{res_seq1:>4}"
                new_line = line[:22] + res_seq_str + line[26:]
                new_info['line'] = new_line.rstrip('\n')
                new_info['original_line'] = new_line
                all_items.append((current_idx, new_line, new_info))
                current_idx += 1
            current_idx = max1 + 1
        elif current_idx == min2:
            for _, line_num, line, atom_info in atoms_res1:
                new_info = atom_info.copy()
                new_info['res_seq'] = res_seq2
                res_seq_str = f"{res_seq2:>4}"
                new_line = line[:22] + res_seq_str + line[26:]
                new_info['line'] = new_line.rstrip('\n')
                new_info['original_line'] = new_line
                all_items.append((current_idx, new_line, new_info))
                current_idx += 1
            current_idx = max2 + 1
        else:
            found = False
            for idx, line_num, line, atom_info in other_atoms:
                if idx == current_idx:
                    all_items.append((current_idx, line, atom_info))
                    found = True
                    break
            if not found:
                pass
            current_idx += 1
    all_items.sort(key=lambda x: x[0])
    updated_atoms = [(line_num, line, atom_info) for _, line, atom_info in all_items]
    updated_atoms = renumber_atoms(updated_atoms)
    return updated_atoms

def change_residue_type(atom_lines, chain_id, res_seq, new_res_name):
    """Change the residue name for a given residue."""
    updated_atoms = []
    changed_count = 0
    for line_num, line, atom_info in atom_lines:
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq:
            new_line = format_pdb_atom_line(
                atom_info.get('atom_num', 0),
                atom_info.get('atom_name', ''),
                atom_info.get('alt_loc', ' '),
                new_res_name,
                chain_id,
                res_seq,
                atom_info.get('icode', ' '),
                atom_info.get('x', 0.0),
                atom_info.get('y', 0.0),
                atom_info.get('z', 0.0),
                atom_info.get('occupancy', '1.00'),
                atom_info.get('temp_factor', '0.00'),
                atom_info.get('element', '')
            )
            atom_info['res_name'] = str(new_res_name).strip().rjust(3)
            atom_info['line'] = new_line.rstrip('\n')
            atom_info['original_line'] = new_line
            changed_count += 1
            updated_atoms.append((line_num, new_line, atom_info))
        else:
            updated_atoms.append((line_num, line, atom_info))
    if changed_count > 0:
        print(f"✓ Changed residue type {chain_id}:{res_seq} to {new_res_name} ({changed_count} atoms)")
    return updated_atoms

# -------------------------------------------------------------------
# ADVANCED OPERATIONS
# -------------------------------------------------------------------
def basic_advanced_change_residue_type(atom_lines, chain_id, res_seq, new_nuc_type):
    """Advanced residue type change with atom set replacement."""
    residue_atoms = []
    for line_num, line, atom_info in atom_lines:
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq:
            residue_atoms.append(atom_info.copy())
    if not residue_atoms:
        print(f"Residue {chain_id}:{res_seq} not found")
        return atom_lines, []
    old_res_name = residue_atoms[0]['res_name'].strip()
    old_nuc_type = RESIDUE_NAME_MAP.get(old_res_name, old_res_name[0] if old_res_name else '?')
    old_atom_list = NUCLEOTIDE_ATOM_ORDER.get(old_nuc_type, [])
    new_atom_list = NUCLEOTIDE_ATOM_ORDER.get(new_nuc_type, [])
    if not new_atom_list:
        print(f"ERROR: No atom list for new type '{new_nuc_type}'")
        return atom_lines, []
    old_atom_set = set(old_atom_list)
    new_atom_set = set(new_atom_list)
    atoms_to_remove_from_old = old_atom_set - new_atom_set
    removed_atoms_info = []
    for atom in residue_atoms:
        if atom['atom_name'] in atoms_to_remove_from_old:
            removed_atoms_info.append({
                'chain_id': chain_id,
                'res_seq': res_seq,
                'atom_name': atom['atom_name']
            })
    atom_mapping = {
        ('C', 'U'): {'N4': 'O4'},
        ('U', 'C'): {'O4': 'N4'},
        ('A', 'G'): {'N6': 'O6'},
        ('G', 'A'): {'O6': 'N6'},
    }
    mapping = atom_mapping.get((old_nuc_type, new_nuc_type), {})
    existing_atoms = {a['atom_name']: a for a in residue_atoms}
    new_residue_atom_infos = []
    for new_atom_name in new_atom_list:
        mapped_from = None
        for k, v in mapping.items():
            if v == new_atom_name:
                mapped_from = k
                break
        if mapped_from and mapped_from in existing_atoms:
            ai = existing_atoms[mapped_from].copy()
            ai['atom_name'] = new_atom_name
            ai['element'] = new_atom_name[0]
            new_residue_atom_infos.append(ai)
        elif new_atom_name in existing_atoms:
            ai = existing_atoms[new_atom_name].copy()
            new_residue_atom_infos.append(ai)
    for ai in residue_atoms:
        if ai['atom_name'] not in old_atom_list and ai['atom_name'] not in [a['atom_name'] for a in new_residue_atom_infos]:
            ai_copy = ai.copy()
            ai_copy['res_name'] = new_nuc_type.strip().rjust(3)
            new_residue_atom_infos.append(ai_copy)
    new_res_name_str = new_nuc_type.strip().rjust(3)
    for ai in new_residue_atom_infos:
        ai['res_name'] = new_res_name_str
    insert_pos = None
    old_count = 0
    for idx, (line_num, line, atom_info) in enumerate(atom_lines):
        if atom_info['chain_id'] == chain_id and atom_info['res_seq'] == res_seq:
            if insert_pos is None:
                insert_pos = idx
            old_count += 1
    if insert_pos is None:
        return atom_lines, removed_atoms_info
    updated_atoms = []
    updated_atoms.extend(atom_lines[:insert_pos])
    for ai in new_residue_atom_infos:
        serial_tmp = ai.get('atom_num', 0) or 0
        new_line = format_pdb_atom_line(
            serial_tmp,
            ai.get('atom_name', ''),
            ai.get('alt_loc', ' '),
            new_nuc_type,
            chain_id,
            res_seq,
            ai.get('icode', ' '),
            ai.get('x', 0.0),
            ai.get('y', 0.0),
            ai.get('z', 0.0),
            ai.get('occupancy', '1.00'),
            ai.get('temp_factor', '0.00'),
            ai.get('element', '')
        )
        ai['line'] = new_line.rstrip('\n')
        ai['original_line'] = new_line
        ai['atom_num'] = 0
        updated_atoms.append((0, new_line, ai))
    updated_atoms.extend(atom_lines[insert_pos + old_count:])
    updated_atoms = renumber_atoms(updated_atoms)
    return updated_atoms, removed_atoms_info

def remove_matching_atoms(atom_lines, atoms_to_remove):
    """Remove atoms matching a list of criteria (chain, res_seq, atom_name)."""
    if not atoms_to_remove:
        return atom_lines
    atoms_to_keep = []
    removed_count = 0
    for line_num, line, atom_info in atom_lines:
        should_remove = False
        for removal_info in atoms_to_remove:
            match = True
            if 'chain_id' in removal_info and atom_info['chain_id'] != removal_info['chain_id']:
                match = False
            if 'res_seq' in removal_info and atom_info['res_seq'] != removal_info['res_seq']:
                match = False
            if 'atom_name' in removal_info and atom_info['atom_name'] != removal_info['atom_name']:
                match = False
            if match:
                should_remove = True
                removed_count += 1
                break
        if not should_remove:
            atoms_to_keep.append((line_num, line, atom_info))
    if removed_count > 0:
        print(f"  Removed {removed_count} matching atoms")
    return atoms_to_keep

def auto_correct_nucleotides_by_residue(solution_atoms, model_atoms):
    """Automatically correct model nucleotides to match solution types by residue number."""
    print("\n" + "="*80)
    print("AUTOMATIC NUCLEOTIDE CORRECTION (SOLUTION → MODEL)")
    print("(Change model nucleotides to solution types based on residue numbering)")
    print("="*80)
    model_residues = {}
    for _, _, atom_info in model_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        if key not in model_residues:
            res_name = atom_info['res_name'].strip()
            nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            model_residues[key] = nuc_type
    solution_residues = {}
    for _, _, atom_info in solution_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        if key not in solution_residues:
            res_name = atom_info['res_name'].strip()
            nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            solution_residues[key] = nuc_type
    print(f"Found {len(model_residues)} residues in model")
    print(f"Found {len(solution_residues)} residues in solution")
    changed_count = 0
    same_count = 0
    not_found_count = 0
    for (chain_id, res_seq), model_nuc_type in sorted(model_residues.items()):
        if (chain_id, res_seq) in solution_residues:
            solution_nuc_type = solution_residues[(chain_id, res_seq)]
            if solution_nuc_type != model_nuc_type:
                model_atoms, removed_atoms = basic_advanced_change_residue_type(
                    model_atoms, chain_id, res_seq, solution_nuc_type
                )
                if removed_atoms:
                    solution_atoms = remove_matching_atoms(solution_atoms, removed_atoms)
                changed_count += 1
            else:
                same_count += 1
        else:
            print(f"  Warning: {chain_id}:{res_seq} does not exist in solution")
            not_found_count += 1
    print(f"\nSummary:")
    print(f"  Changed: {changed_count} nucleotides")
    print(f"  Unchanged: {same_count} nucleotides")
    print(f"  Not found in solution: {not_found_count} residues")
    solution_atoms = renumber_atoms(solution_atoms)
    model_atoms = renumber_atoms(model_atoms)
    print(f"\n✓ Correction finished!")
    print(f"  Solution: {len(solution_atoms)} atoms")
    print(f"  Model: {len(model_atoms)} atoms")
    print("="*80)
    return solution_atoms, model_atoms

def auto_correct_nucleotides_by_residue_reverse(solution_atoms, model_atoms):
    """Automatically correct solution nucleotides to match model types by residue number."""
    print("\n" + "="*80)
    print("AUTOMATIC NUCLEOTIDE CORRECTION (MODEL → SOLUTION)")
    print("(Change solution nucleotides to model types based on residue numbering)")
    print("="*80)
    model_residues = {}
    for _, _, atom_info in model_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        if key not in model_residues:
            res_name = atom_info['res_name'].strip()
            nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            model_residues[key] = nuc_type
    solution_residues = {}
    for _, _, atom_info in solution_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        if key not in solution_residues:
            res_name = atom_info['res_name'].strip()
            nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            solution_residues[key] = nuc_type
    print(f"Found {len(solution_residues)} residues in solution")
    print(f"Found {len(model_residues)} residues in model")
    changed_count = 0
    same_count = 0
    not_found_count = 0
    for (chain_id, res_seq), sol_nuc_type in sorted(solution_residues.items()):
        if (chain_id, res_seq) in model_residues:
            model_nuc_type = model_residues[(chain_id, res_seq)]
            if model_nuc_type != sol_nuc_type:
                solution_atoms, removed_atoms = basic_advanced_change_residue_type(
                    solution_atoms, chain_id, res_seq, model_nuc_type
                )
                if removed_atoms:
                    model_atoms = remove_matching_atoms(model_atoms, removed_atoms)
                changed_count += 1
            else:
                same_count += 1
        else:
            print(f"  Warning: {chain_id}:{res_seq} does not exist in model")
            not_found_count += 1
    print(f"\nSummary:")
    print(f"  Changed: {changed_count} nucleotides")
    print(f"  Unchanged: {same_count} nucleotides")
    print(f"  Not found in model: {not_found_count} residues")
    solution_atoms = renumber_atoms(solution_atoms)
    model_atoms = renumber_atoms(model_atoms)
    print(f"\n✓ Correction finished!")
    print(f"  Solution: {len(solution_atoms)} atoms")
    print(f"  Model: {len(model_atoms)} atoms")
    print("="*80)
    return solution_atoms, model_atoms

def brute_force_sync_atoms(solution_atoms, model_atoms):
    """Keep only atoms that are common to both structures (by chain, residue, atom name)."""
    print("\n" + "="*80)
    print("BRUTE FORCE SYNCHRONIZATION BETWEEN FILES")
    print("(Keep only common residues and common atoms within them)")
    print("="*80)
    solution_res = {}
    for _, _, atom_info in solution_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        solution_res.setdefault(key, set()).add(atom_info['atom_name'])
    model_res = {}
    for _, _, atom_info in model_atoms:
        key = (atom_info['chain_id'], atom_info['res_seq'])
        model_res.setdefault(key, set()).add(atom_info['atom_name'])

    common_residues = set(solution_res.keys()) & set(model_res.keys())
    print(f"Solution residues: {len(solution_res)}, Model residues: {len(model_res)}")
    print(f"Common residues: {len(common_residues)}")

    common_atoms_per_res = {}
    for key in common_residues:
        common = solution_res[key] & model_res[key]
        common_atoms_per_res[key] = common
        # Optional verbose output for small sets
        if len(common_residues) <= 10:
            sol_only = solution_res[key] - common
            mod_only = model_res[key] - common
            if sol_only or mod_only:
                print(f"\nResidue {key[0]}:{key[1]}:")
                if sol_only:
                    print(f"  Only in SOLUTION: {sorted(sol_only)}")
                if mod_only:
                    print(f"  Only in MODEL: {sorted(mod_only)}")
                print(f"  Common atoms: {len(common)}")

    # Build new atom lists
    new_solution = []
    removed_sol_residues = 0
    removed_sol_atoms = 0
    for line_num, line, info in solution_atoms:
        key = (info['chain_id'], info['res_seq'])
        if key in common_residues:
            if info['atom_name'] in common_atoms_per_res[key]:
                new_solution.append((line_num, line, info))
            else:
                removed_sol_atoms += 1
        else:
            removed_sol_residues += 1

    new_model = []
    removed_mod_residues = 0
    removed_mod_atoms = 0
    for line_num, line, info in model_atoms:
        key = (info['chain_id'], info['res_seq'])
        if key in common_residues:
            if info['atom_name'] in common_atoms_per_res[key]:
                new_model.append((line_num, line, info))
            else:
                removed_mod_atoms += 1
        else:
            removed_mod_residues += 1

    print(f"\nSummary:")
    print(f"  Removed from SOLUTION: {removed_sol_residues} whole residues, {removed_sol_atoms} extra atoms")
    print(f"  Removed from MODEL:    {removed_mod_residues} whole residues, {removed_mod_atoms} extra atoms")
    print(f"  Kept in SOLUTION: {len(new_solution)} atoms")
    print(f"  Kept in MODEL:    {len(new_model)} atoms")
    new_solution = renumber_atoms(new_solution)
    new_model = renumber_atoms(new_model)
    print("✓ Synchronization finished!")
    return new_solution, new_model

def remove_nonstandard_atoms(atom_lines, target_name="structure"):
    """
    Remove atoms that are not part of the standard nucleotide atom order.
    For each residue, determine its nucleotide type (A, C, G, U) and keep only
    atoms that appear in NUCLEOTIDE_ATOM_ORDER for that type.
    Non‑nucleotide residues (e.g., water, ligands) are left unchanged.
    """
    print(f"\nAutomatically removing non‑standard atoms from {target_name}...")
    new_atoms = []
    removed = 0
    residues = get_residue_info(atom_lines)
    for (chain_id, res_seq), residue_info in residues.items():
        res_name = residue_info['res_name'].strip()
        nuc_type = RESIDUE_NAME_MAP.get(res_name, None)
        expected_set = set(NUCLEOTIDE_ATOM_ORDER.get(nuc_type, [])) if nuc_type else None

        for atom_info in residue_info['atoms']:
            keep = True
            if expected_set is not None:
                if atom_info['atom_name'] not in expected_set:
                    keep = False
                    removed += 1
            if keep:
                new_atoms.append(atom_info)

    # Re‑create atom lines from the kept info objects (renumbering will happen later)
    result = []
    for i, info in enumerate(new_atoms, 1):
        new_line = format_pdb_atom_line(
            i,
            info['atom_name'],
            info['alt_loc'],
            info['res_name'],
            info['chain_id'],
            info['res_seq'],
            info['icode'],
            info['x'], info['y'], info['z'],
            info['occupancy'],
            info['temp_factor'],
            info['element']
        )
        info['atom_num'] = i
        info['line'] = new_line.rstrip('\n')
        info['original_line'] = new_line
        result.append((0, new_line, info))

    print(f"Removed {removed} non‑standard atoms.")
    return result

# -------------------------------------------------------------------
# RMSD CALCULATION
# -------------------------------------------------------------------
def calculate_rmsd(atoms1, atoms2):
    """Calculate RMSD between two sets of atoms (must be same length)."""
    if len(atoms1) != len(atoms2):
        return None
    if len(atoms1) == 0:
        return None
    coords1 = np.array([(info['x'], info['y'], info['z']) for _, _, info in atoms1])
    coords2 = np.array([(info['x'], info['y'], info['z']) for _, _, info in atoms2])
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

# -------------------------------------------------------------------
# CHANGE HISTORY SYSTEM
# -------------------------------------------------------------------
class ChangeHistory:
    def __init__(self):
        self.history = []
        self.redo_stack = []
    def record_change(self, change_type, target, before_state, description=""):
        before_copy = [(idx, line, info.copy()) for idx, line, info in before_state]
        self.history.append({
            'type': change_type,
            'target': target,
            'before': before_copy,
            'description': description,
            'timestamp': datetime.datetime.now()
        })
        self.redo_stack.clear()
    def undo_last(self, current_solution_atoms, current_model_atoms):
        if not self.history:
            print("No operations to undo!")
            return current_solution_atoms, current_model_atoms, False
        last_change = self.history.pop()
        target = last_change['target']
        before_state = last_change['before']
        print(f"\nUNDOING OPERATION: {last_change['type']} on {target}")
        if target == 'solution':
            current_solution_copy = [(idx, line, info.copy()) for idx, line, info in current_solution_atoms]
            self.redo_stack.append({
                'type': last_change['type'],
                'target': target,
                'before': current_solution_copy,
                'description': f"Redo of {last_change['type']}"
            })
            current_solution_atoms = before_state
        elif target == 'model':
            current_model_copy = [(idx, line, info.copy()) for idx, line, info in current_model_atoms]
            self.redo_stack.append({
                'type': last_change['type'],
                'target': target,
                'before': current_model_copy,
                'description': f"Redo of {last_change['type']}"
            })
            current_model_atoms = before_state
        current_solution_atoms = renumber_atoms(current_solution_atoms)
        current_model_atoms = renumber_atoms(current_model_atoms)
        print(f"✓ Undo completed. {len(self.history)} operations remaining.")
        return current_solution_atoms, current_model_atoms, True
    def display_history(self):
        if not self.history:
            print("History empty.")
            return
        print("\nCHANGE HISTORY:")
        for i, ch in enumerate(reversed(self.history), 1):
            print(f"{i}. {ch['target']}: {ch.get('description', ch['type'])}")
    def clear_history(self):
        self.history.clear()
        self.redo_stack.clear()

# -------------------------------------------------------------------
# OPERATION RECORDER (for batch application)
# -------------------------------------------------------------------
class OperationRecorder:
    def __init__(self):
        self.operations = []
    def add_operation(self, op_type, target, *args):
        self.operations.append({
            'type': op_type,
            'target': target,
            'args': args,
            'timestamp': datetime.datetime.now()
        })
    def clear(self):
        self.operations.clear()
    def get_operations(self):
        return self.operations.copy()
    def display_operations(self):
        if not self.operations:
            print("No operations recorded.")
            return
        print("\n" + "="*80)
        print("RECORDED OPERATIONS:")
        print("="*80)
        for i, op in enumerate(self.operations, 1):
            target_display = "SOLUTION" if op['target'] == 'solution' else "MODEL"
            args_str = str(op['args'])
            print(f"{i}. {op['type']} on {target_display}: {args_str}")
        print("="*80)

# -------------------------------------------------------------------
# REMARK GENERATION
# -------------------------------------------------------------------
def generate_remarks(operations, target):
    """Return a list of REMARK lines for the given target (solution/model)."""
    remarks = []
    target_ops = [op for op in operations if op['target'] == target]
    target_ops.sort(key=lambda op: op['timestamp'])   # preserve order of execution
    for op in target_ops:
        remark = _remark_from_operation(op)
        if remark:
            remarks.append(f"REMARK - {remark}\n")
    return remarks

def _remark_from_operation(op):
    op_type = op['type']
    args = op['args']
    if op_type == 'delete_atoms':
        return f"Deleted atoms: {args[0]}"
    elif op_type == 'change_chain':
        old = args[0]
        new = args[1]
        if len(args) > 2 and args[2] is not None:
            return f"Changed chain {old} -> {new} for residues {args[2]}"
        else:
            return f"Changed chain {old} -> {new}"
    elif op_type == 'auto_renumber':
        return "Automatically renumbered residues"
    elif op_type == 'remove_residue':
        chain = args[0]
        res = args[1]
        return f"Removed residue(s) {chain}:{res}"
    elif op_type == 'swap_residues':
        chain = args[0]
        r1, r2 = args[1], args[2]
        return f"Swapped residues {chain}:{r1} and {chain}:{r2}"
    elif op_type == 'change_residue_type':
        chain = args[0]
        res = args[1]
        new_type = args[2]
        return f"Changed residue type {chain}:{res} to {new_type}"
    elif op_type == 'advanced_change':
        chain = args[0]
        res = args[1]
        new_type = args[2]
        return f"Advanced change of residue {chain}:{res} -> {new_type}"
    elif op_type == 'auto_correct_nucleotides':
        return "Automatically corrected nucleotides (according to solution)"
    elif op_type == 'auto_correct_nucleotides_reverse':
        return "Automatically corrected nucleotides (according to model)"
    elif op_type == 'brute_force_sync':
        return "Synchronized atoms (removed non-matching)"
    elif op_type == 'renumber_atoms':
        return "Renumbered all atoms sequentially"
    elif op_type == 'remove_extra_atoms':
        return f"Removed extra atoms from residues {args[0]}"
    elif op_type == 'change_residue':
        old = args[0]
        new = args[1]
        return f"Changed residue number {old} -> {new}"
    else:
        return f"Performed operation: {op_type}"

# -------------------------------------------------------------------
# HELPER FUNCTIONS FOR MODEL/SOLUTION SELECTION
# -------------------------------------------------------------------
def extract_base_id(filename):
    base_name = pathlib.Path(filename).stem
    patterns_to_remove = [r'_refined$', r'_model$', r'_solution$', r'_clean$']
    for pattern in patterns_to_remove:
        base_name = re.sub(pattern, '', base_name)
    return base_name

def select_models(models, prompt, labels=None):
    """Interactive selection of models from a list."""
    if not models:
        print("No models to choose from.")
        return []
    print("\n" + "="*80)
    print(f"{prompt} - Available models ({len(models)}):")
    print("="*80)
    for i, model in enumerate(models, 1):
        label = ""
        if labels and i <= len(labels):
            label = f" {labels[i-1]}"
        print(f"{i:4d}. {model}{label}")
    print("="*80)
    while True:
        print("\nYou can choose:")
        print("  - 'all' - all models")
        print("  - numbers: 1, 2, 3")
        print("  - ranges: 1-10, 135-148")
        print("  - patterns: TS029, TSR01")
        print("  - multiple patterns separated by commas")
        user_input = input("\nSelect models: ").strip()
        if not user_input:
            continue
        if user_input.lower() == 'all':
            return models.copy()
        parts = [part.strip() for part in user_input.split(',') if part.strip()]
        selected = []
        for part in parts:
            if '-' in part and not any(c.isalpha() for c in part):
                try:
                    start, end = map(int, part.split('-'))
                    start = max(1, min(start, len(models)))
                    end = max(start, min(end, len(models)))
                    selected.extend(models[start-1:end])
                except:
                    print(f"Invalid range: {part}")
            elif any(c.isalpha() for c in part):
                pattern = part.upper()
                matches = [m for m in models if pattern in m.upper()]
                selected.extend(matches)
                if matches:
                    print(f"  Found {len(matches)} models containing '{pattern}'")
            else:
                try:
                    num = int(part)
                    if 1 <= num <= len(models):
                        selected.append(models[num-1])
                    else:
                        print(f"Number {num} out of range")
                except:
                    print(f"Invalid value: {part}")
        if selected:
            unique = []
            seen = set()
            for item in selected:
                if item not in seen:
                    seen.add(item)
                    unique.append(item)
            return unique
        print("No models selected.")

def collect_model_files_from_external(puzzle_path):
    return [f for f in os.listdir(puzzle_path) if pathlib.Path(f).suffix in ['.pdb', '.cif'] and 'solution' not in f.lower()]

def collect_solution_files(puzzle_path):
    return [f for f in os.listdir(puzzle_path) if pathlib.Path(f).suffix in ['.pdb', '.cif'] and 'solution' in f.lower()]

# -------------------------------------------------------------------
# BATCH APPLY OPERATIONS
# -------------------------------------------------------------------
def apply_operations_to_pair(solution_atoms, model_atoms, operations):
    if not operations:
        return solution_atoms, model_atoms
    for op in operations:
        print(f"  [{op['type']}] {op['args']}")
        if op['target'] == 'solution':
            if op['type'] == 'delete_atoms':
                solution_atoms = delete_atoms(solution_atoms, op['args'][0])
                solution_atoms = renumber_atoms(solution_atoms)
            elif op['type'] == 'change_chain':
                old_chain, new_chain = op['args'][:2]
                residue_range = op['args'][2] if len(op['args']) > 2 else None
                solution_atoms = change_chain_id(solution_atoms, old_chain, new_chain, residue_range)
            elif op['type'] == 'auto_renumber':
                solution_atoms = auto_renumber_residues(solution_atoms)
            elif op['type'] == 'remove_residue':
                chain_id, res_seq_input = op['args']
                solution_atoms = remove_residue(solution_atoms, chain_id, res_seq_input)
                solution_atoms = renumber_atoms(solution_atoms)
            elif op['type'] == 'auto_correct_nucleotides':
                solution_atoms, model_atoms = auto_correct_nucleotides_by_residue(solution_atoms, model_atoms)
            elif op['type'] == 'auto_correct_nucleotides_reverse':
                solution_atoms, model_atoms = auto_correct_nucleotides_by_residue_reverse(solution_atoms, model_atoms)
            elif op['type'] == 'brute_force_sync':
                solution_atoms, model_atoms = brute_force_sync_atoms(solution_atoms, model_atoms)
        elif op['target'] == 'model':
            if op['type'] == 'delete_atoms':
                model_atoms = delete_atoms(model_atoms, op['args'][0])
                model_atoms = renumber_atoms(model_atoms)
            elif op['type'] == 'change_chain':
                old_chain, new_chain = op['args'][:2]
                residue_range = op['args'][2] if len(op['args']) > 2 else None
                model_atoms = change_chain_id(model_atoms, old_chain, new_chain, residue_range)
            elif op['type'] == 'auto_renumber':
                model_atoms = auto_renumber_residues(model_atoms)
            elif op['type'] == 'remove_residue':
                chain_id, res_seq_input = op['args']
                model_atoms = remove_residue(model_atoms, chain_id, res_seq_input)
                model_atoms = renumber_atoms(model_atoms)
            elif op['type'] == 'auto_correct_nucleotides':
                solution_atoms, model_atoms = auto_correct_nucleotides_by_residue(solution_atoms, model_atoms)
            elif op['type'] == 'auto_correct_nucleotides_reverse':
                solution_atoms, model_atoms = auto_correct_nucleotides_by_residue_reverse(solution_atoms, model_atoms)
            elif op['type'] == 'brute_force_sync':
                solution_atoms, model_atoms = brute_force_sync_atoms(solution_atoms, model_atoms)
    return solution_atoms, model_atoms

# -------------------------------------------------------------------
# PROCESS A SINGLE MODEL WITH ALL SOLUTIONS (BATCH)
# -------------------------------------------------------------------
def process_model_with_all_solutions(model_file, puzzle_path, output_dir, operations, solutions_filter=None):
    """Process a model against solutions (all if filter is None) and return a status report dict."""
    print(f"\nProcessing model: {model_file}")
    model_path = puzzle_path / model_file
    try:
        model_lines = read_pdb_lines(model_path)
    except Exception as e:
        return {
            'model': model_file,
            'overall_success': False,
            'solutions': [],
            'error': f"Failed to read model file: {e}"
        }

    model_atoms_original = get_atom_lines(model_lines)
    model_atoms_original, _ = standardize_atom_order(model_atoms_original, model_atoms_original)

    # Determine which solution files to process
    if solutions_filter is None:
        solution_files = collect_solution_files(puzzle_path)
    else:
        solution_files = [f for f in solutions_filter if (puzzle_path / f).exists()]

    if not solution_files:
        print(f"  No solution files – skipping.")
        return {
            'model': model_file,
            'overall_success': False,
            'solutions': [],
            'error': "No solution files found"
        }

    saved_pairs = []
    solution_statuses = []

    for sol_file in solution_files:
        print(f"  Comparing with {sol_file}...")
        sol_path = puzzle_path / sol_file
        try:
            sol_lines = read_pdb_lines(sol_path)
        except Exception as e:
            solution_statuses.append({
                'solution': sol_file,
                'status': 'failed',
                'reason': f"Read error: {e}",
                'rmsd': None
            })
            continue

        sol_atoms = get_atom_lines(sol_lines)
        try:
            sol_atoms, model_atoms = standardize_atom_order(sol_atoms, model_atoms_original.copy())
        except Exception as e:
            solution_statuses.append({
                'solution': sol_file,
                'status': 'failed',
                'reason': f"Standardization error: {e}",
                'rmsd': None
            })
            continue

        try:
            sol_atoms, model_atoms = apply_operations_to_pair(sol_atoms, model_atoms, operations)
        except Exception as e:
            solution_statuses.append({
                'solution': sol_file,
                'status': 'failed',
                'reason': f"Operation application error: {e}",
                'rmsd': None
            })
            continue

        if len(sol_atoms) != len(model_atoms):
            solution_statuses.append({
                'solution': sol_file,
                'status': 'failed',
                'reason': f"Atom count mismatch: solution {len(sol_atoms)} vs model {len(model_atoms)}",
                'rmsd': None
            })
            continue

        rmsd = calculate_rmsd(sol_atoms, model_atoms)
        if rmsd is None:
            solution_statuses.append({
                'solution': sol_file,
                'status': 'failed',
                'reason': "RMSD calculation failed",
                'rmsd': None
            })
            continue

        print(f"    RMSD = {rmsd:.4f}")

        # Save files
        model_stem = pathlib.Path(model_file).stem
        sol_stem = pathlib.Path(sol_file).stem
        sol_short = sol_stem.replace('solution', '').replace('SOLUTION', '').strip('_')
        if not sol_short:
            sol_short = sol_stem
        sol_output = output_dir / f"{model_stem}_vs_{sol_short}_solution.pdb"
        model_output = output_dir / f"{model_stem}_vs_{sol_short}_refined.pdb"

        updated_solution_lines = update_full_pdb(sol_lines, sol_atoms)
        updated_model_lines = update_full_pdb(model_lines, model_atoms)

        # Add REMARK lines
        sol_remarks = generate_remarks(operations, 'solution')
        mod_remarks = generate_remarks(operations, 'model')
        updated_solution_lines = sol_remarks + updated_solution_lines
        updated_model_lines = mod_remarks + updated_model_lines

        write_pdb_lines(sol_output, updated_solution_lines)
        write_pdb_lines(model_output, updated_model_lines)

        print(f"    ✓ Saved: {model_output.name} ({len(model_atoms)} atoms)")
        print(f"      Saved: {sol_output.name} ({len(sol_atoms)} atoms)")

        saved_pairs.append((model_file, sol_file, rmsd))
        solution_statuses.append({
            'solution': sol_file,
            'status': 'success',
            'reason': None,
            'rmsd': rmsd
        })

    overall_success = len(saved_pairs) > 0
    if not overall_success:
        print(f"  ✗ Could not save any pair for {model_file}")

    return {
        'model': model_file,
        'overall_success': overall_success,
        'solutions': solution_statuses,
        'error': None
    }

# -------------------------------------------------------------------
# STANDARDIZE ATOM ORDER
# -------------------------------------------------------------------
def standardize_atom_order(atoms1, atoms2):
    def standardize(atoms):
        atoms = remove_nonstandard_atoms(atoms, "standardization")
        residue_map = defaultdict(list)
        for _, _, info in atoms:
            key = (info['chain_id'], info['res_seq'])
            residue_map[key].append(info)
        sorted_residues = sorted(residue_map.keys(), key=lambda x: (x[0], x[1]))
        result = []
        for key in sorted_residues:
            atoms_list = residue_map[key]
            res_name = atoms_list[0]['res_name'].strip()
            nuc_type = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            expected = NUCLEOTIDE_ATOM_ORDER.get(nuc_type, [])
            atom_dict = {a['atom_name']: a for a in atoms_list}
            sorted_atoms = []
            if expected:
                for an in expected:
                    if an in atom_dict:
                        sorted_atoms.append(atom_dict[an])
                for a in atoms_list:
                    if a['atom_name'] not in expected:
                        sorted_atoms.append(a)
            else:
                sorted_atoms = atoms_list
            for info in sorted_atoms:
                result.append(info)
        final = []
        for i, info in enumerate(result, 1):
            info['atom_num'] = i
            new_line = format_pdb_atom_line(
                i,
                info['atom_name'],
                info['alt_loc'],
                info['res_name'],
                info['chain_id'],
                info['res_seq'],
                info['icode'],
                info['x'], info['y'], info['z'],
                info['occupancy'],
                info['temp_factor'],
                info['element']
            )
            info['line'] = new_line.rstrip('\n')
            info['original_line'] = new_line
            final.append((0, new_line, info))
        return final
    return standardize(atoms1), standardize(atoms2)

# -------------------------------------------------------------------
# INTERACTIVE EDITOR
# -------------------------------------------------------------------
def interactive_editor(solution_file, model_file, output_dir, recorder=None, is_first_model=True):
    print(f"\n{'='*80}")
    print("INTERACTIVE PDB EDITOR - RNA RESIDUE MANAGEMENT")
    print(f"Solution: {solution_file.name}")
    print(f"Model:    {model_file.name}")
    print("="*80)

    solution_lines = read_pdb_lines(solution_file)
    model_lines = read_pdb_lines(model_file)

    solution_atoms = get_atom_lines(solution_lines)
    model_atoms = get_atom_lines(model_lines)

    # Initial standardization
    solution_atoms, model_atoms = standardize_atom_order(solution_atoms, model_atoms)

    original_solution_atoms = solution_atoms.copy()
    original_model_atoms = model_atoms.copy()

    change_history = ChangeHistory()

    while True:
        display_comparison(solution_atoms, model_atoms, solution_file.name, model_file.name)
        print("\n" + "="*80)
        print("OPERATION MENU - RESIDUE MANAGEMENT:")
        print("  RESIDUE MANAGEMENT:")
        print("   1. Show residues in SOLUTION")
        print("   2. Show residues in MODEL")
        print("   3. Remove a residue (entire nucleotide) from SOLUTION")
        print("   4. Remove a residue (entire nucleotide) from MODEL")
        print("   5. Swap two residues in SOLUTION")
        print("   6. Swap two residues in MODEL")
        print("   7. Change nucleotide type in SOLUTION")
        print("   8. Change nucleotide type in MODEL")
        print("  ATOM OPERATIONS:")
        print("   9. Delete atom(s) from SOLUTION")
        print("  10. Delete atom(s) from MODEL")
        print("  11. Change residue number")
        print("  12. Automatically renumber residues")
        print("  13. Renumber all atoms sequentially")
        print("  14. Change chain ID")
        print("  TOOLS:")
        print("  15. Display nucleotide sequences")
        print("  16. Validate RNA structure")
        print("  17. Remove extra atoms (e.g., hydrogens)")
        print("  18. Revert to original files (after standardization)")
        print("  SYSTEM:")
        print("  19. Display change history")
        print("  20. Undo last operation")
        print("  21. Advanced nucleotide type change (with atoms)")
        print("  22. Automatic nucleotide correction (according to solution)")
        print("  23. Automatic nucleotide correction (according to model)")
        print("  24. Brute force: synchronize atoms and residues")
        print("  25. SAVE and continue")
        print("   0. Exit without saving")
        print("="*80)

        choice = input("\nSelect operation (0‑25): ").strip()

        if choice == '0':
            return False, recorder

        elif choice == '1':
            display_residues(solution_atoms, "SOLUTION")
        elif choice == '2':
            display_residues(model_atoms, "MODEL")

        elif choice == '3':
            display_residues(solution_atoms, "SOLUTION")
            chain_id = input("Enter chain identifier: ").strip() or " "
            display_residue_ranges(solution_atoms, chain_id, "SOLUTION")
            res_seq_input = input("Enter residue number(s) to delete: ").strip()
            if res_seq_input:
                change_history.record_change('remove_residue', 'solution', solution_atoms,
                                             f"Removed {chain_id}:{res_seq_input}")
                solution_atoms = remove_residue(solution_atoms, chain_id, res_seq_input)
                solution_atoms = renumber_atoms(solution_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('remove_residue', 'solution', chain_id, res_seq_input)

        elif choice == '4':
            display_residues(model_atoms, "MODEL")
            chain_id = input("Enter chain identifier: ").strip() or " "
            display_residue_ranges(model_atoms, chain_id, "MODEL")
            res_seq_input = input("Enter residue number(s) to delete: ").strip()
            if res_seq_input:
                change_history.record_change('remove_residue', 'model', model_atoms,
                                             f"Removed {chain_id}:{res_seq_input}")
                model_atoms = remove_residue(model_atoms, chain_id, res_seq_input)
                model_atoms = renumber_atoms(model_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('remove_residue', 'model', chain_id, res_seq_input)

        elif choice == '5':
            display_residues(solution_atoms, "SOLUTION")
            chain_id = input("Enter chain identifier: ").strip() or " "
            try:
                r1 = int(input("First residue: "))
                r2 = int(input("Second residue: "))
                change_history.record_change('swap_residues', 'solution', solution_atoms,
                                             f"Swapped {chain_id}:{r1} ↔ {chain_id}:{r2}")
                solution_atoms = swap_residues(solution_atoms, chain_id, r1, r2)
                if recorder and is_first_model:
                    recorder.add_operation('swap_residues', 'solution', chain_id, r1, r2)
            except:
                print("Invalid numbers.")

        elif choice == '6':
            display_residues(model_atoms, "MODEL")
            chain_id = input("Enter chain identifier: ").strip() or " "
            try:
                r1 = int(input("First residue: "))
                r2 = int(input("Second residue: "))
                change_history.record_change('swap_residues', 'model', model_atoms,
                                             f"Swapped {chain_id}:{r1} ↔ {chain_id}:{r2}")
                model_atoms = swap_residues(model_atoms, chain_id, r1, r2)
                if recorder and is_first_model:
                    recorder.add_operation('swap_residues', 'model', chain_id, r1, r2)
            except:
                print("Invalid numbers.")

        elif choice == '7':
            display_residues(solution_atoms, "SOLUTION")
            chain_id = input("Chain: ").strip() or " "
            try:
                res = int(input("Residue number: "))
                new_type = input("New type (A/C/G/U): ").strip().upper()
                if new_type in 'ACGU':
                    change_history.record_change('change_residue_type', 'solution', solution_atoms,
                                                 f"{chain_id}:{res} → {new_type}")
                    solution_atoms = change_residue_type(solution_atoms, chain_id, res, new_type)
                    if recorder and is_first_model:
                        recorder.add_operation('change_residue_type', 'solution', chain_id, res, new_type)
            except:
                print("Error.")

        elif choice == '8':
            display_residues(model_atoms, "MODEL")
            chain_id = input("Chain: ").strip() or " "
            try:
                res = int(input("Residue number: "))
                new_type = input("New type (A/C/G/U): ").strip().upper()
                if new_type in 'ACGU':
                    change_history.record_change('change_residue_type', 'model', model_atoms,
                                                 f"{chain_id}:{res} → {new_type}")
                    model_atoms = change_residue_type(model_atoms, chain_id, res, new_type)
                    if recorder and is_first_model:
                        recorder.add_operation('change_residue_type', 'model', chain_id, res, new_type)
            except:
                print("Error.")

        elif choice == '9':
            max_num = len(solution_atoms)
            if max_num == 0:
                print("No atoms.")
                continue
            sel = input("Atom numbers to delete from SOLUTION: ").strip()
            indices = parse_selection(sel, max_num)
            if indices:
                change_history.record_change('delete_atoms', 'solution', solution_atoms,
                                             f"Deleted atoms: {indices}")
                solution_atoms = delete_atoms(solution_atoms, indices)
                solution_atoms = renumber_atoms(solution_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('delete_atoms', 'solution', indices)

        elif choice == '10':
            max_num = len(model_atoms)
            if max_num == 0:
                print("No atoms.")
                continue
            sel = input("Atom numbers to delete from MODEL: ").strip()
            indices = parse_selection(sel, max_num)
            if indices:
                change_history.record_change('delete_atoms', 'model', model_atoms,
                                             f"Deleted atoms: {indices}")
                model_atoms = delete_atoms(model_atoms, indices)
                model_atoms = renumber_atoms(model_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('delete_atoms', 'model', indices)

        elif choice == '11':
            target = input("Change in SOLUTION (s) or MODEL (m)? ").lower()
            try:
                old = int(input("Old residue number: "))
                new = int(input("New residue number: "))
                if target == 's':
                    change_history.record_change('change_residue', 'solution', solution_atoms,
                                                 f"{old}→{new}")
                    solution_atoms = change_residue_number(solution_atoms, old, new)
                    if recorder and is_first_model:
                        recorder.add_operation('change_residue', 'solution', old, new)
                elif target == 'm':
                    change_history.record_change('change_residue', 'model', model_atoms,
                                                 f"{old}→{new}")
                    model_atoms = change_residue_number(model_atoms, old, new)
                    if recorder and is_first_model:
                        recorder.add_operation('change_residue', 'model', old, new)
            except:
                print("Error.")

        elif choice == '12':
            target = input("Renumber in SOLUTION (s) or MODEL (m)? ").lower()
            if target == 's':
                change_history.record_change('auto_renumber', 'solution', solution_atoms,
                                             "Auto renumber")
                solution_atoms = auto_renumber_residues(solution_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('auto_renumber', 'solution')
            elif target == 'm':
                change_history.record_change('auto_renumber', 'model', model_atoms,
                                             "Auto renumber")
                model_atoms = auto_renumber_residues(model_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('auto_renumber', 'model')

        elif choice == '13':
            change_history.record_change('renumber_atoms', 'solution', solution_atoms, "Renumber atoms")
            change_history.record_change('renumber_atoms', 'model', model_atoms, "Renumber atoms")
            solution_atoms = renumber_atoms(solution_atoms)
            model_atoms = renumber_atoms(model_atoms)
            if recorder and is_first_model:
                recorder.add_operation('renumber_atoms', 'solution')
                recorder.add_operation('renumber_atoms', 'model')

        elif choice == '14':
            target = input("Change chain ID in SOLUTION (s) or MODEL (m)? ").lower()
            old = input("Old chain: ").strip()
            new = input("New chain: ").strip()
            if target == 's':
                change_history.record_change('change_chain', 'solution', solution_atoms,
                                             f"{old}→{new}")
                solution_atoms = change_chain_id(solution_atoms, old, new)
                if recorder and is_first_model:
                    recorder.add_operation('change_chain', 'solution', old, new)
            elif target == 'm':
                change_history.record_change('change_chain', 'model', model_atoms,
                                             f"{old}→{new}")
                model_atoms = change_chain_id(model_atoms, old, new)
                if recorder and is_first_model:
                    recorder.add_operation('change_chain', 'model', old, new)

        elif choice == '15':
            display_sequences(solution_atoms, model_atoms, solution_file.name, model_file.name)

        elif choice == '16':
            print("\nValidating SOLUTION:")
            issues = validate_structure(solution_atoms)
            if issues:
                for iss in issues[:5]:
                    print(f"  - {iss}")
            else:
                print("  OK")
            print("\nValidating MODEL:")
            issues = validate_structure(model_atoms)
            if issues:
                for iss in issues[:5]:
                    print(f"  - {iss}")
            else:
                print("  OK")

        elif choice == '17':
            target = input("Remove extra atoms from SOLUTION (s) or MODEL (m)? ").lower()
            range_in = input("Residue range (e.g., 1-END): ").strip()
            if target == 's':
                change_history.record_change('remove_extra_atoms', 'solution', solution_atoms,
                                             f"Range {range_in}")
                solution_atoms = remove_extra_atoms_from_residues(solution_atoms, range_in, 'solution')
                solution_atoms = renumber_atoms(solution_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('remove_extra_atoms', 'solution', range_in)
            elif target == 'm':
                change_history.record_change('remove_extra_atoms', 'model', model_atoms,
                                             f"Range {range_in}")
                model_atoms = remove_extra_atoms_from_residues(model_atoms, range_in, 'model')
                model_atoms = renumber_atoms(model_atoms)
                if recorder and is_first_model:
                    recorder.add_operation('remove_extra_atoms', 'model', range_in)

        elif choice == '18':
            solution_atoms = original_solution_atoms.copy()
            model_atoms = original_model_atoms.copy()
            print("Reverted to original files.")
            change_history.clear_history()
            if recorder and is_first_model:
                recorder.clear()

        elif choice == '19':
            change_history.display_history()

        elif choice == '20':
            solution_atoms, model_atoms, _ = change_history.undo_last(solution_atoms, model_atoms)

        elif choice == '21':
            target = input("Advanced change in SOLUTION (s) or MODEL (m)? ").lower()
            chain_id = input("Chain: ").strip() or " "
            try:
                res = int(input("Residue number: "))
                new_type = input("New type (A/C/G/U): ").strip().upper()
                if new_type not in 'ACGU':
                    print("Invalid type.")
                    continue
                if target == 's':
                    change_history.record_change('advanced_change', 'solution', solution_atoms,
                                                 f"{chain_id}:{res} → {new_type}")
                    solution_atoms, removed = basic_advanced_change_residue_type(solution_atoms, chain_id, res, new_type)
                    if removed:
                        model_atoms = remove_matching_atoms(model_atoms, removed)
                        model_atoms = renumber_atoms(model_atoms)
                    if recorder and is_first_model:
                        recorder.add_operation('advanced_change', 'solution', chain_id, res, new_type)
                elif target == 'm':
                    change_history.record_change('advanced_change', 'model', model_atoms,
                                                 f"{chain_id}:{res} → {new_type}")
                    model_atoms, removed = basic_advanced_change_residue_type(model_atoms, chain_id, res, new_type)
                    if removed:
                        solution_atoms = remove_matching_atoms(solution_atoms, removed)
                        solution_atoms = renumber_atoms(solution_atoms)
                    if recorder and is_first_model:
                        recorder.add_operation('advanced_change', 'model', chain_id, res, new_type)
            except:
                print("Error.")

        elif choice == '22':
            change_history.record_change('auto_correct_nucleotides', 'solution', solution_atoms, "Auto correct")
            change_history.record_change('auto_correct_nucleotides', 'model', model_atoms, "Auto correct")
            solution_atoms, model_atoms = auto_correct_nucleotides_by_residue(solution_atoms, model_atoms)
            if recorder and is_first_model:
                recorder.add_operation('auto_correct_nucleotides', 'solution')
                recorder.add_operation('auto_correct_nucleotides', 'model')

        elif choice == '23':
            change_history.record_change('auto_correct_nucleotides_reverse', 'solution', solution_atoms, "Auto correct reverse")
            change_history.record_change('auto_correct_nucleotides_reverse', 'model', model_atoms, "Auto correct reverse")
            solution_atoms, model_atoms = auto_correct_nucleotides_by_residue_reverse(solution_atoms, model_atoms)
            if recorder and is_first_model:
                recorder.add_operation('auto_correct_nucleotides_reverse', 'solution')
                recorder.add_operation('auto_correct_nucleotides_reverse', 'model')

        elif choice == '24':
            change_history.record_change('brute_force_sync', 'solution', solution_atoms, "Brute force sync")
            change_history.record_change('brute_force_sync', 'model', model_atoms, "Brute force sync")
            solution_atoms, model_atoms = brute_force_sync_atoms(solution_atoms, model_atoms)
            if recorder and is_first_model:
                recorder.add_operation('brute_force_sync', 'solution')
                recorder.add_operation('brute_force_sync', 'model')

        elif choice == '25':
            break

        else:
            print("Invalid choice.")

    # Save files after editing
    print("\nSaving files...")
    updated_solution_lines = update_full_pdb(solution_lines, solution_atoms)
    updated_model_lines = update_full_pdb(model_lines, model_atoms)

    # Add REMARK lines
    if recorder and is_first_model:
        sol_remarks = generate_remarks(recorder.operations, 'solution')
        mod_remarks = generate_remarks(recorder.operations, 'model')
        updated_solution_lines = sol_remarks + updated_solution_lines
        updated_model_lines = mod_remarks + updated_model_lines

    # Construct output filenames with solution identifier
    model_stem = model_file.stem
    sol_stem = solution_file.stem
    sol_short = sol_stem.replace('solution', '').replace('SOLUTION', '').strip('_')
    if not sol_short:
        sol_short = sol_stem
    sol_output = output_dir / f"{model_stem}_vs_{sol_short}_solution.pdb"
    mod_output = output_dir / f"{model_stem}_vs_{sol_short}_refined.pdb"

    output_dir.mkdir(parents=True, exist_ok=True)
    write_pdb_lines(sol_output, updated_solution_lines)
    write_pdb_lines(mod_output, updated_model_lines)

    print(f"Saved: {sol_output.name} ({len(solution_atoms)} atoms)")
    print(f"Saved: {mod_output.name} ({len(model_atoms)} atoms)")
    print(f"\nFull paths:\n  Solution: {sol_output.resolve()}\n  Model:    {mod_output.resolve()}")

    return True, recorder

def display_sequences(solution_atoms, model_atoms, sol_name, mod_name):
    """Show nucleotide sequences of both structures."""
    def extract(atoms):
        seq = []
        nums = []
        residues = {}
        for _, _, info in atoms:
            key = (info['chain_id'], info['res_seq'])
            if key not in residues:
                res_name = info['res_name'].strip()
                nuc = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
                residues[key] = (nuc, info['res_seq'])
        sorted_res = sorted(residues.values(), key=lambda x: x[1])
        for nuc, num in sorted_res:
            if nuc in 'ACGU':
                seq.append(nuc)
                nums.append(f"{nuc}{num}")
        return ''.join(seq), nums
    sol_seq, sol_nums = extract(solution_atoms)
    mod_seq, mod_nums = extract(model_atoms)
    print("\n" + "="*120)
    print("NUCLEOTIDE SEQUENCES")
    print("="*120)
    print(f"{sol_name}: {sol_seq}")
    print(f"{mod_name}:    {mod_seq}")
    print("\nNUMBERING:")
    max_len = max(len(sol_nums), len(mod_nums))
    for i in range(max_len):
        s = sol_nums[i] if i < len(sol_nums) else ""
        m = mod_nums[i] if i < len(mod_nums) else ""
        diff = " |" if i < len(sol_nums) and i < len(mod_nums) and sol_nums[i] != mod_nums[i] else "  "
        print(f"{s:>60s}{diff}{m:<60s}")
    print("="*120)

def validate_structure(atoms):
    """Check each residue for missing or extra atoms."""
    residue_atoms = defaultdict(list)
    issues = []
    for _, _, info in atoms:
        key = (info['chain_id'], info['res_seq'], info['res_name'])
        residue_atoms[key].append(info['atom_name'])
    for (chain, seq, res_name), atoms_list in residue_atoms.items():
        nuc = RESIDUE_NAME_MAP.get(res_name.strip(), res_name[0] if res_name else '?')
        expected = set(NUCLEOTIDE_ATOM_ORDER.get(nuc, []))
        actual = set(atoms_list)
        missing = expected - actual
        extra = actual - expected
        if missing:
            issues.append(f"{chain}:{seq} {nuc} missing {sorted(missing)}")
        if extra:
            issues.append(f"{chain}:{seq} {nuc} extra {sorted(extra)}")
    return issues

def remove_extra_atoms_from_residues(atom_lines, residue_range, target_name):
    """Remove atoms not in the expected list for given residues."""
    if isinstance(residue_range, str):
        if 'END' in residue_range.upper():
            try:
                start = int(residue_range.split('-')[0].strip())
                max_res = max(info['res_seq'] for _, _, info in atom_lines)
                residue_range = (start, max_res)
            except:
                print("Invalid range.")
                return atom_lines
    start, end = residue_range
    print(f"Removing extra atoms from {target_name} for residues {start}-{end}")
    atoms_to_keep = []
    removed = 0
    for line_num, line, info in atom_lines:
        if start <= info['res_seq'] <= end:
            res_name = info['res_name'].strip()
            nuc = RESIDUE_NAME_MAP.get(res_name, res_name[0] if res_name else '?')
            expected = NUCLEOTIDE_ATOM_ORDER.get(nuc, [])
            if expected and info['atom_name'] not in expected:
                removed += 1
                continue
        atoms_to_keep.append((line_num, line, info))
    print(f"Removed {removed} atoms.")
    return atoms_to_keep

def display_residue_ranges(atom_lines, chain_id, target_name):
    """Show residue ranges for a given chain."""
    residues = sorted(set(info['res_seq'] for _, _, info in atom_lines if info['chain_id'] == chain_id))
    if not residues:
        print(f"No residues in chain {chain_id}")
        return
    ranges = []
    start = residues[0]
    end = residues[0]
    for r in residues[1:]:
        if r == end + 1:
            end = r
        else:
            ranges.append((start, end))
            start = r
            end = r
    ranges.append((start, end))
    print(f"Ranges in {target_name} chain {chain_id}:")
    for s, e in ranges:
        if s == e:
            print(f"  {s}")
        else:
            print(f"  {s}-{e}")

# -------------------------------------------------------------------
# MAIN PROCESSING MENU
# -------------------------------------------------------------------
def main_processing_menu(models, puzzle_path, output_dir, results=None):
    if results is None:
        results = []
    recorder = OperationRecorder()
    processed_models = set()

    all_models = collect_model_files_from_external(puzzle_path)
    all_solutions = collect_solution_files(puzzle_path)

    while True:
        print("\n" + "="*80)
        print("MAIN PROCESSING MENU")
        print(f"Selected {len(models)} model(s) for editing. Recorded operations: {len(recorder.get_operations())}")
        print(f"Total models in directory: {len(all_models)} (including {len(all_solutions)} solution files)")
        print("="*80)
        print("1. Edit a model (record operations) – after saving, automatically select best solution by RMSD")
        print("2. Apply recorded operations to other models (automatic solution selection by RMSD)")
        print("3. Clear recorded operations")
        print("4. Display recorded operations")
        print("5. Exit this puzzle set")
        print("="*80)

        choice = input("\nChoose option (1-5): ").strip()

        if choice == '1':
            unprocessed = [m for m in models if m not in processed_models]
            if not unprocessed:
                print("All selected models already processed!")
                continue

            selected = select_models(unprocessed, "Select model to edit")
            if not selected:
                continue
            selected_model = selected[0]

            if not all_solutions:
                print("No solution files! Cannot proceed.")
                continue
            print("\nAvailable solution files:")
            for i, sol in enumerate(all_solutions, 1):
                print(f"{i}. {sol}")
            try:
                idx = int(input("Choose solution file number (or 0 to cancel): ").strip())
                if idx == 0:
                    continue
                if 1 <= idx <= len(all_solutions):
                    chosen_solution = all_solutions[idx-1]
                else:
                    print("Invalid number.")
                    continue
            except ValueError:
                print("Invalid input.")
                continue

            success, recorder = interactive_editor(
                puzzle_path / chosen_solution,
                puzzle_path / selected_model,
                output_dir,
                recorder,
                is_first_model=True
            )

            if success:
                print(f"✓ Model {selected_model} saved after editing.")

                # Determine which solutions are already successful for this model
                report_path = output_dir / "report.csv"
                successful_solutions = get_successful_solutions(selected_model, report_path)

                # Ask user what to compare after editing
                print("\nWhat would you like to compare after editing?")
                print("1. Only the solution I just edited")
                print("2. All solutions")
                print("3. Select specific solutions")
                print("4. Process only solutions not yet successful (the rest)")
                choice_comp = input("Choose (1-4) [default=1]: ").strip()

                if choice_comp == '2':
                    solutions_filter = None
                elif choice_comp == '3':
                    solutions_filter = select_models(all_solutions, "Select solutions to process")
                elif choice_comp == '4':
                    # Process solutions that are not yet successful
                    pending = [sol for sol in all_solutions if sol not in successful_solutions]
                    if not pending:
                        print("All solutions already successful. Nothing to process.")
                        solutions_filter = []   # will skip
                    else:
                        solutions_filter = pending
                        print(f"Will process {len(pending)} solution(s): {', '.join(pending)}")
                else:   # default to option 1
                    solutions_filter = [chosen_solution]

                if solutions_filter is not None and not solutions_filter:
                    print("No solutions selected. Skipping post‑edit processing.")
                else:
                    print("Now comparing with selected solution(s)...")
                    result = process_model_with_all_solutions(
                        selected_model, puzzle_path, output_dir, recorder.get_operations(),
                        solutions_filter=solutions_filter
                    )
                    if result:
                        results.append(result)
                        processed_models.add(selected_model)
                    else:
                        print(f"✗ Could not process {selected_model}")
            else:
                print(f"✗ Edit of {selected_model} cancelled.")

        elif choice == '2':
            if not recorder.get_operations():
                print("No recorded operations!")
                continue

            print(f"\nRecorded operations: {len(recorder.get_operations())}")
            recorder.display_operations()

            available = [m for m in all_models if m not in processed_models]
            if not available:
                print("All models have already been processed.")
                continue

            selected_for_apply = select_models(available, "Select models to apply operations to")
            if not selected_for_apply:
                continue

            # Ask for solution selection once (global for all models)
            print("\nWhich solutions do you want to process?")
            print("1. All solutions")
            print("2. Select specific solutions (same for all models)")
            print("3. Only solutions not yet successful for each model (based on existing report)")
            sol_choice = input("Choose (1-3) [default=1]: ").strip()

            # Prepare global filter if needed
            global_sol_filter = None
            use_pending = False
            if sol_choice == '2':
                # Let user select solutions once
                global_sol_filter = select_models(all_solutions, "Select solutions to apply to all selected models")
                if not global_sol_filter:
                    print("No solutions selected. Exiting batch processing.")
                    continue
            elif sol_choice == '3':
                use_pending = True  # compute per model later
            else:
                # default to all solutions
                pass

            print(f"\nApplying operations to {len(selected_for_apply)} models...")
            for model_file in selected_for_apply:
                # Determine solutions filter for this model
                if sol_choice == '2':
                    sol_filter = global_sol_filter
                elif use_pending:
                    # Get successful solutions from report (if any)
                    report_path = output_dir / "report.csv"
                    successful = get_successful_solutions(model_file, report_path)
                    pending = [sol for sol in all_solutions if sol not in successful]
                    if not pending:
                        print(f"All solutions already successful for {model_file}. Skipping.")
                        continue
                    sol_filter = pending
                    print(f"Will process {len(pending)} pending solution(s) for {model_file}: {', '.join(pending)}")
                else:
                    sol_filter = None   # all solutions

                result = process_model_with_all_solutions(
                    model_file, puzzle_path, output_dir, recorder.get_operations(),
                    solutions_filter=sol_filter
                )
                if result:
                    results.append(result)
                    processed_models.add(model_file)

        elif choice == '3':
            confirm = input("Clear recorded operations? (y/n): ").strip().lower()
            if confirm == 'y':
                recorder.clear()
                print("Cleared.")

        elif choice == '4':
            recorder.display_operations()

        elif choice == '5':
            print("Exiting processing menu.")
            break

        else:
            print("Invalid choice.")

    return results

# -------------------------------------------------------------------
# MAIN FUNCTION
# -------------------------------------------------------------------
def main():
    print("RNA Puzzle - Interactive PDB Editor")
    print("=" * 50)
    print("Compare and edit solution/model files")
    print("Operations are recorded and can be applied to other models")
    print("AUTOMATIC ATOM ORDER STANDARDIZATION")
    print("AUTOMATIC SELECTION OF BEST SOLUTION BY RMSD")
    print(f"External directory: {BASE_EXTERNAL_DIR}")
    print(f"Output directory:   {BASE_PROCESSED_DIR}")

    if not BASE_EXTERNAL_DIR.exists():
        print(f"ERROR: Directory '{BASE_EXTERNAL_DIR}' does not exist!")
        sys.exit(1)

    repositories = [d for d in os.listdir(BASE_EXTERNAL_DIR) if (BASE_EXTERNAL_DIR / d).is_dir()]
    if not repositories:
        print("No repositories found.")
        sys.exit(1)

    print("\nAvailable repositories:")
    for i, repo in enumerate(repositories, 1):
        print(f"{i}. {repo}")
    selected_repos_input = input("\nSelect repositories to process (e.g., 1, 1-3, 1,3,5): ").strip()
    selected_repos = []
    for part in selected_repos_input.replace(' ', '').split(','):
        if '-' in part:
            try:
                s, e = map(int, part.split('-'))
                selected_repos.extend(repositories[s-1:e])
            except:
                print(f"Invalid range: {part}")
        else:
            try:
                idx = int(part)
                if 1 <= idx <= len(repositories):
                    selected_repos.append(repositories[idx-1])
            except:
                print(f"Invalid number: {part}")
    if not selected_repos:
        print("No repositories selected.")
        sys.exit(0)

    print("\nSelected:")
    for r in selected_repos:
        print(f"  - {r}")

    for repo_name in selected_repos:
        repo_path = BASE_EXTERNAL_DIR / repo_name
        puzzles_path = repo_path / "data"

        if not puzzles_path.exists():
            print(f"No 'data' directory in {repo_name}. Skipping.")
            continue

        puzzle_sets = [d for d in os.listdir(puzzles_path) if (puzzles_path / d).is_dir()]
        if not puzzle_sets:
            print(f"No puzzle sets in {repo_name}.")
            continue

        print(f"\nRepository: {repo_name}")
        print("\nAvailable puzzle sets:")
        for i, ps in enumerate(puzzle_sets, 1):
            print(f"{i}. {ps}")
        selected_sets_input = input(f"\nSelect puzzle sets in '{repo_name}' (e.g., 1, 1-3): ").strip()
        selected_sets = []
        for part in selected_sets_input.replace(' ', '').split(','):
            if '-' in part:
                try:
                    s, e = map(int, part.split('-'))
                    selected_sets.extend(puzzle_sets[s-1:e])
                except:
                    print(f"Invalid range: {part}")
            else:
                try:
                    idx = int(part)
                    if 1 <= idx <= len(puzzle_sets):
                        selected_sets.append(puzzle_sets[idx-1])
                except:
                    print(f"Invalid number: {part}")
        if not selected_sets:
            print("No puzzle sets selected.")
            continue

        print(f"\nProcessing {len(selected_sets)} puzzle set(s) in {repo_name}")

        for puzzle_set in selected_sets:
            puzzle_path = puzzles_path / puzzle_set
            print(f"\n{'='*60}")
            print(f"Puzzle set: {puzzle_set}")
            print('-' * 60)

            success = process_puzzle_set_interactive(puzzle_path, repo_name)

            if success:
                print(f"\n  [FINISHED] Processed {puzzle_set}")
            else:
                print(f"\n  [FAILURE] Could not process {puzzle_set}")

    print("\n" + "=" * 60)
    print("Processing finished!")
    print("=" * 60)

def process_puzzle_set_interactive(puzzle_path, repo_name):
    pdb_path = puzzle_path / "pdb"
    if not pdb_path.exists():
        print(f"Missing 'pdb' directory in: {puzzle_path}")
        return False

    puzzle_name = pdb_path.parent.name
    print(f"\nProcessing puzzle: {puzzle_name}")
    print(f"Path: {pdb_path}")

    models = collect_model_files_from_external(pdb_path)
    if not models:
        print(f"No models found in {pdb_path}")
        return False

    print(f"Found {len(models)} model(s)")

    try:
        rel_path = pdb_path.relative_to(BASE_EXTERNAL_DIR)
        output_dir = BASE_PROCESSED_DIR / rel_path.parent / "pdb"
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {output_dir}")
    except ValueError as e:
        print(f"Error creating output path: {e}")
        return False

    print("\n" + "="*80)
    print("MODEL SELECTION")
    print("="*80)
    selected_models = select_models(models, "Select models to process")
    if not selected_models:
        print("No models selected.")
        return False

    print(f"\nSelected {len(selected_models)} model(s).")
    results = main_processing_menu(selected_models, pdb_path, output_dir)

    # Write CSV report
    write_report_csv(results, output_dir)
    return True

# -------------------------------------------------------------------
# ENTRY POINT
# -------------------------------------------------------------------
if __name__ == "__main__":
    main()