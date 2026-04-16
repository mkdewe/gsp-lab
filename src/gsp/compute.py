#!/usr/bin/python

import sys, getopt
import os
from biopandas.pdb import PandasPdb
from biopandas.mmcif import PandasMmcif
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
from bisect import bisect
import csv
from os import listdir
from os.path import isfile, join
import pathlib
from biopandas.pdb.engines import amino3to1dict
from copy import deepcopy

na2to1dict = {
    "DA": "A",
    "DG": "G",
    "DU": "U",
    "DT": "T",
    "DC": "C",
    "A": "A",
    "G": "G",
    "U": "U",
    "T": "T",
    "C": "C"
}

# https://hunterheidenreich.com/posts/kabsch_algorithm/
def superimpose(P, Q):
    assert P.shape == Q.shape, "Matrix dimensions must match"
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    t = centroid_Q - centroid_P
    p = P - centroid_P
    q = Q - centroid_Q
    H = np.dot(p.T, q)
    U, S, Vt = np.linalg.svd(H)
    if np.linalg.det(np.dot(Vt.T, U.T)) < 0.0:
        Vt[-1, :] *= -1.0
    R = np.dot(Vt.T, U.T)
    rmsd = np.sqrt(np.sum(np.square(np.dot(p, R.T) - q)) / P.shape[0])
    return R, t, rmsd
    
def compute_rmsd(tdf, mdf, istpdb, ismpdb):
    target_coords = get_coords(tdf, istpdb)
    model_coords = get_coords(mdf, ismpdb)
    R_opt, t_opt, rmsd = superimpose(target_coords, model_coords)
    return round(rmsd,3)

def floatify(s):
    try:
        return float(s)
    except ValueError:
        return None

def read_config(argv):
    try:
        # Dodano recompute_from_csv do listy długich opcji
        opts, args = getopt.getopt(argv,"sht:m:a:r:d:",["save_structures","help","target_path=","model_path=","central_atoms=","sphere_radii=","rmsd_threshold=","skip_constant_radii","radii_tolerance=","recompute_from_csv="])
    except getopt.GetoptError as error:
        print('{}'.format(error))
        sys.exit(1)
    
    target_path = None
    model_path = None
    recompute_path = None # Nowa zmienna
    central_atoms = "CA,C1'".split(',')
    sphere_radii = '2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0'.split(',')
    sphere_radii[:] = filter(float.__instancecheck__, map(floatify, sphere_radii))
    rmsd_threshold = '5.0'
    save_structures = False
    skip_constant_radii = False
    radii_tolerance = 1e-6

    for opt, arg in opts:
        if opt == '-h':
            print('...')
            sys.exit()
        elif opt == "--recompute_from_csv": 
            recompute_path = arg
        elif opt in ("-t", "--target_path"):
            target_path = arg
        elif opt in ("-m", "--model_path"):
            model_path = arg
        elif opt in ("-a", "--central_atoms"):
            central_atoms = [atom.strip() for atom in arg.split(',')]
        elif opt in ("-r", "--sphere_radii"):
            sphere_radii = [radius.strip() for radius in arg.split(',')]
            sphere_radii[:] = filter(float.__instancecheck__, map(floatify, sphere_radii))
            sphere_radii.sort()
        elif opt in ("-d", "--rmsd_threshold"):
            rmsd_threshold = arg
        elif opt in ("-s", "--save_structures"):
            save_structures = True
        elif opt == "--skip_constant_radii":
            skip_constant_radii = True
        elif opt == "--radii_tolerance":
            radii_tolerance = floatify(arg)

    # Walidacja: jeśli nie recompute, to wymagane -t i -m
    if not recompute_path:
        if target_path is None or model_path is None:
            print('Error: -t and -m are required unless --recompute_from_csv is used.')
            sys.exit(1)

    threshold_list = [float(t.strip()) for t in rmsd_threshold.split(',')]
    return target_path, model_path, central_atoms, sphere_radii, threshold_list, save_structures, skip_constant_radii, radii_tolerance, recompute_path
    
def read_structure(path):
    try:
        return PandasPdb().read_pdb(path).label_models()
    except ValueError:
        return PandasMmcif().read_mmcif(path)
        
def get_unique_column_values(structure, column_name):
    return structure[column_name].unique()
        
def get_list_of_models(structure):
    if isinstance(structure,PandasPdb):
        return get_unique_column_values(structure.df['ATOM'], 'model_id')
    elif isinstance(structure,PandasMmcif):
        return get_unique_column_values(structure.df['ATOM'], 'pdbx_PDB_model_num')
    else:
        return None
        
def get_central_atoms(df, cas, ispdb):
    if ispdb:
        return df[df['atom_name'].isin(cas)]
    else:
        return df[df['auth_atom_id'].isin(cas)]

def get_coord_names(ispdb):
    if ispdb:
        return ['x_coord', 'y_coord', 'z_coord']
    else:
        return ['Cartn_x', 'Cartn_y', 'Cartn_z']
        
def get_coords(df, ispdb):
    coord_names = get_coord_names(ispdb)
    return df[coord_names].to_numpy()

def get_dist_matrix(central_atoms_df, all_atoms_df, ispdb):
    coord_names = get_coord_names(ispdb)
    central_atoms_coords = get_coords(central_atoms_df, ispdb)
    all_atoms_coords = get_coords(all_atoms_df, ispdb)
    return np.sqrt(((central_atoms_coords[:, None] - all_atoms_coords) ** 2).sum(axis=2))
    
def get_atom_id(atom, ispdb):
    if ispdb:
        return get_residue_id(atom,ispdb) + '_' + atom['atom_name']
    else:
        return get_residue_id(atom,ispdb) + '_' + atom['auth_atom_id']
        
def get_residue_id(atom, ispdb):
    if ispdb:
        ins_code = atom['insertion'] if atom['insertion'] != 'None' else ''
        return atom['chain_id'] + str(atom['residue_number']) + ins_code
    else:
        ins_code = atom['pdbx_PDB_ins_code'] if atom['pdbx_PDB_ins_code'] != 'None' else ''
        return atom['auth_asym_id'] + str(atom['auth_seq_id']) + ins_code
        
def get_model_atoms_dict(df, ispdb):
    return {get_atom_id(df.iloc[i],ispdb): i for i in range(len(df))}
    
def get_target_atoms_dict(df, ispdb):
    return {i: get_atom_id(df.iloc[i],ispdb) for i in range(len(df))}

def get_sequence(df, ispdb, selected_central_atom_idxs):
    if ispdb:
        seq_list = extract_sequence_list(df, 'residue_number', 'residue_name', 'insertion')
    else:
        seq_list = extract_sequence_list(df, 'auth_seq_id', 'auth_comp_id', 'pdbx_PDB_ins_code')
    return ''.join([seq_list.iloc[i] for i in selected_central_atom_idxs])
    
def extract_sequence_list(df, residue_no_col_name, residue_name_col_name, icode_col_name):
    cmp = "placeholder"
    indices = []
    residue_number_insertion = df[residue_no_col_name].astype(str) + df[icode_col_name]
    for num, ind in zip(residue_number_insertion, np.arange(df.shape[0])):
        if num != cmp:
            indices.append(ind)
        cmp = num
    dicts = {**amino3to1dict, **na2to1dict}
    return df.iloc[indices][residue_name_col_name].map(dicts)
    
def compute_spheres(sphere_radii, target_central_atoms_df, target_df, istpdb):
    dist_matrix = get_dist_matrix(target_central_atoms_df, target_df,istpdb)
    spheres = {}
    for radius in sphere_radii:
        spheres[radius] = {}
        for i in range(len(target_central_atoms_df)):
            spheres[radius][i] = []
    for i in range(len(target_central_atoms_df)):
        for j in range(len(target_df)):
            current_distance = dist_matrix[i][j]
            selected_radius = sphere_radii[bisect(sphere_radii,current_distance,0,len(sphere_radii)-1)]
            spheres[selected_radius][i].append(j)
    return spheres
    
def save_sphere(struct, ispdb, df, output_file_path):
    if not os.path.isfile(output_file_path):
        local_struct_copy = deepcopy(struct)
        if not ispdb:
            local_struct_copy = local_struct_copy.convert_to_pandas_pdb()
        local_struct_copy.df['ATOM'] = df
        local_struct_copy.to_pdb(path=output_file_path, records=None, gz=False, append_newline=False)

def compute_scores(target_central_atoms_df, sphere_radii, spheres, target_df, current_model_df, istpdb, ismpdb, target_str, target_path, model_str, model_path, save_structures):
    target_atoms_dict = get_target_atoms_dict(target_df, istpdb)
    model_atoms_dict = get_model_atoms_dict(current_model_df, ismpdb)
    scores = {}
    for i in range(len(target_central_atoms_df)):
        scores[i] = []
        for r in range(len(sphere_radii)):
            target_atom_idxs_in_sphere = [x for a in range(len(sphere_radii)) if a <= r for x in spheres[sphere_radii[a]][i]]
            target_atom_idxs_in_sphere.sort()
            selected_target_atoms_df = target_df.loc[target_atom_idxs_in_sphere]
            if save_structures:
                target_file_path = '{}_{}_{}.pdb'.format(target_path, get_atom_id(target_central_atoms_df.iloc[i],istpdb), sphere_radii[r])            
                save_sphere(target_str, istpdb, selected_target_atoms_df, target_file_path)
            model_atom_idxs_in_sphere = [model_atoms_dict[target_atoms_dict[b]] for b in target_atom_idxs_in_sphere if target_atoms_dict[b] in model_atoms_dict]
            selected_model_atoms_df = current_model_df.loc[model_atom_idxs_in_sphere]
            if save_structures:
                model_file_path = '{}_{}_{}.pdb'.format(model_path, get_atom_id(target_central_atoms_df.iloc[i],istpdb), sphere_radii[r])            
                save_sphere(model_str, ismpdb, selected_model_atoms_df, model_file_path)
            if (len(selected_target_atoms_df) == len(selected_model_atoms_df)):
                scores[i].append(compute_rmsd(selected_target_atoms_df, selected_model_atoms_df, istpdb, ismpdb))
            else:
                print('Inconsistent sphere atom sets between target and model built around atom {} (radius {})!'.format(get_atom_id(target_central_atoms_df.iloc[i], istpdb), sphere_radii[r]))
                sys.exit(1)
    return scores
    
def contains_atom(df, ca, ispdb):
    if ispdb:
        return ca in df['atom_name'].values
    else:
        return ca in df['auth_atom_id'].values
        
def get_indices(df, ca, ispdb):
    if ispdb:
        return [i for i in range(len(df)) if df.iloc[i]['atom_name']==ca]
    else:
        return [i for i in range(len(df)) if df.iloc[i]['auth_atom_id']==ca]
        
def save_csv(result_file_path,data):
    with open(result_file_path, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)
        
def compute_ggsp(radii, selected_scores, rmsd_threshold):
    """
    Compute global GSP score using filtered radii and selected_scores.
    """
    ggsp = 0
    for r in range(len(radii)):
        accept_count = 0
        all_count = 0
        for curr_score in selected_scores[r]:
            if not (curr_score > rmsd_threshold):
                accept_count += 1
            all_count += 1
        ggsp += radii[r] * accept_count / all_count
    ggsp /= sum(radii)
    return round(ggsp * 100.0,3)
    
def compute_whole_matrix(radii, selected_scores, target_df, selected_central_atom_idxs, target_central_atoms_df, istpdb):
    """
    Builds the details matrix from already filtered radii and selected_scores.
    """
    target_sequence = get_sequence(target_central_atoms_df, istpdb, selected_central_atom_idxs)
    matrix = []
    row = ['']
    for i in range(len(selected_central_atom_idxs)):
        row.append('"' + get_residue_id(target_central_atoms_df.iloc[selected_central_atom_idxs[i]],istpdb) + '"')
    matrix.append(row)
    matrix.append(['radius'] + [*target_sequence] + ['mean','std dev'])
    
    for r in range(len(radii)):
        row = [str(radii[r])]
        rows_vect = selected_scores[r]
        for val in rows_vect:
            row.append(str(val))
        row.append(str(round(np.mean(rows_vect),3)))
        row.append(str(round(np.std(rows_vect),3)))
        matrix.append(row)
    
    # Column statistics
    col_means_row = ['mean']
    col_std_row = ['std dev']
    for col in range(len(selected_central_atom_idxs)):
        cols_vect = [selected_scores[row][col] for row in range(len(radii))]
        col_means_row.append(str(round(np.mean(cols_vect),3)))
        col_std_row.append(str(round(np.std(cols_vect),3)))
    matrix.append(col_means_row)
    matrix.append(col_std_row)
    return matrix
    
def compute_residue_scores(radii, selected_scores, target_central_atoms_df, selected_central_atom_idxs, rmsd_threshold, istpdb):
    """
    Compute per‑residue GSP scores using filtered radii and selected_scores.
    """
    residue_scores = [['residue','rGSP']]
    for i in range(len(selected_central_atom_idxs)):
        selected_idx = selected_central_atom_idxs[i]
        accept_count = 0
        all_count = 0
        rgsp = 0
        for r in range(len(radii)):
            curr_score = selected_scores[r][i]
            if not (curr_score > rmsd_threshold):
                accept_count += 1
            all_count += 1
            rgsp += radii[r] * accept_count / all_count
        rgsp /= sum(radii)
        row = []
        row.append('"' + get_residue_id(target_central_atoms_df.iloc[selected_central_atom_idxs[i]],istpdb) + '"')
        row.append(str(round(rgsp * 100.0,3)))
        residue_scores.append(row)
    return residue_scores
    
def update_local_scores(local_scores, ca, model_path, residue_scores):
    local_scores[ca][0].append(os.path.splitext(os.path.basename(model_path))[0])
    if len(local_scores[ca])==1:
        for i in range(1,len(residue_scores)):
            local_scores[ca].append(residue_scores[i])
    else:
        for i in range(1,len(local_scores[ca])):
            local_scores[ca][i].append(residue_scores[i][1])
    
def process_model(target_path, model_path, sphere_radii, central_atoms, thresholds_with_prefix, target_df, target_central_atoms_df, spheres, istpdb, per_threshold_accumulators, target_str, save_structures, skip_constant_radii, radii_tolerance):
    target_filename_without_ext = os.path.splitext(os.path.basename(target_path))[0]
    model_filename_without_ext = os.path.splitext(os.path.basename(model_path))[0]
    print('{} vs. {}'.format(target_filename_without_ext, model_filename_without_ext))
    model_str = read_structure(model_path)
    ismpdb = isinstance(model_str,PandasPdb)
    model_model_numbers = get_list_of_models(model_str)
    for model_no in model_model_numbers:
        current_model_str = model_str.get_model(model_no)
        current_model_df = current_model_str.df['ATOM']
        if (len(model_model_numbers)>1):
            scores = compute_scores(target_central_atoms_df, sphere_radii, spheres, target_df, current_model_df, istpdb, ismpdb, target_str, os.path.join(os.path.dirname(target_path), target_filename_without_ext), current_model_str, os.path.join(os.path.dirname(model_path), '{}_{}'.format(model_filename_without_ext,model_no)), save_structures)
        else:
            scores = compute_scores(target_central_atoms_df, sphere_radii, spheres, target_df, current_model_df, istpdb, ismpdb, target_str, os.path.join(os.path.dirname(target_path), target_filename_without_ext), current_model_str, os.path.join(os.path.dirname(model_path), model_filename_without_ext), save_structures)
        
        common_filename = '{}_{}'.format(target_filename_without_ext,model_filename_without_ext)
        if (len(model_model_numbers)>1):
            common_filename = '{}_{}'.format(common_filename,model_no)
        for ca in central_atoms:
            if (contains_atom(target_central_atoms_df,ca,istpdb)):
                selected_central_atom_idxs = get_indices(target_central_atoms_df, ca, istpdb)
                
                # ---- Build full selected_scores matrix (rows = radii, cols = residues) ----
                full_selected_scores = []
                for r in range(len(sphere_radii)):
                    row_vals = [scores[i][r] for i in selected_central_atom_idxs]
                    full_selected_scores.append(row_vals)
                
                # ---- Determine which radii to keep (skip constant radii) ----
                if skip_constant_radii:
                    radii_mask = []
                    for r in range(len(sphere_radii)):
                        row = full_selected_scores[r]
                        if max(row) - min(row) < radii_tolerance:
                            radii_mask.append(False)
                        else:
                            radii_mask.append(True)
                    filtered_radii = [sphere_radii[i] for i, keep in enumerate(radii_mask) if keep]
                    filtered_selected_scores = [row for row, keep in zip(full_selected_scores, radii_mask) if keep]
                    removed = sum(1 for keep in radii_mask if not keep)
                    if removed > 0:
                        print(f"    [INFO] Removed {removed} constant radii (tolerance={radii_tolerance}) for CA {ca}")
                else:
                    filtered_radii = sphere_radii
                    filtered_selected_scores = full_selected_scores
                
                # Build details matrix once (using filtered data)
                matrix = compute_whole_matrix(filtered_radii, filtered_selected_scores, target_df, selected_central_atom_idxs, target_central_atoms_df, istpdb)
                
                # Now for each threshold, write outputs using filtered data
                for thresh, prefix in thresholds_with_prefix:
                    if prefix:
                        details_file = os.path.join(os.path.dirname(model_path), f"{prefix}{common_filename}_{ca}-details.csv")
                        rgsp_file = os.path.join(os.path.dirname(model_path), f"{prefix}{common_filename}_{ca}-rGSP.csv")
                        ggsp_file = os.path.join(os.path.dirname(model_path), f"{prefix}{common_filename}_{ca}-gGSP.csv")
                    else:
                        details_file = os.path.join(os.path.dirname(model_path), f"{common_filename}_{ca}-details.csv")
                        rgsp_file = os.path.join(os.path.dirname(model_path), f"{common_filename}_{ca}-rGSP.csv")
                        ggsp_file = os.path.join(os.path.dirname(model_path), f"{common_filename}_{ca}-gGSP.csv")
                    
                    save_csv(details_file, matrix)
                    
                    residue_scores = compute_residue_scores(filtered_radii, filtered_selected_scores, target_central_atoms_df, selected_central_atom_idxs, thresh, istpdb)
                    save_csv(rgsp_file, residue_scores)
                    
                    ggsp = compute_ggsp(filtered_radii, filtered_selected_scores, thresh)
                    save_csv(ggsp_file, [['gGSP'],[str(ggsp)]])
                    
                    # Update per-threshold accumulators (for summary files)
                    per_threshold_accumulators[thresh][ca]['global'].append([model_filename_without_ext, str(ggsp)])
                    if ca not in per_threshold_accumulators[thresh][ca]['local']:
                        per_threshold_accumulators[thresh][ca]['local'][ca] = [['residue']]
                    local_acc = per_threshold_accumulators[thresh][ca]['local'][ca]
                    local_acc[0].append(os.path.splitext(os.path.basename(model_path))[0])
                    if len(local_acc) == 1:
                        for i in range(1, len(residue_scores)):
                            local_acc.append(residue_scores[i])
                    else:
                        for i in range(1, len(residue_scores)):
                            local_acc[i].append(residue_scores[i][1])
    
def process(target_path, model_path, central_atoms, sphere_radii, thresholds, save_structures, skip_constant_radii, radii_tolerance):
    print('Start processing...')
    print('1. Processing the reference structure...')
    target_str = read_structure(target_path)
    target_model_numbers = get_list_of_models(target_str)
    if len(target_model_numbers)>1:
        target_model_no = target_model_numbers[0]
        target_str = target_str.get_model(target_model_no)
        print('Reference structure includes the following models {}'.format(target_model_numbers))
        print('For further processing only the {} model is used.'.format(str(target_model_no)))
    
    istpdb = isinstance(target_str,PandasPdb)
    target_df = target_str.df['ATOM']
    target_central_atoms_df = get_central_atoms(target_df,central_atoms,istpdb)
    
    spheres = compute_spheres(sphere_radii, target_central_atoms_df, target_df, istpdb)
    print('1. Done.')
    
    print('2. Processing structure predictions...')
    # Prepare thresholds with prefixes
    thresholds_with_prefix = []
    for thresh in thresholds:
        if abs(thresh - 5.0) < 1e-9:
            prefix = ''
        else:
            if thresh.is_integer():
                prefix = f"{int(thresh)}A_"
            else:
                prefix = f"{thresh}A_".replace('.', '_')
        thresholds_with_prefix.append((thresh, prefix))
    
    # Initialize per-threshold accumulators for global and local summary files
    per_threshold_accumulators = {}
    for thresh, _ in thresholds_with_prefix:
        per_threshold_accumulators[thresh] = {}
        for ca in central_atoms:
            per_threshold_accumulators[thresh][ca] = {
                'global': [['model','gGSP']],
                'local': {}
            }
    
    if os.path.isfile(model_path):
        process_model(target_path, model_path, sphere_radii, central_atoms, thresholds_with_prefix, target_df, target_central_atoms_df, spheres, istpdb, per_threshold_accumulators, target_str, save_structures, skip_constant_radii, radii_tolerance) 
    elif os.path.isdir(model_path):
        allfiles = [join(model_path, f) for f in listdir(model_path) if isfile(join(model_path, f))]
        selected_models = [f for f in allfiles if pathlib.Path(f).suffix in ['.pdb','.cif']]
        for selected_model_path in selected_models:
            process_model(target_path, selected_model_path, sphere_radii, central_atoms, thresholds_with_prefix, target_df, target_central_atoms_df, spheres, istpdb, per_threshold_accumulators, target_str, save_structures, skip_constant_radii, radii_tolerance)
    print('2. Done.')
    
    # Write global summary files per threshold and per central atom
    for thresh, prefix in thresholds_with_prefix:
        for ca in central_atoms:
            if (contains_atom(target_central_atoms_df,ca,istpdb)):
                target_filename_without_ext = os.path.splitext(os.path.basename(target_path))[0]
                if prefix:
                    global_summary_file = os.path.join(os.path.dirname(model_path), f"{prefix}{target_filename_without_ext}_{ca}_gGSP.csv")
                    local_summary_file = os.path.join(os.path.dirname(model_path), f"{prefix}{target_filename_without_ext}_{ca}_rGSP.csv")
                else:
                    global_summary_file = os.path.join(os.path.dirname(model_path), f"{target_filename_without_ext}_{ca}_gGSP.csv")
                    local_summary_file = os.path.join(os.path.dirname(model_path), f"{target_filename_without_ext}_{ca}_rGSP.csv")
                save_csv(global_summary_file, per_threshold_accumulators[thresh][ca]['global'])
                local_data = per_threshold_accumulators[thresh][ca]['local'][ca] if ca in per_threshold_accumulators[thresh][ca]['local'] else [['residue']]
                save_csv(local_summary_file, local_data)
    print('Done.')

def run_recompute(csv_path, thresholds):
    print(f'Recomputing from: {csv_path}')
    
    # Utworzenie folderu trimmed
    csv_dir = os.path.dirname(csv_path)
    base_name = os.path.basename(csv_path)
    output_dir = os.path.join(csv_dir, 'trimmed')
    os.makedirs(output_dir, exist_ok=True)

    # Bezpieczny odczyt pliku wiersz po wierszu
    rows = []
    with open(csv_path, 'r', newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            rows.append(row)

    if not rows:
        print("Error: CSV file is empty.")
        return

    # Wyciąganie Residue IDs (pierwszy wiersz, od kolumny 1)
    # Odfiltrowujemy puste wartości, które mogą być na końcu wiersza
    res_ids = [x for x in rows[0][1:] if x.strip()]
    num_residues = len(res_ids)
    
    # Wyciąganie danych (od wiersza 2 w górę)
    data_rows = []
    for i in range(2, len(rows)):
        # Jeśli napotkamy "mean" lub "std dev" w pierwszej kolumnie, kończymy dane
        first_col = str(rows[i][0]).strip().lower()
        if first_col in ['mean', 'std dev', '']:
            break
        data_rows.append(rows[i])
    
    radii = [float(row[0]) for row in data_rows]
    
    # Wyciągamy wartości RMSD (kolumny od 1 do liczby reziduów)
    scores_matrix = [] 
    for row in data_rows:
        # Konwertujemy tylko te kolumny, które odpowiadają reziduom
        scores_matrix.append([float(x) for x in row[1:num_residues+1]])

    # Przeliczanie dla każdego progu
    for thresh in thresholds:
        if abs(thresh - 5.0) < 1e-9: 
            prefix = ''
        else: 
            prefix = f"{int(thresh)}A_" if thresh.is_integer() else f"{thresh}A_".replace('.', '_')
        
        # 1. gGSP
        ggsp_val = compute_ggsp(radii, scores_matrix, thresh)
        ggsp_filename = base_name.replace('-details.csv', '-gGSP.csv')
        save_csv(os.path.join(output_dir, f"{prefix}{ggsp_filename}"), [['gGSP'], [str(ggsp_val)]])
        
        # 2. rGSP
        rgsp_data = [['residue', 'rGSP']]
        for res_idx in range(num_residues):
            accept_count = 0
            rgsp_sum = 0
            for r_idx in range(len(radii)):
                curr_score = scores_matrix[r_idx][res_idx]
                if curr_score <= thresh:
                    accept_count += 1
                # rGSP liczymy jako średnią ważoną promieniami (zgodnie z Twoją logiką w compute_ggsp)
                rgsp_sum += radii[r_idx] * (accept_count / (r_idx + 1))
            
            final_rgsp = (rgsp_sum / sum(radii)) * 100.0
            rgsp_data.append([f'"{res_ids[res_idx]}"', str(round(final_rgsp, 3))])
        
        rgsp_filename = base_name.replace('-details.csv', '-rGSP.csv')
        save_csv(os.path.join(output_dir, f"{prefix}{rgsp_filename}"), rgsp_data)
        
        # 3. Details (kopiujemy oryginał do folderu trimmed)
        details_filename = f"{prefix}{base_name}"
        save_csv(os.path.join(output_dir, details_filename), rows)

    print(f'Done. Results saved in: {output_dir}')

def main(argv):
    target_path, model_path, central_atoms, sphere_radii, thresholds, save_structures, skip_constant_radii, radii_tolerance, recompute_path = read_config(argv)
    
    if recompute_path:
        run_recompute(recompute_path, thresholds)
        return # Kończymy działanie
    
    print('Target path: {}'.format(target_path))
    print('Model path: {}'.format(model_path))
    print('Central atoms: {}'.format(central_atoms))
    print('Sphere radii: {}'.format(sphere_radii))
    print('RMSD thresholds: {}'.format(thresholds))
    print('Save structures: {}'.format(save_structures))
    print('Skip constant radii: {}'.format(skip_constant_radii))
    if skip_constant_radii:
        print('Radii tolerance: {}'.format(radii_tolerance))
    
    if os.path.isfile(target_path):
        process(target_path, model_path, central_atoms, sphere_radii, thresholds, save_structures, skip_constant_radii, radii_tolerance)
    elif os.path.isdir(target_path):
        allfiles = [join(target_path, f) for f in listdir(target_path) if isfile(join(target_path, f))]
        selected_targets = [f for f in allfiles if pathlib.Path(f).suffix in ['.pdb','.cif']]
        for selected_target_path in selected_targets:
            process(selected_target_path, model_path, central_atoms, sphere_radii, thresholds, save_structures, skip_constant_radii, radii_tolerance)

if __name__ == "__main__":
    main(sys.argv[1:])