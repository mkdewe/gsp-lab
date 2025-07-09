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
        opts, args = getopt.getopt(argv,"sht:m:a:r:d:",["save_structures","help","target_path=","model_path=","central_atoms=","sphere_radii=","rmsd_threshold="])
    except getopt.GetoptError as error:
        print('{}'.format(error))
        print('compute.py -t <target_path> (required) -m <model_path> (required) -a <central_atoms> -r <sphere_radii> -d <rmsd_threshold> -s <save_structures>')
        sys.exit(1)
    target_path = None
    model_path = None
    central_atoms = "CA,C1'".split(',')
    sphere_radii = '2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0'.split(',')
    sphere_radii[:] = filter(float.__instancecheck__, map(floatify, sphere_radii))
    rmsd_threshold = 5.0
    required_arguments_count = 0
    target_options = 0
    model_options = 0;
    save_structures = False
    for opt, arg in opts:
        if opt == '-h':
            print('compute.py -t <target_path> (required) -m <model_path> (required) -a <central_atoms> -r <sphere_radii> -d <rmsd_threshold> -s <save_structures>')
            sys.exit()
        elif opt in ("-t", "--target_path"):
            if (os.path.isfile(arg) and os.path.exists(arg)) or (not os.path.isfile(arg) and os.path.exists(arg)):
                target_path = arg
                target_options = 1
            else:
                print('Target path should reference the existing file or directory that includes PDB or CIF files!')
                sys.exit(1)
        elif opt in ("-m", "--model_path"):
            if (os.path.isfile(arg) and os.path.exists(arg)) or (not os.path.isfile(arg) and os.path.exists(arg)):
                model_path = arg
                model_options = 1
            else:
                print('Model path should reference the existing file or directory that includes PDB or CIF files!')
                sys.exit(1)
        elif opt in ("-a", "--central_atoms"):
            central_atoms = [atom.strip() for atom in arg.split(',')]
        elif opt in ("-r", "--sphere_radii"):
            sphere_radii = [radius.strip() for radius in arg.split(',')]
            sphere_radii[:] = filter(float.__instancecheck__, map(floatify, sphere_radii))
            sphere_radii.sort()
        elif opt in ("-d", "--rmsd_threshold"):
            input_value = floatify(arg)
            if input_value != None:
                rmsd_threshold = input_value
        elif opt in ("-s", "--save_structures"):
            save_structures = True
    if target_options == 0:
        print('You have to set a value for the following config parameter: target_path!')
    else:
        required_arguments_count = required_arguments_count + 1
    if model_options == 0:
        print('You have to set a value for the following config parameter: model_path!')
    else:
        required_arguments_count = required_arguments_count + 1
    if required_arguments_count < 2:
        print('compute.py -t <target_path> (required) -m <model_path> (required) -a <central_atoms> -r <sphere_radii> -d <rmsd_threshold> -s <save_structures>')
        sys.exit(1)
    return target_path, model_path, central_atoms, sphere_radii, rmsd_threshold, save_structures
    
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
        
def compute_glrmsd(sphere_radii, selected_scores, rmsd_threshold):
    glrmsd = 0
    for r in range(len(sphere_radii)):
        accept_count = 0
        all_count = 0
        for curr_score in selected_scores[r]:
            if not (curr_score > rmsd_threshold):
                accept_count += 1
            all_count += 1
        glrmsd += sphere_radii[r] * accept_count / all_count
    glrmsd /= sum(sphere_radii)
    return round(glrmsd * 100.0,3)
    
def compute_whole_matrix(sphere_radii, target_df, selected_central_atom_idxs, target_central_atoms_df, istpdb, scores):
    target_sequence = get_sequence(target_central_atoms_df, istpdb, selected_central_atom_idxs)
    matrix = []
    row = ['']
    for i in range(len(selected_central_atom_idxs)):
        row.append('"' + get_residue_id(target_central_atoms_df.iloc[selected_central_atom_idxs[i]],istpdb) + '"')
    matrix.append(row)
    matrix.append(['radius'] + [*target_sequence] + ['mean','std dev'])
    selected_scores = []
    for r in range(len(sphere_radii)):
        row = [str(sphere_radii[r])]
        rows_vect = []
        for i in range(len(selected_central_atom_idxs)):
            selected_idx = selected_central_atom_idxs[i]
            row.append(str(scores[selected_idx][r]))
            rows_vect.append(scores[selected_idx][r])
        row.append(str(round(np.mean(rows_vect),3)))
        row.append(str(round(np.std(rows_vect),3)))
        matrix.append(row)
        selected_scores.append(rows_vect)
    col_means_row = ['mean']
    col_std_row = ['std dev']
    for i in range(len(selected_central_atom_idxs)):
        cols_vect = []
        for r in range(len(sphere_radii)):
            cols_vect.append(scores[selected_central_atom_idxs[i]][r])
        col_means_row.append(str(round(np.mean(cols_vect),3)))
        col_std_row.append(str(round(np.std(cols_vect),3)))
    matrix.append(col_means_row)
    matrix.append(col_std_row)
    return matrix, selected_scores
    
def compute_residue_scores(sphere_radii, target_central_atoms_df, selected_central_atom_idxs, rmsd_threshold, istpdb, scores):
    residue_scores = [['residue','rGSP']]
    for i in range(len(selected_central_atom_idxs)):
        selected_idx = selected_central_atom_idxs[i]
        accept_count = 0
        all_count = 0
        rlrmsd = 0
        row = []
        for r in range(len(sphere_radii)):
            curr_score = scores[selected_idx][r]
            if not (curr_score > rmsd_threshold):
                accept_count += 1
            all_count += 1
            rlrmsd += sphere_radii[r] * accept_count / all_count
        rlrmsd /= sum(sphere_radii)
        row.append('"' + get_residue_id(target_central_atoms_df.iloc[selected_central_atom_idxs[i]],istpdb) + '"')
        row.append(str(round(rlrmsd * 100.0,3)))
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
    
def process_model(target_path, model_path, sphere_radii, central_atoms, rmsd_threshold, target_df, target_central_atoms_df, spheres, istpdb, global_scores, local_scores, target_str, save_structures):
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
                if ca not in global_scores:
                    global_scores[ca] = [['model','gGSP']]
                if ca not in local_scores:
                    local_scores[ca] = [['residue']]
                result_filename = '{}_{}.csv'.format(common_filename,ca)
                result_file_path = os.path.join(os.path.dirname(model_path),result_filename)
                selected_central_atom_idxs = get_indices(target_central_atoms_df, ca, istpdb)
                
                matrix, selected_scores = compute_whole_matrix(sphere_radii, target_df, selected_central_atom_idxs, target_central_atoms_df, istpdb, scores)
                save_csv(result_file_path.replace('.csv','-details.csv'), matrix)
                
                residue_scores = compute_residue_scores(sphere_radii, target_central_atoms_df, selected_central_atom_idxs, rmsd_threshold, istpdb, scores)
                save_csv(result_file_path.replace('.csv','-rGSP.csv'), residue_scores)
                update_local_scores(local_scores, ca, model_path, residue_scores)
                
                glrmsd = compute_glrmsd(sphere_radii, selected_scores, rmsd_threshold)
                save_csv(result_file_path.replace('.csv','-gGSP.csv'), [['gGSP'],[str(glrmsd)]])
                global_scores[ca].append([model_filename_without_ext, str(glrmsd)])
    
def process(target_path, model_path, central_atoms, sphere_radii, rmsd_threshold, save_structures):
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
    global_scores = {}
    local_scores = {}
    
    if os.path.isfile(model_path):
        process_model(target_path, model_path, sphere_radii, central_atoms, rmsd_threshold, target_df, target_central_atoms_df, spheres, istpdb, global_scores, local_scores, target_str, save_structures) 
    elif os.path.isdir(model_path):
        allfiles = [join(model_path, f) for f in listdir(model_path) if isfile(join(model_path, f))]
        selected_models = [f for f in allfiles if pathlib.Path(f).suffix in ['.pdb','.cif']]
        for selected_model_path in selected_models:
            process_model(target_path, selected_model_path, sphere_radii, central_atoms, rmsd_threshold, target_df, target_central_atoms_df, spheres, istpdb, global_scores, local_scores, target_str, save_structures)
    print('2. Done.')
            
    for ca in central_atoms:
        if (contains_atom(target_central_atoms_df,ca,istpdb)):
            target_filename_without_ext = os.path.splitext(os.path.basename(target_path))[0]
            global_scores_summary_file_path = os.path.join(os.path.dirname(model_path),'{}_{}_gGSP.csv'.format(target_filename_without_ext,ca))
            save_csv(global_scores_summary_file_path, global_scores[ca])
            local_scores_summary_file_path = os.path.join(os.path.dirname(model_path),'{}_{}_rGSP.csv'.format(target_filename_without_ext,ca))
            save_csv(local_scores_summary_file_path, local_scores[ca])
    print('Done.')

def main(argv):
    target_path, model_path, central_atoms, sphere_radii, rmsd_threshold, save_structures = read_config(argv)
    print('Target path: {}'.format(target_path))
    print('Model path: {}'.format(model_path))
    print('Central atoms: {}'.format(central_atoms))
    print('Sphere radii: {}'.format(sphere_radii))
    print('RMSD threshold: {}'.format(rmsd_threshold))
    print('Save structures: {}'.format(save_structures))
    
    if os.path.isfile(target_path):
        process(target_path, model_path, central_atoms, sphere_radii, rmsd_threshold, save_structures)
    elif os.path.isdir(target_path):
        allfiles = [join(target_path, f) for f in listdir(target_path) if isfile(join(target_path, f))]
        selected_targets = [f for f in allfiles if pathlib.Path(f).suffix in ['.pdb','.cif']]
        for selected_target_path in selected_targets:
            process(selected_target_path, model_path, central_atoms, sphere_radii, rmsd_threshold, save_structures)

if __name__ == "__main__":
    main(sys.argv[1:])