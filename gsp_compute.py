import os
import sys
import subprocess
import pathlib
import re
import csv
import time
import multiprocessing as mp
import shutil
import tempfile

# Get the absolute path to the compute.py script
BASE_DIR = pathlib.Path(__file__).parent.resolve()
COMPUTE_SCRIPT = BASE_DIR / "src" / "gsp" / "compute.py"
THRESHOLDS = [2.0, 4.0, 5.0, 8.0]
CROSS_COMPARE = False   # ustaw na True, jeśli chcesz porównywać każdy model z każdym rozwiązaniem

def extract_base_id(filename, puzzle_name=None):
    name = pathlib.Path(filename).stem
    for suffix in ['_refined', '_solution', '_rpr']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    
    is_solution = 'solution' in filename.lower()
    
    name = re.sub(r'^\d+_', '', name)
    
    if puzzle_name and name.startswith(puzzle_name + '_'):
        name = name[len(puzzle_name)+1:]
    
    if is_solution:
        if name in ("solution", "solution_"):
            return "solution_0"
    
    if name == "":
        return "solution_0"
    return name

def get_group_base_id(base_id):
    if '_vs_' in base_id:
        return base_id.split('_vs_')[0]
    return base_id

def collect_files(pdb_dir):
    models = []
    solutions = []

    for f in pdb_dir.iterdir():
        if f.is_file() and f.suffix.lower() in ['.pdb', '.cif']:
            if 'solution' in f.name.lower():
                solutions.append(f)
            else:
                models.append(f)

    return models, solutions

def read_gGSP_from_csv(csv_path):
    if csv_path is None or not csv_path.exists():
        return None
    try:
        with open(csv_path, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            col_index = None
            for i, col in enumerate(header):
                if 'ggsp' in col.lower():
                    col_index = i
                    break
            if col_index is None:
                col_index = 0

            first_row = next(reader, None)
            if first_row and col_index < len(first_row):
                return first_row[col_index].strip()
    except Exception as e:
        return f"READ_ERROR: {e}"
    return None

def find_gGSP_files(cwd, output_dir, target_stem, model_stem, thresholds):
    results = {}
    for thresh in thresholds:
        if abs(thresh - 5.0) < 1e-9:
            prefix = ""
        else:
            if thresh.is_integer():
                prefix = f"{int(thresh)}A_"
            else:
                prefix = f"{str(thresh).replace('.', '_')}A_"

        possible_names = [
            f"{prefix}{target_stem}_{model_stem}_C1'-gGSP.csv",
            f"{prefix}{model_stem}_{target_stem}_C1'-gGSP.csv",
            f"{prefix}{target_stem}_{model_stem}_gGSP.csv",
            f"{prefix}{model_stem}_{target_stem}_gGSP.csv",
        ]

        found = None
        for loc in (cwd, output_dir):
            if not loc.exists():
                continue
            for name in possible_names:
                candidate = loc / name
                if candidate.exists():
                    found = candidate
                    break
            if found:
                break

        if not found:
            for loc in (cwd, output_dir):
                if not loc.exists():
                    continue
                for csv_file in loc.glob("*.csv"):
                    if 'gGSP' in csv_file.name and target_stem in csv_file.name and model_stem in csv_file.name:
                        if prefix:
                            if csv_file.name.startswith(prefix):
                                found = csv_file
                                break
                        else:
                            if not any(csv_file.name.startswith(f"{int(t)}A_") for t in thresholds if t != 5.0):
                                found = csv_file
                                break
                if found:
                    break

        if found:
            gGSP = read_gGSP_from_csv(found)
            if gGSP is not None and not gGSP.startswith("READ_ERROR"):
                results[thresh] = (found, gGSP)
            else:
                results[thresh] = (found, f"READ_ERROR: {gGSP}")
        else:
            results[thresh] = (None, "NO_FILE")

    return results

def select_best_for_group(group_results):
    best_key = None
    best_model = None
    best_solution = None
    best_results = None
    for model, sol, res_dict in group_results:
        key = []
        for thresh in THRESHOLDS:
            _, val = res_dict.get(thresh, (None, "N/A"))
            if val in ("NO_FILE", "READ_ERROR: None") or val is None:
                key.append(-1e9)
            else:
                try:
                    key.append(float(val))
                except ValueError:
                    key.append(-1e9)
        if best_key is None or key > best_key:
            best_key = key
            best_model = model
            best_solution = sol
            best_results = res_dict
    return best_model, best_solution, best_results

def run_compute_single(args):
    model_file, group_base_id, solution_file, output_dir, cwd, puzzle_name = args
    model_name = model_file.name
    solution_name = solution_file.name
    target_stem = solution_file.stem
    model_stem = model_file.stem

    thresholds_str = ','.join(str(t) for t in THRESHOLDS)

    cmd = [
        "python",
        str(COMPUTE_SCRIPT),
        "-t", solution_name,
        "-m", model_name,
        "-d", thresholds_str,
    ]

    env = os.environ.copy()
    env["OUTPUT_DIR"] = os.path.abspath(output_dir)

    try:
        result = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=6000)
    except subprocess.TimeoutExpired:
        return (model_name, group_base_id, solution_name, "TIMEOUT", "Timeout after 600s", None)

    if result.returncode != 0:
        details = f"exit {result.returncode}: {result.stderr[:200]}"
        return (model_name, group_base_id, solution_name, "COMPUTE_ERROR", details, None)

    results_dict = find_gGSP_files(cwd, output_dir, target_stem, model_stem, THRESHOLDS)
    if not results_dict:
        return (model_name, group_base_id, solution_name, "NO_CSV", "No CSV file found after computation", None)

    any_found = any(path is not None for path, _ in results_dict.values())
    if not any_found:
        return (model_name, group_base_id, solution_name, "NO_CSV", "No CSV file found after computation", None)

    details = ', '.join(f"{thresh}A: {val}" for thresh, (_, val) in results_dict.items() if val not in ("NO_FILE", "READ_ERROR: None"))
    return (model_name, group_base_id, solution_name, "SUCCESS", details, results_dict)

def safe_move(src, dst, max_retries=3, delay=0.5):
    for attempt in range(max_retries):
        try:
            shutil.copy2(src, dst)
            os.unlink(src)
            return True
        except (PermissionError, OSError) as e:
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                print(f"  [WARNING] Could not move {src} to {dst}: {e}")
                return False
    return False

def safe_unlink(path, max_retries=3, delay=0.5):
    for attempt in range(max_retries):
        try:
            path.unlink()
            return True
        except PermissionError:
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                print(f"  [WARNING] Could not delete {path} (still in use)")
                return False
        except Exception as e:
            print(f"  [WARNING] Failed to delete {path}: {e}")
            return False
    return False

def parse_details_csv(csv_path):
    """
    Parse a -details.csv file generated by compute.py.
    Returns a dict with keys: path, header_rows, radii, scores_matrix,
    data_rows, num_residues, stat_rows.
    stat_rows holds the trailing 'mean' / 'std dev' lines (if any) so
    they can be reproduced (and recalculated) in the trimmed output.
    Returns None on any parse error.
    """
    try:
        rows = []
        with open(csv_path, 'r', newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                rows.append(row)

        if len(rows) < 3:
            return None

        res_ids = [x for x in rows[0][1:] if x.strip()]
        num_residues = len(res_ids)
        if num_residues == 0:
            return None

        header_rows = rows[:2]  # residue-ID row + sequence row

        # Data rows: stop at 'mean', 'std dev', or empty first cell
        data_rows = []
        stat_rows = []
        in_stats = False
        for i in range(2, len(rows)):
            first_col = str(rows[i][0]).strip().lower()
            if first_col in ('mean', 'std dev') or in_stats:
                in_stats = True
                stat_rows.append(rows[i])
            elif not rows[i][0].strip():
                # empty separator – skip
                pass
            else:
                data_rows.append(rows[i])

        if not data_rows:
            return None

        radii = []
        scores_matrix = []
        for row in data_rows:
            try:
                r = float(row[0])
            except (ValueError, IndexError):
                continue
            radii.append(r)
            vals = []
            for x in row[1:num_residues + 1]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(0.0)
            scores_matrix.append(vals)

        if not radii:
            return None

        return {
            'path': pathlib.Path(csv_path),
            'header_rows': header_rows,
            'radii': radii,
            'scores_matrix': scores_matrix,
            'data_rows': data_rows,
            'num_residues': num_residues,
            'stat_rows': stat_rows,   # original mean/std dev rows (may be empty)
        }
    except Exception as e:
        print(f"  [WARNING] Could not parse {csv_path}: {e}")
        return None


def find_global_max_useful_radius_idx(all_parsed, constant_tolerance=1e-9):
    """
    Find the last radius index where AT LEAST ONE model (across ALL parsed
    details files passed in) still has non-constant per-residue row values.

    Call this with ALL parsed files for the whole puzzle (not grouped by target)
    so that the worst-performing model in the entire puzzle set drives the
    common trim point.

    Returns the inclusive upper index to keep, or -1 if all rows are constant
    in every file (degenerate case → callers should keep everything).
    """
    if not all_parsed:
        return -1

    num_radii = max(len(p['radii']) for p in all_parsed)
    last_useful = -1

    for r_idx in range(num_radii):
        for parsed in all_parsed:
            if r_idx >= len(parsed['radii']):
                continue
            row_vals = parsed['scores_matrix'][r_idx]
            if len(row_vals) > 1 and (max(row_vals) - min(row_vals)) > constant_tolerance:
                last_useful = r_idx
                break   # at least one model varies at this radius – keep it

    return last_useful


def _compute_mean_std(scores_matrix, num_residues):
    """
    Compute column-wise mean and std-dev over all data rows in scores_matrix.
    Returns (mean_row, std_row) as lists of floats (length = num_residues).
    std_row uses population std-dev (ddof=0) to match the original compute.py
    behaviour.
    """
    import math

    if not scores_matrix or num_residues == 0:
        return [], []

    n = len(scores_matrix)
    means = []
    stds = []

    for col in range(num_residues):
        vals = [row[col] for row in scores_matrix if col < len(row)]
        if not vals:
            means.append(0.0)
            stds.append(0.0)
            continue
        m = sum(vals) / len(vals)
        variance = sum((v - m) ** 2 for v in vals) / len(vals)
        means.append(m)
        stds.append(math.sqrt(variance))

    return means, stds


def write_trimmed_details_csv(parsed, max_radius_idx, output_path):
    """
    Write a trimmed details CSV keeping only data rows up to max_radius_idx
    (inclusive).

    The trailing 'mean' and 'std dev' rows are RECOMPUTED from the trimmed
    data so the statistics reflect the shorter radius range, matching the
    format expected by downstream tools.
    """
    trimmed_data_rows = parsed['data_rows'][:max_radius_idx + 1]
    trimmed_scores    = parsed['scores_matrix'][:max_radius_idx + 1]
    num_residues      = parsed['num_residues']

    means, stds = _compute_mean_std(trimmed_scores, num_residues)

    mean_row   = ['mean']   + [f"{v:.3f}" for v in means]
    std_row    = ['std dev'] + [f"{v:.3f}" for v in stds]

    all_rows = parsed['header_rows'] + trimmed_data_rows + [mean_row, std_row]

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(all_rows)


def get_base_details_files(directory):
    """
    Return all *-details.csv files in *directory* that do NOT carry a
    threshold prefix (e.g. 2A_, 4A_, 8A_).
    """
    directory = pathlib.Path(directory)
    if not directory.exists():
        return []

    exclude_prefixes = set()
    for thresh in THRESHOLDS:
        if thresh != 5.0:
            if float(thresh).is_integer():
                exclude_prefixes.add(f"{int(thresh)}A_")
            else:
                exclude_prefixes.add(f"{str(thresh).replace('.', '_')}A_")

    result = []
    for f in directory.glob("*-details.csv"):
        if not any(f.name.startswith(pfx) for pfx in exclude_prefixes):
            result.append(f)
    return result



def _get_nonfive_details_prefixes():
    """Return the set of threshold prefixes that are NOT 5 Å (e.g. '2A_', '4A_', '8A_')."""
    prefixes = set()
    for thresh in THRESHOLDS:
        if abs(thresh - 5.0) < 1e-9:
            continue
        if float(thresh).is_integer():
            prefixes.add(f"{int(thresh)}A_")
        else:
            prefixes.add(f"{str(thresh).replace('.', '_')}A_")
    return prefixes


def _remove_nonfive_details(directory):
    """
    Delete all *-details.csv files that carry a non-5 Å threshold prefix
    (e.g. 2A_…-details.csv, 4A_…-details.csv, 8A_…-details.csv).
    """
    directory = pathlib.Path(directory)
    if not directory.exists():
        return
    prefixes = _get_nonfive_details_prefixes()
    removed = 0
    for f in list(directory.glob("*-details.csv")):
        if any(f.name.startswith(pfx) for pfx in prefixes):
            if safe_unlink(f):
                removed += 1
    if removed:
        print(f"  [CLEAN] Removed {removed} redundant non-5Å details file(s) from {directory}")


def _apply_trimming_inplace(finished_dir, puzzle_name):
    """
    Read all base (5 Å) details CSVs in *finished_dir*, determine the
    puzzle-wide trim point (worst model across the whole set), then
    overwrite each details file in place with the trimmed + recalculated
    mean/std-dev version.  Also reruns compute.py --recompute_from_csv so
    that all gGSP/rGSP CSV files are updated to match the trimmed data.
    Outputs go directly to *finished_dir* (no trimmed/ subdirectory).
    """
    finished_dir = pathlib.Path(finished_dir)
    details_files = get_base_details_files(finished_dir)
    if not details_files:
        print(f"  [TRIM] No base details files found in {finished_dir} – skipping")
        return

    all_parsed = []
    for dp in details_files:
        parsed = parse_details_csv(dp)
        if parsed:
            all_parsed.append(parsed)
        else:
            print(f"  [WARNING] Could not parse {dp.name} – skipping in trim pass")

    if not all_parsed:
        return

    max_radius_idx = find_global_max_useful_radius_idx(all_parsed)
    num_radii_orig = max(len(p['radii']) for p in all_parsed)

    if max_radius_idx < 0:
        print(f"  [TRIM] All rows constant – no trimming needed")
        max_radius_idx = num_radii_orig - 1
    else:
        sample_radii = next(
            (p['radii'] for p in all_parsed if max_radius_idx < len(p['radii'])), None
        )
        max_r     = sample_radii[max_radius_idx] if sample_radii else "?"
        n_trimmed = num_radii_orig - (max_radius_idx + 1)
        if n_trimmed > 0:
            print(f"  [TRIM] Puzzle-wide cut: keeping radii up to {max_r} Å "
                  f"(removing {n_trimmed} constant trailing row(s) across all {len(all_parsed)} files)")
        else:
            print(f"  [TRIM] No constant trailing rows – last useful radius: {max_r} Å")

    thresholds_str = ','.join(str(t) for t in THRESHOLDS)
    total_ok  = 0
    total_err = 0

    for parsed in all_parsed:
        original_name = parsed['path'].name

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = pathlib.Path(tmpdir)

            # Write trimmed details to a temp copy (keeps filename for compute.py)
            temp_input = tmpdir_path / original_name
            write_trimmed_details_csv(parsed, max_radius_idx, temp_input)

            # Recompute gGSP/rGSP from the trimmed details
            cmd = [
                "python", str(COMPUTE_SCRIPT),
                "--recompute_from_csv", str(temp_input),
                "-d", thresholds_str,
            ]
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                if proc.returncode != 0:
                    print(f"  [ERROR] compute.py failed for {original_name}:\n"
                          f"          {proc.stderr[:300]}")
                    total_err += 1
                    continue
            except subprocess.TimeoutExpired:
                print(f"  [ERROR] Timeout while processing {original_name}")
                total_err += 1
                continue
            except Exception as e:
                print(f"  [ERROR] {original_name}: {e}")
                total_err += 1
                continue

            # compute.py writes results to <tmpdir>/trimmed/
            tmp_trimmed = tmpdir_path / "trimmed"
            if not tmp_trimmed.exists():
                print(f"  [WARNING] No trimmed output generated for {original_name}")
                total_err += 1
                continue

            # Overwrite original details file with trimmed version
            shutil.copy2(temp_input, parsed['path'])

            # Overwrite gGSP/rGSP files in finished_dir in-place
            n_moved = 0
            for out_f in tmp_trimmed.iterdir():
                # Skip details files produced here – we already wrote the trimmed one above
                if out_f.name.endswith("-details.csv"):
                    continue
                dest = finished_dir / out_f.name
                if dest.exists():
                    dest.unlink()
                shutil.copy2(out_f, dest)
                n_moved += 1

            total_ok += 1

    print(f"  [TRIM] Done – {total_ok} file(s) trimmed in-place, {total_err} error(s)")


def process_puzzle_set_parallel(puzzle_path, puzzle_name, global_results):
    pdb_dir = pathlib.Path(puzzle_path)
    if not pdb_dir.exists():
        print(f"  [ERROR] Directory {pdb_dir} does not exist.")
        return False

    print(f"  [INFO] Processing puzzle: {puzzle_name}")
    print(f"  [INFO] Looking in: {pdb_dir}")

    models, solutions = collect_files(pdb_dir)

    if not models:
        print(f"  [ERROR] No model files found in {pdb_dir}")
        return False

    if not solutions:
        print(f"  [ERROR] No solution files found in {pdb_dir}")
        return False

    print(f"  [INFO] Found {len(models)} model(s) and {len(solutions)} solution(s)")

    temp_output_dir = pathlib.Path("./results/gsp") / puzzle_name
    temp_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [INFO] Temporary results will be saved to: {temp_output_dir}")

    tasks = []
    if CROSS_COMPARE:
        for model_file in models:
            base_id = extract_base_id(model_file.name, puzzle_name)
            group_id = get_group_base_id(base_id)
            for sol_file in solutions:
                tasks.append((model_file, group_id, sol_file, temp_output_dir, pdb_dir, puzzle_name))
    else:
        solutions_by_exact = {}
        for sol in solutions:
            exact_base = extract_base_id(sol.name, puzzle_name)
            solutions_by_exact[exact_base] = sol

        missing_solution_count = 0
        for model_file in models:
            exact_base = extract_base_id(model_file.name, puzzle_name)
            solution_file = solutions_by_exact.get(exact_base)
            if not solution_file:
                generic_base = "solution_0"
                solution_file = solutions_by_exact.get(generic_base)
                if not solution_file:
                    print(f"  [WARNING] No solution for model {model_file.name}, skipping.")
                    missing_solution_count += 1
                    continue
            group_id = get_group_base_id(exact_base)
            tasks.append((model_file, group_id, solution_file, temp_output_dir, pdb_dir, puzzle_name))

    if not tasks:
        print(f"  [INFO] No tasks to process.")
        return False

    print(f"  [INFO] Submitting {len(tasks)} tasks to process pool...")

    pool = mp.Pool(processes=mp.cpu_count())
    results = pool.map(run_compute_single, tasks)
    pool.close()
    pool.join()

    groups = {}
    for res in results:
        model_name, group_base_id, sol_name, status, details, res_dict = res
        if status != "SUCCESS" or res_dict is None:
            continue
        groups.setdefault(group_base_id, []).append((model_name, sol_name, res_dict))

    best_per_group = {}
    for group_id, group_results in groups.items():
        best_model, best_sol, best_res_dict = select_best_for_group(group_results)
        best_per_group[group_id] = (best_model, best_sol, best_res_dict)

    all_group_ids = set()
    for model_file in models:
        exact_base = extract_base_id(model_file.name, puzzle_name)
        group_id = get_group_base_id(exact_base)
        all_group_ids.add(group_id)
    total_groups = len(all_group_ids)
    success = len(best_per_group)
    errors = total_groups - success
    print(f"  [INFO] All tasks completed. Best results obtained for {success}/{total_groups} model groups (errors: {errors}).")

    successful_models = {best_model for _, (best_model, _, _) in best_per_group.items()}
    failed_models_files = {model_file.name for model_file in models} - successful_models
    failure_details = {}
    for res in results:
        model_name, _, sol_name, status, details, _ = res
        if model_name not in successful_models:
            failure_details[model_name] = (sol_name, status, details)

    if best_per_group or failed_models_files:
        report_path = temp_output_dir / "report.txt"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(f"Report for set: {puzzle_name}\n")
            f.write(f"Processed {success} / {total_groups} model groups (best score per group)\n")
            f.write(f"Failed: {len(failed_models_files)} models\n")
            f.write("\n")

            if best_per_group:
                f.write("=== SUCCESSFUL MODEL GROUPS (BEST SCORES) ===\n")
                header = f"{'ModelGroup':<50} {'BestModel':<50} {'BestSolution':<50}"
                for thresh in THRESHOLDS:
                    header += f" {thresh:<12}"
                f.write(header + "\n")
                f.write("-" * (150 + 12 * len(THRESHOLDS)) + "\n")
                for group_id, (best_model, best_sol, res_dict) in best_per_group.items():
                    row = f"{group_id:<50} {best_model:<50} {best_sol:<50}"
                    for thresh in THRESHOLDS:
                        _, val = res_dict.get(thresh, (None, "N/A"))
                        if val in ("NO_FILE", "READ_ERROR: None") or val is None or val == "N/A":
                            row += " N/A        "
                        else:
                            try:
                                num_val = float(val)
                                row += f" {num_val:<12.3f}"
                            except ValueError:
                                row += " N/A        "
                    f.write(row + "\n")
                f.write("\n")

            if failed_models_files:
                f.write("=== FAILED MODELS ===\n")
                f.write(f"{'Model':<50} {'Solution':<50} {'Status':<20} {'Details':<30}\n")
                f.write("-" * 150 + "\n")
                for model in sorted(failed_models_files):
                    sol, status, details = failure_details.get(model, ("", "UNKNOWN", ""))
                    f.write(f"{model:<50} {sol:<50} {status:<20} {details:<30}\n")
                f.write("\n")

        print(f"  [REPORT] Saved to {report_path}")

        global_results[puzzle_name] = {
            "total_groups": total_groups,
            "success": success,
            "best_per_group": best_per_group,
            "failed_models": list(failed_models_files)
        }

    moved_csv = 0
    tracked_csv_paths = set()
    for res in results:
        if res[5] is not None:
            for thresh, (csv_path, _) in res[5].items():
                if csv_path is None:
                    continue
                src = pathlib.Path(csv_path)
                tracked_csv_paths.add(src.resolve())
                if src.exists():
                    dest = temp_output_dir / src.name
                    if dest.exists():
                        dest.unlink()
                    if safe_move(src, dest):
                        moved_csv += 1
                    else:
                        print(f"  [WARNING] Failed to move {src.name}")

    for csv_file in pdb_dir.glob("*.csv"):
        if csv_file.resolve() in tracked_csv_paths:
            continue
        dest = temp_output_dir / csv_file.name
        if dest.exists():
            dest.unlink()
        if safe_move(csv_file, dest):
            moved_csv += 1
        else:
            print(f"  [WARNING] Failed to move {csv_file.name}")

    final_dir = pathlib.Path("results/gsp/finished") / puzzle_name
    final_dir.mkdir(parents=True, exist_ok=True)
    final_moved = 0
    for csv_file in temp_output_dir.glob("*.csv"):
        if csv_file.name == "report.txt":
            continue
        dest = final_dir / csv_file.name
        if dest.exists():
            dest.unlink()
        for attempt in range(5):
            try:
                shutil.move(str(csv_file), str(dest))
                final_moved += 1
                break
            except PermissionError:
                if attempt == 4:
                    print(f"  [ERROR] Failed to move {csv_file.name} after 5 attempts")
                else:
                    time.sleep(0.5)
    print(f"  [MOVE] Moved {final_moved} CSV files to {final_dir}")

    report_src = temp_output_dir / "report.txt"
    if report_src.exists():
        report_dest = final_dir / "report.txt"
        if report_dest.exists():
            report_dest.unlink()
        report_src.rename(report_dest)
        print(f"  [MOVE] Moved report to {report_dest}")

    # ------------------------------------------------------------------ #
    # Auto-step: apply global row trimming in-place and remove redundant  #
    # details files for non-5A thresholds.                                #
    # ------------------------------------------------------------------ #
    print(f"  [INFO] Running automatic trimming pass on {final_dir} ...")
    _apply_trimming_inplace(final_dir, puzzle_name)
    _remove_nonfive_details(final_dir)

    return success > 0

# ---------------------------------------------------------------------------
# RECOMPUTE WITH GLOBAL ROW TRIMMING  –  per-puzzle grouping
# ---------------------------------------------------------------------------

def process_recompute_trimmed_puzzle(puzzle_name):
    """
    Recompute GSP scores from existing -details.csv files applying a
    globally consistent row-trimming strategy for the WHOLE PUZZLE SET.

    All files are updated in-place inside results/gsp/finished/<puzzle>/.
    Non-5 Å details files (2A_, 4A_, 8A_ prefixes) are removed afterwards.
    """
    finished_dir = pathlib.Path("results/gsp/finished") / puzzle_name

    print(f"  [INFO] Recomputing trimmed scores for: {puzzle_name}")
    print(f"  [INFO] Source dir : {finished_dir}")

    if not finished_dir.exists():
        print(f"  [ERROR] Results directory not found: {finished_dir}")
        return False

    _apply_trimming_inplace(finished_dir, puzzle_name)
    _remove_nonfive_details(finished_dir)
    return True


# ---------------------------------------------------------------------------

def parse_selection(input_str, max_value):
    selected = set()
    range_pattern = r'(\d+)-(\d+)'
    for match in re.finditer(range_pattern, input_str):
        start = int(match.group(1))
        end = int(match.group(2))
        if start > end:
            start, end = end, start
        selected.update(range(start, end + 1))
    number_pattern = r'\d+'
    for match in re.findall(number_pattern, input_str):
        num = int(match)
        selected.add(num)
    return [idx for idx in sorted(selected) if 1 <= idx <= max_value]

def select_from_list(items, prompt):
    if not items:
        print("No items available for selection.")
        return []
    print("\nAvailable options:")
    for i, item in enumerate(items, 1):
        print(f"{i}. {item}")
    while True:
        user_input = input(f"\n{prompt} (e.g., 1, 1-3, 1,3,5): ").strip()
        if not user_input:
            print("Selection cannot be empty. Please try again.")
            continue
        selected_indices = parse_selection(user_input, len(items))
        if not selected_indices:
            print("No valid selections. Please try again.")
            continue
        print("\nYou selected:")
        for idx in selected_indices:
            print(f"  - {items[idx-1]}")
        return [items[idx-1] for idx in selected_indices]

def main():
    print("RNA Puzzle Processing – FINISHED DATA FLOW (PARALLEL)")
    print("=" * 60)

    print("\nSelect processing mode:")
    print("  1. Normal GSP computation (from PDB files)")
    print("  2. Recompute from existing details CSV (with global row trimming)")

    while True:
        mode_input = input("\nMode [1/2]: ").strip()
        if mode_input in ("1", "2"):
            break
        print("Please enter 1 or 2.")

    # ------------------------------------------------------------------
    # MODE 1 – original parallel GSP computation
    # ------------------------------------------------------------------
    if mode_input == "1":
        print("\nProcessing files from data/finished/<puzzle>/pdb/")

        data_root = "data/finished"
        if not os.path.exists(data_root):
            print(f"Error: Data directory '{data_root}' not found!")
            sys.exit(1)

        puzzles = [d for d in os.listdir(data_root)
                   if os.path.isdir(os.path.join(data_root, d)) and
                   os.path.isdir(os.path.join(data_root, d, "pdb"))]

        if not puzzles:
            print("No puzzle directories with 'pdb' subfolder found in data/finished.")
            sys.exit(1)

        selected_puzzles = select_from_list(puzzles, "Select puzzles to process")
        if not selected_puzzles:
            print("No puzzles selected. Exiting.")
            sys.exit(0)

        manager = mp.Manager()
        global_results = manager.dict()

        for puzzle_name in selected_puzzles:
            puzzle_path = os.path.join(data_root, puzzle_name, "pdb")
            print(f"\n{'='*60}")
            print(f"Puzzle: {puzzle_name}")
            print('-' * 60)

            success = process_puzzle_set_parallel(puzzle_path, puzzle_name, global_results)

            if success:
                print(f"\n  [COMPLETED] Finished processing {puzzle_name}")
            else:
                print(f"\n  [FAILED] Failed to process {puzzle_name}")

        if global_results:
            global_report_path = pathlib.Path("results/gsp/finished/global_report.txt")
            global_report_path.parent.mkdir(parents=True, exist_ok=True)
            with open(global_report_path, 'w', encoding='utf-8') as f:
                f.write("Global report for all processed sets (best per model group)\n")
                f.write("=" * 60 + "\n")
                for puzzle_name, data in global_results.items():
                    f.write(f"\n{puzzle_name}\n")
                    f.write(f"Processed {data['success']} / {data['total_groups']} model groups\n")
                    if data['best_per_group']:
                        f.write("  Successful model groups:\n")
                        header = f"{'ModelGroup':<50} {'BestModel':<50} {'BestSolution':<50}"
                        for thresh in THRESHOLDS:
                            header += f" {thresh:<12}"
                        f.write(header + "\n")
                        f.write("-" * (150 + 12 * len(THRESHOLDS)) + "\n")
                        for group_id, (best_model, best_sol, res_dict) in data['best_per_group'].items():
                            row = f"{group_id:<50} {best_model:<50} {best_sol:<50}"
                            for thresh in THRESHOLDS:
                                _, val = res_dict.get(thresh, (None, "N/A"))
                                if val in ("NO_FILE", "READ_ERROR: None") or val is None or val == "N/A":
                                    row += " N/A        "
                                else:
                                    try:
                                        num_val = float(val)
                                        row += f" {num_val:<12.3f}"
                                    except ValueError:
                                        row += " N/A        "
                            f.write(row + "\n")
                    if data['failed_models']:
                        f.write("  Failed models:\n")
                        for model in data['failed_models']:
                            f.write(f"    {model}\n")
                    f.write("\n")
            print(f"\nGlobal report saved to: {global_report_path}")

    # ------------------------------------------------------------------
    # MODE 2 – recompute from existing details CSVs with row trimming
    # ------------------------------------------------------------------
    else:
        print("\nRecomputing from existing details CSV files in results/gsp/finished/")

        results_root = pathlib.Path("results/gsp/finished")
        if not results_root.exists():
            print(f"Error: Results directory '{results_root}' not found!")
            sys.exit(1)

        puzzles = sorted(
            d.name for d in results_root.iterdir()
            if d.is_dir() and any(get_base_details_files(d))
        )

        if not puzzles:
            print("No puzzle result directories with base details CSV files found.")
            sys.exit(1)

        selected_puzzles = select_from_list(puzzles, "Select puzzles to recompute")
        if not selected_puzzles:
            print("No puzzles selected. Exiting.")
            sys.exit(0)

        for puzzle_name in selected_puzzles:
            print(f"\n{'='*60}")
            print(f"Puzzle: {puzzle_name}")
            print('-' * 60)

            success = process_recompute_trimmed_puzzle(puzzle_name)

            if success:
                print(f"\n  [COMPLETED] Finished recomputing {puzzle_name}")
            else:
                print(f"\n  [FAILED] Failed to recompute {puzzle_name}")

    print("\n" + "="*60)
    print("Processing completed!")

if __name__ == "__main__":
    mp.freeze_support()
    main()