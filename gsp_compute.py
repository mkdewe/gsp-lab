import os
import sys
import subprocess
import pathlib
import re
import csv
import time
import multiprocessing as mp
import shutil

# Get the absolute path to the compute.py script
BASE_DIR = pathlib.Path(__file__).parent.resolve()
COMPUTE_SCRIPT = BASE_DIR / "src" / "gsp" / "compute.py"
THRESHOLDS = [2.0, 4.0, 5.0, 8.0]
CROSS_COMPARE = False   # ustaw na True, jeśli chcesz porównywać każdy model z każdym rozwiązaniem

def extract_base_id(filename, puzzle_name=None):
    """
    Extracts base identifier from a filename.
    Removes typical suffixes (_refined, _solution, _rpr),
    then removes an optional leading number (e.g., "1_")
    and an optional puzzle prefix (e.g., "PZ16a_").
    For solutions: if after removing the prefix we get "solution" or "solution_", return "solution_0".
    """
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
    """
    Returns common base for variants (e.g., CR1108_TS029_1_vs_CR1108__0 -> CR1108_TS029_1).
    """
    if '_vs_' in base_id:
        return base_id.split('_vs_')[0]
    return base_id

def collect_files(pdb_dir):
    """
    Returns a tuple (models, solutions) for .pdb/.cif files in the given directory.
    Models: files without '_solution' in the name.
    Solutions: files with '_solution' in the name.
    """
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
    """
    Reads the gGSP value from a CSV file.
    Assumes the file has a header and a column containing 'gGSP' (case-insensitive).
    """
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
    """
    For each threshold, finds the gGSP CSV file.
    Returns dict: {threshold: (path_to_csv, gGSP_value)}.
    """
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
    """
    group_results: list of (model_filename, solution_filename, results_dict)
    results_dict: {threshold: (csv_path, gGSP_value)}
    Selects the best result according to: highest gGSP for 2A, then 5A, then 8A.
    Returns (best_model, best_solution, best_results_dict).
    """
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
    """
    Function executed in a child process.
    args: (model_file, group_base_id, solution_file, output_dir, cwd, puzzle_name)
    Returns tuple (model_filename, group_base_id, solution_filename, status, details, results_dict).
    """
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
        "--skip_constant_radii",
        "--radii_tolerance", "1e-9"
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
    """Safely move a file with retries."""
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
    """Try to delete a file, retrying on error."""
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

def process_puzzle_set_parallel(puzzle_path, puzzle_name, global_results):
    """
    Processes a single puzzle set in parallel.
    """
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

    # Prepare temporary directory for results
    temp_output_dir = pathlib.Path("./results/gsp") / puzzle_name
    temp_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [INFO] Temporary results will be saved to: {temp_output_dir}")

    # Prepare tasks: one task per model, using exact matching (full base ID)
    tasks = []
    if CROSS_COMPARE:
        # Cross-comparison mode: each model vs every solution
        for model_file in models:
            base_id = extract_base_id(model_file.name, puzzle_name)
            group_id = get_group_base_id(base_id)
            for sol_file in solutions:
                tasks.append((model_file, group_id, sol_file, temp_output_dir, pdb_dir, puzzle_name))
    else:
        # One-to-one mode: each model matched to its exact solution
        solutions_by_exact = {}
        for sol in solutions:
            exact_base = extract_base_id(sol.name, puzzle_name)
            solutions_by_exact[exact_base] = sol

        missing_solution_count = 0
        for model_file in models:
            exact_base = extract_base_id(model_file.name, puzzle_name)
            solution_file = solutions_by_exact.get(exact_base)
            if not solution_file:
                # Fallback to generic solution_0
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

    # Run in parallel
    pool = mp.Pool(processes=mp.cpu_count())
    results = pool.map(run_compute_single, tasks)
    pool.close()
    pool.join()

    # Group results by group_base_id
    groups = {}
    for res in results:
        model_name, group_base_id, sol_name, status, details, res_dict = res
        if status != "SUCCESS" or res_dict is None:
            continue
        groups.setdefault(group_base_id, []).append((model_name, sol_name, res_dict))

    # For each group, select the best solution
    best_per_group = {}
    for group_id, group_results in groups.items():
        best_model, best_sol, best_res_dict = select_best_for_group(group_results)
        best_per_group[group_id] = (best_model, best_sol, best_res_dict)

    # Compute total number of model groups (unique group_base_id among models)
    all_group_ids = set()
    for model_file in models:
        exact_base = extract_base_id(model_file.name, puzzle_name)
        group_id = get_group_base_id(exact_base)
        all_group_ids.add(group_id)
    total_groups = len(all_group_ids)
    success = len(best_per_group)
    errors = total_groups - success
    print(f"  [INFO] All tasks completed. Best results obtained for {success}/{total_groups} model groups (errors: {errors}).")

    # Prepare detailed failure list
    successful_models = {best_model for _, (best_model, _, _) in best_per_group.items()}
    failed_models_files = {model_file.name for model_file in models} - successful_models
    failure_details = {}
    for res in results:
        model_name, _, sol_name, status, details, _ = res
        if model_name not in successful_models:
            failure_details[model_name] = (sol_name, status, details)

    # Save report for the set
    if best_per_group or failed_models_files:
        report_path = temp_output_dir / "report.txt"
        report_path.parent.mkdir(parents=True, exist_ok=True)   # <-- ZABEZPIECZENIE
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

        # Store in global_results
        global_results[puzzle_name] = {
            "total_groups": total_groups,
            "success": success,
            "best_per_group": best_per_group,
            "failed_models": list(failed_models_files)
        }

    # Move all CSV files from results_dict to temp_output_dir
    moved_csv = 0
    for res in results:
        if res[5] is not None:
            for thresh, (csv_path, _) in res[5].items():
                if csv_path is None:
                    continue
                src = pathlib.Path(csv_path)
                if src.exists():
                    dest = temp_output_dir / src.name
                    if dest.exists():
                        dest.unlink()
                    if safe_move(src, dest):
                        moved_csv += 1
                    else:
                        print(f"  [WARNING] Failed to move {src.name}")

    # Move all CSV files from temp_output_dir to the final directory
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

    # Move report
    report_src = temp_output_dir / "report.txt"
    if report_src.exists():
        report_dest = final_dir / "report.txt"
        if report_dest.exists():
            report_dest.unlink()
        report_src.rename(report_dest)
        print(f"  [MOVE] Moved report to {report_dest}")

    # Clean up: delete all CSV files from the source directory (pdb_dir)
    cleaned = 0
    for csv_file in pdb_dir.glob("*.csv"):
        if safe_unlink(csv_file):
            cleaned += 1
    if cleaned:
        print(f"  [CLEAN] Removed {cleaned} CSV files from {pdb_dir}")

    return success > 0

def parse_selection(input_str, max_value):
    """Parses user input like '1', '1-5', '1,3,5' into list of indices."""
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
    """Displays a numbered list and lets user select items."""
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
    print("Processing files from data/finished/<puzzle>/pdb/")

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

    selected_puzzles = select_from_list(
        puzzles,
        "Select puzzles to process"
    )
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

    # Generate global report
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

    print("\n" + "="*60)
    print("Processing completed!")

if __name__ == "__main__":
    mp.freeze_support()
    main()