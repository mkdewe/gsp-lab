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

def extract_base_id(filename, puzzle_name=None):
    """
    Extracts base identifier from a filename.
    Removes typical suffixes (_refined, _solution, _rpr),
    then removes an optional leading number (e.g., "1_")
    and an optional puzzle prefix (e.g., "PZ16a_").
    For solutions: if after removing the prefix we get "solution" or "solution_", return "solution_0".
    """
    name = pathlib.Path(filename).stem
    # Remove typical suffixes
    for suffix in ['_refined', '_solution', '_rpr']:
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    
    # Check if it is a solution (contains "solution" in original name)
    is_solution = 'solution' in filename.lower()
    
    # Remove leading number and underscore (if present)
    name = re.sub(r'^\d+_', '', name)
    
    # Remove puzzle prefix (if provided and present)
    if puzzle_name and name.startswith(puzzle_name + '_'):
        name = name[len(puzzle_name)+1:]
    
    # For solutions: if it became "solution" or "solution_", treat as generic
    if is_solution:
        if name in ("solution", "solution_"):
            return "solution_0"
    
    # If after all removals it is empty, it's a generic solution
    if name == "":
        return "solution_0"
    return name

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
            # Find column with 'gGSP'
            col_index = None
            for i, col in enumerate(header):
                if 'ggsp' in col.lower():
                    col_index = i
                    break
            if col_index is None:
                # If not found, take the first column
                col_index = 0

            first_row = next(reader, None)
            if first_row and col_index < len(first_row):
                return first_row[col_index].strip()
    except Exception as e:
        return f"READ_ERROR: {e}"
    return None

def run_compute_single(args):
    """
    Function executed in a child process.
    args: (model_file, solution_file, output_dir, cwd, puzzle_name)
    Returns a tuple (model_name, solution_name, status, details, csv_path).
    """
    model_file, solution_file, output_dir, cwd, puzzle_name = args
    model_name = model_file.name
    solution_name = solution_file.name
    model_stem = model_file.stem
    target_stem = solution_file.stem

    cmd = [
        "python",
        str(COMPUTE_SCRIPT),
        "-t", solution_name,
        "-m", model_name
    ]

    env = os.environ.copy()
    env["OUTPUT_DIR"] = os.path.abspath(output_dir)

    try:
        result = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        return (model_name, solution_name, "TIMEOUT", "Timeout after 600s", None)

    if result.returncode != 0:
        details = f"exit {result.returncode}: {result.stderr[:200]}"
        return (model_name, solution_name, "COMPUTE_ERROR", details, None)

    # Expected CSV file names (various combinations)
    possible_names = [
        f"{target_stem}_{model_stem}_C1'-gGSP.csv",
        f"{model_stem}_{target_stem}_C1'-gGSP.csv",
        f"{target_stem}_{model_stem}_gGSP.csv",
        f"{model_stem}_{target_stem}_gGSP.csv",
    ]

    # Search in cwd and output_dir
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
        # Fallback: find any CSV file containing both stems
        for loc in (cwd, output_dir):
            if not loc.exists():
                continue
            for f in loc.glob("*.csv"):
                if model_stem in f.name and target_stem in f.name:
                    found = f
                    break
            if found:
                break

    if found:
        gGSP = read_gGSP_from_csv(found)
        if gGSP is None or gGSP.startswith("READ_ERROR"):
            return (model_name, solution_name, "CSV_READ_ERROR", str(gGSP), found)
        # Success
        print(f"  [OK] {model_name} -> {gGSP}")
        return (model_name, solution_name, "SUCCESS", gGSP, found)
    else:
        return (model_name, solution_name, "NO_CSV", "No CSV file found after computation", None)

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

    # Index solutions by base_id
    solution_by_baseid = {}
    for sol in solutions:
        base_id = extract_base_id(sol.name, puzzle_name)
        solution_by_baseid[base_id] = sol
        # print(f"  [DEBUG] Solution: {sol.name} -> base_id: {base_id}")

    # Prepare list of tasks for each model
    tasks = []
    missing_solution_count = 0
    for model_file in models:
        base_id = extract_base_id(model_file.name, puzzle_name)
        solution_file = solution_by_baseid.get(base_id)
        if not solution_file:
            # Try generic solution
            solution_file = solution_by_baseid.get('solution_0')
        if not solution_file:
            print(f"  [WARNING] No solution for model {model_file.name}, skipping.")
            missing_solution_count += 1
            continue
        tasks.append((model_file, solution_file, temp_output_dir, pdb_dir, puzzle_name))

    if not tasks:
        print(f"  [INFO] No tasks to process (missing solution for {missing_solution_count} models).")
        return False

    print(f"  [INFO] Submitting {len(tasks)} tasks to process pool...")

    # Run in parallel
    pool = mp.Pool(processes=mp.cpu_count())
    results = pool.map(run_compute_single, tasks)
    pool.close()
    pool.join()

    # Summary
    total = len(tasks)
    success = sum(1 for r in results if r[2] == "SUCCESS")
    errors = total - success
    print(f"  [INFO] All tasks completed. Successfully processed {success}/{total} models (errors: {errors}).")

    # Move CSV files to temp_output_dir
    moved_csv = 0
    for res in results:
        if res[4] is not None:  # csv_path
            src = pathlib.Path(res[4])
            if src.exists():
                dest = temp_output_dir / src.name
                if dest.exists():
                    dest.unlink()
                if safe_move(src, dest):
                    moved_csv += 1
                else:
                    print(f"  [WARNING] Failed to move {src.name}")

    # Save report for the set
    if results:
        report_path = temp_output_dir / "report.txt"
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(f"Report for set: {puzzle_name}\n")
            f.write(f"Processed {success} / {total} models (missing solution for {missing_solution_count} models)\n")
            f.write("=" * 120 + "\n")
            f.write(f"{'Model':<50} {'Solution':<50} {'Status':<15} {'Details':<30}\n")
            f.write("-" * 120 + "\n")
            for model, sol, status, details, _ in results:
                f.write(f"{model:<50} {sol:<50} {status:<15} {details:<30}\n")
        print(f"  [REPORT] Saved to {report_path}")

        # Add to global results
        global_results[puzzle_name] = {
            "total": total,
            "success": success,
            "missing": missing_solution_count,
            "results": [(model, sol, status, details) for model, sol, status, details, _ in results]
        }

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
        csv_file.rename(dest)
        final_moved += 1
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

    # Find all puzzles (directories in data_root)
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

    # Dictionary for global results
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
            f.write("Global report for all processed sets\n")
            f.write("=" * 60 + "\n")
            for puzzle_name, data in global_results.items():
                f.write(f"\n{puzzle_name}\n")
                f.write(f"Processed {data['success']} / {data['total']} models (missing solution for {data['missing']})\n")
                f.write("-" * 120 + "\n")
                f.write(f"{'Model':<50} {'Solution':<50} {'Status':<15} {'Details':<30}\n")
                f.write("-" * 120 + "\n")
                for model, sol, status, details in data['results']:
                    f.write(f"{model:<50} {sol:<50} {status:<15} {details:<30}\n")
        print(f"\nGlobal report saved to: {global_report_path}")

    print("\n" + "="*60)
    print("Processing completed!")

if __name__ == "__main__":
    # Required for multiprocessing on Windows
    mp.freeze_support()
    main()