import os
import sys
import subprocess
import pathlib
import re
import csv
import time
import multiprocessing as mp
import shutil

BASE_DIR = pathlib.Path(__file__).parent.resolve()
LDDT_SCRIPT = BASE_DIR / "src" / "utils" / "lddt.py"

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def docker_is_running():
    try:
        result = subprocess.run(['docker', 'info'], capture_output=True, text=True, timeout=5)
        return result.returncode == 0
    except:
        return False

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

def run_lddt_single(args):
    model_file, solution_file, temp_output_dir, puzzle_name, debug_flag = args
    model_name = model_file.name
    solution_name = solution_file.name

    temp_output_dir = temp_output_dir.resolve()
    task_dir = temp_output_dir / f"task_{os.getpid()}_{hash((model_name, solution_name))}"
    task_dir.mkdir(parents=True, exist_ok=True)
    output_json = task_dir / "output.json"

    cmd = [
        "python",
        str(LDDT_SCRIPT),
        "-t", str(solution_file),
        "-m", str(model_file),
        "-o", str(output_json),
        "--root", str(BASE_DIR)
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, cwd=task_dir)
        if debug_flag:
            print(f"    [DEBUG] Running: {' '.join(cmd)}")
            print(f"    [DEBUG] Return code: {result.returncode}")
            if result.stderr:
                print(f"    [DEBUG] STDERR: {result.stderr[:500]}")
        if result.returncode != 0:
            details = f"exit {result.returncode}: {result.stderr[:200]}"
            return (model_name, solution_name, "ERROR", None, details)
        lddt = result.stdout.strip()
        try:
            lddt_float = float(lddt)
            return (model_name, solution_name, "SUCCESS", lddt_float, "")
        except ValueError:
            return (model_name, solution_name, "INVALID", None, f"Output not a number: {lddt}")
    except subprocess.TimeoutExpired:
        return (model_name, solution_name, "TIMEOUT", None, "Timeout after 300s")
    except Exception as e:
        return (model_name, solution_name, "EXCEPTION", None, str(e))
    finally:
        try:
            shutil.rmtree(task_dir)
        except:
            pass

def process_puzzle_set(puzzle_path, puzzle_name, global_results, debug=False):
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

    temp_output_dir = pathlib.Path("./results/lddt") / puzzle_name
    temp_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [INFO] Temporary results will be saved to: {temp_output_dir}")

    solutions_by_exact = {}
    for sol in solutions:
        exact_base = extract_base_id(sol.name, puzzle_name)
        solutions_by_exact[exact_base] = sol
        if debug:
            print(f"    [DEBUG] Solution: {sol.name} -> exact base: '{exact_base}'")

    tasks = []
    missing_solution_count = 0
    for model_file in models:
        exact_base = extract_base_id(model_file.name, puzzle_name)
        if debug:
            print(f"    [DEBUG] Model: {model_file.name} -> exact base: '{exact_base}'")
        solution_file = solutions_by_exact.get(exact_base)
        if not solution_file:
            if debug:
                print(f"    [DEBUG] No exact match, trying fallback 'solution_0'")
            solution_file = solutions_by_exact.get("solution_0")
            if not solution_file:
                print(f"  [WARNING] No solution for model {model_file.name}, skipping.")
                missing_solution_count += 1
                continue
        # Use absolute paths
        model_abs = model_file.resolve()
        sol_abs = solution_file.resolve()
        tasks.append((model_abs, sol_abs, temp_output_dir, puzzle_name, debug))

    if not tasks:
        print(f"  [INFO] No tasks to process.")
        return False

    print(f"  [INFO] Submitting {len(tasks)} tasks to process pool...")

    pool = mp.Pool(processes=mp.cpu_count())
    results = pool.map(run_lddt_single, tasks)
    pool.close()
    pool.join()

    debug_log = temp_output_dir / "debug.log"
    with open(debug_log, 'w', encoding='utf-8') as f:
        for r in results:
            f.write(f"{r}\n")
    print(f"  [DEBUG] Detailed results saved to {debug_log}")

    groups = {}
    for model_name, sol_name, status, lddt, details in results:
        if status != "SUCCESS" or lddt is None:
            if debug:
                print(f"    [DEBUG] Failed: {model_name} vs {sol_name} -> {status}: {details}")
            continue
        base = extract_base_id(model_name, puzzle_name)
        group_id = get_group_base_id(base)
        groups.setdefault(group_id, []).append((model_name, sol_name, lddt))

    best_per_group = {}
    for group_id, group_results in groups.items():
        best = max(group_results, key=lambda x: x[2])
        best_per_group[group_id] = best

    all_group_ids = set()
    for model_file in models:
        base = extract_base_id(model_file.name, puzzle_name)
        group_id = get_group_base_id(base)
        all_group_ids.add(group_id)
    total_groups = len(all_group_ids)
    success = len(best_per_group)
    errors = total_groups - success
    print(f"  [INFO] All tasks completed. Best results obtained for {success}/{total_groups} model groups (errors: {errors}).")

    successful_models = {best[0] for best in best_per_group.values()}
    failed_models = {model_file.name for model_file in models} - successful_models

    report_path = temp_output_dir / "report.csv"
    with open(report_path, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["ModelGroup", "BestModel", "BestSolution", "lDDT"])
        for group_id, (model, sol, lddt) in sorted(best_per_group.items()):
            writer.writerow([group_id, model, sol, f"{lddt:.3f}"])
        if failed_models:
            writer.writerow([])
            writer.writerow(["Failed models:"])
            for model in sorted(failed_models):
                writer.writerow([model])

    print(f"  [REPORT] Saved to {report_path}")

    global_results[puzzle_name] = {
        "total_groups": total_groups,
        "success": success,
        "best_per_group": best_per_group,
        "failed_models": list(failed_models)
    }

    final_dir = pathlib.Path("results/lddt/finished") / puzzle_name
    final_dir.mkdir(parents=True, exist_ok=True)
    final_report = final_dir / "report.csv"
    shutil.move(str(report_path), str(final_report))
    print(f"  [MOVE] Moved report to {final_report}")

    return success > 0

# -----------------------------------------------------------------------------
# Selection and main
# -----------------------------------------------------------------------------
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
    print("RNA Puzzle Processing – lDDT BATCH (PARALLEL)")
    print("=" * 60)
    print("Processing files from data/finished/<puzzle>/pdb/")

    if not docker_is_running():
        print("\n[ERROR] Docker is not running! Please start Docker Desktop and try again.\n")
        sys.exit(1)

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

    debug_input = input("\nEnable debug mode? (y/n) [n]: ").strip().lower()
    debug = debug_input == 'y'

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

        success = process_puzzle_set(puzzle_path, puzzle_name, global_results, debug)

        if success:
            print(f"\n  [COMPLETED] Finished processing {puzzle_name}")
        else:
            print(f"\n  [FAILED] Failed to process {puzzle_name}")

    if global_results:
        global_report_path = pathlib.Path("results/lddt/finished/global_report.csv")
        global_report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(global_report_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Puzzle", "ModelGroup", "BestModel", "BestSolution", "lDDT"])
            for puzzle_name, data in global_results.items():
                for group_id, (model, sol, lddt) in sorted(data["best_per_group"].items()):
                    writer.writerow([puzzle_name, group_id, model, sol, f"{lddt:.3f}"])
        print(f"\nGlobal report saved to: {global_report_path}")

    print("\n" + "="*60)
    print("Processing completed!")

if __name__ == "__main__":
    mp.freeze_support()
    main()