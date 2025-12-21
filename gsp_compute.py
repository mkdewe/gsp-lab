import os
import sys
import subprocess
import pathlib
import re
import shutil

# Get the absolute path to the compute.py script
BASE_DIR = pathlib.Path(__file__).parent.resolve()
COMPUTE_SCRIPT = BASE_DIR / "src" / "gsp" / "compute.py"

# Base directories
BASE_PROCESSED_DIR = pathlib.Path("data/processed")
BASE_EXTERNAL_DIR = pathlib.Path("data/external")

def extract_base_id(filename):
    """Extracts the base ID from filename (e.g., Bujnicki_1 from Bujnicki_1.pdb)"""
    base_name = pathlib.Path(filename).stem
    
    # Lista patternów do usuwania typowych sufiksów
    patterns_to_remove = [
        r'_refined$', r'_model$', r'_solution$', r'_clean$', r'_processed$',
        r'_predicted$', r'_optimized$', r'_final$'
    ]
    
    for pattern in patterns_to_remove:
        base_name = re.sub(pattern, '', base_name)
    
    # Usuwanie prefixu solution_ jeśli istnieje
    if base_name.startswith('solution_'):
        base_name = base_name.replace('solution_', '', 1)
    
    # Usuwanie prefixu puzzle (np. PZ3_) jeśli istnieje
    puzzle_prefix_pattern = r'^(PZ\d+_|RP\d+_)'
    base_name = re.sub(puzzle_prefix_pattern, '', base_name)
    
    return base_name

def find_processed_files_by_base_id(repo_name, puzzle_name, base_id):
    """
    Znajduje przetworzone pliki dla danego base_id w katalogu processed.
    Zwraca krotkę (processed_model, processed_solution) - każdy może być None.
    """
    processed_dir = BASE_PROCESSED_DIR / repo_name / "data" / puzzle_name / "pdb"
    if not processed_dir.exists():
        return None, None
    
    processed_model = None
    processed_solution = None
    
    for candidate in processed_dir.iterdir():
        if candidate.is_file() and candidate.suffix in ['.pdb', '.cif']:
            candidate_base_name = candidate.stem
            candidate_base_id = extract_base_id(candidate.name)
            
            # Debug info
            # print(f"  [DEBUG] Checking {candidate.name}: base_id={candidate_base_id}, target={base_id}")
            
            if candidate_base_id == base_id:
                # Sprawdzamy czy to jest solution czy model
                if 'solution' in candidate.name.lower():
                    processed_solution = candidate
                    # print(f"  [DEBUG] Found solution: {candidate.name}")
                else:
                    # To jest model
                    processed_model = candidate
                    # print(f"  [DEBUG] Found model: {candidate.name}")
    
    return processed_model, processed_solution

def find_default_solution_in_processed(repo_name, puzzle_name, model_base_id=None):
    """
    Znajduje domyślny solution w processed.
    Jeśli podano model_base_id, szuka solution z tym samym base_id.
    W przeciwnym razie szuka dowolnego solution.
    """
    processed_dir = BASE_PROCESSED_DIR / repo_name / "data" / puzzle_name / "pdb"
    if not processed_dir.exists():
        return None
    
    solutions = []
    
    for candidate in processed_dir.iterdir():
        if candidate.is_file() and candidate.suffix in ['.pdb', '.cif']:
            if 'solution' in candidate.name.lower():
                candidate_base_id = extract_base_id(candidate.name)
                
                if model_base_id:
                    # Szukamy solution z tym samym base_id co model
                    if candidate_base_id == model_base_id:
                        return candidate
                else:
                    solutions.append(candidate)
    
    # Jeśli nie znaleziono solution z tym samym base_id, zwróć pierwszy znaleziony
    if solutions and not model_base_id:
        return solutions[0]
    
    return None

def find_default_solution_in_external(puzzle_path, model_base_id=None):
    """
    Znajduje domyślny solution w katalogu external.
    """
    solutions = []
    
    for fname in os.listdir(puzzle_path):
        if 'solution' in fname.lower() and pathlib.Path(fname).suffix in ['.pdb', '.cif']:
            if model_base_id:
                # Sprawdzamy czy solution ma ten sam base_id co model
                solution_base_id = extract_base_id(fname)
                if solution_base_id == model_base_id:
                    return fname
            else:
                solutions.append(fname)
    
    # Jeśli nie znaleziono solution z tym samym base_id, zwróć pierwszy znaleziony
    if solutions and not model_base_id:
        return solutions[0]
    
    return None

def collect_model_files_from_external(puzzle_path):
    """
    Zbiera wszystkie pliki modeli z katalogu external.
    """
    models = []
    for fname in os.listdir(puzzle_path):
        if pathlib.Path(fname).suffix in ['.pdb', '.cif']:
            # Pomijamy pliki z 'solution' w nazwie
            if 'solution' in fname.lower():
                continue
            models.append(fname)
    return models

def copy_file_to_working_dir(source_path, dest_dir, new_name=None):
    """
    Kopiuje plik do katalogu roboczego.
    """
    source_path = pathlib.Path(source_path)
    if new_name:
        dest_path = dest_dir / new_name
    else:
        dest_path = dest_dir / source_path.name
    
    try:
        shutil.copy2(source_path, dest_path)
        return dest_path.name
    except Exception as e:
        print(f"  [ERROR] Failed to copy {source_path}: {e}")
        return None

def run_compute_simple(target_file, model_file, output_dir, cwd):
    """
    Uruchamia compute.py bez fallbacków (pliki są już przygotowane).
    """
    model_stem = pathlib.Path(model_file).stem
    target_stem = pathlib.Path(target_file).stem
    output_filename = f"{model_stem}_vs_{target_stem}.csv"
    output_path = pathlib.Path(output_dir) / output_filename

    if output_path.exists():
        print(f"  [SKIP] Output file {output_filename} already exists. Overwriting...")
        output_path.unlink()

    cmd = [
        "python",
        str(COMPUTE_SCRIPT),
        "-t", target_file,
        "-m", model_file
    ]
    
    env = os.environ.copy()
    env["OUTPUT_DIR"] = os.path.abspath(output_dir)
    
    result = subprocess.run(cmd, cwd=cwd, env=env)
    
    if result.returncode == 0:
        print(f"  [SUCCESS] Created {output_filename}")
        return True
    else:
        print(f"  [ERROR] Processing {model_file} vs {target_file} (exit {result.returncode})")
        return False

def process_single_model(model_file, puzzle_path, repo_name, puzzle_name, output_dir):
    """
    Przetwarza pojedynczy model zgodnie z nowym flow.
    """
    print(f"\n  [PROCESSING] Model: {model_file}")
    
    # 1. Pobierz base_id z pliku w external
    base_id = extract_base_id(model_file)
    print(f"  [BASE_ID] {base_id}")
    
    # 2. Szukaj przetworzonego modelu i predefiniowanego solution w processed
    processed_model, processed_solution = find_processed_files_by_base_id(
        repo_name, puzzle_name, base_id
    )
    
    # 3. Przygotuj pliki do użycia
    files_to_cleanup = []
    
    # Przygotuj model
    model_to_use = None
    if processed_model:
        print(f"  [FOUND] Processed model: {processed_model.name}")
        model_to_use = copy_file_to_working_dir(
            processed_model, puzzle_path, f"processed_{model_file}"
        )
        if model_to_use:
            files_to_cleanup.append(model_to_use)
    else:
        print(f"  [INFO] Using original model from external")
        model_to_use = model_file
    
    if not model_to_use:
        print(f"  [ERROR] Failed to prepare model file")
        return False
    
    # Przygotuj solution
    solution_to_use = None
    
    # Najpierw spróbuj predefiniowany solution z processed
    if processed_solution:
        print(f"  [FOUND] Predefined solution for base_id {base_id}: {processed_solution.name}")
        solution_to_use = copy_file_to_working_dir(
            processed_solution, puzzle_path, f"solution_{base_id}.pdb"
        )
        if solution_to_use:
            files_to_cleanup.append(solution_to_use)
    
    # Jeśli nie ma predefiniowanego, spróbuj domyślny solution z processed z tym samym base_id
    if not solution_to_use:
        default_solution_processed = find_default_solution_in_processed(repo_name, puzzle_name, base_id)
        if default_solution_processed:
            print(f"  [FOUND] Default solution in processed with matching base_id: {default_solution_processed.name}")
            solution_to_use = copy_file_to_working_dir(
                default_solution_processed, puzzle_path, f"default_solution_{base_id}.pdb"
            )
            if solution_to_use:
                files_to_cleanup.append(solution_to_use)
    
    # Jeśli nadal nie ma, spróbuj jakikolwiek solution z processed
    if not solution_to_use:
        any_solution_processed = find_default_solution_in_processed(repo_name, puzzle_name)
        if any_solution_processed:
            print(f"  [FOUND] Any solution in processed: {any_solution_processed.name}")
            solution_to_use = copy_file_to_working_dir(
                any_solution_processed, puzzle_path, f"any_solution.pdb"
            )
            if solution_to_use:
                files_to_cleanup.append(solution_to_use)
    
    # Jeśli nadal nie ma, spróbuj solution z external z tym samym base_id
    if not solution_to_use:
        default_solution_external = find_default_solution_in_external(puzzle_path, base_id)
        if default_solution_external:
            print(f"  [FOUND] Default solution in external with matching base_id: {default_solution_external}")
            solution_to_use = default_solution_external
    
    # Jeśli nadal nie ma, spróbuj jakikolwiek solution z external
    if not solution_to_use:
        any_solution_external = find_default_solution_in_external(puzzle_path)
        if any_solution_external:
            print(f"  [FOUND] Any solution in external: {any_solution_external}")
            solution_to_use = any_solution_external
    
    if not solution_to_use:
        print(f"  [ERROR] No solution found for model {model_file} (base_id: {base_id})")
        # Wyczyść tymczasowe pliki
        for temp_file in files_to_cleanup:
            try:
                temp_path = puzzle_path / temp_file
                if temp_path.exists():
                    temp_path.unlink()
            except:
                pass
        return False
    
    # 4. Uruchom compute.py
    success = run_compute_simple(
        solution_to_use, model_to_use, output_dir, puzzle_path
    )
    
    # 5. Wyczyść tymczasowe pliki
    for temp_file in files_to_cleanup:
        try:
            temp_path = puzzle_path / temp_file
            if temp_path.exists():
                temp_path.unlink()
                # print(f"  [CLEANUP] Removed temporary file: {temp_file}")
        except Exception as e:
            print(f"  [WARNING] Failed to remove {temp_file}: {e}")
    
    return success

def process_puzzle_set_new(puzzle_path, repo_name):
    """
    Nowy flow przetwarzania zestawu puzzli.
    """
    pdb_path = pathlib.Path(puzzle_path) / "pdb"
    if not pdb_path.exists():
        print(f"  [WARNING] No pdb directory found in {puzzle_path}")
        return False
    
    puzzle_name = pdb_path.parent.name
    
    print(f"  [INFO] Processing puzzle: {puzzle_name}")
    print(f"  [INFO] Looking in: {pdb_path}")
    
    # Zbierz modele z external
    models = collect_model_files_from_external(pdb_path)
    
    if not models:
        print(f"  [ERROR] No models found in {pdb_path}")
        return False
    
    print(f"  [INFO] Found {len(models)} model(s)")
    
    # Przygotuj katalog wyników
    output_dir = pathlib.Path("./results/gsp") / repo_name / puzzle_name
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"  [INFO] Results will be saved to: {output_dir}")
    
    # Przetwórz każdy model
    success_count = 0
    for i, model_file in enumerate(models, 1):
        print(f"\n  [{i}/{len(models)}] Processing model: {model_file}")
        
        try:
            success = process_single_model(
                model_file=model_file,
                puzzle_path=pdb_path,
                repo_name=repo_name,
                puzzle_name=puzzle_name,
                output_dir=output_dir
            )
            
            if success:
                success_count += 1
            else:
                print(f"  [FAILED] Failed to process {model_file}")
                
        except Exception as e:
            print(f"  [EXCEPTION] Error processing {model_file}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    # Przenieś pozostałe CSV (jeśli jakieś zostały w pdb_path)
    move_remaining_csv_files(pdb_path, output_dir)
    
    print(f"\n  [SUMMARY] Successfully processed {success_count}/{len(models)} models")
    return success_count > 0

def move_remaining_csv_files(source_dir, dest_dir):
    """Moves any remaining CSV files from source to destination directory."""
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv"):
            source_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(dest_dir, filename)
            if os.path.exists(dest_path):
                print(f"  [OVERWRITE] Overwriting existing file: {filename}")
                os.remove(dest_path)
            try:
                shutil.move(source_path, dest_path)
                print(f"  [MOVE] Moved {filename} to results directory")
            except Exception as e:
                print(f"  [ERROR] Failed to move {filename}: {str(e)}")

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
    print("RNA Puzzle Processing - NEW FLOW")
    print("=" * 50)
    print("Flow: external base_id -> processed lookup -> predef. solution or default")
    
    data_root = "data/external"
    if not os.path.exists(data_root):
        print(f"Error: Data directory '{data_root}' not found!")
        sys.exit(1)
    
    repositories = [d for d in os.listdir(data_root) 
                  if os.path.isdir(os.path.join(data_root, d))]
    if not repositories:
        print("No repositories found in data directory.")
        sys.exit(1)
    
    selected_repos = select_from_list(
        repositories, 
        "Select repositories to process (enter numbers/ranges)"
    )
    if not selected_repos:
        print("No repositories selected. Exiting.")
        sys.exit(0)
    
    for repo_name in selected_repos:
        repo_path = os.path.join(data_root, repo_name)
        puzzles_path = os.path.join(repo_path, "data")
        if not os.path.exists(puzzles_path):
            print(f"\n[WARNING] No 'data' directory found in {repo_name}. Skipping.")
            continue
        
        puzzle_sets = [d for d in os.listdir(puzzles_path) 
                      if os.path.isdir(os.path.join(puzzles_path, d))]
        if not puzzle_sets:
            print(f"\n[WARNING] No puzzle sets found in {repo_name}/data. Skipping.")
            continue
        
        print(f"\nRepository: {repo_name}")
        selected_sets = select_from_list(
            puzzle_sets, 
            f"Select puzzle sets in '{repo_name}' to process"
        )
        if not selected_sets:
            print(f"  [SKIP] No puzzle sets selected for {repo_name}. Skipping.")
            continue
        
        print(f"\nProcessing {len(selected_sets)} puzzle sets in {repo_name}")
        for puzzle_set in selected_sets:
            puzzle_path = os.path.join(puzzles_path, puzzle_set)
            print(f"\n{'='*60}")
            print(f"Puzzle set: {puzzle_set}")
            print('-'*60)
            
            success = process_puzzle_set_new(puzzle_path, repo_name)
            
            if success:
                print(f"\n  [COMPLETED] Finished processing {puzzle_set}")
            else:
                print(f"\n  [FAILED] Failed to process {puzzle_set}")
    
    print("\n" + "="*60)
    print("Processing completed!")

if __name__ == "__main__":
    main()