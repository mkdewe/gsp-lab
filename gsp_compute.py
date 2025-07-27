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

def extract_base_id(filename):
    """Extracts the base ID from filename (e.g., Chen_1 from Chen_1.pdb)"""
    # Remove common prefixes and suffixes
    base_name = pathlib.Path(filename).stem
    patterns = [
        r'^(.*?)(?:\.pdb|\.cif|_refined|_solution|_model|_clean|_processed)?$',
        r'^(.*?)_(?:chain|split|processed|refined|model|solution)'
    ]
    
    for pattern in patterns:
        match = re.match(pattern, base_name)
        if match:
            return match.group(1)
    
    # Fallback to first part before underscore
    parts = base_name.split('_')
    if len(parts) > 1:
        return '_'.join(parts[:-1])
    return base_name

def find_processed_file(repo_name, puzzle_name, base_id):
    """Finds a processed file with the same base ID in the processed directory"""
    # Search in the corresponding processed directory
    processed_dir = BASE_PROCESSED_DIR / repo_name / "data" / puzzle_name / "pdb"
    
    if not processed_dir.exists():
        return None
    
    # Look for files with the same base ID
    for candidate in processed_dir.iterdir():
        if candidate.is_file() and candidate.suffix in ['.pdb', '.cif']:
            candidate_id = extract_base_id(candidate.name)
            if candidate_id == base_id:
                return candidate
    
    return None

def find_solution_file(folder_path):
    """Finds a file containing 'solution' in its name within the folder."""
    for fname in os.listdir(folder_path):
        if 'solution' in fname.lower() and pathlib.Path(fname).suffix in ['.pdb', '.cif']:
            return fname
    return None

def collect_model_files(folder_path, solution_file):
    """Collects all model files except the solution file."""
    models = []
    for fname in os.listdir(folder_path):
        if (fname != solution_file and 
            pathlib.Path(fname).suffix in ['.pdb', '.cif'] and
            'solution' not in fname.lower()):
            models.append(fname)
    return models

def run_compute(target_file, model_file, output_dir, cwd):
    """Runs the compute.py script with error handling and fallback"""
    # Determine expected output file name
    model_stem = pathlib.Path(model_file).stem
    target_stem = pathlib.Path(target_file).stem
    output_filename = f"{model_stem}_vs_{target_stem}.csv"
    output_path = pathlib.Path(output_dir) / output_filename
    
    # Check if output file already exists
    if output_path.exists():
        print(f"  [SKIP] Output file {output_filename} already exists. Overwriting...")
        # Remove existing file to ensure fresh results
        output_path.unlink()

    cmd = [
        "python",
        str(COMPUTE_SCRIPT),
        "-t", target_file,
        "-m", model_file
    ]

    env = os.environ.copy()
    # Set OUTPUT_DIR for compute.py to write directly to the results
    env["OUTPUT_DIR"] = os.path.abspath(output_dir)

    # First try with the original file
    result = subprocess.run(cmd, cwd=cwd, env=env)
    
    if result.returncode == 0:
        return True
    
    print(f"  [ERROR] Processing {model_file} (exit code: {result.returncode})")
    
    # Extract base ID for fallback
    base_id = extract_base_id(model_file)
    print(f"  [SEARCH] Looking for processed file with ID: {base_id}")
    
    # Try to find a processed version
    puzzle_name = pathlib.Path(cwd).parent.name
    repo_name = pathlib.Path(cwd).parent.parent.parent.name
    
    processed_file = find_processed_file(repo_name, puzzle_name, base_id)
    
    if processed_file:
        print(f"  [FOUND] Processed file: {processed_file}")
        print(f"  [RETRY] Retrying with processed file...")
        
        # Temporarily copy processed file to working directory
        temp_file = pathlib.Path(cwd) / f"processed_{model_file}"
        shutil.copy(processed_file, temp_file)
        
        # Run compute.py with the processed file
        result = subprocess.run(cmd + ["-m", temp_file.name], cwd=cwd, env=env)
        
        # Clean up temporary file
        temp_file.unlink()
        
        if result.returncode == 0:
            print(f"  [SUCCESS] Successfully processed with fallback file")
            return True
        else:
            print(f"  [ERROR] Error processing fallback file (exit code: {result.returncode})")
    else:
        print(f"  [WARNING] No processed file found for {base_id}")
    
    return False

def process_puzzle_set(puzzle_path, repo_name):
    """Processes a single puzzle set directory with fallback to processed files"""
    pdb_path = os.path.join(puzzle_path, "pdb")
    if not os.path.exists(pdb_path):
        print(f"  [WARNING] No pdb directory found in {puzzle_path}")
        return False

    solution_file = find_solution_file(pdb_path)
    if not solution_file:
        print(f"  [ERROR] No solution file found in {pdb_path}")
        return False

    models = collect_model_files(pdb_path, solution_file)
    if not models:
        print(f"  [ERROR] No models found in {pdb_path}")
        return False

    # Extract puzzle set name
    puzzle_name = os.path.basename(puzzle_path)
    
    # Create output directory
    output_dir = os.path.join("./results/gsp", repo_name, puzzle_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"  [INFO] Results will be saved to: {output_dir}")

    # Process each model
    for i, model_file in enumerate(models, 1):
        print(f"  [MODEL] Processing model {i}/{len(models)}: {model_file}")
        
        success = run_compute(solution_file, model_file, output_dir, pdb_path)
        
        if success:
            print(f"  [SUCCESS] Successfully processed {model_file}")
        else:
            print(f"  [FAILED] Failed to process {model_file} with fallback")
    
    # After processing all models, move any remaining CSV files
    move_remaining_csv_files(pdb_path, output_dir)
    
    return True

def move_remaining_csv_files(source_dir, dest_dir):
    """Moves any remaining CSV files from source to destination directory."""
    for filename in os.listdir(source_dir):
        if filename.endswith(".csv"):
            source_path = os.path.join(source_dir, filename)
            dest_path = os.path.join(dest_dir, filename)
            
            # Nadpisz istniejący plik
            if os.path.exists(dest_path):
                print(f"  [OVERWRITE] Overwriting existing file: {filename}")
                os.remove(dest_path)  # Usuń istniejący plik przed przeniesieniem
                
            try:
                shutil.move(source_path, dest_path)
                print(f"  [MOVE] Moved {filename} to results directory")
            except Exception as e:
                print(f"  [ERROR] Failed to move {filename}: {str(e)}")

def parse_selection(input_str, max_value):
    """Parses user input like '1', '1-5', '1,3,5' into list of indices."""
    selected = set()
    
    # Handle ranges like 1-5
    range_pattern = r'(\d+)-(\d+)'
    for match in re.finditer(range_pattern, input_str):
        start = int(match.group(1))
        end = int(match.group(2))
        if start > end:
            start, end = end, start
        selected.update(range(start, end + 1))
    
    # Handle individual numbers like 1,3,5
    number_pattern = r'\d+'
    for match in re.findall(number_pattern, input_str):
        num = int(match)
        selected.add(num)
    
    # Filter out invalid selections
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
    print("Interactive RNA Puzzle Processing")
    print("=" * 50)
    
    # Step 1: Select repositories
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
    
    # Step 2: Process each selected repository
    for repo_name in selected_repos:
        repo_path = os.path.join(data_root, repo_name)
        puzzles_path = os.path.join(repo_path, "data")
        
        if not os.path.exists(puzzles_path):
            print(f"\n[WARNING] No 'data' directory found in {repo_name}. Skipping.")
            continue
            
        # Get puzzle sets in this repository
        puzzle_sets = [d for d in os.listdir(puzzles_path) 
                      if os.path.isdir(os.path.join(puzzles_path, d))]
        
        if not puzzle_sets:
            print(f"\n[WARNING] No puzzle sets found in {repo_name}/data. Skipping.")
            continue
            
        # Select puzzle sets for this repository
        print(f"\nRepository: {repo_name}")
        selected_sets = select_from_list(
            puzzle_sets, 
            f"Select puzzle sets in '{repo_name}' to process"
        )
        
        if not selected_sets:
            print(f"  [SKIP] No puzzle sets selected for {repo_name}. Skipping.")
            continue
            
        # Process selected puzzle sets
        print(f"\nProcessing {len(selected_sets)} puzzle sets in {repo_name}")
        for puzzle_set in selected_sets:
            puzzle_path = os.path.join(puzzles_path, puzzle_set)
            print(f"\nPuzzle set: {puzzle_set}")
            success = process_puzzle_set(puzzle_path, repo_name)
            
            if success:
                print(f"  [COMPLETED] Finished processing {puzzle_set}")
            else:
                print(f"  [FAILED] Failed to process {puzzle_set}")
    
    print("\nProcessing completed!")

if __name__ == "__main__":
    main()