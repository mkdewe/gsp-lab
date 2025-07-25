from pathlib import Path

def print_directory_listing(dirs):
    for i, d in enumerate(dirs, 1):
        print(f"{i}. {d.name}")

def parse_selection(selection, max_index):
    selection = selection.replace(" ", "")
    chosen = set()
    for part in selection.split(","):
        if "-" in part:
            start, end = part.split("-")
            chosen.update(range(int(start), int(end) + 1))
        else:
            chosen.add(int(part))
    return [i - 1 for i in chosen if 0 < i <= max_index]

def navigate_and_select_pdbs(start_path=None):
    if start_path is None:
        current_dir = Path.cwd()
    else:
        current_dir = Path(start_path)
    
    current_dir = current_dir.resolve()
    
    while True:
        subdirs = sorted([d for d in current_dir.iterdir() if d.is_dir()])
        pdb_files = sorted([f for f in current_dir.iterdir() if f.is_file() and f.suffix.lower() == ".pdb"])
        
        if pdb_files:
            print(f"\nFound {len(pdb_files)} PDB files in: {current_dir}")
            for i, f in enumerate(pdb_files, 1):
                print(f"{i}. {f.name}")
                
            selection = input("Enter file numbers (e.g. 1,3,5-7) or 'b' to go back: ").strip()
            if selection.lower() == 'b':
                current_dir = current_dir.parent
                continue
            if selection.lower() == 'q':
                return []
            
            try:
                indices = parse_selection(selection, len(pdb_files))
                selected_files = [pdb_files[i] for i in indices]
                return selected_files
            except ValueError:
                print("Invalid selection. Please try again.")
                continue
        
        print(f"\nCurrent directory: {current_dir}")
        print("\nSubdirectories:")
        for i, d in enumerate(subdirs, 1):
            print(f"{i}. {d.name}")
        print("B. Go back")
        print("Q. Quit")
        
        choice = input("Select directory (number), 'B' to go back, or 'Q' to quit: ").strip().lower()
        
        if choice == 'q':
            return []
        if choice == 'b':
            if current_dir.parent != current_dir:  # Prevent going above root
                current_dir = current_dir.parent
            continue
        
        try:
            choice_index = int(choice) - 1
            if 0 <= choice_index < len(subdirs):
                current_dir = subdirs[choice_index]
            else:
                print("Invalid selection. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number, 'B', or 'Q'.")