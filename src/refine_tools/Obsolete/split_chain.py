import os
import sys
import re
from pathlib import Path
import json

# Get current script location
current_dir = Path(__file__).resolve().parent

# Calculate project root (two levels up)
project_root = current_dir.parent.parent

# Add project root to Python path
sys.path.insert(0, str(project_root))

# Now import the navigator module correctly
from src.utils.pdb_navigator import navigate_and_select_pdbs

# Configuration
DEFAULT_OUTPUT_BASE = project_root / "data/processed"
BASE_INPUT_DIR = project_root / "data/external"  # Base directory for input files

# Config file for storing last location
CONFIG_FILE = project_root / ".split_chain_config.json"

def load_config():
    """Load configuration from file if exists"""
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r') as f:
                return json.load(f)
        except json.JSONDecodeError:
            return {}
    return {}

def save_config(config):
    """Save configuration to file"""
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f)

def split_chain(input_file, output_file, chain_id='A', new_chain_id='B', split_point=23):
    """
    Splits a chain in a PDB file at a specified residue number
    and assigns the second part to a new chain. Automatically converts
    undefined chains (non-A/B) to A/B.
    
    Args:
        input_file: Path to input PDB file
        output_file: Path to output PDB file
        chain_id: Chain ID to split (default 'A')
        new_chain_id: New chain ID for the split part (default 'B')
        split_point: Residue number to split at (default 23)
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    output_lines = []
    found_chainB = False
    chain_A_done = False
    auto_convert = chain_id not in ['A', 'B']  # Check if we need auto-conversion
    
    for line in lines:
        if line.startswith('ATOM'):
            current_chain = line[21]
            try:
                resSeq = int(line[22:26].strip())
            except ValueError:
                resSeq = 0
            
            # Handle all cases: auto-convert and standard chains
            if current_chain == chain_id:
                if auto_convert:
                    # Auto-convert undefined chain to A/B
                    if resSeq <= split_point:
                        new_line = line[:21] + 'A' + line[22:]
                        output_lines.append(new_line)
                    else:
                        if not chain_A_done:
                            output_lines.append("TER\n")
                            chain_A_done = True
                        new_resSeq = resSeq - split_point
                        new_resSeq_str = str(new_resSeq).rjust(4)
                        new_line = line[:21] + 'B' + new_resSeq_str + line[26:]
                        output_lines.append(new_line)
                        found_chainB = True
                else:
                    # Standard chain processing (A or B)
                    if resSeq <= split_point:
                        output_lines.append(line)
                    else:
                        if not chain_A_done:
                            output_lines.append("TER\n")
                            chain_A_done = True
                        new_resSeq = resSeq - split_point
                        new_resSeq_str = str(new_resSeq).rjust(4)
                        new_line = line[:21] + new_chain_id + new_resSeq_str + line[26:]
                        output_lines.append(new_line)
                        found_chainB = True
            else:
                output_lines.append(line)
        elif line.startswith('TER'):
            # Preserve TER records for other chains
            output_lines.append(line)
        else:
            output_lines.append(line)
    
    if found_chainB:
        output_lines.append("TER\n")
    
    with open(output_file, 'w') as f:
        f.writelines(output_lines)

def process_files(file_paths, chain_id='A', new_chain_id='B', split_point=23):
    """
    Processes multiple PDB files and saves them to the output directory
    while preserving the full relative path structure.
    
    Args:
        file_paths: List of input file paths (Path objects)
        chain_id: Chain ID to split
        new_chain_id: New chain ID for split part
        split_point: Residue number to split at
    """
    for input_path in file_paths:
        # Calculate relative path from base input directory
        try:
            relative_path = input_path.relative_to(BASE_INPUT_DIR)
        except ValueError:
            print(f"Error: Input file {input_path} is not within the base input directory ({BASE_INPUT_DIR}).")
            continue
        
        # Construct output path preserving full relative structure
        output_path = DEFAULT_OUTPUT_BASE / relative_path
        output_path = output_path.with_name(f"{input_path.stem}_refined.pdb")
        
        # Ensure parent directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing: {input_path} -> {output_path}")
        split_chain(str(input_path), str(output_path), chain_id, new_chain_id, split_point)
    print("\nProcessing completed successfully!")

def main():
    config = load_config()
    
    if len(sys.argv) == 1:  # Interactive mode
        print("=== INTERACTIVE MODE ===")
        
        # Use last location if available
        start_path = config.get('last_location')
        selected_files = navigate_and_select_pdbs(start_path)
        
        if not selected_files:
            print("No files selected. Exiting.")
            return
        
        # Save last location for next time
        if selected_files:
            config['last_location'] = str(selected_files[0].parent)
            save_config(config)
        
        # Get parameters from user
        chain_id = input("Enter chain ID to split (default: A): ").strip() or 'A'
        new_chain_id = input("Enter new chain ID (default: B): ").strip() or 'B'
        split_input = input("Enter split point residue number (default: 23): ").strip()
        split_point = int(split_input) if split_input else 23
        
        process_files(selected_files, chain_id, new_chain_id, split_point)
    
    elif len(sys.argv) >= 3:  # Direct mode
        print("=== DIRECT MODE ===")
        input_file = Path(sys.argv[1])
        output_file = Path(sys.argv[2])
        
        if not input_file.exists():
            print(f"Error: Input file {input_file} does not exist!")
            sys.exit(1)
        
        # Handle optional parameters
        chain_id = sys.argv[3] if len(sys.argv) > 3 else 'A'
        new_chain_id = sys.argv[4] if len(sys.argv) > 4 else 'B'
        split_point = int(sys.argv[5]) if len(sys.argv) > 5 else 23
        
        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        print(f"Processing: {input_file} -> {output_file}")
        split_chain(str(input_file), str(output_file), chain_id, new_chain_id, split_point)
        print("Processing completed!")
    
    else:
        print("Usage:")
        print("  Interactive mode: python split_chain.py")
        print("  Direct mode:      python split_chain.py <input.pdb> <output.pdb> [chain_id] [new_chain_id] [split_point]")
        sys.exit(1)

if __name__ == "__main__":
    main()