import argparse

def replace_star_in_pdb(input_file, output_file=None):
    with open(input_file, "r") as f:
        lines = f.readlines()

    output_lines = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")):
            atom_name = line[12:16]
            if "*" in atom_name:
                atom_name = atom_name.replace("*", "'")
                line = line[:12] + atom_name.ljust(4) + line[16:]
        output_lines.append(line)

    if output_file is None:
        output_file = input_file

    with open(output_file, "w") as f:
        f.writelines(output_lines)

    print(f"Saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Replace '*' with \"'\" in PDB atom names")
    parser.add_argument("input", help="Path to input PDB file")
    parser.add_argument("-o", help="Path to output PDB file (optional)")
    args = parser.parse_args()

    replace_star_in_pdb(args.input, args.o)
