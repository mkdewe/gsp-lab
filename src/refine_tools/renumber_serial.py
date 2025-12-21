#!/usr/bin/env python3
import argparse
import os

def renumber_pdb(start, infile, outfile):
    counter = start
    for line in infile:
        if line.startswith(("ATOM  ","HETATM", "TER  ")):
            outfile.write(f"{line[:6]}{counter:5d}{line[11:]}")
            counter += 1
        else:
            outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--offset", type=int, default=1, help="Starting atom index")
    parser.add_argument("-p", "--postfix", default="_renum", help="Postfix for output filename")
    parser.add_argument("pdbfile", help="Input PDB file")
    args = parser.parse_args()

    base, ext = os.path.splitext(args.pdbfile)
    new_filename = f"{base}{args.postfix}{ext}"

    with open(args.pdbfile, "r") as f_in, open(new_filename, "w") as f_out:
        renumber_pdb(args.offset, f_in, f_out)

    print(f"Saved to: {new_filename}")
