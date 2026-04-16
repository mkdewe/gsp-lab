#!/usr/bin/env python3
import sys
import os
import json
import subprocess
import argparse
import pathlib

def run_lddt(target_path, model_path, output_json, project_root):
    target_path = pathlib.Path(target_path).resolve()
    model_path = pathlib.Path(model_path).resolve()
    output_json = pathlib.Path(output_json).resolve()
    root = pathlib.Path(project_root).resolve()

    # Upewnij się, że pliki znajdują się pod root
    try:
        rel_target = target_path.relative_to(root).as_posix()
        rel_model = model_path.relative_to(root).as_posix()
        rel_output = output_json.relative_to(root).as_posix()
    except ValueError as e:
        raise Exception(f"Plik nie znajduje się pod root: {e}")

    container_target = f'/data/{rel_target}'
    container_model = f'/data/{rel_model}'
    container_json = f'/data/{rel_output}'

    output_json.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        'docker', 'run', '--rm',
        '-v', f'{root}:/data',
        'registry.scicore.unibas.ch/schwede/openstructure:latest',
        'compare-structures',
        '-r', container_target,
        '-m', container_model,
        '--lddt',
        '--lddt-no-stereochecks',
        '-o', container_json
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Command: {' '.join(cmd)}", file=sys.stderr)
            print(f"Return code: {result.returncode}", file=sys.stderr)
            print("STDOUT:", file=sys.stderr)
            print(result.stdout, file=sys.stderr)
            print("STDERR:", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            return None
        with open(output_json, 'r') as f:
            data = json.load(f)
        lddt = data.get('lddt')
        if lddt is None:
            print(f"Warning: No lddt field in JSON: {data}", file=sys.stderr)
            return None
        return lddt
    except Exception as e:
        print(f"Exception: {e}", file=sys.stderr)
        return None

def main():
    parser = argparse.ArgumentParser(description='Compute lDDT for a model vs target')
    parser.add_argument('-t', '--target', required=True, help='Target/solution PDB file')
    parser.add_argument('-m', '--model', required=True, help='Model PDB file')
    parser.add_argument('-o', '--output', required=True, help='Output JSON file path')
    parser.add_argument('--root', required=True, help='Project root directory')
    args = parser.parse_args()

    lddt = run_lddt(args.target, args.model, args.output, args.root)
    if lddt is not None:
        print(lddt)
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()