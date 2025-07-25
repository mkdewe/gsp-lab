# GSP-Lab Project: Biomolecular Structure Validation

The **GSP-Lab** project utilizes the **Geometric Spherical Precision (GSP)** tool to evaluate the structural validity of proteins, RNA, and DNA by comparing them against reference structures. The primary goals are to validate molecular models and compare GSP metrics with traditional structural quality measures (RMSD, TM-Score, etc.).

> **Note**: The GSP tool is a fork of Prof. Maciej Antczak's original implementation: [github.com/mantczak/gsp](https://github.com/mantczak/gsp)

---

## Folder Structure

```
gsp-lab/
├── src/                    # Source code
│   ├── gsp/                # Core GSP tool (submodule)
│   ├── refine-tools/       # Structure refinement utilities
│   └── utils/              # Data processing utilities
│
├── data/                   # Data collections
│   ├── external/           # Raw data (Git submodules)
│   └── processed/          # Processed data for analysis
│
├── results/                # Computational results
│   ├── correlations/       # Correlation analysis results
│   ├── gsp/                # GSP calculation outputs
│   └── other-metrics/      # External metrics for comparison
│
├── notebook/               # Analytical notebooks
│   └── gsp-jupyter.ipynb   # Main analysis notebook
│
├── .vscode/                # IDE configuration
│   └── settings.json
│
├── gsp_compute.py          # Main computational script
├── gsp_correlation.py      # Correlation analysis script
├── gsp-env_VSCode.bat      # VS Code environment setup
├── requirements.txt        # Python dependencies
├── .gitignore              # Ignored files
├── .gitmodules             # Submodule configuration
└── README.md               # This file
```

---

## Installation

Clone the repository with submodules:

```bash
git clone --recurse-submodules https://github.com/yourusername/gsp-lab.git
cd gsp-lab
```

Then create and activate a conda environment:

```bash
conda create -n gsp-lab-env python=3.8
conda activate gsp-lab-env
pip install -r requirements.txt
```

---

## Usage

### Example GSP computation:

```bash
python gsp_compute.py \
    --target_path data/processed/reference_structures/1ABC.pdb \
    --model_path data/processed/models/1ABC_model.pdb
```

### Example correlation analysis:

```bash
python gsp_correlation.py \
    --gsp_dir results/gsp/ \
    --metrics_dir results/other-metrics/ \
    --output_dir results/correlations/
```

---

## Output

Output files are saved under `results/`:

- `results/gsp/` — GSP scores and detailed matrices
- `results/other-metrics/` — external validation scores (e.g., RMSD, TM-Score)
- `results/correlations/` — correlation values and plots comparing GSP to other metrics

---

## License

MIT License © 2025