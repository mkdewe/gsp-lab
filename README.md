--------------------------------------------------------------------------------
GSP-Lab Project: Biomolecular Structure Validation
The GSP-Lab project utilizes the Geometric Spherical Precision (GSP) tool to evaluate the structural validity of proteins, RNA, and DNA by comparing them against reference structures
. The primary goals are to validate molecular models and compare GSP metrics with traditional structural quality measures, such as RMSD, TM-Score, lDDT, and Interaction Network Fidelity (INF_all)
.
Note: The GSP tool is a fork of Prof. Maciej Antczak's original implementation: github.com/mantczak/gsp.

--------------------------------------------------------------------------------
Folder Structure
gsp-lab/
├── src/                    # Source code
│   ├── gsp/                # Core GSP tool (submodule)
│   ├── refine-tools/       # Structure refinement utilities
│   └── utils/              # Data processing utilities
│
├── data/                   # Data collections
│   └── finished/           # Finalized structural data for analysis
│
├── results/                # Computational results
│   ├── correlations/       # Correlation analysis (Spearman, Kendall tau)
│   ├── gsp/                # GSP calculation outputs (gGSP, rGSP, details)
│   ├── lddt/               # Local Distance Difference Test reports
│   ├── inf_all/            # Interaction Network Fidelity metrics
│   └── other-metrics/      # Additional benchmarks (RMSD, TM-Score)
│
├── .vscode/                # IDE configuration
├── gsp_compute.py          # Main computational script
├── gsp_correlation.py      # Correlation analysis and plot generation script
├── requirements.txt        # Python dependencies
├── .gitignore              # Ignored files
├── .gitmodules             # Submodule configuration
└── README.md               # This file

--------------------------------------------------------------------------------

Installation
Clone the repository with submodules:
git clone --recurse-submodules https://github.com/mkdewe/gsp-lab.git
cd gsp-lab
Then create and activate a conda environment (Python 3.8 recommended):
conda create -n gsp-lab-env python=3.8
conda activate gsp-lab-env
pip install -r requirements.txt

--------------------------------------------------------------------------------

Usage
1. Multi-scale GSP Computation
The gsp_compute.py script serves as a high-level wrapper for the core GSP engine. It automates the calculation of quality scores for large structural datasets (e.g., RNA-Puzzles or CASP) using parallel processing.
By default, the tool performs the following:
Calculates scores for multiple RMSD thresholds: 2.0, 4.0, 5.0, and 8.0 Å.
Utilizes a global row-trimming strategy to ensure consistent matrix dimensions across entire puzzle sets.
Generates global scores (gGSP), per-residue scores (rGSP), and detailed residue-radius RMSD matrices.
# Basic execution to process standard puzzle/benchmark directories
python gsp_compute.py
2. Statistical and Correlation Analysis
The gsp_correlation.py script evaluates the relationship between GSP and established metrics like RMSD, TM-score, INF_all, and lDDT. 
It automatically handles metric orientations (where higher is better for GSP, lDDT, and INF, but lower is better for RMSD).
Key features:
Interactive Selection: Allows the user to select specific targets or process "all" available data via a command-line menu.
Data Integration: Automatically merges gGSP results with global lDDT reports and other benchmark CSV files.
Comprehensive Reporting: Produces both simple summaries and detailed reports including p-values, standard deviations, and min/max ranges for Pearson, Spearman, and Kendall correlations.
# Run interactive correlation analysis
python gsp_correlation.py

--------------------------------------------------------------------------------

Output
Output files are systematically organized under results/:
results/gsp/ — GSP scores and detailed residue-radius RMSD matrices.
results/lddt/ — Superposition-free local accuracy estimates.
results/inf_all/ — Fidelity of base-pair and stacking interaction recovery.
results/correlations/ — Statistical analysis (Pearson, Spearman ρ, Kendall τ) and publication-ready plots.


--------------------------------------------------------------------------------

External Metric Calculation Tools
To ensure a consistent and comprehensive comparison across the entire dataset, supplementary calculations for targets with missing metrics were performed using the following specialized tools:
Interaction Network Fidelity (INF): Calculated using MC-Annotate for base-pairing and stacking interaction extraction.
Local Distance Difference Test (lDDT): Computed using the OpenStructure framework (with the no stereochecks option).
TM-score: Calculated using tmtools.

--------------------------------------------------------------------------------

Benchmark Datasets
This repository includes a finalized benchmark set of 4,599 models across 61 targets (tRNA excluded). The data includes:
RNA-Puzzles: Targets PZ1 through PZ39.
CASP15-RNA: Official blind prediction submissions.

--------------------------------------------------------------------------------

License
MIT License © 2025