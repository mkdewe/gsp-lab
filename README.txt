GSP-Lab Project: Biomolecular Structure Validation
Project Overview
The GSP-Lab project utilizes the Geometric Soundness Potential (GSP) tool to evaluate the structural validity of proteins, RNA, and DNA by comparing them against reference structures. The primary goals are to validate molecular models and compare GSP metrics with traditional structural quality measures (RMSD, TM-Score, etc.).

Note: The GSP tool is a fork of Prof. Maciej Antczak's original implementation:
github.com/mantczak/gsp


Folder Structure

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

Installation
Clone the repository with submodules:

bash
git clone --recurse-submodules https://github.com/your-username/gsp-lab.git
cd gsp-lab
Create and activate Conda environment:

bash
conda create -n gsp-env python=3.9
conda activate gsp-env
Install dependencies:

bash
pip install -r requirements.txt
pip install -r src/gsp/requirements.txt
Configure VS Code environment (optional):

bash
gsp-env_VSCode.bat

Usage
1. Running GSP Calculations
bash
python gsp_compute.py
This interactive script processes structural data and calculates GSP scores.

2. Performing Correlation Analysis
bash
python gsp_correlation.py
Analyze correlations between GSP metrics and other quality measures.

3. Using Refinement Tools
bash
python src/refine-tools/fix_pdb_numbers.py input.pdb
Prepare structural files for analysis with various refinement tools.

4. Exploring Results in Jupyter
bash
jupyter notebook notebook/gsp-jupyter.ipynb
Visualize and explore your analysis results.

Workflow Overview
Prepare Data:

Download datasets: python src/utils/fetch_PZs.py

Preprocess files using tools in src/refine-tools/

Calculate GSP Scores:

bash
python gsp_compute.py
Results saved to results/gsp/

Analyze Correlations:

bash
python gsp_correlation.py
Results saved to results/correlations/

Visualize Results:

Use notebook/gsp-jupyter.ipynb for visualization

Supported Metrics
Metric	Orientation	Description
gGSP	+	Global GSP score
RMSD	-	Root Mean Square Deviation
TM-score	+	Template Modeling score
Clash_Score	-	Atomic clash score
P-value	-	Statistical significance
DI_all	-	Deformation index
INF_all	+	Interaction network fidelity
mcq	+	Model quality estimate
License
This project is licensed under the MIT License - see the LICENSE file for details.

The GSP core tool is a fork of mantczak/gsp, which retains its original MIT license.

Acknowledgments
Prof. Maciej Antczak for the original GSP implementation

RNA Puzzles community for structural data resources

Developers of scientific Python stack (NumPy, SciPy, Pandas)