#!/bin/bash
find . -name "*.csv" -exec rm {} \;
# RNA
python compute.py -t examples/RNA/R1107_D_1292119758_model-annotate_P1human.pdb -m examples/RNA/models
# protein
python compute.py -t examples/protein/T1104-D1.pdb -m examples/protein/models -d 2.0
# RNP
python compute.py -t examples/RNAP/7yr6.pdb -m examples/RNAP/models
