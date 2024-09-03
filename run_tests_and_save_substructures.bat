del /s /q *.csv
:: delete all 3D substructure PDB files
del /s /q *.0.pdb
:: execute computations and save all 3D substructure PDB files
:: RNA
python -W ignore compute.py -t examples\RNA\R1107_D_1292119758_model-annotate_P1human.pdb -m examples\RNA\models -s
:: protein
python -W ignore compute.py -t examples\protein\T1104-D1.pdb -m examples\protein\models -d 2.0 -s
:: RNP
python -W ignore compute.py -t examples\RNAP\7yr6.pdb -m examples\RNAP\models -s
