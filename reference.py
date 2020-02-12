import numpy as np
from feature import *
from rdkit import Chem


with open('NR-AR-LBD_wholetraining.smiles') as f:
    smi_list = f.readlines()

smi_list = [smi.split(' ')[0] for smi in smi_list]
out_list = []

for smi in smi_list:
    mol = Chem.MolFromSmiles(smi)
    int_vec = mol_to_feature(mol, -1, 400)
    vec = [float(value) for value in int_vec]
    out_list.append(str(vec))

with open('reference_output.csv', 'w') as f:
    f.write('\n'.join(out_list))

