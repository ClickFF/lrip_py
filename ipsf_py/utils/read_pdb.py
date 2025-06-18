import os
import sys
import re
import argparse
import numpy as np

def read_pdb(pdb_path):
    if os.path.exists(pdb_path):
        pass
    else:
        print('ERROR: PDB file does not exist: %s'%pdb_path)
        sys.exit(1)
    mol_atom = {
        'num_mol': 0,
        'num_atom':[],
        'atm_idx':[],
        'atm_name':[],
        'res_name':[],
        'res_idx':[],
        'atm_crd':np.empty((0,3)) # this one is for temptory use
    }
    crystal=np.empty((0,6))
    n_mol = 0
    pre_mol_idx = -1
    with open(pdb_path, 'r') as pdb_file:
        read_flag = False
        for line in pdb_file:
            tmp_line = line.split()
            if len(tmp_line) > 0:
                if tmp_line[0] == 'CRYST1' and len(tmp_line) == 8:
                    crystal = np.array(tmp_line[1:8])
                if tmp_line[0] == 'ATOM' and len(tmp_line) == 11:
                    tmp_crd = np.array(tmp_line[5:8])
                    mol_atom['atm_idx'].append(tmp_line[4])
                    mol_atom['atm_name'].append(tmp_line[2])
                    mol_atom['res_name'].append(tmp_line[3])
                    mol_atom['res_idx'].append(tmp_line[1])
                    mol_atom['atm_crd'] = np.vstack([mol_atom['atm_crd'], tmp_crd])
                    if int(tmp_line[4]) != pre_mol_idx:
                        n_mol += 1
                        pre_mol_idx = int(tmp_line[4])
                        mol_atom['num_atom'].append(n_atm)
                        n_atm = 0
                    else:
                        n_atm += 1
        mol_atom['num_mol'] = n_mol

    return crystal, mol_atom
