import os
import sys
import re
import numpy as np

def read_mol2_atom(mol2_path):
    if os.path.exists(mol2_path):
        pass
    else:
        print('ERROR: mol2 file does not exist: %s'%mol2_path)
        raise SystemExit
    mol_atom = {
        'num_atom':0,
        'atm_idx':[],
        'atm_name':[],
        'atm_crd':np.empty((0,3)),
        'atm_type':[],
        'sub_id':[],
        'sub_name':[],
        'atm_charge':[],
        'stat_bit':[]
    }
    with open(mol2_path, 'r') as mol2_file:
        read_flag = False
        for line in mol2_file:
            tmp_mol = line.split()
            if len(tmp_mol) > 0:
                if tmp_mol[0] == '@<TRIPOS>ATOM':
                    read_flag = True
                    continue
                elif re.search('@<TRIPOS>',tmp_mol[0]) and read_flag:
                    read_flag = False
                    break
                if read_flag:
                    tmp_crd = np.array(tmp_mol[2:5])
                    mol_atom['num_atom'] += 1
                    mol_atom['atm_idx'].append(tmp_mol[0])
                    mol_atom['atm_name'].append(tmp_mol[1])
                    mol_atom['atm_crd'] = np.vstack([mol_atom['atm_crd'], tmp_crd])
                    mol_atom['atm_type'].append(tmp_mol[5])
                    mol_atom['sub_id'].append(tmp_mol[6])
                    mol_atom['sub_name'].append(tmp_mol[7])
                    mol_atom['atm_charge'].append(tmp_mol[8])
                    if len(tmp_mol) > 9:
                        mol_atom['stat_bit'].append(tmp_mol[9])
        return(mol_atom)
    
def read_mol2_bond(mol2_path):
    if os.path.exists(mol2_path):
        pass
    else:
        print('ERROR: mol2 file does not exist: %s'%mol2_path)
        raise SystemExit
    mol_bond = {
        'num_bond':0,
        'bond_idx':[],
        'bond':np.empty((0,2)),
        'bond_type':[],
        'stat_bit':[]
    }
    with open(mol2_path, 'r') as mol2_file:
        read_flag = False
        for line in mol2_file:
            tmp_mol = line.split()
            if len(tmp_mol) > 0:
                if tmp_mol[0] == '@<TRIPOS>BOND':
                    read_flag = True
                    continue
                elif re.search('@<TRIPOS>',tmp_mol[0]) and read_flag:
                    read_flag = False
                    break
                if read_flag:
                    tmp_bond = np.array(tmp_mol[1:3])
                    mol_bond['num_bond'] += 1
                    mol_bond['bond_idx'].append(tmp_mol[0])
                    mol_bond['bond'] = np.vstack([mol_bond['bond'], tmp_bond])
                    mol_bond['bond_type'].append(tmp_mol[3])
                    if len(tmp_mol) > 4:
                        mol_bond['stat_bit'].append(tmp_mol[4])
        return(mol_atom)
    
## below need to wirte more !! ##
if __name__ == '__main__':
    # if len(sys.argv) < 5:
    #   print("Usage: %s -i <mol2_file> -j <job_type>"%sys.argv[0])
    #   raise SystemExit
    inp_file = sys.argv[1]
    mol_atom = read_mol2_atom(inp_file)
    mol_bond = read_mol2_bond(inp_file)


            
