import numpy as np
from .read_prmcrd import read_prmcrd
from .read_mol2 import read_mol2_atom
import argparse

def mod_prmcrd_from_mol2(ref_prmcrd_path, mol2_for_crd_path):
    prmcrd = read_prmcrd(ref_prmcrd_path)
    mol2 = read_mol2_atom(mol2_for_crd_path)
    num_atom = mol2['num_atom']
    for i in range(num_atom):
        for j in range(3):
            prmcrd['crd'][-(i+1),j] = mol2['atm_crd'][-(i+1),j]
    return(prmcrd)
    
def write_prmcrd(prmcrd, new_prmcrd_path):
    with open(new_prmcrd_path, 'w') as new_prmcrd:
        new_prmcrd.write('%s\n'%prmcrd['name'])
        new_prmcrd.write('  %d\n'%prmcrd['num_atom'])
        for i in range(0, len(prmcrd['crd']), 2):
            line = []
            for j in range(2):
                if i + j < len(prmcrd['crd']):
                    line += ['%11.7f'%float(value) for value in prmcrd['crd'][i+j]]
            new_prmcrd.write(' '.join(line) + '\n')

def parse_arguments():
    parser = argparse.ArgumentParser(description='Modify ligand coordinates in a protein-ligands .pamcrd file based on a reference .prmcrd file and a ligand .mol2 file')
    parser.add_argument('-rp', '--ref_prmcrd', required=True, help='Reference complex prmtop file')
    parser.add_argument('-rm', '--ref_mol2', required=True, help='Reference ligand mol2 file.')
    parser.add_argument('-o', '--output', required=True, help='Output prmcrd file')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    # ref_prmcrd_path = './test_d4r_10.prmcrd'
    # mol2_for_crd_path = './test_d4r_10.mol2'
    # new_prmcrd_path = './test_mod_d4r_10.prmcrd'
    args = parse_arguments()
    ref_prmcrd_path = args.ref_prmcrd
    mol2_for_crd_path = args.ref_mol2
    new_prmcrd_path = args.output
    prmcrd = mod_prmcrd_from_mol2(ref_prmcrd_path, mol2_for_crd_path)
    write_prmcrd(prmcrd, new_prmcrd_path)