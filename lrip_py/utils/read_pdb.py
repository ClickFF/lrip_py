import os
import sys
import re
import argparse
import numpy as np
from lrip_py.utils.monitor import _write_log
import subprocess

def read_pdb(pdb_path, log_file=None):
    if os.path.exists(pdb_path):
        pass
    else:
        _write_log(log_file, 'ERROR: PDB file does not exist: %s'%pdb_path)
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


def brf_chk_pdb(pdb_path, log_file=None):
    pdb_path = os.path.abspath(pdb_path)
    if not os.path.exists(pdb_path):
        _write_log(log_file, 'ERROR: PDB file does not exist: %s'%pdb_path)
        sys.exit(1)
    
    has_atom = False
    last_line = ''
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            stripped = line.strip()
            if stripped.startswith('ATOM') or stripped.startswith('HETATM'):
                has_atom = True
            if stripped[0:3] == 'END':
                _write_log(log_file, 'ERROR: PDB file contains END record. It is not suitable for LRIP-SF run: %s'%pdb_path)
                sys.exit(1)
            if stripped:
                last_line = stripped

    if not has_atom:
        _write_log(log_file, 'ERROR: PDB file is empty or does not contain ATOM/HETATM records: %s'%pdb_path)
        sys.exit(1)
    if last_line != 'TER':
        _write_log(log_file, 'ERROR: PDB file does not end with TER record: %s'%pdb_path)
        sys.exit(1)

def clean_pdb(pdb_path, pdb_clean_command, log_file=None):
    import pdbclean
    current_path = os.getcwd()

    if not os.path.exists(pdb_path):
        _write_log(log_file, 'ERROR: PDB file does not exist: %s'%pdb_path)
        sys.exit(1)
    
    # try:
    #     result = subprocess.run(pdb_clean_command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    #     _write_log(log_file, f'pdbclean output:\n{result.stdout}')
    #     if result.stderr:
    #         _write_log(log_file, f'pdbclean error output:\n{result.stderr}')
    # except subprocess.CalledProcessError as e:
    #     _write_log(log_file, f'ERROR: Failed to clean PDB file {pdb_path}: {e}\n{e.stderr}')
    #     sys.exit(1)
    # _write_log(log_file, f'Successfully cleaned PDB file: {pdb_path}. Cleaned file saved as {current_path}/receptor_clean.pdb')

    try:
        pdbclean.pdbclean(pdb_clean_command)
        _write_log(log_file, f'Successfully cleaned PDB file: {pdb_path}. Cleaned file saved as {current_path}/receptor_clean.pdb')
    except Exception as e:
        _write_log(log_file, f'ERROR: Failed to clean PDB file {pdb_path}: {e}')
        sys.exit(1)


def re_write_pdb(pdb_path, out_file, log_file=None):
    """
    Re-write the PDB file to ensure it is in the correct format.
    """
    if os.path.exists(pdb_path):
        pass
    else:
        _write_log(log_file, 'ERROR: PDB file does not exist: %s'%pdb_path)
        sys.exit(1)

    out_file = os.path.abspath(out_file)
    if os.path.exists(out_file):
        _write_log(log_file, 'WARNING: Output PDB file already exists, it will be overwritten: %s'%out_file)
        os.remove(out_file)
    
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            tmp_line = line.strip()
            if len(tmp_line) > 0:
                if  tmp_line[0:4] == 'ATOM' or \
                    tmp_line[0:3] == 'TER':
                    with open(out_file, 'a') as out_pdb:
                        out_pdb.write(tmp_line + '\n')
    with open(out_file, 'a') as out_pdb:
        out_pdb.write('TER\n')
    _write_log(log_file, 'Re-written PDB file saved as: %s'%out_file)


def _get_resi(pdb_path, log_file=None):
    """
    Get the last residue index from the PDB file.
    """
    pdb_path = os.path.abspath(pdb_path)
    if not os.path.exists(pdb_path):
        _write_log(log_file, 'ERROR: PDB file does not exist: %s'%pdb_path)
        sys.exit(1)

    resi_l = []
    with open(pdb_path, 'r') as pdb_file:
        for line in pdb_file:
            if line[0:4] == 'ATOM':
                resi = int(line[20:26].strip())
                resi_l.append(resi)
                try:
                    last_resi = int(resi)
                except ValueError:
                    _write_log(log_file, f'ERROR: Invalid residue index found in PDB file: {resi}')
                    sys.exit(1)
    if len(resi_l) == 0:
        _write_log(log_file, 'ERROR: No ATOM records found in PDB file: %s'%pdb_path)
        sys.exit(1)
    elif len(resi_l) > 0:
        if resi_l[-1] - resi_l[0] < len(resi_l) - 1:
            _write_log(log_file, 'WARNING: Residue indices are not continuous in PDB file: %s'%pdb_path)
    return resi_l[0], resi_l[-1]
                        

