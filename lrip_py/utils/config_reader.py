import os
import sys
from lrip_py.utils.monitor import _write_log

def read_dock_parm(glide_in_file, log_file=None):
    read_parm = {}
    #parm_pattern = r'\s*([a-zA-Z_]\w*)\s*=\s*([\w]+)\s*'
    #parm_pattern = r'\s*([\w\d\W]+?)\s*=\s*([\w\d\W]+)\s*'
    if os.path.exists(glide_in_file):
        pass
    else:
        _write_log(log_file, '\nERROR: Docking control file doesn\'t exist, exit: %s'%glide_in_file)
        sys.exit(1)
    with open(glide_in_file, 'r') as list:
        content = list.read()
        entries = content.split(',')
        for entry in entries:
            if '=' in entry:
                parm, value = map(str.strip, entry.split('=', 1))
                read_parm[parm] = value
    return read_parm

def read_list(list_file, log_file=None):
    read_list = {
        'id'        :[],
        'net_charge':[]
    }
    if os.path.exists(list_file):
        pass
    else:
        _write_log(log_file, '\nERROR: Ligands list doesn\'t exist, exit: %s'%list_file)
        sys.exit(1)
    with open(list_file, 'r') as list:
        for line in list:
            tmp_line = line.split(',')
            read_list['id'].append(tmp_line[0])
            if len(tmp_line) < 2:
                read_list['net_charge'].append(str('NaN'))
            else:
                read_list['net_charge'].append(float(tmp_line[1]))

    return read_list

def update_ligs(ligs, miss_list_file=None, miss_lig_list=None, log_file=None):
    """
    Update the ligands list with the ligands that failed in each stage.
    """
    if miss_lig_list is not None and isinstance(miss_lig_list, list):
        if len(miss_lig_list) > 0:
            for lig in miss_lig_list:
                if lig in ligs['id']:
                    i = ligs['id'].index(lig)
                    ligs['id'].pop(i)
                    ligs['net_charge'].pop(i)
                else:
                    _write_log(log_file, "WARNING: Ligand %s not found in the ligands list. It will not be removed."% lig)
    elif miss_list_file is not None and os.path.exists(miss_list_file):
        with open(miss_list_file, 'r') as f:
            for line in f:
                if len(line.strip()) > 0:
                    tmp_line = line.strip().split(',')
                    if tmp_line[0].startswith('#') or len(tmp_line) < 1:
                        continue
                    if tmp_line[0] in ligs['id']:
                        i = ligs['id'].index(tmp_line[0])
                        ligs['id'].pop(i)
                        ligs['net_charge'].pop(i)
                    else:
                        _write_log(log_file, "WARNING: Ligand %s not found in the ligands list. It will not be removed."%tmp_line[0])
                else:
                    continue
    else:
        _write_log(log_file, "ERROR: Missing list file or variable is not provided.")
        sys.exit(1)
    return ligs
