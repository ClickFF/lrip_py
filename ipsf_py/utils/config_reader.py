import os
import sys

def read_dock_parm(glide_in_file):
    read_parm = {}
    #parm_pattern = r'\s*([a-zA-Z_]\w*)\s*=\s*([\w]+)\s*'
    #parm_pattern = r'\s*([\w\d\W]+?)\s*=\s*([\w\d\W]+)\s*'
    if os.path.exists(glide_in_file):
        pass
    else:
        print('ERROR: Docking control file doesn\'t exist, exit: %s'%glide_in_file)
        sys.exit(1)
    with open(glide_in_file, 'r') as list:
        content = list.read()
        entries = content.split(',')
        for entry in entries:
            if '=' in entry:
                parm, value = map(str.strip, entry.split('=', 1))
                read_parm[parm] = value
    return read_parm

def read_list(list_file):
    read_list = {
        'id'        :[],
        'net_charge':[]
    }
    if os.path.exists(list_file):
        pass
    else:
        print('ERROR: Ligands list doesn\'t exist, exit: %s'%list_file)
        sys.exit(1)
    with open(list_file, 'r') as list:
        for line in list:
            tmp_line = line.split(',')
            read_list['id'].append(tmp_line[0])
            read_list['net_charge'].append(float(tmp_line[1]))

    return read_list
