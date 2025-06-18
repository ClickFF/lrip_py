import os
import sys
import numpy as np

def read_prmcrd(crd_path):
    if os.path.exists(crd_path):
        pass
    else:
        print('ERROR: the prmcrd file does not exist: %s'%crd_path)
        raise SystemExit
    prmcrd = {
        'name':'',
        'num_atom':0,
        'crd':np.empty((0,3))
    }
    with open(crd_path, 'r') as crd_file:
        num_atom = 0
        name_flag = False
        num_atom1_flag = False
        for line in crd_file:
            crd_tmp = line.split()
            if len(crd_tmp) == 1:
                if not name_flag and not num_atom1_flag:
                    prmcrd['name'] = crd_tmp[0]
                    name_flag = True
                    continue
                if name_flag and not num_atom1_flag:
                    num_atom = int(crd_tmp[0])
                    num_atom1_flag = True
                    continue
            elif name_flag and num_atom1_flag and (len(crd_tmp) == 6 or len(crd_tmp) == 3):
                if len(crd_tmp) == 6:
                    tmp_coord = np.array(crd_tmp[0:3])
                    prmcrd['crd'] = np.vstack([prmcrd['crd'], tmp_coord])
                    tmp_coord = np.array(crd_tmp[3:6])
                    prmcrd['crd'] = np.vstack([prmcrd['crd'], tmp_coord])
                    prmcrd['num_atom'] += 2
                if len(crd_tmp) == 3:
                    tmp_coord = np.array(crd_tmp[0:3])
                    prmcrd['crd'] = np.vstack([prmcrd['crd'], tmp_coord])
                    prmcrd['num_atom'] += 1
        if num_atom == prmcrd['num_atom']:
            pass
        else:
            print('ERROR: prmcrd file inconsistant, number of atoms: %s'%crd_path)
            raise SystemExit
    return(prmcrd)



if __name__ == '__main__':
    if len(sys.argv) < 2:
      print("Usage: %s <parmcrd>"%sys.argv[0])
      raise SystemExit
    inp_file = sys.argv[1]
    prmcrd = read_prmcrd(inp_file)
    print('Number of atoms of this prmcrd file: %s'%prmcrd['num_atom'])
    print('Shape of output x,y,z array: %s'%str(np.shape(prmcrd['crd'])))
    # add what u want