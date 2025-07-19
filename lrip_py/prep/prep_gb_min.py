import shutil
import re
import os
from lrip_py.utils.monitor import _write_log

def prep_gb_min(ligs, top_path, input_path, min_file_list, comp_info, num_mpi, ncyc=100, log_file=None):
    current_path = os.getcwd()
    input_dir_name = 'input_ipsf'
    run_template_name = 'RUN.template'
    min_file_name_list = []
    #min_end_idx += 1
    recsrt = comp_info['recsrt']
    recend = comp_info['recend']
    ligsrt = comp_info['ligsrt']
    ligend = comp_info['ligend']
    mpi = num_mpi


    if mpi != 'cuda' and mpi > 1:
        mpi_command = 'mpirun -np %s'%mpi + ' ' + 'pmemd.MPI '
        mpi_flag = True
    elif mpi == 'cuda':
        mpi_command = 'pmemd.cuda'
        mpi_flag = False
    else:
        mpi_command = 'pmemd'
        mpi_flag = False

    if os.path.exists('%s/%s'%(current_path, input_dir_name)):
        shutil.rmtree('%s/%s'%(current_path, input_dir_name))
    os.mkdir('%s/%s'%(current_path, input_dir_name))

    shutil.copy('%s/RUN.template' % input_path, '%s/%s/%s' % (current_path, input_dir_name, run_template_name))
    for min_file in min_file_list:
        min_file_name = os.path.basename(min_file)
        min_file_name_list.append(min_file_name)
        if not os.path.exists(min_file):
            _write_log(log_file, 'ERROR: %s does not exist, exit.' % min_file)
            exit(1)
        shutil.copy('%s' % min_file, '%s/%s/%s' % (current_path, input_dir_name, min_file_name))

    with open('%s/%s/%s'%(current_path, input_dir_name, run_template_name), 'r') as run_tmp:
        run_tmp_contents = run_tmp.read()
        run_tmp_contents = re.sub('%MPI%', mpi_command, run_tmp_contents)
    with open('%s/%s/RUN'%(current_path, input_dir_name), 'w') as run_in:
        run_in.write(run_tmp_contents)

    num_min = 1
    for min_file in min_file_name_list:
        # Prepare min input
        with open('%s/%s/%s'%(current_path, input_dir_name, min_file), 'r') as min_tmp:
            min_tmp_contents = min_tmp.read()
        min_tmp_contents = re.sub('%RECIDX%', '%s  %s'%(recsrt,recend), min_tmp_contents)
        min_tmp_contents = re.sub('%LIGIDX%', '%s  %s'%(ligsrt,ligend), min_tmp_contents)
        min_tmp_contents = re.sub('%NCYC%'  , '%s'%ncyc,                min_tmp_contents)
        if os.path.exists('%s/%s/min%s.in'%(current_path, input_dir_name, num_min)):
            print('ERROR: min%s.in is repeated with internal file name, please use other name, such as min<x>.inp.'% num_min)
            exit(1)
        with open('%s/%s/min%s.in'%(current_path, input_dir_name, num_min), 'w') as min_in:
            min_in.write(min_tmp_contents)
        num_min += 1


    for lig in ligs:
        if os.path.isdir(lig):
            pass
        else:
            os.mkdir(lig)
        shutil.copy('%s/%s.prmtop'%(top_path, lig), '%s/prmtop'%lig)
        shutil.copy('%s/%s.prmcrd'%(top_path, lig), '%s/prmcrd'%lig)
        shutil.copy('./input_ipsf/RUN', lig)

    return num_min - 1



        

