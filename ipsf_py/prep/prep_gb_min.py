import shutil
import re
import os

def prep_gb_min(ligs, top_path, input_path, min_end_idx, comp_info, num_mpi):
    current_path = os.getcwd()
    #top_path = current_path + "/TOP"
    #input_path = './input_ipsf'
    #min_end_idx = 5
    min_end_idx += 1
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

    if os.path.exists('%s/input_ipsf'%current_path):
        shutil.rmtree('%s/input_ipsf'%current_path)
        shutil.copytree('%s/'%input_path, 'input_ipsf')
    else:
        shutil.copytree('%s/'%input_path, 'input_ipsf')

    with open('./input_ipsf/RUN.template', 'r') as run_tmp:
        run_tmp_contents = run_tmp.read()
        run_tmp_contents = re.sub('%MPI%', mpi_command, run_tmp_contents)
    with open('./input_ipsf/RUN', 'w') as run_in:
        run_in.write(run_tmp_contents)

    for num_min in range(1,min_end_idx):
        # Prepare min input
        with open('./input_ipsf/min%s.in.template'%num_min, 'r') as min_tmp:
            min_tmp_contents = min_tmp.read()
        min_tmp_contents = re.sub('%RECIDX%', '%s  %s'%(recsrt,recend), min_tmp_contents)
        min_tmp_contents = re.sub('%LIGIDX%', '%s  %s'%(ligsrt,ligend), min_tmp_contents)
        with open('./input_ipsf/min%s.in'%num_min, 'w') as min_in:
            min_in.write(min_tmp_contents)


    for lig in ligs:
        if os.path.isdir(lig):
            pass
        else:
            os.mkdir(lig)
        shutil.copy('%s/%s.prmtop'%(top_path, lig), '%s/prmtop'%lig)
        shutil.copy('%s/%s.prmcrd'%(top_path, lig), '%s/prmcrd'%lig)
        shutil.copy('./input_ipsf/RUN', lig)



        

