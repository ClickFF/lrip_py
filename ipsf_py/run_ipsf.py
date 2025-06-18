#!/usr/bin/env python3
import sys
import os
import argparse
import subprocess
import json
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
# from docking.dock import glide_docking
# from docking.export_str import export_str
from prep.prep_lig_abcg2 import assign_abcg2, abcg2_pre_chk, assign_abcg2_parallel
from prep.prep_top import prep_top, top_pre_chk, prep_top_parallel
from prep.prep_gb_min import prep_gb_min
from utils.submit_all import submit_jobs
from prep.prep_gb_decomp import prep_gb_decomp
from utils.get_keyres import get_keyres
from utils.config_reader import read_dock_parm, read_list
from utils.restart import save_state
from prep.prep_ie import extract_ie
from ml.ML_pred import ml_save, ml_load, ml_bs

def parse_arguments():
    parser = argparse.ArgumentParser(description="This script will execute Glide docking, molecular dynamic simulations using AMBER, \
                                     and energy decompositions to generate ligand-residue interaction profiles. \
                                     Then a machine learning model will be trained or loaded to predict the ligand-receptor binding free energies.\
                                     (Haven't done 1/1/2025)This script also provides an option to train the model first, and iteratively perform docking pose screening.")
    parser.add_argument('-i',   '--input',              required=True, help="Path to the ligands list.")
    parser.add_argument('-c',   '--config',             required=True, help="Path to the config file.")
    parser.add_argument('-p',   '--path2lig',           required=True, help="Path to the directory of ligand files.")
    parser.add_argument('-g',   '--maestro_grid',       required=True, help="Path to the maestro grid files.")
    parser.add_argument('-j',   '--job_type',           required=True, \
                                                        choices=['train', 'pred', 'bs', 'ite', ],\
                                                                       help="Job type to be executed.")
    parser.add_argument('-x',   '--expt_data',                         help="(Optional)Path to experimental data file. This file must have rigorous format (refer to README or example).")
    parser.add_argument('-n',   '--num_ite',                           help="(Optional) Number of iterations. Must be provided if job type is ite.")
    parser.add_argument('-r',   '--restart_file',                      help="(Optional) restart_file file to execute the job from a specific point.")
    parser.add_argument('-v',   '--verbose',            action='store_true', \
                                                                       help="Enable verbose mode to display additional information.")
    args = parser.parse_args()
    return args



def main():
    ### Set variables ###
    current_path        = os.getcwd()
    ipsf_py_home        = os.environ.get('IPSF_PY_HOME')
    schrodinger_home    = os.environ.get('SCHRODINGER')
    cuda_dev            = os.environ.get('CUDA_VISIBLE_DEVICES')
    ligs                = {}
    lig_list            = []
    lig_list_full_path  = []
    glide_parm          = {}
    restart_point       = -1
    restart_flag        = False
    state               = -1
    abcg_flag           = True
    top_flag            = True

    args                = parse_arguments()
    list_file           = args.input
    config_file         = args.config
    path2lig            = args.path2lig
    grid_file           = args.maestro_grid
    job_type            = args.job_type
    expt_file           = args.expt_data
    num_ite             = args.num_ite
    restart_file        = args.restart_file
    verbose_mode        = args.verbose

    if expt_file:
        expt_file = os.path.abspath(expt_file)
    if not ipsf_py_home:
        print('ERROR: IPSF_PY_HOME not exist, exit.')
        sys.exit(1)
    if not schrodinger_home:
        print('ERROR: SCHRODINGER not exist, exit.')
        sys.exit(1)
    if job_type == 'train' or job_type == 'pred' or job_type == 'bs':
        if num_ite:
            print('\nnumber of iteration will be set = 0 due to the job type: %s'%job_type)
        num_ite = 0
    elif job_type == 'ite' and (not num_ite or num_ite == 0):
        print('ERROR: number of iteration (num_ite) doesnt set or <= 0, num_ite = %d. Exit'%num_ite)
        sys.exit(1)
    if (job_type == 'train' or job_type == 'bs' or job_type == 'ite') and not expt_file:
        print('ERROR: No experimental value file was provided for this job type.')
        sys.exit(1)
    elif expt_file:
        if not os.path.exists(expt_file):
            print('ERROR: Experimental value file not exist, exit.')
            sys.exit(1)
    if restart_file:
        restart_flag = True
    else:
        restart_flag = False


    # Read config file
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
    else:
        print('ERROR: Configuration file does not exist, exit.')
        sys.exit(1)
    job_root            = config.get('job_root', 'ITE')
    job_name            = config.get('job_name', 'TST') # job name for glide docking
    njobs               = config.get('njobs', 1) # number of jobs for glide docking
    out_property        = config.get('out_property', 'r_i_glide_gscore') # Output property of from Glide docking
    lig_format          = config.get('lig_format', 'mol2') # Ligand format for abcg2 calculations
    receptor_path       = config.get('receptor_path', None)
    tleap_template_path = config.get('tleap_template_path', '%s/ipsf_py/docs/tleap.template'%ipsf_py_home)
    gb_min_in_path      = config.get('gb_min_in_path', '%s/ipsf_py/docs/md'%ipsf_py_home)
    min_end_idx         = config.get('min_end_idx', 5)
    comp_info           = config.get('comp_info', {
                            'recsrt':1,
                            'recend':123,
                            'ligsrt':124,
                            'ligend':125
                        }) # must read from control file
    num_mpi             = config.get('num_mpi', 4) # num of threads for minimization
    num_parts           = config.get('num_parts', 3) # the set of ligands will be devided into <num_parts> groups to run minization.
    gb_min_run_command  = config.get('gb_min_run_command', ['bash', './RUN']) 
    gb_ene_dec_in_path  = config.get('gb_ene_dec_in_path', '%s/ipsf_py/docs/ene_dec/min.in'%ipsf_py_home)
    gb_dec_command      = config.get('gb_dec_command',  [
                            #sander -O -i min.in -o minout/minout_1.out -p prmtop -c inpcrd
                            'sander',
                            '-O',
                            '-i', 'min.in',
                            '-o', 'minout/minout_1.out',
                            '-p', 'prmtop',
                            '-c', 'inpcrd'
                        ])
    ene_dec_command     = config.get('ene_dec_command', [
                            #ie_ave -d minout -i minout -o ie_ave.dat -s 1 -e 1 -st 1
                            'ie_ave', 
                            '-d', 'minout',
                            '-i', 'minout',
                            '-o', 'ie_ave.dat',
                            '-s', '1',
                            '-e', '1',
                            '-st', '1'
                        ])
    keyres_file_name    = config.get('keyres_file_name', 'keyres.csv') 
    ml_model_list       = config.get('ml_model_list', [1, 2, 3, 4, 5, 6, 7, 8]) # must be a list of integers
    glide_in_file       = config.get('glide_in_file', '%s/ipsf_py/docs/glide_docking_control'%ipsf_py_home)
    glide_out_sum_file  = config.get('glide_out_sum_file', 'glide_out.csv') # out file is in csv format
    glide_count_file    = config.get('glide_count_file', 'glide_count.csv')
    split_ratio         = config.get('split_ratio', 0.2)
    n_bs                = config.get('n_bs', 10)
    model_root_path     = config.get('model_root_path', None)
    top_pose_for_ml     = config.get('top_pose_for_ml', 1)
    restart_out_file    = config.get('restart_out_file', '%s/%s.json'%(current_path, job_root))
    pre_kres_file       = config.get('pre_kres_file', keyres_file_name)
    


    # Pre chk and process of input file existance
    njobs               = int(njobs)
    lig_format          = lig_format.lstrip('.')
    min_end_idx         = int(min_end_idx)
    grid_file           = os.path.abspath(grid_file)
    path2lig            = os.path.abspath(path2lig)
    model_root_path     = os.path.abspath(model_root_path)
    if num_mpi == 'cuda' or num_mpi == 'gpu' or num_mpi == 'GPU' or num_mpi == 'Cuda' or num_mpi == 'CUDA': 
        num_mpi         = 'cuda'
    else:
        num_mpi         = int(num_mpi)
    num_parts           = int(num_parts)
    if restart_flag:
        restart_file    = os.path.abspath(restart_file)
    for i in comp_info:
        comp_info[i]    = int(comp_info[i])
    c=0
    for i in ml_model_list:
        ml_model_list[c]= int(ml_model_list[c])
        c+=1
    if not glide_out_sum_file.endswith(".csv"):
        glide_out_sum_file += ".csv"
    if not glide_count_file.endswith(".csv"):
        glide_count_file += ".csv"
    if not receptor_path:
        print('ERROR: No receptor pdb file was defined, exit.')
        sys.exit(1)
    else:
        receptor_path = os.path.abspath(receptor_path)
    if job_type == 'pred' and model_root_path is None:
        print('ERROR: No model was given for pred job type, exit.')
        sys.exit(1)
    if num_mpi == 'cuda':
        if not cuda_dev:
            print('ERROR: CUDA_VISIBLE_DEVICES not set, exit.')
            sys.exit(1)
    if pre_kres_file:
        pre_kres_file = os.path.abspath(pre_kres_file)
    else:
        print('ERROR: pre_kres_file not defined, exit.')
        sys.exit(1)
    

    ### Load restart file ###
    if restart_flag:
        if os.path.exists(restart_file):
            with open(restart_file, 'r') as f:
                restart = json.load(f)
        else:
            print('ERROR: restart_file file does not exist, exit.')
            sys.exit(1)
        restart_point = restart
    else:
        restart_point = -1
        
    

    # Read control
    print('\nReading ligands list and net charge...')
    ligs = read_list(list_file) # id: str, net_charge: float. read ligands' ID and net charge in a dictionary
    i = 0
    for lig in ligs['id']:
        lig_list.append(str(lig))
        lig_list_full_path.append('%s/%s.%s'%(path2lig, lig, lig_format))
        i += 1
    print('\nReading parameters for GLIDE docking...')
    glide_parm = read_dock_parm(glide_in_file) # read parameters to control GLIDE docking
    

    time_file = open(f'{current_path}/time_record.txt', 'w')

    # Setup and run the job
    print('#######################################################\n')
    print('                   IPSF Job start                      \n')
    print('       Start time: %s\n'%datetime.now()                   )
    print('#######################################################\n')

    time_file.write(f'Start time: {datetime.now()}\n')

    # Docking and ABCG2 charge assignment are two stages must be performed for all job types.
    if job_type == 'bs' or job_type == 'ite' or job_type == 'train' or job_type == 'pred':
        num_ite = 0
        if os.path.isdir("%s/%s_%s" % (current_path, job_root, num_ite)):
            print('\nWork directory exists: %s/%s_%s'%(current_path, job_root, num_ite))
            pass
        else:
            print('\nMaking work directory: %s/%s_%s'%(current_path, job_root, num_ite))
            os.mkdir("%s/%s_%s" % (current_path, job_root, num_ite))




        ##### Glide docking stage #####
        if restart_point <= 1:
            state=1
            print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            print('\n\n\nStart ligand docking...')
            if os.path.isdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            print('\nCurrent directory: %s/%s_%s/DOCKING'%(current_path, job_root, num_ite))

            ### working dir: ./ITE_0/DOCKING ###
            glide_parm['LIGANDFILE'] = lig_list_full_path
            if os.path.exists(grid_file):
                glide_parm['GRIDFILE'] = grid_file
            else:
                print('ERROR: Glide grid file not exist, exit: %s'%grid_file)
                sys.exit(1)
            glide = '%s/glide'%schrodinger_home
            glide_run = [glide, '%s.in' % job_name, "-OVERWRITE", "-NJOBS", str(njobs)]
            json_glide_run = json.dumps(glide_run)
            json_glide_parm = json.dumps(glide_parm)
            docking_py = '%s/ipsf_py/docking/dock.py'%ipsf_py_home
            subprocess.run([
                '%s/run'%schrodinger_home, 'python', docking_py,
                '-j', job_name,
                '-r', json_glide_run,
                '-p', json_glide_parm
            ])
            save_state(state, restart_out_file)

            ### Function for integrated python run ###
            # glide_docking(job_name, glide_run, glide_parm) # could be used directly but if the schrodinger version is old, some commonly used module 

        ### Extract structure from docking output ### 
        if restart_point < 2:
            state=2
            if os.path.isdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            json_lig_list = json.dumps(ligs['id'])
            export_str_py = '%s/ipsf_py/docking/export_str.py'%ipsf_py_home
            if os.path.exists('%s_pv.maegz'%job_name):
                if out_property != 'r_i_glide_gscore':
                    print('\nWARNING: The output ligand order is inconsistent with the order of their properties!')
                subprocess.run([
                    '%s/run'%schrodinger_home, 'python', export_str_py,
                    '-i', '%s_pv.maegz'%job_name, 
                    '-o', glide_out_sum_file,
                    '-p', out_property,
                    '-c', glide_count_file,
                    '-l', json_lig_list,
                    '-f', lig_format
                ])
                # export_str('%s_pv.maegz'%job_name, glide_out_sum_file, out_property, glide_count_file, ligs['id'], lig_format)
                save_state(state, restart_out_file)
                time_file.write(f'End of docking: {datetime.now()}\n')
            else:
                print('\nERROR: Docking job failed!')
                sys.exit(1)
            

        ### Assign abcg2 charge ###
        if restart_point < 3:
            state=3
            print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            print('\n\n\nAssigning ABCG2 charge and preparing MD inputs...')
            if os.path.isdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            print('\nCurrent directory: %s/%s_%s/LIGFILES'%(current_path, job_root, num_ite))
            ### working dir: ./ITE_0/LIGFILES ###
            lig_dir_path_dock = "%s/%s_%s/DOCKING" % (current_path, job_root, num_ite)
            abcg2_pre_chk()
            print('\nAssign ABCG2 charges start...\n')

            ### Single job ###
            # assign_abcg2(ligs, lig_format, lig_dir_path_dock, top_pose_for_ml) # single job

            ### Parallel job ###
            with ThreadPoolExecutor() as executor:
                results = executor.map(
                    lambda lig, i: assign_abcg2_parallel(lig, lig_format, lig_dir_path_dock, ligs, i, top_pose_for_ml),
                    ligs["id"], range(len(ligs["id"]))
                )
            abcg_flag = True
            for lig, status in results:
                if not status:
                    abcg_flag = False
                print(f"Ligand {lig} charge assignment {'succeeded' if status else 'failed'}.")
            if abcg_flag:
                save_state(state, restart_out_file)
            time_file.write(f'End of ABCG2 charge assignment: {datetime.now()}\n')

        ### Prep top and crd for md ###
        if restart_point < 4:
            state=4
            print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            print('\n\n\nPreparing topology and coordinate files...')
            if os.path.isdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            lig_dir_path_abcg = "%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)
            top_pre_chk()
            with ThreadPoolExecutor() as executor:
                results = executor.map(
                    lambda lig: prep_top_parallel(lig, lig_dir_path_abcg, receptor_path, tleap_template_path),
                    ligs["id"]
                )
            top_flag = True
            for lig, status in results:
                if not status:
                    top_flag = False
                print(f"Ligand {lig} top and crd {'succeeded' if status else 'failed'}.")
            #prep_top(lig_list, lig_dir_path_abcg, receptor_path, tleap_template_path)
            if top_flag:
                save_state(state, restart_out_file)
            time_file.write(f'End of topology and coordinate files preparation: {datetime.now()}\n')
        #exit(1)   
        ### GB min ###
        if restart_point < 5:
            state=5
            print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            print('\n\n\nMD stage...')
            if os.path.isdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite))
            #### working dir: ./ITE_0/GB_MIN ###
            lig_dir_path_abcg = "%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)
            prep_gb_min(lig_list, lig_dir_path_abcg, gb_min_in_path, min_end_idx, comp_info, num_mpi)
            #exit(1)
            if num_mpi != 'cuda':
                submit_jobs(ligs['id'], "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), num_parts, gb_min_run_command)
            elif num_mpi == 'cuda':
                submit_jobs(ligs['id'], "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), len(ligs['id']), gb_min_run_command)
            #!! check md success or not !!#
            save_state(state, restart_out_file)
            time_file.write(f'End of MD: {datetime.now()}\n')


        ### Prep ene decomposition ###
        if restart_point < 6:
            state=6
            print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            print('\n\n\nEnergy decomposition...')
            if os.path.isdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite))
            ### working dir: ./ITE_0/GB_DEC ###
            prep_gb_decomp(ligs, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), gb_ene_dec_in_path, comp_info, min_end_idx) # error: bad atom type: i, solution: change to amber16
            submit_jobs(ligs['id'], "%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), num_parts, gb_dec_command)
            ## add job check function and re-run function ##
            submit_jobs(ligs['id'], "%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), num_parts, ene_dec_command)
            for lig in ligs['id']:
                extract_ie('%s/%s_%s/GB_DEC/%s/ie_ave.dat'%(current_path, job_root, num_ite, lig), comp_info, '%s/%s_%s/GB_DEC/%s/%s.ie'%(current_path, job_root, num_ite, lig, lig))
            save_state(state, restart_out_file)
            time_file.write(f'End of energy decomposition: {datetime.now()}\n')
        #exit(1)

        ### ML scorning function ###
        if job_type == 'train':
            if restart_point <= 7:
                state=7
                print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                print('\n\n\nJob type: Train and save ML models.')
                if os.path.isdir("%s/%s_%s/TRN" % (current_path, job_root, num_ite)):
                    pass
                else:
                    os.mkdir("%s/%s_%s/TRN" % (current_path, job_root, num_ite))
                os.chdir("%s/%s_%s/TRN" % (current_path, job_root, num_ite))
                ### working dir: ./ITE_0/TRN ###
                tmp_file_list = []
                for lig in ligs['id']:
                    tmp_file_list.append("%s/%s_%s/GB_DEC/%s/%s.ie" % (current_path, job_root, num_ite, lig, lig))
                get_keyres(ligs, tmp_file_list, keyres_file_name)
                ml_save("%s/%s_%s/TRN/%s" % (current_path, job_root, num_ite, keyres_file_name), expt_file, ml_model_list) # model_list elements must be integers
                save_state(state, restart_out_file)

        elif job_type == 'bs':
            if restart_point <= 7:
                state=7
                print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                print('\n\n\nJob type: Models validation using bootstrap.')
                if os.path.isdir("%s/%s_%s/BS" % (current_path, job_root, num_ite)):
                    pass
                else:
                    os.mkdir("%s/%s_%s/BS" % (current_path, job_root, num_ite))
                os.chdir("%s/%s_%s/BS" % (current_path, job_root, num_ite))   
                ### working dir: ./ITE_0/BS ###
                tmp_file_list = []
                for lig in ligs['id']:
                    tmp_file_list.append("%s/%s_%s/GB_DEC/%s/%s.ie" % (current_path, job_root, num_ite, lig, lig))
                get_keyres(ligs, tmp_file_list, keyres_file_name)
                ml_bs("%s/%s_%s/BS/%s" % (current_path, job_root, num_ite, keyres_file_name), expt_file, ml_model_list, split_ratio, n_bs) 
                # model_list elements must be integers, split_ratio must be [0.0,1.0], n_bs must be integer
                save_state(state, restart_out_file)
                time_file.write(f'End of bootstrap: {datetime.now()}\n')
        elif job_type == 'pred':
            if restart_point <= 7:
                state=7
                print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                print('\n\n\nJob type: Load pre_trained ML model for predictions.')
                if os.path.isdir("%s/%s_%s/PRE" % (current_path, job_root, num_ite)):
                    pass
                else:
                    os.mkdir("%s/%s_%s/PRE" % (current_path, job_root, num_ite))
                os.chdir("%s/%s_%s/PRE" % (current_path, job_root, num_ite)) 
                ### working dir: ./ITE_0/PRE ###
                tmp_file_list = []
                for lig in ligs['id']:
                    tmp_file_list.append("%s/%s_%s/GB_DEC/%s/%s.ie" % (current_path, job_root, num_ite, lig, lig))
                get_keyres(ligs, tmp_file_list, keyres_file_name, pre_kres_file=pre_kres_file)
                ml_load("%s/%s_%s/PRE/%s" % (current_path, job_root, num_ite, keyres_file_name), ml_model_list, model_root_path, expt_path=expt_file)
                # expt_path could be the path to experimental value file for validation, and it could also be None or not provided to mute the validation stage.
                save_state(state, restart_out_file)
        elif job_type == 'ite':
            if restart_point <= 7:
                state=7
                print('\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                print('\n\n\nThis function has not been well established, exit.')
                save_state(state, restart_out_file)
                sys.exit(1)
    
    print('#######################################################\n')
    print('                   IPSF Job finished                   \n')
    print('      Finish time: %s\n'%datetime.now()                   )
    print('#######################################################\n')
    time_file.write(f'End time: {datetime.now()}\n')
    time_file.close()
##### Start of Main Program #####
if __name__ == '__main__':
    main()


        
        

    
    
