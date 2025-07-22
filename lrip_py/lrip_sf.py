#!/usr/bin/env python3
import sys
import os
import math
import argparse
import subprocess
import json
import ie_ave
import pdbclean
from datetime                          import datetime
from concurrent.futures                import ThreadPoolExecutor
from lrip_py.prep.prep_lig_abcg2       import abcg2_pre_chk, assign_abcg2_parallel
from lrip_py.prep.prep_top             import top_pre_chk, prep_top_parallel
from lrip_py.prep.prep_gb_min          import prep_gb_min
from lrip_py.utils.submit_all          import submit_jobs, ie_parallel
from lrip_py.prep.prep_gb_decomp       import prep_gb_decomp
from lrip_py.utils.get_keyres          import get_keyres
from lrip_py.utils.config_reader       import read_dock_parm, read_list, update_ligs
from lrip_py.utils.restart             import save_state
from lrip_py.utils.resource            import _get_resource
from lrip_py.utils.monitor             import monitor_cpu_usage, md_check, decomp_check, _write_log
from lrip_py.prep.prep_ie              import extract_ie
from lrip_py.ml.ML_pred                import ml_save, ml_load, ml_bs
from lrip_py.utils.read_pdb            import brf_chk_pdb, clean_pdb, re_write_pdb, _get_resi
from lrip_py.utils.calc_mpi            import _calc_mpi

# from docking.dock                      import glide_docking
# from docking.export_str                import export_str


def parse_arguments():
    parser = argparse.ArgumentParser(description="This script will execute Glide docking, molecular dynamic simulations using AMBER, \
                                     and energy decompositions to generate ligand-residue interaction profiles. \
                                     Then a machine learning model will be trained or loaded to predict the ligand-receptor binding free energies.\
                                     ")
    parser.add_argument('-i',   '--input',              required=True, help="Path to the ligands list.")
    parser.add_argument('-c',   '--config',             required=True, help="Path to the config file.")
    parser.add_argument('-p',   '--path2lig',           required=True, help="Path to the directory of ligand files.")
    parser.add_argument('-g',   '--maestro_grid',       required=True, help="Path to the maestro grid files.")
    parser.add_argument('-j',   '--job_type',           required=True, \
                                                        choices=['train', 'pred', 'bs', 'ite', ],\
                                                                       help="Job type to be executed.")
    parser.add_argument('-x',   '--expt_data',                         help="(Optional) Path to experimental data file. This file must have rigorous format (refer to README or example).")
    parser.add_argument('-n',   '--num_ite',                           help="(Optional) Number of iterations. Must be provided if job type is ite.")
    parser.add_argument('-r',   '--restart_file',                      help="(Optional) restart_file file to execute the job from a specific point.")
    parser.add_argument('-l',   '--log_file',                          help="(Optional) Path to the log file to save the output of the job. Default is <jobname>.log in the current directory.")
    parser.add_argument('-v',   '--verbose',            action='store_true', \
                                                                       help="Enable verbose mode to display additional information.")
    args = parser.parse_args()
    return args

def run_lrip_command_line():
    """
    This function is used to run the LRIP-SF workflow from command line.
    It will parse the command line arguments and execute the run_lrip function.
    """
    args                = parse_arguments()
    list_file           = args.input
    config_file         = args.config
    path2lig            = args.path2lig
    grid_file           = args.maestro_grid
    job_type            = args.job_type
    expt_file           = args.expt_data
    num_ite             = args.num_ite
    restart_file        = args.restart_file
    log_file            = args.log_file
    verbose_mode        = args.verbose

    run_lrip(list_file, config_file, path2lig, grid_file, job_type, \
             expt_file, num_ite, restart_file, verbose_mode)


def run_lrip(list_file: str, config_file: str, path2lig: str, grid_file: str, job_type: str, \
             expt_file: str=None, num_ite: int=0, restart_file: str=None, log_file: str=None, \
             over_write: bool=False, verbose_mode: bool=False):
    """
    Description
    ----------
        This function is used to run the LRIP-SF workflow.
        It will read the input file, config file, and other parameters, and execute the workflow.
    Parameters
    ----------
        list_file (str): Path to the ligands list file.
        config_file (str): Path to the config json file.
        path2lig (str): Path to the directory of ligand files.
        grid_file (str): Path to the maestro grid .zip file.
        job_type (str): Job type to be executed. Options: 'train', 'pred', 'bs', 'ite'.
        expt_file (str, optional): Path to experimental data file. Default is None.
        num_ite (int, optional): Number of iterations. Default is 0.
        restart_file (str, optional): Path to the restart file. Default is None.
        log_file (str, optional): Path to the log file. If it is not provided, log info will be print on the screen. Default is None.
        over_write (bool, optional): Whether to overwrite the log file if it exists. Default is False.
        verbose_mode (bool, optional): Enable verbose mode to display additional information. Default is False.
    Returns
    ----------
        0 (if success)
    """
    class SET:
        def __init__(self, list_file, config_file, path2lig, grid_file, job_type, expt_file=None, num_ite=0, restart_file=None, log_file=None, over_write=False, verbose_mode=False):
            self.list_file = list_file
            self.config_file = config_file
            self.path2lig = path2lig
            self.grid_file = grid_file
            self.job_type = job_type
            self.expt_file = expt_file
            self.num_ite = num_ite
            self.restart_file = restart_file
            self.log_file = log_file
            self.over_write = over_write
            self.verbose_mode = verbose_mode
            
            self.current_path = os.getcwd()
            self.schrodinger_home = os.environ.get('SCHRODINGER')
            self.cuda_dev = os.environ.get('CUDA_VISIBLE_DEVICES')
            self.ligs = {}
            self.lig_list = []
            self.lig_list_full_path = []
            self.glide_parm = {}
            self.restart_point = -1
            self.restart_flag = False
            self.state = -1
            self.abcg_flag = True
            self.top_flag = True
            self.failed_ligand = []

    ### Set variables ###
    current_path        = os.getcwd()
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
    failed_ligand       = []

    if log_file:
        log_file = os.path.abspath(log_file)
        
    if over_write and log_file:
        if os.path.exists(log_file):
            os.remove(log_file)

    # Setup and run the job
    _write_log(log_file,'#######################################################\n')
    _write_log(log_file,'                LRIP-SF Job start                      \n')
    _write_log(log_file,'       Start time: %s\n'%datetime.now()                   )
    _write_log(log_file,'#######################################################\n')

    # Read config file
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            config = json.load(f)
    else:
        _write_log(log_file, 'ERROR: Config file does not exist, exit.')
        sys.exit(1)
    
    job_root            = config.get('job_root', 'ITE')
    job_name            = config.get('job_name', 'TST') # job name for glide docking
    njobs               = config.get('njobs', 1) # number of jobs for glide docking
    out_property        = config.get('out_property', 'r_i_glide_gscore') # Output property of from Glide docking
    lig_format          = config.get('lig_format', 'mol2') # Ligand format for abcg2 calculations
    receptor_path       = config.get('receptor_path', None)
    tleap_template_path = config.get('tleap_template_path', '{}/docs/tleap.template'.format(_get_resource()))
    gb_min_in_path      = config.get('gb_min_in_path', '{}/docs/md'.format(_get_resource()))
    gb_min_in_file_list = config.get('gb_min_in_file_list', [
                            'min1.in.template',
                            'min2.in.template',
                            'min3.in.template',
                            'min4.in.template',
                            'min5.in.template',
                        ])
    comp_info           = config.get('comp_info', {
                            'recsrt':1,
                            'recend':123,
                            'ligsrt':124,
                            'ligend':125
                        }) # must read from control file
    num_mpi             = config.get('num_mpi', 1) # num of threads for minimization
    num_parts           = config.get('num_parts', 9999) # the set of ligands will be devided into <num_parts> groups to run minization.
    mpi_threads         = config.get('mpi_threads', 4) # number of threads for MPI
    gb_min_run_command  = config.get('gb_min_run_command', ['bash', './RUN'])
    gb_ene_dec_in_path  = config.get('gb_ene_dec_in_path', '{}/docs/ene_dec/min.in'.format(_get_resource()))
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
    glide_in_file       = config.get('glide_in_file', '{}/docs/glide_docking_control'.format(_get_resource()))
    glide_out_sum_file  = config.get('glide_out_sum_file', 'glide_out.csv') # out file is in csv format
    glide_count_file    = config.get('glide_count_file', 'glide_count.csv')
    split_ratio         = config.get('bs_test_ratio', 0.2)
    n_bs                = config.get('n_bs', 10)
    model_root_path     = config.get('model_root_path', current_path)
    top_pose_for_ml     = config.get('top_pose_for_ml', 1)
    restart_out_file    = config.get('restart_out_file', '%s/%s.json'%(current_path, job_root))
    pre_kres_file       = config.get('pre_kres_file', keyres_file_name)
    usr_def_tleap       = config.get('usr_def_tleap', False)
    pect_used_cpu       = config.get('pect_used_cpu', 0.8)
    usage_threshold     = config.get('cpu_usage_det', 90) # CPU usage threshold to monitor the CPU usage per thread
    miss_list_file      = config.get('miss_list_file', 'dock_fail_ligands.txt') # file to save the ligands that failed to dock
    pdb_clean_command   = config.get('pdb_clean_command', [
                            'pdbclean',
                            '-i', os.path.abspath(receptor_path),
                            '-o', 'receptor_clean.pdb'
                        ])
    clean_pdb_f         = config.get('if_clean_pdb', True)
    
    


    # Pre chk and process of input file existance
    
    if not schrodinger_home:
        _write_log(log_file, 'ERROR: SCHRODINGER not exist, exit.\n')
        sys.exit(1)
    if expt_file:
        expt_file = os.path.abspath(expt_file)
    if job_type == 'train' or job_type == 'pred' or job_type == 'bs':
        if num_ite:
            _write_log(log_file, 'number of iteration will be set = 0 due to the job type: %s\n'%job_type)
        num_ite = 0
    elif job_type == 'ite' and (not num_ite or num_ite == 0):
        _write_log(log_file, 'ERROR: number of iteration (num_ite) doesnt set or <= 0, num_ite = %d. Exit\n'%num_ite)
        sys.exit(1)
    if (job_type == 'train' or job_type == 'bs' or job_type == 'ite') and not expt_file:
        _write_log(log_file, 'ERROR: No experimental value file was provided for this job type.')
        sys.exit(1)
    elif expt_file:
        if not os.path.exists(expt_file):
            _write_log(log_file, 'ERROR: Experimental value file not exist, exit.\n')
            sys.exit(1)
    for gb_min_in_file in gb_min_in_file_list:
        if not os.path.exists('%s/%s'%(gb_min_in_path, gb_min_in_file)):
            _write_log(log_file, 'ERROR: %s/%s does not exist, exit. The <gb_min_in_path> now is: %s. If it is not provided in input json, \
                  please provide this variable in json, pointing to the directory of MD input files.'%(gb_min_in_path, gb_min_in_file, gb_min_in_path))
            sys.exit(1)
        else:
            gb_min_in_file_list[gb_min_in_file_list.index(gb_min_in_file)] = os.path.abspath('%s/%s'%(gb_min_in_path, gb_min_in_file))
    if restart_file:
        restart_flag = True
    else:
        restart_flag = False

    njobs               = int(njobs)
    lig_format          = lig_format.lstrip('.')
    grid_file           = os.path.abspath(grid_file)
    path2lig            = os.path.abspath(path2lig)
    model_root_path     = os.path.abspath(model_root_path)
    if  num_mpi == 'cuda'   or \
        num_mpi == 'gpu'    or \
        num_mpi == 'GPU'    or \
        num_mpi == 'Cuda'   or \
        num_mpi == 'CUDA':
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
        _write_log(log_file, 'ERROR: No receptor pdb file was defined, exit.')
        sys.exit(1)
    else:
        receptor_path = os.path.abspath(receptor_path)
        if clean_pdb_f:
            ### clean receptor pdb file ###
            clean_pdb(receptor_path, pdb_clean_command, log_file=log_file)
            re_write_pdb(current_path + '/receptor_clean.pdb', current_path + '/receptor_clean_re.pdb', log_file=log_file)
            brf_chk_pdb(current_path + '/receptor_clean_re.pdb', log_file=log_file)
            comp_info['recsrt'], comp_info['recend'] = _get_resi(current_path + '/receptor_clean_re.pdb', log_file=log_file)
            comp_info['ligsrt'] = comp_info['recend'] + 1
            comp_info['ligend'] = comp_info['ligsrt'] + 1
            receptor_path = os.path.abspath(current_path + '/receptor_clean_re.pdb')
        else:
            brf_chk_pdb(receptor_path, log_file=log_file)
            comp_info['recsrt'], comp_info['recend'] = _get_resi(receptor_path, log_file=log_file)
            comp_info['ligsrt'] = comp_info['recend'] + 1
            comp_info['ligend'] = comp_info['ligsrt'] + 1

            


    if job_type == 'pred' and model_root_path is None:
        _write_log(log_file, 'ERROR: No model was given for pred job type, exit.')
        sys.exit(1)
    if num_mpi == 'cuda':
        if not cuda_dev:
            _write_log(log_file, 'ERROR: CUDA_VISIBLE_DEVICES not set, exit.')
            sys.exit(1)
    if pre_kres_file:
        pre_kres_file = os.path.abspath(pre_kres_file)
    else:
        _write_log(log_file, 'ERROR: pre_kres_file not defined, exit.')
        sys.exit(1)
    

    ### Load restart file ###
    if restart_flag:
        if os.path.exists(restart_file):
            with open(restart_file, 'r') as f:
                restart = json.load(f)
        else:
            _write_log(log_file, 'ERROR: restart_file file does not exist, exit.')
            sys.exit(1)
        restart_point = restart
    else:
        restart_point = -1


    # Read control
    _write_log(log_file, '\nReading ligands list and net charge...')
    ligs_static = read_list(list_file) # id: str, net_charge: float. read ligands' ID and net charge in a dictionary
    ligs = ligs_static
    i = 0
    for lig in ligs['id']:
        lig_list.append(str(lig))
        lig_list_full_path.append('%s/%s.%s'%(path2lig, lig, lig_format))
        i += 1

    # Calculate optimal number of parallel jobs
    num_ava_cpu = monitor_cpu_usage(usage_threshold)
    if num_ava_cpu <= 0:
        _write_log(log_file, 'Number of availiable CPU cores is %d, not enough for this job, exit.'%num_ava_cpu)
        sys.exit(1)
    num_mpi, num_parts = _calc_mpi(ligs['id'], num_ava_cpu, num_mpi, num_parts, pect_used_cpu, mpi_threads, log_file=log_file)
        
    # Read parameters for GLIDE docking
    _write_log(log_file, '\nReading parameters for GLIDE docking...')
    glide_parm = read_dock_parm(glide_in_file) # read parameters to control GLIDE docking

    # Docking and ABCG2 charge assignment are two stages must be performed for all job types.
    if job_type == 'bs' or job_type == 'ite' or job_type == 'train' or job_type == 'pred':
        num_ite = 0
        if os.path.isdir("%s/%s_%s" % (current_path, job_root, num_ite)):
            _write_log(log_file, '\nWork directory exists: %s/%s_%s'%(current_path, job_root, num_ite))
            pass
        else:
            _write_log(log_file, '\nMaking work directory: %s/%s_%s'%(current_path, job_root, num_ite))
            os.mkdir("%s/%s_%s" % (current_path, job_root, num_ite))


        ##### Glide docking stage #####
        if restart_point <= 1:
            state=1
            _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            _write_log(log_file, '\n\n\nStart ligand docking...')
            if os.path.isdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/DOCKING" % (current_path, job_root, num_ite))
            _write_log(log_file, '\nCurrent directory: %s/%s_%s/DOCKING'%(current_path, job_root, num_ite))

            ### working dir: ./ITE_0/DOCKING ###
            glide_parm['LIGANDFILE'] = lig_list_full_path
            if os.path.exists(grid_file):
                glide_parm['GRIDFILE'] = grid_file
            else:
                _write_log(log_file, 'ERROR: Glide grid file not exist, exit: %s'%grid_file)
                sys.exit(1)
            glide = '%s/glide'%schrodinger_home
            glide_run = [glide, '%s.in' % job_name, "-OVERWRITE", "-NJOBS", str(njobs)]
            json_glide_run = json.dumps(glide_run)
            json_glide_parm = json.dumps(glide_parm)
            with open(f'{job_name}_glide_parm.json', 'w') as f:
                f.write(json_glide_parm)
            # docking_py = '%s/ipsf_py/docking/dock.py'%ipsf_py_home
            docking_py = '{}/docking/dock.py'.format(_get_resource())
            subprocess.run([
                '%s/run'%schrodinger_home, 'python', docking_py,
                '-j', job_name,
                '-r', json_glide_run,
                '-p', f'{job_name}_glide_parm.json'
                #'-p', json_glide_parm
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
            with open(f'{job_name}_lig_list.json', 'w') as f:
                f.write(json_lig_list)
            ll_path = os.path.abspath("%s/%s_%s/DOCKING/%s_lig_list.json" % (current_path, job_root, num_ite, job_name))
            # export_str_py = '%s/ipsf_py/docking/export_str.py'%ipsf_py_home
            export_str_py = '{}/docking/export_str.py'.format(_get_resource())
            if os.path.exists('%s_pv.maegz'%job_name):
                if out_property != 'r_i_glide_gscore':
                    _write_log(log_file, '\nWARNING: The output ligand order is inconsistent with the order of their properties!')
                subprocess.run([
                    '%s/run'%schrodinger_home, 'python', export_str_py,
                    '-i', '%s_pv.maegz'%job_name, 
                    '-o', glide_out_sum_file,
                    '-p', out_property,
                    '-c', glide_count_file,
                    '-ll', ll_path,
                    '-f', lig_format
                ])
                # export_str('%s_pv.maegz'%job_name, glide_out_sum_file, out_property, glide_count_file, ligs['id'], lig_format)
                save_state(state, restart_out_file)
            else:
                _write_log(log_file, '\nERROR: Docking job failed!')
                sys.exit(1)

            ligs = update_ligs(ligs, miss_list_file='../DOCKING/{}'.format(miss_list_file), log_file=log_file)


        ### Assign abcg2 charge ###
        if restart_point < 3:
            state=3
            _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            _write_log(log_file, '\n\n\nAssigning ABCG2 charge and preparing MD inputs...')
            if os.path.isdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            _write_log(log_file, '\nCurrent directory: %s/%s_%s/LIGFILES'%(current_path, job_root, num_ite))
            ### working dir: ./ITE_0/LIGFILES ###
            lig_dir_path_dock = "%s/%s_%s/DOCKING" % (current_path, job_root, num_ite)
            abcg2_pre_chk(log_file=log_file)
            _write_log(log_file, '\nAssign ABCG2 charges start...\n')

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
                    failed_ligand.append(lig)
                #print(f"Ligand {lig} charge assignment {'succeeded' if status else 'failed'}.")
            if abcg_flag:
                save_state(state, restart_out_file)
            elif len(failed_ligand) > 0 and not status:
                _write_log(log_file, f"\nWARNING: The following ligands failed to assign ABCG2 charge: {', '.join(failed_ligand)}\n")
                ligs = update_ligs(ligs, miss_lig_list=failed_ligand)
            failed_ligand = []

        ### Prep top and crd for md ###
        if restart_point < 4:
            state=4
            _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            _write_log(log_file, '\n\n\nPreparing topology and coordinate files...')
            if os.path.isdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite))
            lig_dir_path_abcg = "%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)
            top_pre_chk(log_file=log_file)
            with ThreadPoolExecutor() as executor:
                results = executor.map(
                    lambda lig: prep_top_parallel(lig, lig_dir_path_abcg, receptor_path, tleap_template_path),
                    ligs["id"]
                )
            top_flag = True
            for lig, status in results:
                if not status:
                    top_flag = False
                    failed_ligand.append(lig)
                #print(f"Ligand {lig} top and crd {'succeeded' if status else 'failed'}.")
            if top_flag:
                save_state(state, restart_out_file)
                _write_log(log_file, '\nTopology and coordinate files prepared successfully.\n')
            else:
                _write_log(log_file, f"WARNING: The following ligands failed to prepare topology and coordinate files: {', '.join(failed_ligand)}\n")
                ligs = update_ligs(ligs, miss_lig_list=failed_ligand)
            failed_ligand = []


        ### GB min ###
        if restart_point < 5:
            state=5
            _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            _write_log(log_file, '\n\n\nMD stage...')
            if os.path.isdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite))
            #### working dir: ./ITE_0/GB_MIN ###
            lig_dir_path_abcg = "%s/%s_%s/LIGFILES" % (current_path, job_root, num_ite)
            num_min = prep_gb_min(ligs['id'], lig_dir_path_abcg, gb_min_in_path, gb_min_in_file_list, comp_info, num_mpi, ncyc=100, log_file=log_file)

            if num_mpi != 'cuda':
                submit_jobs(ligs['id'], "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), num_parts, gb_min_run_command, log_file=log_file)
            elif num_mpi == 'cuda':
                submit_jobs(ligs['id'], "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), len(ligs['id']), gb_min_run_command, log_file=log_file)

            #!! check md success or not !!#
            failed_ligand = []
            failed_ligand = md_check(lig_list=ligs['id'], f_idx=len(gb_min_in_file_list), \
                                     dir="%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), \
                                    log_file=log_file)
            
            # re-run failed jobs x2
            if len(failed_ligand) > 0:
                num_ava_cpu = monitor_cpu_usage(usage_threshold)
                if num_ava_cpu <= 0:
                    _write_log(log_file, 'Number of availiable CPU cores is %d, not enough for this job, exit.'%num_ava_cpu)
                    sys.exit(1)
                num_mpi, num_parts = _calc_mpi(ligs['id'], num_ava_cpu, num_mpi, num_parts, pect_used_cpu, mpi_threads, log_file=log_file)
                prep_gb_min(failed_ligand, lig_dir_path_abcg, gb_min_in_path, gb_min_in_file_list, comp_info, num_mpi, ncyc=200, log_file=log_file)
                if num_mpi != 'cuda':
                    submit_jobs(failed_ligand, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), num_parts, gb_min_run_command, log_file=log_file)
                elif num_mpi == 'cuda':
                    submit_jobs(failed_ligand, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), len(failed_ligand), gb_min_run_command, log_file=log_file)

            failed_ligand = []
            failed_ligand = md_check(lig_list=failed_ligand, f_idx=len(gb_min_in_file_list), \
                                     dir="%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), \
                                    log_file=log_file)
            
            # re-run failed jobs x3
            if len(failed_ligand) > 0:
                num_ava_cpu = monitor_cpu_usage(usage_threshold)
                if num_ava_cpu <= 0:
                    _write_log(log_file, 'Number of availiable CPU cores is %d, not enough for this job, exit.'%num_ava_cpu)
                    sys.exit(1)
                num_mpi, num_parts = _calc_mpi(ligs['id'], num_ava_cpu, num_mpi, num_parts, pect_used_cpu, mpi_threads, log_file=log_file)
                prep_gb_min(failed_ligand, lig_dir_path_abcg, gb_min_in_path, gb_min_in_file_list, comp_info, num_mpi, ncyc=300, log_file=log_file)
                if num_mpi != 'cuda':
                    submit_jobs(failed_ligand, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), num_parts, gb_min_run_command, log_file=log_file)
                elif num_mpi == 'cuda':
                    submit_jobs(failed_ligand, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), len(failed_ligand), gb_min_run_command, log_file=log_file)

            failed_ligand = []
            failed_ligand = md_check(lig_list=failed_ligand, f_idx=len(gb_min_in_file_list), \
                                     dir="%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), \
                                    log_file=log_file)

            ligs = update_ligs(ligs, miss_lig_list=failed_ligand)

            if len(failed_ligand) > 0:
                _write_log(log_file, "WARNING: The following ligands failed to run MD: %s\n" % ', '.join(failed_ligand))
            else:
                _write_log(log_file, 'MD run successfully.\n')

            if len(ligs['id']) == 0:
                _write_log(log_file, 'ERROR: No ligands left after MD stage, exit.\n')
                sys.exit(1)

            failed_ligand = []

            save_state(state, restart_out_file)


        ### Prep ene decomposition ###
        if restart_point < 6:
            state=6
            _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
            _write_log(log_file, '\n\n\nEnergy decomposition...')
            if os.path.isdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite)):
                pass
            else:
                os.mkdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite))
            os.chdir("%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite))
            ### working dir: ./ITE_0/GB_DEC ###

            num_ava_cpu = monitor_cpu_usage(usage_threshold)
            if num_ava_cpu <= 0:
                _write_log(log_file, 'Number of availiable CPU cores is %d, not enough for this job, exit.'%num_ava_cpu)
                sys.exit(1)
            num_mpi, num_parts = _calc_mpi(ligs['id'], num_ava_cpu, num_mpi, num_parts, pect_used_cpu, mpi_threads, log_file=log_file)

            prep_gb_decomp(ligs, "%s/%s_%s/GB_MIN" % (current_path, job_root, num_ite), gb_ene_dec_in_path, comp_info, num_min) # error: bad atom type: i, solution: change to amber16
            submit_jobs(ligs['id'], "%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), num_parts, gb_dec_command, log_file=log_file)
            ## add job check function and re-run function ##
            submit_jobs(ligs['id'], "%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), num_parts, ene_dec_command, log_file=log_file)

            with ThreadPoolExecutor() as executor:
                results = executor.map(
                    lambda lig, i: ie_parallel(lig, "%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), ene_dec_command),
                    ligs["id"], range(len(ligs["id"]))
                )
            
            #!! check decomp success or not !!#
            failed_ligand = []
            failed_ligand = decomp_check(lig_list=ligs['id'], dir="%s/%s_%s/GB_DEC" % (current_path, job_root, num_ite), log_file=log_file)
            ligs = update_ligs(ligs, miss_lig_list=failed_ligand, log_file=log_file)

            if len(failed_ligand) > 0:
                _write_log(log_file, "WARNING: The following ligands failed to run energy decomposition: %s\n" % ', '.join(failed_ligand))

            if len(ligs['id']) == 0:
                _write_log(log_file, 'ERROR: No ligands left after energy decomposition stage, exit.\n')
                sys.exit(1)

            for lig in ligs['id']:
                extract_ie('%s/%s_%s/GB_DEC/%s/ie_ave.dat'%(current_path, job_root, num_ite, lig), comp_info, '%s/%s_%s/GB_DEC/%s/%s.ie'%(current_path, job_root, num_ite, lig, lig))

            
            save_state(state, restart_out_file)


        ### ML scorning function ###
        if job_type == 'train':
            if restart_point <= 7:
                state=7
                _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                _write_log(log_file, '\n\n\nJob type: Train and save ML models.')
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
                _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                _write_log(log_file, '\n\n\nJob type: Models validation using bootstrap.')
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

        elif job_type == 'pred':
            if restart_point <= 7:
                state=7
                _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                _write_log(log_file, '\n\n\nJob type: Load pre_trained ML model for predictions.')
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
                _write_log(log_file, '\n\n\nGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad')
                _write_log(log_file, '\n\n\nThis function has not been well established, exit.')
                save_state(state, restart_out_file)
                sys.exit(1)
    
    _write_log(log_file, '#######################################################\n')
    _write_log(log_file, '                   IPSF Job finished                   \n')
    _write_log(log_file, '      Finish time: %s\n'%datetime.now()                   )
    _write_log(log_file, '#######################################################\n')
    return 0
##### Start of the program from command line #####
if __name__ == '__main__':
    run_lrip_command_line()