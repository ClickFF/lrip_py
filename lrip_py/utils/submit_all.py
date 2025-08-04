import os
import subprocess
import shutil
from datetime import datetime
from lrip_py.utils.monitor import _write_log

def submit_jobs(dir_list, work_dir, num_parts, run_command, log_file=None):
    #cwd = os.getcwd()
    arr = [d for d in dir_list if os.path.isdir('%s/%s'%(work_dir,d))]
    elements_per_part = len(arr) // num_parts
    remainder = len(arr) % num_parts

    start = 0
    end = elements_per_part - 1
    new_arrays = {}

    for i in range(1, num_parts + 1):
        if i <= remainder:
            end += 1
        
        new_arrays[i] = arr[start:end + 1]

        #_write_log(log_file, 'Group %s: %s'%(i, new_arrays[i]))

        processes = []
        for mol in new_arrays[i]:
            os.chdir(work_dir)
            os.chdir(mol)
            #print("Working in: %s"%os.getcwd())
            log_file_ = open('log', "w")
            process = subprocess.Popen(
                run_command,                  # Command to execute
                stdout=log_file_,             # Redirect stdout to a log file
                stderr=subprocess.STDOUT      # Redirect stderr to the same log file
            )
            processes.append(process)

        for process in processes:
            process.wait()

        start = end + 1
        end = start + elements_per_part - 1

def ie_parallel(lig, lig_dir_path, commands):
    import ie_ave
    import tempfile
    from glob import glob
    current_path = os.getcwd()
    with tempfile.TemporaryDirectory() as temp_dir:
        os.mkdir('minout')
        files = glob(lig_dir_path + '/' + lig + '/' + 'minout' + '/' + '*')
        #file = lig_dir_path + '/' + lig + '/' + 'minout' + '/' + '*'
        for file in files:
            shutil.copy(file, './minout')

        try:
            # Execute commands
            ie_ave.ie_ave(commands)
            
            # Copy files
            shutil.copy('./ie_ave.dat', f"{current_path}/{lig}/ie_ave.dat")

        except subprocess.CalledProcessError as e:
            return lig, False  # Indicate failure for this ligand
        
    return lig, True  # Indicate success

def submit_abcg2_jobs(ligs, num_parts, lig_format, lig_dir_path, top_pose_for_ml, log_file=None):
    #cwd = os.getcwd()
    current_path = os.getcwd()
    arr = ligs['id']
    elements_per_part = len(arr) // num_parts
    remainder = len(arr) % num_parts

    start = 0
    end = elements_per_part - 1
    new_arrays = {}
    n = 0
    for i in range(1, num_parts + 1):
        if i <= remainder:
            end += 1
        
        new_arrays[i] = arr[start:end + 1]

        _write_log(log_file, "Group %s: %s"%(i, new_arrays[i]))
        _write_log(log_file, 'Start time: %s'%datetime.now())

        processes = []
        for lig in new_arrays[i]:
            lig_file = lig + '_' + str(top_pose_for_ml) + '.' + lig_format # ligands structure file must have .xxx extension for read.
            lig_sybyl_name = "%s_sybyl.%s" % (lig, lig_format)
            lig_abcg2_name = "%s_gaff2.%s" % (lig, 'mol2')
            lig_pdb_name   = "%s.%s"       % (lig, 'pdb')
            lig_frcmod_name= "%s.%s"       % (lig, 'frcmod')

            shutil.copy("%s/%s" % (lig_dir_path, lig_file), "%s/%s" % (current_path, lig_sybyl_name))

            abcg2_command1 = [
                            'antechamber',\
                            '-fi',lig_format,\
                            '-fo','mol2',\
                            '-i','./%s'%lig_sybyl_name,\
                            '-o','./%s'%lig_abcg2_name,\
                            '-c','bcc',\
                            '-at','gaff2',\
                            '-nc',str(ligs['net_charge'][n]),\
                            '-rn','LIG',\
                            '-dr','no'
                            ]
            abcg2_command2 = [
                            'antechamber',\
                            '-fi','mol2',\
                            '-fo','pdb',\
                            '-i','./%s'%lig_abcg2_name,\
                            '-o','./%s'%lig_pdb_name,\
                            '-dr','no'
                            ]
            parmchk2_command = [
                                'parmchk2',\
                                '-f','mol2',\
                                '-i','./%s'%lig_abcg2_name,\
                                '-o','./%s'%lig_frcmod_name,\
                                '-a','Y',\
                                '-s',str(2)
                            ]


            log_file = open('log', "w")
            process = subprocess.Popen(
                abcg2_command1,                  # Command to execute
                stdout=log_file,              # Redirect stdout to a log file
                stderr=subprocess.STDOUT      # Redirect stderr to the same log file
            )
            processes.append(process)
            n += 1

        for process in processes:
            process.wait()

        start = end + 1
        end = start + elements_per_part - 1

        _write_log(log_file, 'End time: %s'%datetime.now())


if __name__ == "__main__":
    dir_str = input("Enter a list of directories working in:")
    dir_list = dir_str.split()
    #dir_list = ['TOP', 'LIGFILES']
    #work_dir = input("Enter the absolute path of the working directory: ")
    work_dir = '/data6/tan77/gpcr/docking/test'
    num_parts = int(input("Enter the number of parallel jobs: "))  # Number of job groups
    run_str = input("Enter the script name to run: ")  # Job script name
    run_command = run_str.split()
    submit_jobs(dir_list, work_dir, num_parts, run_command)
