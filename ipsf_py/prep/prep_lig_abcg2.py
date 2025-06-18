import shutil
import subprocess
import os
import sys
import tempfile

def abcg2_pre_chk():
    antechamber_path = shutil.which('antechamber')
    parmchk_path = shutil.which('parmchk2')
    amberhome = os.environ.get('AMBERHOME')
    if amberhome:
        print('\nAMBERHOME: %s' % amberhome)
    else:
        print('\nERROR: AMBERHOME env veriable not found, exit.')
        sys.exit(1)
    if antechamber_path:
        print('\nANTECHAMBER path: %s'%antechamber_path)
    else:
        print('\nERROR: System executable antechamber not found, exit.')
        sys.exit(1)
    if parmchk_path:
        print('\nParmchk2 path: %s'%parmchk_path)
    else:
        print('\nERROR: System executable parmchk2 not found, exit.')
        sys.exit(1)

def assign_abcg2(ligs, lig_format, lig_dir_path, top_pose_for_ml):

    current_path = os.getcwd()
    print('\nJOB: Assign ABCG2 charges start...\n')
    
    # could be paralleled to save time
    i = 0
    for lig in ligs['id']:
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
                        '-c','abcg2',\
                        '-at','gaff2',\
                        '-nc',str(ligs['net_charge'][i]),\
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
        
        # For python3.5+

        # abcg2_log1 = subprocess.run(abcg2_command1, capture_output=True, text=True)
        # with open('./LIGFILES/%s'%abcg2_log_file, 'a') as log_file:
        #     log_file.write('-- %s --'%lig)
        #     log_file.write(abcg2_log1.stdout)

        # abcg2_log2 = subprocess.run(abcg2_command2, capture_output=True, text=True)
        # with open('./LIGFILES/%s'%abcg2_log_file, 'a') as log_file:
        #     log_file.write(abcg2_log2.stdout)
        
        # abcg2_log3 = subprocess.run(parmchk2_command, capture_output=True, text=True)
        # with open('./LIGFILES/%s'%abcg2_log_file, 'a') as log_file:
        #     log_file.write(abcg2_log3.stdout)
        #     log_file.write('\n\n')


        # for python2.7 released with schrodinger2017
        abcg2_log1 = subprocess.Popen(abcg2_command1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = abcg2_log1.communicate()
        # with open('./%s'%abcg2_log_file, 'a') as log_file:
        # log_file.write('-- %s --\n'%lig)
        # log_file.write(stdout)
        # log_file.write(stderr)
        print('\n-- %s --'%lig)
        print('\nABCG2 STDOUT: \n%s'%stdout)
        print('\nABCG2 STDERR: \n%s'%stderr)
        if abcg2_log1.returncode == 1:
            print("\nERROR: Assign charge for %s failed. Inconsistant net charge would be the most possible issue."%lig)
            sys.exit(1)

        abcg2_log2 = subprocess.Popen(abcg2_command2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = abcg2_log2.communicate()
        # with open('./%s'%abcg2_log_file, 'a') as log_file:
        # log_file.write(stdout)
        # log_file.write(stderr)
        print('\nmol2 to pdb STDOUT: \n%s'%stdout)
        print('\nmol2 to pdb STDERR: \n%s'%stderr)
        if abcg2_log2.returncode == 1:
            print("\nERROR: Convert %s from mol2 to pdb failed."%lig)
            sys.exit(1)
        
        abcg2_log3 = subprocess.Popen(parmchk2_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = abcg2_log3.communicate()
        # with open('./%s'%abcg2_log_file, 'a') as log_file:
        # log_file.write(stdout)
        # log_file.write(stderr)
        print('\nGenerate GAFF2 STDOUT: \n%s'%stdout)
        print('\nGenerate GAFF2 STDERR: \n%s'%stderr)
        if abcg2_log2.returncode == 1:
            print("\nERROR: Generate GAFF2 for %s failed."%lig)
            sys.exit(1)
        i += 1


def assign_abcg2_parallel(lig, lig_format, lig_dir_path, ligs, i, top_pose_for_ml):
    current_path = os.getcwd()
    with tempfile.TemporaryDirectory() as temp_dir:
        #########################################################################################
        #print(os.getcwd()) # ??!!!!????(#U&$)(#!@)*# without this line the parallel job will fail
        #########################################################################################
        lig_file = lig + '_' + str(top_pose_for_ml) + '.' + lig_format  # Ligand structure file must have .xxx extension for read.
        lig_sybyl_name = os.path.join(temp_dir, "%s_sybyl.%s"%(lig,lig_format))
        lig_abcg2_name = os.path.join(temp_dir, "%s_gaff2.mol2"%lig)
        lig_pdb_name = os.path.join(temp_dir, "%s.pdb"%lig)
        lig_frcmod_name = os.path.join(temp_dir, "%s.frcmod"%lig)

        shutil.copy(f"{lig_dir_path}/{lig_file}", lig_sybyl_name)

        abcg2_command1 = [
            'antechamber',
            '-fi', lig_format,
            '-fo', 'mol2',
            '-i', lig_sybyl_name,
            '-o', lig_abcg2_name,
            '-c', 'abcg2',
            '-at', 'gaff2',
            '-nc', str(ligs['net_charge'][i]),
            '-rn', 'LIG',
            '-dr', 'no',
            '-pf', 'yes'
        ]
        abcg2_command2 = [
            'antechamber',
            '-fi', 'mol2',
            '-fo', 'pdb',
            '-i', lig_abcg2_name,
            '-o', lig_pdb_name,
            '-dr', 'no',
            '-pf', 'yes'
        ]
        parmchk2_command = [
            'parmchk2',
            '-f', 'mol2',
            '-i', lig_abcg2_name,
            '-o', lig_frcmod_name,
            '-a', 'Y',
            '-s', str(2)
        ]

        try:
            # Execute commands
            subprocess.run(abcg2_command1, cwd=temp_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subprocess.run(abcg2_command2, cwd=temp_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            subprocess.run(parmchk2_command, cwd=temp_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Copy files
            shutil.copy(lig_abcg2_name, f"{current_path}/{os.path.basename(lig_abcg2_name)}")
            shutil.copy(lig_pdb_name, f"{current_path}/{os.path.basename(lig_pdb_name)}")
            shutil.copy(lig_frcmod_name, f"{current_path}/{os.path.basename(lig_frcmod_name)}")

            #print(f"Processing of ligand {lig} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error processing ligand {lig}: {e}")
            return lig, False  # Indicate failure for this ligand
    return lig, True  # Indicate success

##### Unfinished part for assigning ABCG2 charge parallelly.
# Parallel execution
# with ThreadPoolExecutor() as executor:
#     results = executor.map(
#         lambda lig, i: process_ligand(lig, lig_format, lig_dir_path, current_path, ligs, i),
#         ligs["id"], range(len(ligs["id"]))
#     )

# Print results
# for lig, status in results:
#     print(f"Ligand {lig} processing {'succeeded' if status else 'failed'}.")

if __name__ == '__main__':
    ligs={'d4r_3': 0} # Key: ligands ID, value: net charges
    lig_format = "mol2" # should be mol2, otherwise schrodinger can not recongize.
    lig_dir_path = "/data6/tan77/gpcr/docking/test"
    abcg2_log_file = 'prepare_prepi.log'
    assign_abcg2(ligs, lig_format, lig_dir_path, abcg2_log_file)
    



    




    
    




