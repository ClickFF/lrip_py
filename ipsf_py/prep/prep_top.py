import shutil
import re
import subprocess
import os
import sys
import tempfile

def top_pre_chk():
    tleap_path = shutil.which('tleap')
    ambpdb_path = shutil.which('ambpdb')
    amberhome = os.environ.get('AMBERHOME')
    if amberhome:
        print('\nAMBERHOME: %s' % amberhome)
    else:
        print('\nERROR: AMBERHOME env veriable not found, exit.')
        sys.exit(1)
    if tleap_path:
        print('\ntleap path: %s'%tleap_path)
    else:
        print('\nERROR: System executable tleap not found, exit.')
        sys.exit(1)
    if ambpdb_path:
        print('\nambpdb path: %s'%ambpdb_path)
    else:
        print('\nERROR: System executable ambpdb not found, exit.')
        sys.exit(1)


def prep_top(ligs, lig_dir_path, receptor_path, tleap_template_path):
############!!!!!!!!!!! this function need a mechanism to detect the tleap error !!!!!!!!!!!!!!###############

    current_path = os.getcwd()
    #ligs={'d4r_3': 0} # Key: ligands ID, value: net charges
    #lig_dir_path = "/data6/tan77/gpcr/docking/test"
    #lig_abcg2_path = current_path + "/LIGFILES"
    #top_path = current_path# + "/TOP"
    #top_log_file = 'prepare_top.log'
    #receptor_path = './D4R_5WIU.pdb' # should be in pdb format
    #tleap_template_path = './tleap.template'

    print('\nGenerate AMBER topology and coordinate start...\n')

    for lig in ligs:
        lig_gaff_path = '%s/%s_gaff2.mol2'%(lig_dir_path, lig)
        lig_frcmod_path = '%s/%s.frcmod'%(lig_dir_path, lig)
        prmtop_path = './%s.prmtop'%lig
        prmcrd_path = './%s.prmcrd'%lig
        tleap_command = [
                        'tleap',\
                        '-s',\
                        '-f',\
                        'tleap_%s.in'%lig
                        ]
        ambpdb_command = [
                        'ambpdb',\
                        '-p',\
                        '%s.prmtop'%lig,\
                        '-c',\
                        '%s.prmcrd',\
                        '>','%s.pdb'%lig
                        ]
        # combine receptor and ligands pdb files
        comp_pdb_path = './comp_%s.pdb'%lig
        with open('comp_%s.pdb'%lig, 'w') as comp:
            with open(receptor_path, 'r') as receptor:
                comp.write(receptor.read())
            with open('%s/%s.pdb'%(lig_dir_path, lig), 'r') as abcg2:
                comp.write(abcg2.read())
        
        #prepare tleap input files for each complex
        shutil.copy(tleap_template_path, './tleap.template')
        with open('tleap.template', 'r') as tleap_template:
            tleap_template_contents = tleap_template.read()
        
        tleap_template_contents = re.sub('%LIGANDGAFF%', lig_gaff_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDFRCMOD%', lig_frcmod_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPDB%', comp_pdb_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPRMTOP%', prmtop_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPRMCRD%', prmcrd_path, tleap_template_contents)
        with open('tleap_%s.in'%lig, 'w') as tleap_in:
            tleap_in.write(tleap_template_contents)
        # run tleap
        top_log1 = subprocess.Popen(tleap_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = top_log1.communicate()
        print('\n-- %s --'%lig)
        if top_log1.returncode == 1: # this error detection mechanism doesn't work for tleap
            print("ERROR: tleap for %s failed."%lig)
            sys.exit(1)

def prep_top_parallel(lig, lig_dir_path, receptor_path, tleap_template_path):

############!!!!!!!!!!! this function need a mechanism to detect the tleap error !!!!!!!!!!!!!!###############
    current_path = os.getcwd()
    with tempfile.TemporaryDirectory() as temp_dir: 
        lig_gaff_path   = os.path.join(lig_dir_path, '%s_gaff2.mol2'%lig)
        lig_frcmod_path = os.path.join(lig_dir_path, '%s.frcmod'%lig)
        prmtop_path     = os.path.join(current_path, '%s.prmtop'%lig)  
        prmcrd_path     = os.path.join(current_path, '%s.prmcrd'%lig)
        comp_pdb_path   = os.path.join(temp_dir, 'comp_%s.pdb'%lig)
        tleap_in_path   = os.path.join(temp_dir, 'tleap_%s.in'%lig)
        tleap_command = [
                        'tleap',
                        '-s',
                        '-f',
                        tleap_in_path
                        ]
        ambpdb_command = [
                        'ambpdb',\
                        '-p',\
                        '%s.prmtop'%lig,\
                        '-c',\
                        '%s.prmcrd',\
                        '>','%s.pdb'%lig
                        ]
        # combine receptor and ligands pdb files
        with open(comp_pdb_path, 'w') as comp:
            with open(receptor_path, 'r') as receptor:
                comp.write(receptor.read())
            with open('%s/%s.pdb'%(lig_dir_path, lig), 'r') as abcg2:
                comp.write(abcg2.read())
        
        #prepare tleap input files for each complex
        shutil.copy(tleap_template_path, '%s/tleap.template'%temp_dir)
        with open(tleap_template_path, 'r') as tleap_template:
            tleap_template_contents = tleap_template.read()
        tleap_template_contents = re.sub('%LIGANDGAFF%', lig_gaff_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDFRCMOD%', lig_frcmod_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPDB%', comp_pdb_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPRMTOP%', prmtop_path, tleap_template_contents)
        tleap_template_contents = re.sub('%LIGANDPRMCRD%', prmcrd_path, tleap_template_contents)
        with open(tleap_in_path, 'w') as tleap_in:
            tleap_in.write(tleap_template_contents)
        # run tleap
        try:
            # Execute commands
            subprocess.run(tleap_command, cwd=temp_dir, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Error processing ligand {lig}: {e}")
            print(f"Command output: {e.stdout.decode()}")
            print(f"Command error: {e.stderr.decode()}")
            return lig, False  # Indicate failure for this ligand
    return lig, True  # Indicate success

        

    
