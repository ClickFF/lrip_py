from schrodinger.structure import StructureReader
from schrodinger.application.glide.glide import Dock
from schrodinger.job.jobcontrol import launch_job

# set job name
job_name = "test_docking"

# set input ligand file name for Glide docking
ligand_file = ["d4r_3.mol2", "d4r_5.mol2", "d4r_6.mol2"]  # Ligand files name list

# set glide job commands
"""
command from Glide job log:
/opt/schrodinger/glide glide-dock_SP_1.in 
-OVERWRITE 
-NJOBS 1 
-HOST localhost:128 
-PROJ /data1/tan77/.schrodinger/tmp/tproj35959a856534 
-DISP append 
-VIEWNAME glide_docking_gui.DockingPanel 
-TMPLAUNCHDIR 
-ATTACHED
"""
run_command = ["/opt/schrodinger/glide", '%s.in' % job_name, "-OVERWRITE", "-NJOBS", "1", "-HOST", "localhost:128"]

dock_params = {
    'CALC_INPUT_RMS': 'True',
    'EPIK_PENALTIES': 'False',
    'FORCEPLANAR': 'True',
    'GRIDFILE': '/data6/tan77/gpcr/docking/D4R_5WIU/glide-grid_1/glide-grid_1.zip',  # Replace with the path to your grid file
    'LIGANDFILE': ligand_file, # Input ligand file
    'POSES_PER_LIG':   10,
    'POSTDOCK_NPOSE': 50,
    'PRECISION': 'SP',  # Docking precision: SP or XP
    'REWARD_INTRA_HBONDS': 'True',
}


glide_job = Dock(dock_params)
#print(glide_job)
glide_job.writeSimplified(filename='%s.in' % job_name)
job = launch_job(run_command)
job.wait()
