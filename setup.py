import os
import sys
import shutil
from pathlib import Path
from sysconfig import get_paths
from setuptools import setup, Extension

def set_environment_variable(key, value):
  """Set an environment variable."""
  os.environ[key] = value
  print(f"Environment variable {key} set to {value}.")

def check_environment_variable(key):
  """Check if an environment variable exists and return its value."""
  global err
  value = os.environ.get(key)
  if value:
      print(f"Environment variable {key} is set to {value}.")
  else:
      print(f"\033[31mEnvironment variable {key} is not set.\033[0m")
      err += 1

def check_required_executables(executables):
    global err
    missing_executables = []
    for executable in executables:
            whi = shutil.which(executable)
            if whi != None:
                print(f"Executable '{executable}' is found.")
            elif whi == None:
                print(f"\033[31mExecutable '{executable}' is missing.\033[0m")
                missing_executables.append(executable)
                err += 1
    
    if missing_executables:
        print("\n\033[31mMissing executable:")
        for executable in missing_executables:
            print(f" - {executable}")
        print("\033[31mPlease check the above executables.\033[0m\n")

current_path    = os.getcwd()
set_env_var     = {'IPSF_PY_HOME' : current_path}
chk_env_var     = ['SCHRODINGER', 'AMBERHOME']
chk_executables = ['pmemd', 'sander', 'antechamber', 'parmchk2', 'pmemd.MPI', 'pmemd.cuda']
include_root    = get_paths()['include']
err = 0

sr1 = Extension(name='pdbclean', 
                sources=['lrip_py/csrc/pdbclean.c'],
                language='c',
                libraries=['m'],
                library_dirs=[include_root],
                )
sr2 = Extension(name='ie_ave', 
                sources=['lrip_py/csrc/ie_ave.c'],
                language='c',
                libraries=['m'],
                library_dirs=[include_root],
                )


try:
  setup(
      name='lrip-py',
      version='1.0.1',
      description='LRIP-SF: An automated workflow for machine learning-molecular mechanics scoring function. Just for internal use, please do not share.',
      author='Taoyu Niu',
      author_email='tan77@pitt.edu, niutaoyu@gmail.com',
      #license='MIT',
      packages=['lrip_py', 
                'lrip_py.ml', 
                'lrip_py.prep', 
                'lrip_py.utils'],
      package_dir={'lrip_py'             : 'lrip_py',
                  'lrip_py.ml'           : 'lrip_py/ml',
                  'lrip_py.prep'         : 'lrip_py/prep',
                  'lrip_py.utils'        : 'lrip_py/utils',},
      include_package_data=True,
      python_requires='>=3.12',
      install_requires=[  'concurrently',
                          'pandas',
                          'numpy',
                          'scikit-learn==1.7.0',
                          'scipy',
                          'pickledb',
                          'datetime',
                          'psutil==5.9.0',
                        ],
      ext_modules=[sr1, sr2],
      entry_points={
        'console_scripts': [
            'pdbclean   = pdbclean:pdbclean',
            'ie_ave     = ie_ave:ie_ave',
        ]}
  )
except Exception as e:
  print(f"\033[31mError during setup: {e}\033[0m")
  err += 1
  sys.exit(1)


for sev in set_env_var:
  set_environment_variable(sev, set_env_var[sev])

for cev in chk_env_var:
  check_environment_variable(cev)

check_required_executables(chk_executables)

if err > 0:
  print("\033[31mSetup failed. Please check above info.\033[0m")
  sys.exit(1)