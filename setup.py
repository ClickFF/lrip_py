#!/usr/bin/env python3

import os
import sys
import importlib
import shutil
from pathlib import Path

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

def check_required_packages(packages):
    """Check if required Python packages are installed."""
    global err
    missing_packages = []
    for package in packages:
        try:
            importlib.import_module(package)
            print(f"Package '{package}' is installed.")
        except ImportError:
            print(f"\033[31mPackage '{package}' is missing.\033[0m")
            missing_packages.append(package)
            err += 1

    if missing_packages:
        print("\n\033[31mMissing packages:\033[0m")
        for package in missing_packages:
            print(f" - {package}")
        print("Please install them using pip:\n")
        print(f"pip install {' '.join(missing_packages)}")

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

def make_script_executable():
    """Ensure the script is executable from the command line."""
    script_path = Path(__file__).resolve()
    print(f"\033[32mTo execute this program from anywhere, you can:\033[0m")
    print(f"1. Add the directory '{script_path.parent}/ipsf_py' to your PATH.")
    print(f"2. Add the environment veriable export IPSF_PY_HOME={current_path} to your bashrc.")
    print(f"3. Run run_ipsf.py -h for help message.")

current_path = os.getcwd()
set_env_var = {'IPSF_PY_HOME' : current_path}
err = 0

if __name__ == "__main__":
    chk_env_var = ['SCHRODINGER', 'AMBERHOME']
    chk_executables = ['pmemd', 'sander', 'antechamber', 'parmchk2', 'pmemd.MPI']
    chk_py_pac  = ['concurrent',
                   'pandas',
                   'numpy',
                   'sklearn',
                   'math',
                   'scipy',
                   'pickle',
                   'datetime']
    for sev in set_env_var:
        set_environment_variable(sev, set_env_var[sev])
    for cev in chk_env_var:
        check_environment_variable(cev)
    
    check_required_executables(chk_executables)

    check_required_packages(chk_py_pac)
    
    if err == 0:
        print("\033[32mSetup success!\033[0m")

        make_script_executable()
    elif err > 0:
        print("\033[31mSetup failed. Please check above info.\033[0m")

    
