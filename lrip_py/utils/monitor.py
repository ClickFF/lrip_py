import psutil
import time
import math
import os


def monitor_cpu_usage(usage_threshold=70, log_file=None):
    """
    Monitor CPU usage per thread and print the usage.
    """
    if usage_threshold > 0 and usage_threshold <= 1:
        _write_log(log_file, "Warning: cpu_usage_det %f is too small, set to fraction by times 10." % usage_threshold)
        usage_threshold = math.floor(100 * usage_threshold)
        
    # Wait a moment to get accurate readings
    time.sleep(1)

    cpu_usages = psutil.cpu_percent(percpu=True)

    count = 0
    for i, usage in enumerate(cpu_usages):
        if usage < usage_threshold:
            count += 1
        else:
            pass
    #     print(f"CPU Thread {i}: {usage}% used")
    # print(f"Total CPU Threads with usage > 20%: {count}")

    return count

def md_check(lig_list=None,f_idx=None, dir=None, log_file=None):
    """
    Check if the MD jobs were successful.
    Returns a list of ligands that failed.
    """
    failed_ligand = []
    if lig_list is None:
        _write_log(log_file, "No ligand list provided for MD check.")
        return failed_ligand
    if f_idx is None:
        _write_log(log_file, "No file index provided for MD check.")
        return failed_ligand
    dir = os.path.abspath(dir)
    for lig in lig_list:
        # Check if the MD output file exists
        if isinstance(f_idx, int):
            md_output_file = "%s/%s/min%s.out"%(dir, lig, f_idx)
            md_output_file = os.path.abspath(md_output_file)
            if not os.path.exists(md_output_file):
                _write_log(log_file, "MD output file for %s does not exist. Adding to failed list."%lig)
                failed_ligand.append(lig)
            else:
                with open(md_output_file, 'r') as f:
                    lines = f.readlines()
                if not lines or \
                    ("|  Total wall time:"           not in lines[-1] and \
                    '|  Master Total wall time:'    not in lines[-1]):
                    _write_log(log_file, "MD output file for %s does not end with '|  Total wall time:'. Adding to failed list."%lig)
                    failed_ligand.append(lig)
    return failed_ligand

def decomp_check(lig_list=None, dir=None, log_file=None):
    """
    Check if the decomp jobs were successful.
    """
    failed_ligand = []
    if lig_list is None:
        _write_log(log_file, "No ligand list provided for decomp check.")
        return failed_ligand
    dir = os.path.abspath(dir)
    for lig in lig_list:
        # Check if the decomp output file exists
        decomp_output_file = f"{dir}/{lig}/minout/minout_1.out"
        if not os.path.exists(decomp_output_file):
            _write_log(log_file, f"Decomp output file for {decomp_output_file} does not exist. Adding to failed list.")
            failed_ligand.append(lig)
        else:
            with open(decomp_output_file, 'r') as f:
                lines = f.readlines()
            if not lines or "|Total" not in lines[-1]:
                _write_log(log_file, f"Decomp output file for {decomp_output_file} does not end with 'Decomposition completed successfully'. Adding to failed list.")
                failed_ligand.append(lig)
    return failed_ligand

def _write_log(log_file, message):
    """
    Write a message to the log file.
    """
    if log_file:
        with open(log_file, 'a') as f:
            f.write(f"{message}\n")
    else:
        print(message)  # Fallback to print if no log file is provided