import os 
import sys
import math
from lrip_py.utils.monitor import _write_log


def _calc_mpi(elements, num_ava_cpu, num_mpi, num_parts, pect_used_cpu, mpi_threads, log_file=None):
    act_num_cpu = math.floor(num_ava_cpu * pect_used_cpu)
    i = len(elements)
    if i > 0:
        if num_mpi == 'cuda':
            num_parts = len(elements)
        else:
            if num_mpi == 1 and num_parts == 9999:
                if act_num_cpu > mpi_threads:
                    num_mpi = mpi_threads
                    num_parts = math.ceil(i / (act_num_cpu / num_mpi))
                    if num_parts < 1:
                        num_parts = 1
            elif math.ceil(i / num_parts) * num_mpi > act_num_cpu:
                _write_log(log_file, 'WARNING: The mpi settings are too large for the available CPU cores, re-assign the number of jobs.')
                num_mpi = mpi_threads
                num_parts = math.ceil(i / (act_num_cpu / num_mpi))
                if num_parts < 1:
                        num_parts = 1

    return num_mpi, num_parts