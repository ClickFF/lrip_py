"""
Utility functions for the LRIP-SF Python package.
This module provides various utility functions for the IPSF package, including environment variable management and executable checks
"""

from lrip_py.utils.config_reader import read_dock_parm, read_list
from lrip_py.utils.get_keyres import get_keyres
from lrip_py.utils.mod_prmcrd import mod_prmcrd_from_mol2, write_prmcrd
from lrip_py.utils.read_mol2 import read_mol2_atom, read_mol2_bond
from lrip_py.utils.read_pdb import read_pdb
from lrip_py.utils.read_prmcrd import read_prmcrd
from lrip_py.utils.resource import _get_resource, _get_root
from lrip_py.utils.submit_all import submit_jobs, submit_abcg2_jobs

__version__ = '1.0.1'
__author__ = 'Taoyu Niu'
__copyright__ = 'Copyright (c) 2025 Taoyu Niu'  

__all__ = [
    'read_dock_parm',
    'read_list',
    'get_keyres',
    'mod_prmcrd_from_mol2',
    'write_prmcrd',
    'read_mol2_atom',
    'read_mol2_bond',
    'read_pdb',
    'read_prmcrd',
    '_get_resource',
    '_get_root',
    'submit_jobs',
    'submit_abcg2_jobs'
]