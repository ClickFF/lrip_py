import os
from lrip_py import __file__

def _get_resource():
    base = os.path.dirname(os.path.realpath(__file__))

    return base

def _get_root():
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))