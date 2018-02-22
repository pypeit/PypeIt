""" Module to deal with meta level routines for PYPIT"""
from __future__ import absolute_import, division, print_function, unicode_literals

from pkg_resources import resource_filename
import glob


def instr_list():
    settings_path = resource_filename('pypit', 'data/settings')
    settings_files = glob.glob(settings_path+'/settings.*')
    instruments = [ifile.split('.')[1] for ifile in settings_files]  # This includes base
    # Trim
    for kpop in ['armed', 'armlsd', 'baseargflag', 'basespect', 'py']:
        instruments.remove(kpop)
    # Return
    return instruments
