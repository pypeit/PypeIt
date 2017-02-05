""" Module for examining/archiving/etc. settings files.
This cannot be located in the pypit/ folder (import issues)
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
from glob import glob
import numpy as np
from os.path import dirname

from pypit import pyputils
msgs = pyputils.get_dummy_logger()#develop=True)
from pypit import arparse

# CANNOT LOAD DEBUGGER AS THIS MODULE IS CALLED BY ARDEBUG
import pdb as debugger

try:
    basestring
except NameError:
    basestring = str

# Logging
#msgs = pyputils.get_dummy_logger(develop=True)


def compare_dicts(top_key, dict1, dict2, skip_keys=()):
    for key in dict1.keys():
        if isinstance(dict1[key], dict):
            if key in dict2.keys():
                compare_dicts(top_key+','+key, dict1[key], dict2[key])
        else:
            try:
                test = dict1[key] == dict2[key]
            except KeyError:
                pass
            else:
                if test:
                    if key not in skip_keys:
                        msgs.info("{:s},{:s} is a duplicate".format(top_key,key))
                else:
                    msgs.warn("{:s},{:s} is different".format(top_key,key))


def argf_diff_and_dup():
    """ Compares default argf values against those in the ARMLSD and AMRED files
    Returns
    -------

    """
    # Load default argf file
    baseargf = arparse.get_argflag_class(('BaseArgFlag', '.tmp'))
    base_lines = baseargf.load_file()
    baseargf.set_paramlist(base_lines)

    # ARMLSD
    msgs.info("===============================================")
    msgs.info("Working on ARMLSD vs. Base")
    armlsd = arparse.get_argflag_class(('ARMLSD', '.tmp'))
    armlsd_lines = armlsd.load_file()
    armlsd.set_paramlist(armlsd_lines)
    # Look for duplicates and diffs
    for key in baseargf._argflag.keys():
        if key in armlsd._argflag.keys():
            compare_dicts(key, baseargf._argflag[key], armlsd._argflag[key])

    # ARMED
    msgs.info("===============================================")
    msgs.info("Working on ARMED vs. Base")
    armed = arparse.get_argflag_class(('ARMED', '.tmp'))
    armed_lines = armed.load_file()
    armed.set_paramlist(armed_lines)
    # Look for duplicates and diffs
    for key in baseargf._argflag.keys():
        if key in armed._argflag.keys():
            compare_dicts(key, baseargf._argflag[key], armed._argflag[key])


def spect_diff_and_dup():
    import pypit
    # Load default spect file
    path = pypit.__path__[0]+'/data/settings/'
    basespect = arparse.BaseSpect(path+'settings.basespect', '.tmp')
    base_lines = basespect.load_file()
    basespect.set_paramlist(base_lines)

    # ARMLSD instruments
    for specname in ['kast_blue', 'kast_red', 'lris_blue', 'lris_red', 'isis_blue']:
        msgs.info("===============================================")
        msgs.info("Working on {:s}".format(specname))
        spect = arparse.get_spect_class(('ARMLSD', specname, ".tmp"))
        spect_lines = spect.load_file()
        spect.set_paramlist(spect_lines)
        msgs.info("===============================================")
        for key in basespect._spect.keys():
            if key in ['set']:
                continue
            if key in spect._spect.keys():
                compare_dicts(key, basespect._spect[key], spect._spect[key], skip_keys=('index'))

if __name__ == '__main__':

    flg_sett = 0
    #flg_sett += 2**0  # argflag checking
    flg_sett += 2**1 # spect checking

    if flg_sett & (2**0):
        argf_diff_and_dup()
    if flg_sett & (2**1):
        spect_diff_and_dup()
