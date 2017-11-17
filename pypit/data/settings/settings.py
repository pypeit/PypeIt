""" Module for examining/archiving/etc. settings files.
This cannot be located in the pypit/ folder (import issues)
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
from glob import glob
import numpy as np
from os.path import dirname
import filecmp
import time
from shutil import copyfile

from astropy.time import Time

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


def archive():
    """ Generate archival file for the baseargf file or spect file
    and instrument setting files
    """
    import pypit
    settings_path = pypit.__path__[0]+'/data/settings/'
    archive_path = pypit.__path__[0]+'/data/settings/archive/'
    # Files
    sett_files = glob(settings_path+'settings.*')
    for sfile in sett_files:
        # Extension
        ext = sfile.split('.')[-1]
        if ext in ['py', 'pyc']:
            continue
        msgs.info("===============================================")
        sroot = sfile.split('/')[-1]
        msgs.info("Working on settings file {:s}".format(sroot))
        msgs.info("===============================================")
        # Archive
        arch_file = current_sett_file(archive_path, sfile)
        if arch_file is None:
            match = False
            arch_root = None
        else: # Compare
            match = filecmp.cmp(sfile, arch_file)
            arch_root = arch_file.split('/')[-1]
        if not match:
            msgs.warn("Current archive {:s} does not match {:s}".format(arch_root, sfile))
            new_arch = archive_path+'/settings.{:s}.{:s}'.format(time.strftime("%Y-%m-%d"), ext)
            msgs.warn("Generating a new archive file: {:s}".format(new_arch.split('/')[-1]))
            copyfile(sfile, new_arch)
            msgs.warn("Add it to the repository!")
        else:
            msgs.info("Current archive file {:s} matches current settings file".format(arch_root))


def archive_instr():
    pass


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


def current_sett_file(apath, sfile):
    # Load archive files
    ftype = sfile.split('.')[-1]
    archive_files = glob(apath+'settings.*.{:s}'.format(ftype))
    if len(archive_files) == 0:
        msgs.warn("No archival files found for {:s}".format(sfile))
        return None
    # Find the most recent
    dates = []
    for afile in archive_files:
        dates.append(afile.split('.')[-2])
    times = Time(dates)
    imax = np.argmax(times)
    arch_file = archive_files[imax]
    # Return
    return arch_file


def spect_diff_and_dup():
    import pypit
    # Load default spect file
    path = pypit.__path__[0]+'/data/settings/'
    basespect = arparse.BaseSpect(path+'settings.basespect', '.tmp')
    base_lines = basespect.load_file()
    basespect.set_paramlist(base_lines)

    # ARMLSD instruments
    for specname in ['shane_kast_blue', 'shane_kast_red', 'keck_lris_blue', 'keck_lris_red', 'wht_isis_blue']:
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
    flg_sett += 2**0  # argflag checking
    flg_sett += 2**1 # spect checking
    flg_sett += 2**2 # archive baseargflag

    if flg_sett & (2**0):
        argf_diff_and_dup()
    if flg_sett & (2**1):
        spect_diff_and_dup()
    if flg_sett & (2**2):
        archive()
