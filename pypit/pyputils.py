""" Some PYPIT utilities for the package
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

import re
from os.path import realpath, join

this_file = realpath(__file__)
this_path = this_file[:this_file.rfind('/')]


def get_version():
    """Get the value of ``__version__`` without having to import the module.
    Copied from DESI desiutils

    Parameters
    ----------

    Returns
    -------
    :class:`str`
        The value of ``__version__``.
    """
    ver = 'unknown'
    version_file = join(this_path, '_version.py')
    with open(version_file, "r") as f:
        for line in f.readlines():
            mo = re.match("__version__ = '(.*)'", line)
            lu = re.match("__lastupdate__ = '(.*)'", line)
            if mo:
                ver = mo.group(1)
            if lu:
                upd = lu.group(1)
    return ver, upd


def get_dummy_logger():
    """ Useful for testing

    Returns
    -------

    """
    from pypit import ardebug
    from pypit import armsgs as pyparm
    debug = ardebug.init()

    pyparm.pypit_logger = pyparm.Messages(None, debug, 0)
    return pyparm.pypit_logger


def make_pypit_file(pyp_file, spectrograph, dfnames, parlines=None,
                    setup_script=False, calcheck=False, spclines=None,
                    skip_files=[]):
    """ Generate a default PYPIT settings file

    Parameters
    ----------
    pyp_file : str
      Name of PYPIT file to be generated
    spectrograph : str
    dfnames : list
      Path + file root of datafiles
    parlines : list, optional
      Standard parameter calls
    spclines : list, optional
      Lines related to filetype and calibrations
    setup_script : bool, optional
      Running setup script?
    calcheck : bool, optional
      Run calcheck?
    skip_files : list, optional
      Files to skip

    Returns
    -------
    Creates a PYPIT File

    """
    # Error checking
    if not isinstance(dfnames, list):
        raise IOError("files_root needs to be a list")

    # Defaults
    if parlines is None:
        parlines = ['run ncpus 1\n',
                    'run calcheck True\n',
                    'output overwrite True\n']
    if spclines is None:
        if setup_script:
            spclines = [' pixelflat number 0\n', ' arc number 1\n',
                        ' slitflat number 0\n',
                        ' bias number 0\n', ' standard number 0\n']
        else:
            spclines = [' pixelflat number 1\n', ' arc number 1\n',
                        ' slitflat number 1\n',
                        ' bias number 1\n', ' standard number 1\n']

    # Here we go
    with open(pyp_file, 'w') as f:
        f.write("# This is a comment line\n")
        f.write("\n")
        f.write("# Change the default settings\n")
        for parline in parlines:
            f.write(parline)
        f.write("run ncpus 1\n")
        if calcheck:
            f.write("run calcheck True\n")
        f.write("run spectrograph {:s}\n".format(spectrograph))
        f.write("\n")
        f.write("# Reduce\n")
        f.write("\n")
        f.write("# Read in the data\n")
        f.write("data read\n")
        # Data file paths
        for dfname in dfnames:
            f.write(dfname)
        # Skip files
        for skip_file in skip_files:
            f.write(' skip '+skip_file+'\n')
        f.write("data end\n")
        f.write("\n")
        f.write("spect read\n")
        for spcline in spclines:
            f.write(spcline)
        f.write("spect end\n")
