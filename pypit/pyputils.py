""" Some PYPIT utilities for the package
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import re
import os

from pypit import armsgs
from pypit import ardebug

#from pypit import __version__, __last_updated__

# this_file = realpath(__file__)
# this_path = this_file[:this_file.rfind('/')]

#def get_version():
#    """Get the value of ``__version__`` without having to import the module.
#    Copied from DESI desiutils
#
#    Parameters
#    ----------
#
#    Returns
#    -------
#    :class:`str`
#        The value of ``__version__``.
#    """
#    ver = 'unknown'
#
#    this_file = os.path.realpath(__file__)
#    this_path = this_file[:this_file.rfind('/')]
#
#    version_file = os.path.join(this_path, '_version.py')
#    with open(version_file, "r") as f:
#        for line in f.readlines():
#            mo = re.match("__version__ = '(.*)'", line)
#            lu = re.match("__lastupdate__ = '(.*)'", line)
#            if mo:
#                ver = mo.group(1)
#            if lu:
#                upd = lu.group(1)
#    return ver, upd


def get_dummy_logger(develop=False):
    """ Useful for testing

    Returns
    -------

    """
    debug = ardebug.init()
    debug['develop'] = develop
    return armsgs.Messages(log=None, debug=debug, verbosity=0)

#    from pypit import ardebug
#    from pypit import armsgs as pyparm
#    debug = ardebug.init()
#    debug['develop'] = develop
#
#    pyparm.pypit_logger = pyparm.Messages(None, debug, 0)
#    return pyparm.pypit_logger


def make_pypit_file(pyp_file, spectrograph, dfnames, parlines=None,
                    setup_script=False, calcheck=False, spclines=None,
                    setuplines=None, setupfiles=None, paths=None):
    """ Generate a default PYPIT settings file

    Parameters
    ----------
    pyp_file : str
      Name of PYPIT file to be generated
    spectrograph : str
    dfnames : list
      Path + file root of datafiles
      Includes skip files
    parlines : list, optional
      Standard parameter calls
    spclines : list, optional
      Lines related to filetype and calibrations
    setup_script : bool, optional
      Running setup script?
    calcheck : bool, optional
      Run calcheck?

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
                    'output overwrite True\n']
        parlines += ["run spectrograph {:s}\n".format(spectrograph)]
    if spclines is None:
        if setup_script:
            spclines = [' pixelflat number 0\n', ' arc number 1\n',
                        ' pinhole number 0\n',
                        ' trace number 0\n',
                        ' bias number 0\n', ' standard number 0\n']
        else:
            #spclines = [' pixelflat number -3\n', ' arc number 1\n',
            #            ' slitflat number -3\n',
            #            ' bias number 0\n', ' standard number 1\n']
            spclines = [' arc number 1\n']

    # Here we go
    with open(pyp_file, 'w') as f:
        f.write("# This is a comment line\n")
        f.write("\n")
        f.write("# Change the default settings\n")
        for parline in parlines:
            f.write(parline)
        if setup_script:
            f.write("run setup True\n")
        if calcheck:
            f.write("run calcheck True\n")
        f.write("\n")
        f.write("# Reduce\n")
        f.write("\n")
        # Setup stuff
        f.write("# Spect\n")
        f.write("spect read\n")
        for spcline in spclines:
            f.write(spcline)
        f.write("spect end\n")
        f.write("\n")
        #
        if setuplines is not None:
            f.write("# Setup\n")
            f.write("setup read\n")
            for sline in setuplines:
                f.write(' '+sline)
            f.write("setup end\n")
            f.write("\n")
        # Data
        f.write("# Read in the data\n")
        f.write("data read\n")
        # Old school
        for dfname in dfnames:
            f.write(' '+dfname+'\n')
        # paths and Setupfiles
        if paths is not None:
            for path in paths:
                f.write(' path '+path+'\n')
        if setupfiles is not None:
            for sfile in setupfiles:
                f.write(sfile)
        f.write("data end\n")
        f.write("\n")

