"""
Functions that create/check the PypeIt directory paths.
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
from pypeit import msgs

def make_dirs(spectrograph, caldir, scidir, qadir, redux_path=None, overwrite=False):
    """
    Make the directories for the pypeit output.

    .. todo::
        I think this should just fault if the directories exist and
        `overwrite` is False.

    Args:
        spectrograph (str):
            The name of the spectrograph that provided the data to be
            reduced.
        caldir (str):
            The directory to use for saving the master calibration
            frames.
        scidir (str):
            The directory to use for the main reduction output files.
        qadir (str):
            The directory to use for the quality assessment output.
        redux_dir (str):
            The directory for the reductions
        overwrite(:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
    """

    # First, get the current working directory
    if redux_path is None:
        redux_path = os.getcwd()
    msgs.info("Creating Science directory")
    newdir = os.path.join(redux_path, scidir)
    if os.path.exists(newdir):
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
        # JFH I am disabling this for now. If the directory exists it does not get re-created. If it does not exist
        # it is created.
#        if not overwrite:
#            rmdir = ''
#            while os.path.exists(newdir):
#                while rmdir != 'n' and rmdir != 'y' and rmdir != 'r':
#                    rmdir = input(msgs.input() + 'Remove this directory and its contents?'
#                                  + '([y]es, [n]o, [r]ename) - ')
#                if rmdir == 'n':
#                    msgs.warn("Any previous calibration files may be overwritten")
#                    break
#                elif rmdir == 'r':
#                    newdir = input(msgs.input()+"Enter a new directory name: ")
#                elif rmdir == 'y':
#                    shutil.rmtree(newdir)
#                    os.mkdir(newdir)
#                    break
#            if rmdir == 'r':
#                os.mkdir(newdir)
    else:
        os.mkdir(newdir)

    # Create a directory for each object in the Science directory
    msgs.info("Creating Object directories")

    # Create a directory where all of the master calibration frames are stored.
    msgs.info("Creating Master Calibrations directory")
    newdir = os.path.join(redux_path, caldir+'_'+spectrograph)
    if os.path.exists(newdir):
        # JFH Disabling this overwrite stuff
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
#        if not overwrite:
#            rmdir = ''
#            while rmdir != 'n' and rmdir != 'y':
#                rmdir = input(msgs.input() + 'Remove this directory and its contents?'
#                              '([y]es, [n]o) - ')
#            if rmdir == 'n':
#                msgs.warn("Any previous calibration files will be overwritten")
#            else:
#                shutil.rmtree(newdir)
#                os.mkdir(newdir)
    else:
        os.mkdir(newdir)

    # Create a directory where all of the QA is stored
    # TODO: I'd rather that this still consider overwrite and fault
    # instead of just proceeding
    msgs.info("Creating QA directory")
    newdir = os.path.join(redux_path, qadir)
    if os.path.exists(newdir):
        msgs.warn("Pre-existing QA plots will be overwritten")
        if not os.path.exists(newdir+'/PNGs'):
            os.mkdir(newdir+'/PNGs')
    else:
        os.mkdir(newdir)
        os.mkdir(newdir+'/PNGs')


