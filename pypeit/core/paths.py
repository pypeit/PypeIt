"""
Functions that create/check the PypeIt directory paths.
"""
import os
from pypeit import msgs

# Moved to pypeit.pypeit.PypeIt
# TODO: These should be created by the classes that actually need to use
# them
def make_dirs(spectrograph, caldir, scidir, qadir, redux_path=None)
    """
    Make the directories for the pypeit output.

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
    """
    _redux_path = os.getcwd() if redux_path is None else redux_path

    # Science directory
    # TODO: Defined by the Reduce class?
    newdir = os.path.join(_redux_path, scidir)
    msgs.info('Science data output to: {0}'.format(newdir))
    if not os.path.exists(newdir):
        os.makedirs(newdir)

#    # Create a directory for each object in the Science directory
#    msgs.info("Creating Object directories")

    # Calibrations/Masters directory
    # TODO: This should be defined by the Calibrations class
    newdir = os.path.join(_redux_path, caldir+'_'+spectrograph)
    msgs.info('Master calibration data output to: {0}'.format(newdir))
    if not os.path.exists(newdir):
        os.makedirs(newdir)

    # QA directory

    # TODO: Defined by the PypeIt class?  Is anything written to the qa
    # dir or only to qa/PNGs?  Should we have separate calibration and
    # science QA directories?
    newdir = os.path.join(_redux_path, qadir)
    msgs.info('Quality assessment plots output to: {0}'.format(newdir))
    newdir = os.path.join(_redux_path, qadir, 'PNGs')
    if not os.path.exists(newdir):
        os.makedirs(newdir)

