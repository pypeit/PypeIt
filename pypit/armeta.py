""" Module to deal with meta level routines for PYPIT"""
from __future__ import absolute_import, division, print_function, unicode_literals

from pkg_resources import resource_filename
import glob


def allowed_file_types():
    """ List of file types used in PYPIT

    Returns
    -------
    ftype_list : list
    """
    # Define
    ftype_list = [     # NOTE:  arc must be listed first!
        'arc',         # Exposure of one or more arc calibration lamps for wavelength calibration
        'bias',        # Exposure for assessing detector bias (usually 0s)
        'dark',        # Exposure for assessing detector dark current
        'pinhole',     # Exposure for tracing the orders or slits
        'pixelflat',   # Exposure for assessing pixel-to-pixel variations
        'science',   # Exposure on one or more science targets
        'standard',    # Exposure on a 'standard star' used for flux calibration
        'trace',       # Exposure for tracing slits or echelle orders (usually twilight sky or flat lamp)
        'unknown',     # Unknown..
                  ]
    return ftype_list


def instr_list():
    """ List of allowed Instruments in PYPIT

    Returns
    -------
    instruments : list
    """

    settings_path = resource_filename('pypit', 'data/settings')
    settings_files = glob.glob(settings_path+'/settings.*')
    instruments = [ifile.split('.')[1] for ifile in settings_files]  # This includes base
    # Trim
    for kpop in ['armed', 'armlsd', 'baseargflag', 'basespect', 'py']:
        instruments.remove(kpop)
    # Return
    return instruments

