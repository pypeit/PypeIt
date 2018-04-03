""" Module for setting up Development-suite tests
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np
import scipy
import glob

#from pypit import armsgs
from pypit import msgs
from pypit import arparse as settings

from pypit import ardebug as debugger

# Logging
#msgs = armsgs.get_logger()


def set_param(argf, specname):
    """ Instrument specific parameters for Development-suite testing
      (may eventually need to be grating specific)
    Parameters
    ----------
    argf :
    specname : str
    """
    if specname == 'keck_lris_red':
        argf.set_param('arc calibrate method arclines')
    # Return
    return
