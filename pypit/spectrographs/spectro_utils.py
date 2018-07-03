""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from pypit import msgs
from pypit.spectrographs import keck_deimos
from pypit.spectrographs import keck_lris
from pypit.spectrographs import keck_nirspec
from pypit.spectrographs import spectroclass

from pypit import ardebug as debugger


def load_spec_class(spectrograph=None, data_file=None):
    """
    Instantiate a Spectrograph class based on the given input

    Args:
        spectrograph: str, optional
        data_file: str, optional
          NOT YET IMPLEMENTED

    Returns:
        spec_class: Spectrograph class (or one of its children)

    """
    spec_class = None
    if spectrograph is not None:
        if 'keck_lris_blue' in spectrograph:
            spec_class = keck_lris.KeckLRISBSpectrograph()
        elif 'keck_lris_red' in spectrograph:
            spec_class = keck_lris.KeckLRISRSpectrograph()
        elif 'keck_deimos' in spectrograph:
            spec_class = keck_deimos.KeckDEIMOSSpectrograph()
        elif 'keck_nirspec' in spectrograph:
            spec_class = keck_nirspec.KeckNIRSPECSpectrograph()
        else:
            spec_class = spectroclass.Spectrograph()
    #
    return spec_class

