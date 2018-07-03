""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from pypit import msgs
from pypit import spectrographs

from pypit import ardebug as debugger

def load_spec_class(spectrograph=None) #, data_file=None):
    """
    Instantiate a Spectrograph class based on the given input

    Args:
        spectrograph: str, optional

        data_file: str, optional
            KBW: Should not be instantiated from a file.  Will need to
            provide an example.

    Returns:
        spec_class: Spectrograph class (or one of its children)

    """

    if spectrograph is None:
        return None

    # TODO: Do we want "in" here, or "=="
    if 'keck_lris_blue' in spectrograph:
        return spectrographs.keck_lris.KeckLRISBSpectrograph()

    elif 'keck_lris_red' in spectrograph:
        return spectrographs.keck_lris.KeckLRISRSpectrograph()

    elif 'keck_deimos' in spectrograph:
        return spectrographs.keck_deimos.KeckDEIMOSSpectrograph()

    elif 'keck_nirspec' in spectrograph:
        return spectrographs.keck_nirspec.KeckNIRSPECSpectrograph()

    elif 'shane_kast_blue' in spectrograph:
        return spectrographs.shane_kast.ShaneKastBlueSpectrograph()

    elif 'shane_kast_red' in spectrograph:
        return spectrographs.shane_kast.ShaneKastRedSpectrograph()

    elif 'shane_kast_red_ret' in spectrograph:
        return spectrographs.shane_kast.ShaneKastRedRetSpectrograph()

    elif 'wht_isis_blue' in spectrograph:
        return spectrographs.wht_isis.WhtIsisBlueSpectrograph()
    
    elif 'tng_dolores' in spectrograph:
        return spectrographs.tng_dolores.TngDoloresSpectrograph()
    
    msgs.error("Spectrograph not supported")
    # TODO: Should never get here, right?
    return None

