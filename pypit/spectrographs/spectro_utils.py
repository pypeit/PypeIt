""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from pypit import msgs
from pypit.spectrographs import gemini_gmos
from pypit.spectrographs import keck_deimos
from pypit.spectrographs import keck_lris
from pypit.spectrographs import keck_nirspec
from pypit.spectrographs import keck_nires
from pypit.spectrographs import shane_kast
from pypit.spectrographs import wht_isis
from pypit.spectrographs import tng_dolores

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
        elif 'shane_kast_blue' in spectrograph:
            spec_class = shane_kast.ShaneKastBlueSpectrograph()
        elif 'shane_kast_red' in spectrograph:
            spec_class = shane_kast.ShaneKastRedSpectrograph()
        elif 'shane_kast_red_ret' in spectrograph:
            spec_class = shane_kast.ShaneKastRedRetSpectrograph()
        elif 'wht_isis_blue' in spectrograph:
            spec_class = wht_isis.WhtIsisBlueSpectrograph()
        elif 'tng_dolores' in spectrograph:
            spec_class = tng_dolores.TngDoloresSpectrograph()
        elif 'keck_nires' in spectrograph:
            spec_class = keck_nires.KeckNIRESpectrograph()
        elif 'gemini_gmos_south' in spectrograph:
            spec_class = gemini_gmos.GeminiGMOSSSpectrograph()
        else:
            msgs.error("Spectrograph not supported")
    #
    return spec_class

