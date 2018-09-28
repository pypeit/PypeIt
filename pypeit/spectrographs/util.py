""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from pypeit import msgs
from pypeit import spectrographs

from pypeit import debugger

# TODO: Allow the spectrographs to be identified by their camera?  Won't
# work for 'shane_kast_red' and 'shane_kast_red_ret'.

def valid_spectrographs():
    # TODO: Is there a more clever way to do this?
    return ['gemini_gnirs','keck_deimos', 'keck_lris_blue', 'keck_lris_red', 'keck_nires', 'keck_nirspec',
            'shane_kast_blue', 'shane_kast_red', 'shane_kast_red_ret', 'tng_dolores',
            'wht_isis_blue', 'vlt_xshooter_uvb', 'vlt_xshooter_vis', 'vlt_xshooter_nir']

def load_spectrograph(spectrograph=None):
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

    if spectrograph == 'gemini_gnirs':
        return spectrographs.gemini_gnirs.GeminiGNIRSSpectrograph()

    if spectrograph == 'keck_deimos':
        return spectrographs.keck_deimos.KeckDEIMOSSpectrograph()

    if spectrograph == 'keck_lris_blue':
        return spectrographs.keck_lris.KeckLRISBSpectrograph()

    if spectrograph == 'keck_lris_red':
        return spectrographs.keck_lris.KeckLRISRSpectrograph()

    if spectrograph == 'keck_nires':
        return spectrographs.keck_nires.KeckNIRESSpectrograph()

    if spectrograph == 'keck_nirspec':
        return spectrographs.keck_nirspec.KeckNIRSPECSpectrograph()

    if spectrograph == 'shane_kast_blue':
        return spectrographs.shane_kast.ShaneKastBlueSpectrograph()

    if spectrograph == 'shane_kast_red':
        return spectrographs.shane_kast.ShaneKastRedSpectrograph()

    if spectrograph == 'shane_kast_red_ret':
        return spectrographs.shane_kast.ShaneKastRedRetSpectrograph()

    if spectrograph == 'wht_isis_blue':
        return spectrographs.wht_isis.WhtIsisBlueSpectrograph()
    
    if spectrograph == 'tng_dolores':
        return spectrographs.tng_dolores.TNGDoloresSpectrograph()

    if spectrograph == 'vlt_xshooter_uvb':
        return spectrographs.vlt_xshooter.VLTXShooterUVBSpectrograph()

    if spectrograph == 'vlt_xshooter_vis':
        return spectrographs.vlt_xshooter.VLTXShooterVISSpectrograph()

    if spectrograph == 'vlt_xshooter_nir':
        return spectrographs.vlt_xshooter.VLTXShooterNIRSpectrograph()

    debugger.set_trace()
    msgs.error("Spectrograph not supported")


def checkme(chk_dict, headarr):
    #
    skip = False
    for head_idx in chk_dict.keys():
        for ch, value in chk_dict[head_idx].items():
            if (value in str(headarr[head_idx][ch]).strip()) is False:
                msgs.info(ch, head_idx, value)
                msgs.warn("Skipping the file..")
                skip = True
    # Return
    return skip
