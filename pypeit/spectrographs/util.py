""" Utilities for spectograph codes
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)
import numpy as np

from pypeit import msgs
from pypeit import spectrographs

# TODO: Allow the spectrographs to be identified by their camera?  Won't
# work for 'shane_kast_red' and 'shane_kast_red_ret'.

def valid_spectrographs():
    # TODO: Is there a more clever way to do this?  If we change these
    # names, we could do something like what's done in
    # pypeit.instantiate_me.
    return ['gemini_gnirs','keck_deimos', 'keck_lris_blue', 'keck_lris_red', 'keck_nires', 'keck_nirspec',
            'shane_kast_blue', 'shane_kast_red', 'shane_kast_red_ret', 'tng_dolores',
            'wht_isis_blue', 'vlt_xshooter_uvb', 'vlt_xshooter_vis', 'vlt_xshooter_nir',
            'gemini_gmos_south', 'gemini_gmos_north_e2v', 'gemini_gmos_north_ham']

def load_spectrograph(spectrograph):
    """
    Instantiate a :class:`spectrographs.spectrograph.Spectrograph`, if
    possible.

    Args:

        spectrograph (:obj:`str`,
            :class:`spectrographs.spectrograph.Spectrograph`): The
            spectrograph to instantiate.  If the input is a spectrograph
            instance, the instance is simply returned.  If a string, the
            string is used to select the spectrograph to instantiate.
            If None, None is returned.

    Returns:
        :class:`spectrographs.spectrograph.Spectrograph`: The
        spectrograph used to obtain the data to be reduced.
    """

    # TODO: I'm not crazy about the idea that that spectrograph can be
    # undefined, even when one is working with master files.  It means
    # that objects are not fully defined.
    if spectrograph is None:
        return None

    if isinstance(spectrograph, spectrographs.spectrograph.Spectrograph):
        return spectrograph

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

    if spectrograph == 'keck_nirspec_low':
        return spectrographs.keck_nirspec.KeckNIRSPECLowSpectrograph()

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

    if spectrograph == 'gemini_gmos_south':
        return spectrographs.gemini_gmos.GeminiGMOSSSpectrograph()

    if spectrograph == 'gemini_gmos_north_e2v':
        return spectrographs.gemini_gmos.GeminiGMOSNE2VSpectrograph()

    if spectrograph == 'gemini_gmos_north_ham':
        return spectrographs.gemini_gmos.GeminiGMOSNHamSpectrograph()

    msgs.error('{0} is not a supported spectrograph.'.format(spectrograph))


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

