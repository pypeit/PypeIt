"""
Generate the wavelength templates for the MMTO's Binospec
"""

import os
from pypeit.core.wavecal import templates


def mmt_binospec_270(overwrite=False):
    """
    Template for Binospec's 270 l/mm grating

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_binospec_270.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_Binospec', '270')

    basefiles = ['mmto_binospec_270_6500.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0]
    slits = [2055]
    wv_cuts = None
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


if __name__ == '__main__':
    mmt_binospec_270(overwrite=True)
