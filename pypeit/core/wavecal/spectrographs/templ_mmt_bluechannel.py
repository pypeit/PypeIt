"""
Generate the wavelength templates for the MMTO's Blue Channel Spectrograph
"""

import os
from pypeit.core.wavecal import templates


def mmt_bluechannel_300(overwrite=False):
    """
    Template for Blue Channel's 300 l/mm grating

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_300GPM.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '300')

    basefiles = ['mmt_bluechannel_300GPM_5700.fits', 'mmt_bluechannel_300GPM_6600.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1]
    slits = [126, 126]
    wv_cuts = [6000]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


def mmt_bluechannel_500(overwrite=False):
    """
    Template for Blue Channel's 500 l/mm grating

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_500GPM.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '500')

    basefiles = ['mmt_bluechannel_500GPM_4500.fits', 'mmt_bluechannel_500GPM_5600.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1]
    slits = [126, 126]
    wv_cuts = [5500]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


def mmt_bluechannel_800(overwrite=False):
    """
    Template for Blue Channel's 800 l/mm grating

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_800GPM.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '800')

    basefiles = ['mmt_bluechannel_800GPM_4100.fits', 'mmt_bluechannel_800GPM_5000.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1]
    slits = [126, 126]
    wv_cuts = [4500]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


def mmt_bluechannel_832_order1(overwrite=False):
    """
    Template for Blue Channel's 832 l/mm grating in 1st order (red)

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_832GPM_order1.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '832')

    basefiles = ['mmt_bluechannel_832GPM_7500.fits', 'mmt_bluechannel_832GPM_7700.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1]
    slits = [203, 126]
    wv_cuts = [7600]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


def mmt_bluechannel_832_order2(overwrite=False):
    """
    Template for Blue Channel's 832 l/mm grating in 2nd order (blue)

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_832GPM_order2.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '832')

    basefiles = ['mmt_bluechannel_832GPM_3700.fits', 'mmt_bluechannel_832GPM_4200.fits', 'mmt_bluechannel_832GPM_4500.fits', 'mmt_bluechannel_832GPM_4800.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1, 2, 3]
    slits = [126, 126, 126, 126]
    wv_cuts = [4000, 4400, 4800]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


def mmt_bluechannel_1200(overwrite=False):
    """
    Template for Blue Channel's 1200 l/mm grating

    Args:
        overwrite: bool, optional
            Overwrite the existing file? [Default: False]
    """
    binspec = 1
    outroot = 'mmt_bluechannel_1200GPM.fits'
    wpath = os.path.join(templates.template_path, 'MMTO_BlueChannel', '1200')

    basefiles = ['mmt_bluechannel_1200GPM_4200.fits', 'mmt_bluechannel_1200GPM_5300.fits', 'mmt_bluechannel_1200GPM_6350.fits', 'mmt_bluechannel_1200GPM_7000.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    ifiles = [0, 1, 2, 3]
    slits = [126, 126, 126, 126]
    wv_cuts = [4700, 5800, 6700]
    det_cut = None

    templates.build_template(
        wfiles, slits, wv_cuts, binspec, outroot,
        ifiles=ifiles, det_cut=det_cut, chk=True,
        normalize=False, lowredux=False,
        subtract_conti=True, overwrite=overwrite,
        shift_wave=True
    )


if __name__ == '__main__':
    mmt_bluechannel_300(overwrite=True)
    mmt_bluechannel_500(overwrite=True)
    mmt_bluechannel_800(overwrite=True)
    mmt_bluechannel_832_order1(overwrite=True)
    mmt_bluechannel_832_order2(overwrite=True)
    mmt_bluechannel_1200(overwrite=True)
