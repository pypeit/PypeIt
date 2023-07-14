""" Generate the wavelength templates for LDT/DeVeny"""
import os

from pypeit.core.wavecal import templates


# LDT/DeVeny Wavelength Templates stitched together from multiple grating angles
def ldt_deveny_150(overwrite=False):
    """ldt_deveny_150 Template from the DV1 grating (150 g/mm)

    Args:
        overwrite: bool, optional
          Overwrite the existing file? [Default: False]
    """

    binspec = 1
    outroot = 'ldt_deveny_150_HgCdAr.fits'

    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'LDT_DeVeny', '150')

    basefiles = ['MasterWaveCalib_A_DV1_7200.fits'] 
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    # Snippets
    ifiles = [0]
    slits = [248]
    wv_cuts = []
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def ldt_deveny_300(overwrite=False):
    """ldt_deveny_300 Template from the DV2 & DV3 gratings (300 g/mm)

    Args:
        overwrite: bool, optional
          Overwrite the existing file? [Default: False]
    """

    binspec = 1
    outroot = 'ldt_deveny_300_HgCdAr.fits'

    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'LDT_DeVeny', '300')

    basefiles = ['MasterWaveCalib_A_DV2_5200.fits',
                 'MasterWaveCalib_B_DV3_8000.fits',
                 'MasterWaveCalib_C_DV3_9000.fits'] 
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    # Snippets
    ifiles = [0, 1, 2]
    slits = [248, 248, 248]
    wv_cuts = [6600., 9000.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def ldt_deveny_600(overwrite=False):
    """ldt_deveny_600 Template from the DV6 & DV7 gratings (600 g/mm)

    Args:
        overwrite: bool, optional
          Overwrite the existing file? [Default: False]
    """

    binspec = 1
    outroot = 'ldt_deveny_600_HgCdAr.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'LDT_DeVeny', '600')

    basefiles = ['MasterWaveCalib_D_DV6_4000.fits',
                 'MasterWaveCalib_E_DV6_5200.fits',
                 'MasterWaveCalib_F_DV7_7000.fits',
                 'MasterWaveCalib_G_DV7_8500.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    # Snippets
    ifiles = [0, 1, 2, 3]
    slits = [248, 248, 248, 248]
    wv_cuts = [4600., 6100., 7800.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def ldt_deveny_1200(overwrite=False):
    """ldt_deveny_1200 Template from the DV9 grating (1200 g/mm)

    Args:
        overwrite: bool, optional
          Overwrite the existing file? [Default: False]
    """

    binspec = 1
    outroot = 'ldt_deveny_1200_HgCdAr.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'LDT_DeVeny', '1200')

    basefiles = ['MasterWaveCalib_H_DV9_3800.fits',
                 'MasterWaveCalib_I_DV9_4800.fits',
                 'MasterWaveCalib_J_DV9_5800.fits',
                 'MasterWaveCalib_K_DV9_6800.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    # Snippets
    ifiles = [0, 1, 2, 3]
    slits = [248, 248, 248, 248]
    wv_cuts = [4300., 5300., 6300.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


if __name__ == '__main__':
    ldt_deveny_150(overwrite=True)
    ldt_deveny_300(overwrite=True)
    ldt_deveny_600(overwrite=True)
    ldt_deveny_1200(overwrite=True)
