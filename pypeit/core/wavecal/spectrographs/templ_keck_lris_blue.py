""" Generate the wavelength templates for Keck/LRIS Blue"""
import os

from pypeit.core.wavecal import templates


def keck_lris_blue_B300_5000(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_blue_B300_5000_d680_ArCdHgKrNeXeZnFeAr.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_blue', 'B300_5000')

    basefiles = ['WaveCalib_A_0_DET01_S1676.fits', 'WaveCalib_A_0_DET01_S0947.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [1676, 947]
    wv_cuts = [6830.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_blue_B400_3400(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_blue_B400_3400_d560_ArCdHgNeZnFeAr.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_blue', 'B400_3400')

    basefiles = ['WaveCalib_A_0_DET02_S1642.fits', 'WaveCalib_A_0_DET01_S1753.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [1642, 1753]
    wv_cuts = [4970.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_blue_B1200_3400(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_blue_B1200_3400_d560_ArCdHgNeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_blue', 'B1200_3400')

    basefiles = ['WaveCalib_A_0_DET01_S1252.fits', 'WaveCalib_A_0_DET02_S0985.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [1252, 985]
    wv_cuts = [3800.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)

# Run em
if __name__ == '__main__':
    # keck_lris_blue_B300_5000(overwrite=False)
    # keck_lris_blue_B400_3400(overwrite=False)
    keck_lris_blue_B1200_3400(overwrite=False)

