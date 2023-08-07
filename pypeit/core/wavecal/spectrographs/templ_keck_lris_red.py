""" Generate the wavelength templates for Keck/LRIS"""
import os

from pypeit.core.wavecal import templates


# Keck/DEIMOS
def keck_lris_red_mark4_R400(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_mark4_R400.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'Mark4', 'R400')

    basefiles = ['MasterWaveCalib_A_1_01_long.fits',
                 'MasterWaveCalib_A_1_01_sunil.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [2048, 2045]
    wv_cuts = [9320.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_red_R150_7500(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R150_7500_ArCdHgNeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R150_7500')

    basefiles = ['WaveCalib_A_0_DET02_S0303.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0]
    slits = [303]
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


def keck_lris_red_R300_5000(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R300_5000_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R300_5000')

    basefiles = ['WaveCalib_A_0_DET02_S0309.fits', 'WaveCalib_A_0_DET02_S1045.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [309, 1045]
    wv_cuts = [5680.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)

def keck_lris_red_R400_8500(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R400_8500_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R400_8500')

    basefiles = ['WaveCalib_A_0_DET01_S1549.fits', 'WaveCalib_A_0_DET02_S0876.fits', 'WaveCalib_A_0_DET01_S0694.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1,2]
    slits = [1549, 876, 694]
    wv_cuts = [5510., 6800.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)



# Run em
if __name__ == '__main__':
    # keck_lris_red_mark4_R400()#overwrite=True)
    # keck_lris_red_R150_7500(overwrite=False)
    keck_lris_red_R300_5000(overwrite=False)
    # keck_lris_red_R400_8500(overwrite=False)


