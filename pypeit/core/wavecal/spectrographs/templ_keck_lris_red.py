""" Generate the wavelength templates for Keck/LRIS RED"""
import os

from pypeit.core.wavecal import templates


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

    basefiles = ['WaveCalib_A_0_DET02_S0302.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0]
    slits = [302]
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


def keck_lris_red_R600_5000(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R600_5000_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R600_5000')

    basefiles = ['WaveCalib_A_0_DET01_1783.fits','WaveCalib_A_0_DET01_S1781.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [1783,1781]
    wv_cuts = [6800.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_red_R600_7500(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R600_7500_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R600_7500')

    basefiles = ['WaveCalib_A_0_DET02_S0302.fits', 'WaveCalib_A_0_DET02_S1077.fits', 'WaveCalib_A_0_DET01_S0757.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1,2]
    slits = [302, 1077, 757]
    wv_cuts = [7416., 8560.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_red_R600_10000(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R600_10000_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R600_10000')

    basefiles = ['WaveCalib_A_0_DET01_S0305.fits', 'WaveCalib_A_0_DET01_S0577.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [305, 577]
    wv_cuts = [8170.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_red_R831_8200(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R831_8200_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R831_8200')

    basefiles = ['WaveCalib_A_0_DET02_S0328.fits', 'WaveCalib_A_0_DET02_S1090.fits', 'WaveCalib_A_0_DET01_S0862.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1,2]
    slits = [328, 1090, 862]
    wv_cuts = [7780., 9080.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)

def keck_lris_red_R900_5500(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R900_5500_ArCdHgNeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R900_5500')

    basefiles = ['WaveCalib_A_0_DET02_S1723.fits', 'WaveCalib_A_0_DET02_S0285.fits', 'WaveCalib_A_0_DET01_S1240.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1,2]
    slits = [1723,285,1240]
    wv_cuts = [5700., 6800.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_lris_red_R1200_7500(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_R1200_7500_ArCdHgKrNeXeZn.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_red', 'R1200_7500')

    basefiles = ['WaveCalib_A_0_DET02_S0046.fits', 'WaveCalib_A_0_DET01_S0290.fits', 'WaveCalib_A_0_DET01_S0899.fits',
                 'WaveCalib_A_0_DET01_S1134.fits', 'WaveCalib_A_0_DET02_S0125.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1,2,3,4]
    slits = [46,290,899,1134,125]
    wv_cuts = [5700., 6800., 7800., 8700.]
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
    # keck_lris_red_R300_5000(overwrite=False)
    # keck_lris_red_R400_8500(overwrite=False)
    # keck_lris_red_R600_5000(overwrite=False)
    # keck_lris_red_R600_7500(overwrite=False)
    # keck_lris_red_R600_10000(overwrite=False)
    # keck_lris_red_R831_8200(overwrite=False)
    # keck_lris_red_R900_5500(overwrite=False)
    keck_lris_red_R1200_7500(overwrite=False)

