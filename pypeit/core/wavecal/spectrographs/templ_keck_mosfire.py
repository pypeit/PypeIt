""" Generate the wavelength templates for Keck/MOSFIRE"""
import os

from pypeit.core.wavecal import templates


def keck_mosfire_OH_H(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_OH_H.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'OH_lines', 'H')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS0434.fits', 'MasterWaveCalib_A_1_DET01_useS1585.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [434, 1585]
    wv_cuts = [17150.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_OH_J(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_OH_J.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'OH_lines', 'J')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1014.fits', 'MasterWaveCalib_A_1_DET01_useS1632.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1014, 1632]
    wv_cuts = [13190.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_OH_J2(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_OH_J2.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'OH_lines', 'J2')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1102.fits', 'MasterWaveCalib_A_1_DET01_useS1565.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1102, 1565]
    wv_cuts = [11740.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_OH_K(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_OH_K.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'OH_lines', 'K')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1032.fits', 'MasterWaveCalib_A_1_DET01_useS01275.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1032, 1275]
    wv_cuts = [20970.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_OH_Y(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_OH_Y.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'OH_lines', 'Y')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1258.fits', 'MasterWaveCalib_A_1_DET01_useS1099.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1258, 1099]
    wv_cuts = [9930.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


# Below are to create templates using the arc frames (Ne, Ar)
def keck_mosfire_arcs_K(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_arcs_K.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'arclines', 'K')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1370.fits', 'MasterWaveCalib_A_1_DET01_useS1278.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1370, 1278]
    wv_cuts = [21630.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_arcs_H(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_arcs_H.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'arclines', 'H')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1058.fits', 'MasterWaveCalib_A_1_DET01_useS1211.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1058, 1211]
    wv_cuts = [17090.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_arcs_J(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_arcs_J.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'arclines', 'J')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1324.fits', 'MasterWaveCalib_A_1_DET01_useS0545.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1324, 545]
    wv_cuts = [11850.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_arcs_J2(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_arcs_J2.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'arclines', 'J2')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS0838.fits', 'MasterWaveCalib_A_1_DET01_useS0126.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [838, 126]
    wv_cuts = [12530.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_mosfire_arcs_Y(overwrite=False):
    binspec = 1
    outroot = 'keck_mosfire_arcs_Y.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_MOSFIRE', 'arclines', 'Y')

    basefiles = ['MasterWaveCalib_A_1_DET01_useS1478.fits', 'MasterWaveCalib_A_1_DET01_useS1543.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1]
    slits = [1478, 1543]
    wv_cuts = [10200.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


if __name__ == '__main__':
    # templates using OH lines
    # keck_mosfire_OH_H(overwrite=False)
    # keck_mosfire_OH_J(overwrite=False)
    # keck_mosfire_OH_J2(overwrite=False)
    # keck_mosfire_OH_K(overwrite=False)
    # keck_mosfire_OH_Y(overwrite=False)

    # templates using arc lines
    # keck_mosfire_arcs_K(overwrite=False)
    # keck_mosfire_arcs_H(overwrite=False)
    # keck_mosfire_arcs_J(overwrite=False)
    # keck_mosfire_arcs_J2(overwrite=False)
    # keck_mosfire_arcs_Y(overwrite=False)
    pass
