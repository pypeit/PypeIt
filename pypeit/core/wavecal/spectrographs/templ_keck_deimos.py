""" Generate the wavelength templates for Keck/DEIMOS"""
import os

from pypeit.core.wavecal import templates


# Keck/DEIMOS

def keck_deimos_600ZD(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_600ZD.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_DEIMOS', '600ZD')

    basefiles = ['MasterWaveCalib_A_1_02_useS0896.fits', 'MasterWaveCalib_A_1_02_useS0477.fits',
                 'MasterWaveCalib_A_1_08_useS1096.fits', 'MasterWaveCalib_A_1_07_useS0209.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1, 2, 3]
    slits = [896, 477, 1096, 209]
    wv_cuts = [5500., 7560., 9404.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_deimos_830G(overwrite=False):

    binspec = 1
    outroot = 'keck_deimos_830G.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_DEIMOS', '830G')

    basefiles = ['MasterWaveCalib_A_1_04_useS1460.fits', 'MasterWaveCalib_A_1_04_useS0933.fits',
                 'MasterWaveCalib_A_1_05_useS1682.fits', 'MasterWaveCalib_A_1_08_useS0844.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1, 2, 3]
    slits = [1460, 933, 1682, 844]
    wv_cuts = [6477., 8342., 9336.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


def keck_deimos_1200G(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_1200G.fits'
    # 3-3 = blue  6268.23 -- 7540
    # 3-14 = red   6508 -- 7730
    # 7-3 = blue  7589 -- 8821
    # 7-17 = red  8000 - 9230
    # 7c-0 = red  9120 -- 9950
    ifiles = [3, 5, 4, 0, 0, 1, 1, 2]
    slits = [1261, 1652, 132, 3, 14, 3, 17, 0]
    lcut = [5200., 5580., 6800., 7450., 7730., 8170, 9120]
    wfile1 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_03.json')
    wfile2 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07.json')
    wfile3 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07c.json')
    wfile4 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
                          'MasterWaveCalib_B_1_02_useS1261.fits')
    wfile5 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
                          'MasterWaveCalib_B_1_06_useS0132.fits')
    wfile6 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
                          'MasterWaveCalib_B_1_02_useS1652.fits')
    #wfile7 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
    #                      'MasterWaveCalib_B_1_06_useS1649.fits')
    files = [wfile1, wfile2, wfile3, wfile4, wfile5, wfile6] #, wfile7]

    # det_dict
    det_cut = None
    # det_cut = {}
    # det_cut['dets'] = [[1,2,3,4], [5,6,7,8]]
    # det_cut['wcuts'] = [[0,9000.], [8200,1e9]]  # Significant overlap is fine
    #
    templates.build_template(files, slits, lcut, binspec, outroot, lowredux=False,
                   ifiles=ifiles, det_cut=det_cut, chk=True, subtract_conti=True,
                             overwrite=overwrite, shift_wave=True)


def keck_deimos_1200B(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_1200B.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200B')
    basefiles = ['MasterWaveCalib_A_1_02_useS0106.fits', 'MasterWaveCalib_A_1_02_useS0291.fits',
                 'MasterWaveCalib_A_1_06_useS0106.fits', 'MasterWaveCalib_A_1_06_useS0287.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [1, 0, 1, 0, 3, 2]
    slits = [291, 106, 291, 106, 287, 106]
    wv_cuts = [4493., 4870., 5100., 5260., 5810.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                   ifiles=ifiles, det_cut=det_cut, chk=True, normalize=False, lowredux=False,
                   subtract_conti=True, overwrite=overwrite, shift_wave=True)

def keck_deimos_900ZD(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_900ZD.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_DEIMOS', '900ZD')

    basefiles = ['MasterWaveCalib_A_1_01_useS1046.fits', 'MasterWaveCalib_A_1_03_useS0600.fits',
                 'MasterWaveCalib_A_1_06_useS0054.fits', 'MasterWaveCalib_A_1_02_useS0066.fits',
                 'MasterWaveCalib_A_1_06_useS0193.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 1, 2, 3, 4, 5]
    slits = [1046, 600, 54, 66, 193]
    wv_cuts = [5250., 5878., 7100., 8245.]
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
    #keck_deimos_600ZD(overwrite=False)
    #keck_deimos_830G(overwrite=False) # False for Testing; True for real
    #keck_deimos_1200G(overwrite=False)
    #keck_deimos_1200B()
    #keck_deimos_900ZD(overwrite=False)
    pass
