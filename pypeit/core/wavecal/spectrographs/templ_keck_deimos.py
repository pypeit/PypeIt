""" Generate the wavelength templates for Keck/DEIMOS"""
import os

from pypeit.core.wavecal import templates


# Keck/DEIMOS

def keck_deimos_600ZD():
    binspec = 1
    slits = [0, 1]
    lcut = [7192.]
    xidl_file = os.path.join(templates.template_path, 'Keck_DEIMOS', '600ZD', 'deimos_600.sav')
    outroot = 'keck_deimos_600.fits'
    templates.build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

def keck_deimos_830G(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_830G.fits'
    # 3-12 = blue  6508 -- 8410
    # 7-24 = blue  8497 -- 9925 (no lines after XeI)
    ifiles = [0, 0, 1]
    slits = [12, 14, 24]
    lcut = [8400., 8480]
    wfile1 = os.path.join(templates.template_path, 'Keck_DEIMOS', '830G_M_8600', 'MasterWaveCalib_A_1_03.json')
    wfile2 = os.path.join(templates.template_path, 'Keck_DEIMOS', '830G_M_8600', 'MasterWaveCalib_A_1_07.json')
    # det_dict
    det_cut = {}
    det_cut['dets'] = [[1, 2, 3, 4], [5, 6, 7, 8]]
    det_cut['wcuts'] = [[0, 9000.], [8200, 1e9]]  # Significant overlap is fine
    #
    templates.build_template([wfile1, wfile2], slits, lcut, binspec, outroot, lowredux=False,
                   ifiles=ifiles, det_cut=det_cut, chk=True, overwrite=overwrite)


def keck_deimos_1200G(overwrite=False):
    binspec = 1
    outroot = 'keck_deimos_1200G.fits'
    # 3-3 = blue  6268.23 -- 7540
    # 3-14 = red   6508 -- 7730
    # 7-3 = blue  7589 -- 8821
    # 7-17 = red  8000 - 9230
    # 7c-0 = red  9120 -- 9950
    ifiles = [3, 4, 0, 0, 1, 1, 2]
    slits = [1261, 132, 3, 14, 3, 17, 0]
    lcut = [5500., 6800., 7450., 7730., 8170, 9120]
    wfile1 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_03.json')
    wfile2 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07.json')
    wfile3 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07c.json')
    wfile4 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
                          'MasterWaveCalib_B_1_02_useS1261.fits')
    wfile5 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200G', '1200G_bluetilt',
                          'MasterWaveCalib_B_1_06_useS0132.fits')
    files = [wfile1, wfile2, wfile3, wfile4, wfile5]

    # det_dict
    det_cut = None
    # det_cut = {}
    # det_cut['dets'] = [[1,2,3,4], [5,6,7,8]]
    # det_cut['wcuts'] = [[0,9000.], [8200,1e9]]  # Significant overlap is fine
    #
    templates.build_template(files, slits, lcut, binspec, outroot, lowredux=False,
                   ifiles=ifiles, det_cut=det_cut, chk=True, subtract_conti=True)


def keck_deimos_1200B():
    binspec = 1
    outroot = 'keck_deimos_1200B.fits'
    # file1 = blue  4063 - 5382.8
    # file2 = blue  ??   - 5425.
    # file3 = red   5394.4 - 6709.2
    ifiles = [0, 1]  # , 2]
    slits = [0, 0]  # , 0]  # Not used
    # wv_cuts = [5100., 5400]
    wv_cuts = [5400]
    # Outputs from IRAF by Carlos
    # wfile1 = os.path.join(template_path, 'Keck_DEIMOS', '1200B', 'deimos_calibrated_arc_bluechip_1200B_tilt5200.dat')
    wfile2 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200B',
                          'deimos_calibrated_arc_bluechip_1200B_tilt5200_slit2.dat')
    wfile3 = os.path.join(templates.template_path, 'Keck_DEIMOS', '1200B', 'deimos_calibrated_arc_redchip_1200B_tilt5200.dat')
    # det_dict
    det_cut = None
    #
    templates.build_template([wfile2, wfile3], slits, wv_cuts, binspec, outroot,
                   ascii_tbl=True, ifiles=ifiles, det_cut=det_cut,
                   chk=True, normalize=True, lowredux=False, in_vac=False,
                   subtract_conti=True)

if __name__ == '__main__':
    #keck_deimos_600ZD()
    #keck_deimos_830G(overwrite=False) # False for Testing; True for real
    keck_deimos_1200G(overwrite=True)
    #keck_deimos_1200B()
    pass
