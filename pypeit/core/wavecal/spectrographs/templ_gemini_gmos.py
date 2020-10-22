""" Generate the wavelength templates for Keck/DEIMOS"""
import os

from pypeit.core.wavecal import templates

from IPython import embed

# ##############################
def gemini_gmos_r400_hama(overwrite=False):  # GMOS R400 Hamamatsu

    binspec = 2
    outroot = 'gemini_gmos_r400_ham.fits'
    #
    ifiles = [0, 1, 2, 3, 4]
    slits = [0, 2, 3, 0, 0]  # Be careful with the order..
    lcut = [5400., 6620., 8100., 9000.]
    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_01_aa.json')
    embed(header='the file below is missing...')
    wfile5 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_05_aa.json')  # 5190 -- 6679
    # wfile2 = os.path.join(template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_02_aa.json')
    wfile3 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_04_aa.json')
    wfile4 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_03_aa.json')
    wfile6 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_06_aa.json')
    #
    templates.build_template([wfile1, wfile5, wfile3, wfile4, wfile6], slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=True, subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)

# ##############################
def gemini_gmos_r400_e2v(overwrite=False):  # GMOS R400 E2V
    binspec = 2
    outroot = 'gemini_gmos_r400_e2v.fits'
    #
    ifiles = [0, 1, 2]
    slits = [0, 0, 0]
    lcut = [6000., 7450]
    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_1_01.json')
    wfile2 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_1_02.json')
    wfile3 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_1_03.json')
    #
    templates.build_template([wfile1, wfile2, wfile3], slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=True, overwrite=overwrite, subtract_conti=True)

# ##############################
def gemini_gmos_b600_ham(overwrite=False):
    binspec = 2
    outroot = 'gemini_gmos_b600_ham.fits'
    #
    wfile1 = os.path.join(templates.template_path, 'GMOS', 'B600', 'MasterWaveCalib_C_1_01.json')
    wfile5 = os.path.join(templates.template_path, 'GMOS', 'B600', 'MasterWaveCalib_D_1_01.json')  # - 4547
    wfile2 = os.path.join(templates.template_path, 'GMOS', 'B600', 'MasterWaveCalib_C_1_02.json')
    wfile4 = os.path.join(templates.template_path, 'GMOS', 'B600', 'MasterWaveCalib_D_1_02.json')  # 4610-5608
    wfile3 = os.path.join(templates.template_path, 'GMOS', 'B600', 'MasterWaveCalib_C_1_03.json') # xx-6615
    # 1x1 binning from Shenli
    wfile6 = os.path.join(templates.template_path, 'GMOS', 'B600', '1x1',
                          'B600_0.660', 'chip3', 'wvcalib.fits')
    files = [wfile1, wfile5, wfile2, wfile4, wfile3, wfile6]

    ifiles = [0, 1, 2, 3, 4, 5]
    slits = [0, 0, 0, 0, 0, 1]
    lcut = [4250., 4547., 5250., 5615., 6600.]
    #
    templates.build_template(files,
        slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=True, subtract_conti=True, miny=-100., overwrite=overwrite,
                   shift_wave=True)

if __name__ == '__main__':
    #gemini_gmos_r400_e2v(overwrite=True)
    gemini_gmos_b600_ham(overwrite=True)
