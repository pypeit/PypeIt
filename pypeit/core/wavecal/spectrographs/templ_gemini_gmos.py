""" Generate the wavelength templates for Gemini/GMOS"""
import os

from pypeit.core.wavecal import templates

from IPython import embed

# ##############################
def gemini_gmos_r400_hama(overwrite=False):  # GMOS R400 Hamamatsu

    binspec = 2
    outroot = 'gemini_gmos_r400_ham.fits'
    #
    ifiles = [0, 1, 2, 3, 4, 5, 6]
    slits = [1, 0, 2, 0, 3, 0, 0]  # Be careful with the order..
    lcut = [5270., 5990., 6635., 7580., 8170., 9000.]
    wfile0 = os.path.join(templates.template_path, 'GMOS', 
                          'R400', 'wvcalib_r400_470.fits')
    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 
                          'wvcalib_r400_520.fits')
    #wfile1b = os.path.join(templates.template_path, 'GMOS', 'R400', 
    #                      'MasterWaveCalib_A_01_aa.json')
    wfile5 = os.path.join(templates.template_path, 'GMOS', 
                          'R400', 'MasterWaveCalib_A_05_aa.json')  # 5190 -- 6679
    wfile2 = os.path.join(templates.template_path, 'GMOS', 'R400', 
                          'MasterWaveCalib_A_02_aa.json')
    wfile3 = os.path.join(templates.template_path, 'GMOS', 
                          'R400', 'MasterWaveCalib_A_04_aa.json')
    wfile4 = os.path.join(templates.template_path, 'GMOS', 
                          'R400', 'MasterWaveCalib_A_03_aa.json')
    wfile6 = os.path.join(templates.template_path, 'GMOS', 
                          'R400', 'MasterWaveCalib_A_06_aa.json')
    #
    templates.build_template([wfile0, wfile1, 
                              wfile5, wfile2, wfile3, wfile4, wfile6], 
                             slits, lcut, binspec,
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
def gemini_gmos_r400_e2v_mosaic(overwrite=False):  # GMOS R400 E2V
    # Turns off normalize
    binspec = 2
    outroot = 'gemini_gmos_r400_e2v_mosaic.fits'
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
                   normalize=False, overwrite=overwrite, subtract_conti=True)

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
                          'B600_0.660', 'chip3', 'wvcalib.fits')  # 5 6386 - 7386
    wfile7 = os.path.join(templates.template_path, 'GMOS', 'B600', '1x1',
                          'B600_0.580', 'chip3', 'wvcalib.fits')  # 6 6873 - 7725
    # 2x2 from Sunil
    wfile8 = os.path.join(templates.template_path, 'GMOS', 'B600', 
                          'wvcalib_chip3_8500.fits')  # - 8500
    files = [wfile1, wfile5, wfile2, wfile4, wfile3, wfile6, wfile7, wfile8]

    ifiles = [0, 1, 2, 3, 4, 5, 6, 7]
    slits = [0, 0, 0, 0, 0, 1, 1, 1724]
    lcut = [4250., 4547., 5250., 5615., 6600., 6900., 7500.]
    binning = [2,2,2,2,2,2,1,2]
    # Run
    templates.build_template(files,
        slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=True, subtract_conti=True, miny=-100., overwrite=overwrite,
                   shift_wave=True, binning=binning)

def gemini_gmos_r831_ham(overwrite=False):
    binspec = 1
    outroot = 'gemini_gmos_r831_ham.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'GMOS', 'R831')
    basefiles = ['R831_0.740/chip1/wvcalib.fits',  # 0 6170-6970
                 'R831_0.740/chip2/wvcalib.fits',  # 1 7032-7725
                 'R831_0.740/chip3/wvcalib.fits',  # 2 8000-8523
                 'R831_0.830/chip1/wvcalib.fits',  # 3 7200-7725
                 'R831_0.830/chip2/wvcalib.fits',  # 4 7950 - 8670
        'R831_0.830/chip3/wvcalib.fits',  # 5 9125 - 9356
        'R831_0.860/chip1/wvcalib.fits',  # 6 7437 - 8117
        'R831_0.860/chip2/wvcalib.fits',  # 7 8266 - 8670. (9125)
        'R831_0.860/chip3/wvcalib.fits',  # 8 9125 - 9800
                 'R831_6679/wvcalib.fits',  # 9 6679 -7274
                 ]
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0, 9, 1, 6, 2, 7, 5, 8]
    slits = [1, 1, 1, 1, 1, 1, 1, 2]
    wv_cuts = [6800., 7200., 7500., 8050., 8300., 9125, 9300.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True, normalize=False, lowredux=False,
                             subtract_conti=True, overwrite=overwrite, shift_wave=True)

# ##############################
def gemini_gmos_r400_nham_mosaic(overwrite=False):  
    """This was a one-off for a wide slit with GMOS
    Not for general usage and it won't be archived

    Args:
        overwrite (bool, optional): _description_. Defaults to False.
    """
    # Turns off normalize
    binspec = 2
    outroot = 'gemini_gmos_r400_Nham_wide_mosaic.fits'
    #
    ifiles = [0]
    slits = [388]
    lcut = []
    wfile1 = os.path.join(templates.template_path, 'GMOS', 'R400', 'MasterWaveCalib_A_NHam.fits')
    #
    templates.build_template([wfile1], slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=False, overwrite=overwrite, subtract_conti=True)


if __name__ == '__main__':
    #gemini_gmos_r400_hama()#overwrite=True)
    #gemini_gmos_r400_e2v(overwrite=True)
    #gemini_gmos_r400_e2v_mosaic(overwrite=True)
    gemini_gmos_b600_ham(overwrite=False)
    #gemini_gmos_r831_ham(overwrite=False)
    #gemini_gmos_r400_nham_mosaic(overwrite=True)
