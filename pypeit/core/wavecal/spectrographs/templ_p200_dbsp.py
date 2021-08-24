""" Generate the wavelength templates for P200/DBSP"""
import os

from pypeit.core.wavecal import templates

##
# DBSP red arm
##

# ##############################
def p200_dbsp_red_1200_7100_d68(overwrite=False):  # DBSPr 1200/7100 D68

    binspec = 1
    outroot = 'p200_dbsp_red_1200_7100_d68.fits'
    #
    ifiles = [0, 1]
    slits = [222, 0]  # Be careful with the order..
    lcut = [7740.]
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'R1200_7100_D68', 'MasterWaveCalib_A_1_01_7600.fits')
    wfile1 = os.path.join(templates.template_path, 'P200_DBSP',
        'R1200_7100_D68', 'wvcalib_8200.fits')

    #
    templates.build_template([wfile0, wfile1], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_red_1200_9400_d55(overwrite=False):  # DBSPr 1200/9400 D55

    binspec = 1
    outroot = 'p200_dbsp_1200_9400_d55.fits'
    #
    ifiles = [0]
    slits = [0]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'R1200_9400_D55', 'wvcalib_8800.fits')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_red_600_10000_d55(overwrite=False):  # DBSPr 600/10000 D55

    binspec = 1
    outroot = 'p200_dbsp_red_600_10000_d55.fits'
    #
    ifiles = [0]
    slits = [221]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'R600_10000_D55', 'MasterWaveCalib_A_1_01.fits')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_red_316_7500_d55(overwrite=False):  # DBSPr 316/7500 D55

    binspec = 1
    outroot = 'p200_dbsp_red_316_7500_d55.fits'
    #
    ifiles = [0]
    slits = [221]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'R316_7500_D55', 'P200_DBSP_Red.json')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

##
# DBSP blue arm
##

# ##############################
def p200_dbsp_blue_300_3990_d55(overwrite=False):  # DBSPb 300/3990 D55

    binspec = 1
    outroot = 'p200_dbsp_blue_300_3990_d55.fits'
    #
    ifiles = [0]
    slits = [0]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'B300_3990_D55', 'wvcalib.fits')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_blue_600_4000_d55(overwrite=False):  # DBSPb 600/4000 D55

    binspec = 1
    outroot = 'p200_dbsp_blue_600_4000_d55.fits'
    #
    ifiles = [0]
    slits = [231]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'B600_4000_D55', 'P200_DBSP_Blue.json')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_blue_1200_5000_d55(overwrite=False):  # DBSPb 1200/5000 D55

    binspec = 1
    outroot = 'p200_dbsp_blue_1200_5000_d55_4700.fits'
    #
    ifiles = [0]
    slits = [0]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'B1200_5000_D55', 'wvcalib_4700.fits')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

# ##############################
def p200_dbsp_blue_1200_5000_d68(overwrite=False):  # DBSPb 1200/5000 D68

    binspec = 1
    outroot = 'p200_dbsp_blue_1200_5000_d68_6000.fits'
    #
    ifiles = [0]
    slits = [180]  # Be careful with the order..
    lcut = None
    wfile0 = os.path.join(templates.template_path, 'P200_DBSP',
        'B1200_5000_D68', 'MasterWaveCalib_A_1_01_6000.fits')

    #
    templates.build_template([wfile0], slits, lcut, binspec, outroot,
        lowredux=False, ifiles=ifiles, normalize=True, overwrite=overwrite)

if __name__ == '__main__':
    # p200_dbsp_red_316_7500_d55(overwrite=True)
    # p200_dbsp_red_600_10000_d55(overwrite=True)
    p200_dbsp_red_1200_7100_d68(overwrite=True)
