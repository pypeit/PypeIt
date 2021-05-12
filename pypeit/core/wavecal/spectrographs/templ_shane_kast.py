""" Generate the wavelength templates for Shane Kast"""
import os

from pypeit.core.wavecal import templates

from IPython import embed

# Shane Kastb
def shane_kastb_452(): #    if flg & (2**4):  # 452/3306
    binspec = 1
    slits = [0]
    xidl_file = os.path.join(templates.template_path, 'Shane_Kast', '452_3306', 'kast_452_3306.sav')
    outroot = 'shane_kast_blue_452.fits'
    templates.build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

def shane_kastb_600():  # if flg & (2**5):  # 600/4310
    binspec = 1
    slits = [0,3]
    lcut = [4550.]
    xidl_file = os.path.join(templates.template_path, 'Shane_Kast', '600_4310', 'kast_600_4310.sav')
    outroot = 'shane_kast_blue_600.fits'
    templates.build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

def shane_kastb_830():  # if flg & (2**6):  # 830/3460
    binspec = 1
    slits = [0]
    xidl_file = os.path.join(templates.template_path, 'Shane_Kast', '830_3460', 'kast_830_3460.sav')
    outroot = 'shane_kast_blue_830.fits'
    templates.build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

# ##############################
def shane_kastr_300_7500_Ne(overwrite=False):  # 300/7500
    """ Build archive for 300/7500
    
    Now uses an un-windowed detector
    """

    binspec = 1
    outroot = 'shane_kast_red_300_7500.fits'
    #
    ifiles = [0]
    slits = [222]
    lcut = []
    #slits = [0]
    #wfile1 = os.path.join(templates.template_path, 'Shane_Kast', '300_7500', 'shane_kast_red_300_7500_ArICdIHeIHgINeI.fits')
    wfile1 = os.path.join(templates.template_path, 'Shane_Kast', '300_7500', 'shane_kast_red_300_7500_full_det.fits')
    #
    templates.build_template([wfile1], slits, lcut, binspec,
                   outroot, lowredux=False, ifiles=ifiles, chk=True,
                   normalize=True, subtract_conti=True, overwrite=overwrite, shift_wave=False)

# ##############################
def shane_kastr_300_7500_NoNe(overwrite=False):  # 300/7500
    """ Warning:  This is *not* the full detector """

    binspec = 1
    outroot = 'shane_kast_red_300_7500_NoNe.fits'
    #
    ifiles = [0]
    slits = [0]
    lcut = []
    wfile1 = os.path.join(templates.template_path, 'Shane_Kast', '300_7500', 'shane_kast_red_300_7500_HeICdIHgIArI.fits')
    #
    templates.build_template([wfile1], slits, lcut, binspec,
                             outroot, lowredux=False, ifiles=ifiles, chk=True,
                             normalize=True, subtract_conti=True, overwrite=overwrite, shift_wave=False)

if __name__ == '__main__':
    shane_kastr_300_7500_Ne()#overwrite=True)
    #shane_kastr_300_7500_NoNe()
