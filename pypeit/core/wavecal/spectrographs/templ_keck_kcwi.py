""" Generate the wavelength templates for Keck KCWI"""
import os

from pypeit.core.wavecal import templates


def keck_kcwi_BL(overwrite=False):
    # FeAr BL
    wfile1 = os.path.join(templates.template_path, 'KCWI', 'BL', 'Keck_KCWI_BL_4500.fits')
    outroot = 'keck_kcwi_BL.fits'
    binspec = 2
    slits = [45]
    lcut = [3000.0, 8000.0]
    templates.build_template([wfile1], slits, lcut, binspec, outroot,
                   lowredux=False, overwrite=overwrite, normalize=True)

def keck_kcwi_BM(overwrite=False):
    # FeAr BM
    wfile1 = os.path.join(templates.template_path, 'KCWI', 'BM', 'Keck_KCWI_BM_4060.json')
    wfile2 = os.path.join(templates.template_path, 'KCWI', 'BM', 'Keck_KCWI_BM_4670.json')
    outroot = 'keck_kcwi_BM.fits'
    binspec = 1
    slits = [1026, 1021]
    lcut = [4350.0, 8000.0]
    templates.build_template([wfile1, wfile2], slits, lcut, binspec, outroot,
                   lowredux=False, overwrite=overwrite, normalize=True)

def keck_kcwi_BH2(overwrite=False):
    # FeAr BH2
    wfile1 = os.path.join(templates.template_path, 'KCWI', 'BH2', 'Keck_KCWI_BH2_4200.json')
    outroot = 'keck_kcwi_BH2.fits'
    binspec = 1
    slits = [1015]
    lcut = [4350.0, 8000.0]
    templates.build_template([wfile1], slits, lcut, binspec, outroot,
                             lowredux=False, overwrite=overwrite, normalize=True)


if __name__ == '__main__':
    keck_kcwi_BL(overwrite=True)
    keck_kcwi_BM(overwrite=True)
    keck_kcwi_BH2(overwrite=True)
