""" Generate the wavelength templates for P200/NGPS"""
import os

from pypeit.core.wavecal import templates


def p200_ngps_R(overwrite=False):  # NGPS_R

    binspec = 2 
    outroot = 'p200_ngps_R.fits'
    
    slits = [98, 254, 424] 
    lcut = None

    wpath = os.path.join(templates.template_path, 'P200_NGPS', 'wvcalib_R.fits')
    basefiles = ['wvarxiv_p200_ngps_20250117T1639.fits', 'wvarxiv_p200_ngps_20250117T1639.fits', 'wvarxiv_p200_ngps_20250117T1639.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    templates.build_template(wfiles, slits, lcut, binspec, outroot,
        lowredux=False, normalize=True, overwrite=overwrite)

if __name__ == '__main__':
     p200_ngps_R(overwrite=True)