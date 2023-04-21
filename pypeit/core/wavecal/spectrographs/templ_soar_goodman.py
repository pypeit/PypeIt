""" Generate the wavelength templates for SOAR Goodman"""
import os

from pypeit.core.wavecal import templates

from IPython import embed

def soar_goodman_400(overwrite=False):
    binspec = 2
    outroot = 'soar_goodman_red_400_SYZY.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'SOAR_Goodman', '400_SYZY')

    basefiles = ['MasterWaveCalib_A_1_01_M1.fits', 'MasterWaveCalib_A_1_01_M2.fits'] 
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]

    # Snippets
    ifiles = [0,1]
    slits = [495, 496]
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


if __name__ == '__main__':
    soar_goodman_400(overwrite=True)