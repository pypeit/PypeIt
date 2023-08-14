""" Generate the wavelength templates for Keck/LRIS Blue"""
import os

from pypeit.core.wavecal import templates



def keck_lris_blue_B300_5000(overwrite=False):
    binspec = 1
    outroot = 'keck_lris_red_B300_5000_d680_ArCdHgKrNeXeZnFeAr.fits'
    # PypeIt fits
    wpath = os.path.join(templates.template_path, 'Keck_LRIS', 'keck_lris_blue', 'B300_5000')

    basefiles = ['WaveCalib_A_0_DET01_S1676.fits', 'WaveCalib_A_0_DET01_S0947.fits']
    wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    # Snippets
    ifiles = [0,1]
    slits = [1676, 947]
    wv_cuts = [6830.]
    assert len(wv_cuts) == len(slits)-1
    # det_dict
    det_cut = None
    #
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             ifiles=ifiles, det_cut=det_cut, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)


# Run em
if __name__ == '__main__':
    keck_lris_blue_B300_5000(overwrite=False)
