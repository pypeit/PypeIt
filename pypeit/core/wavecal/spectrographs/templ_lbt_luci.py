""" Generate the wavelength templates for LBT LUCI"""
import os

from pypeit.core.wavecal import templates



# LBT LUCI1

def lbt_luci1_g200_zJ(overwrite=False, dev_path=None):
    """Build the wavelength template for LBT LUCI.

    :param overwrite: Boolean to indicate whether the reid_arxiv solution
     should be overwritten or not.
    :type overwrite: bool
    :param dev_path: Path to the development suite
    :type dev_path: string
    :return:
    """

    binspec = 1
    outroot = 'lbt_luci1_g200_zJ.fits'

    if dev_path is not None:
        dev_path = os.path.join(dev_path, 'dev_algorithms/wavelengths/template_files')
    else:
        dev_path = templates.template_path

    wfile1 = os.path.join(dev_path,
                          'LBT', 'LUCI1', 'G200_zJ',
                          'LBT_LUCI1_G200_1p1699999571_date20211112.fits')

    slits = [1020]  # assuming this is spat_id
    wvcuts = [9500, 13750]  # lower and upper wavelength boundary for each slit


    templates.build_template([wfile1], slits, wvcuts, binspec, outroot,
                             lowredux=False, normalize=True,
                             subtract_conti=True,
                             overwrite=overwrite, chk=True)


def lbt_luci1_g200_HK(overwrite=False, dev_path=None):
    """Build the wavelength template for LBT LUCI.

    :param overwrite: Boolean to indicate whether the reid_arxiv solution
     should be overwritten or not.
    :type overwrite: bool
    :param dev_path: Path to the development suite
    :type dev_path: string
    :return:
    """

    binspec = 1
    outroot = 'lbt_luci1_g200_HK.fits'

    if dev_path is not None:
        dev_path = os.path.join(dev_path, 'dev_algorithms/wavelengths/template_files')
    else:
        dev_path = templates.template_path

    wfile1 = os.path.join(dev_path,
                          'LBT', 'LUCI1', 'G200_HK',
                          'LBT_LUCI1_G200_1p9299999475_date20211112.fits')

    slits = [1026]  # assuming this is spat_id
    wvcuts = [14700, 23500]  # lower and upper wavelength boundary for each slit

    templates.build_template([wfile1], slits, wvcuts, binspec, outroot,
                             lowredux=False, normalize=True,
                             subtract_conti=True,
                             overwrite=overwrite, chk=True)


# LBT LUCI2

def lbt_luci2_g200_zJ(overwrite=False, dev_path=None):
    """Build the wavelength template for LBT LUCI.

    :param overwrite: Boolean to indicate whether the reid_arxiv solution
     should be overwritten or not.
    :type overwrite: bool
    :param dev_path: Path to the development suite
    :type dev_path: string
    :return:
    """

    binspec = 1
    outroot = 'lbt_luci2_g200_zJ.fits'

    if dev_path is not None:
        dev_path = os.path.join(dev_path, 'dev_algorithms/wavelengths/template_files')
    else:
        dev_path = templates.template_path

    wfile1 = os.path.join(dev_path,
                          'LBT', 'LUCI2', 'G200_zJ',
                          'LBT_LUCI2_G200_1p1699999571_date20211112.fits')

    slits = [1071]  # assuming this is spat_id
    wvcuts = [9550, 13850]  # lower and upper wavelength boundary for each slit


    templates.build_template([wfile1], slits, wvcuts, binspec, outroot,
                             lowredux=False, normalize=True,
                             subtract_conti=True,
                             overwrite=overwrite, chk=True)


def lbt_luci2_g200_HK(overwrite=False, dev_path=None):
    """Build the wavelength template for LBT LUCI.

    :param overwrite: Boolean to indicate whether the reid_arxiv solution
     should be overwritten or not.
    :type overwrite: bool
    :param dev_path: Path to the development suite
    :type dev_path: string
    :return:
    """

    binspec = 1
    outroot = 'lbt_luci2_g200_HK.fits'

    if dev_path is not None:
        dev_path = os.path.join(dev_path, 'dev_algorithms/wavelengths/template_files')
    else:
        dev_path = templates.template_path

    wfile1 = os.path.join(dev_path,
                          'LBT', 'LUCI2', 'G200_HK',
                          'LBT_LUCI2_G200_1p9299999475_date20211112.fits')

    slits = [1069]  # assuming this is spat_id
    wvcuts = [14950, 23750]  # lower and upper wavelength boundary for each slit

    templates.build_template([wfile1], slits, wvcuts, binspec, outroot,
                             lowredux=False, normalize=True,
                             subtract_conti=True,
                             overwrite=overwrite, chk=True)

if __name__ == '__main__':

    # Insert your Pypeit development suite path here
    dev_path = '/Users/schindler/Software/PypeIt-development-suite/'

    lbt_luci1_g200_zJ(overwrite=True,
                      dev_path=dev_path)

    lbt_luci1_g200_HK(overwrite=True,
                      dev_path=dev_path)

    lbt_luci2_g200_zJ(overwrite=True,
                      dev_path=dev_path)

    lbt_luci2_g200_HK(overwrite=True,
                      dev_path=dev_path)