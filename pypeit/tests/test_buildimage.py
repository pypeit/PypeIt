"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os
import pathlib

from IPython import embed

import numpy as np

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph

@dev_suite_required
def test_combine_deimos_flats():

    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                      '830G_L_8400', ifile) 
                            for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    assert len(deimos_flat_files) == 2, 'Incorrect number of files.'

    spectrograph = load_spectrograph('keck_deimos')
    par = spectrograph.default_pypeit_par()
    par['calibrations']['pixelflatframe']['process']['use_biasimage'] = False
    # DEIMOS
    deimos_flat = buildimage.buildimage_fromlist(spectrograph, 3,
                                                  par['calibrations']['pixelflatframe'],
                                                  deimos_flat_files)
    # Process steps
    # Test
    assert deimos_flat.image.shape == (4096,2048)


@dev_suite_required
def test_lris_red_biases():

    path = pathlib.Path(os.getenv('PYPEIT_DEV')).absolute() / 'RAW_DATA' / 'keck_lris_red' \
                / 'multi_600_5000_d560'
    files = [str(path / f) for f in ['LR.20170324.06731.fits.gz', 'LR.20170324.06791.fits.gz',
                                     'LR.20170324.06849.fits.gz', 'LR.20170324.06908.fits.gz',
                                     'LR.20170324.06967.fits.gz']]

    spectrograph = load_spectrograph('keck_lris_red')
    par = spectrograph.default_pypeit_par() 

    bias = buildimage.buildimage_fromlist(spectrograph, 1, par['calibrations']['biasframe'], files)


