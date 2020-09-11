"""
Module to run tests on TraceImage class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '830G_L_8400', ifile)
                for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]


@dev_suite_required
def test_process(deimos_flat_files):
    keck_deimos = load_spectrograph('keck_deimos')
    par = keck_deimos.default_pypeit_par()
    par['calibrations']['traceframe']['process']['use_biasimage'] = False
    # Instantiate
    traceImage = buildimage.buildimage_fromlist(keck_deimos, 1, par['calibrations']['traceframe'],
                                                deimos_flat_files)
    # Run
    assert isinstance(traceImage.image, np.ndarray)
    for key in ['subtract_overscan', 'apply_gain']:
        assert key in traceImage.process_steps

