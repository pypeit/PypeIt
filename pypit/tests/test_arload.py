# Module to run tests on arload

### TEST_UNICODE_LITERALS

import numpy as np
import os
import pytest

from astropy import units as u

from pypit import pyputils
import pypit
msgs = pyputils.get_dummy_logger()


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_load_1dspec():
    from pypit import arload as arl
    from linetools.spectra.xspectrum1d import XSpectrum1D

    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    spec = arl.load_1dspec(spec_file)
    # Test
    assert isinstance(spec, XSpectrum1D)
    # Boxcar
    spec = arl.load_1dspec(spec_file, extract='box')
    assert isinstance(spec, XSpectrum1D)


