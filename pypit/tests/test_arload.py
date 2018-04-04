# Module to run tests on arload
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import os
import pytest

from pypit import arutils
from pypit import arload
from pypit import pyputils

# TODO: Not used...
#import numpy as np
#from astropy import units as u
#import pypit

msgs = pyputils.get_dummy_logger()

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_load_headers():
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    kast_files = [data_path('b1.fits.gz'), data_path('b27.fits.gz')]
    fitsdict, updates = arload.load_headers(kast_files)
    # Test
    headers = fitsdict['headers']
    assert len(headers) == 2
    assert headers[0][0]['OBJECT'] == 'Arcs'

def test_load_specobj():
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    specobjs = arload.load_specobj(spec_file)
    # Test
    assert isinstance(specobjs, list)
    assert len(specobjs[0].boxcar['counts']) == 1199

def test_load_1dspec():
    from linetools.spectra.xspectrum1d import XSpectrum1D

    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    spec = arload.load_1dspec(spec_file)
    # Test
    assert isinstance(spec, XSpectrum1D)
    # Boxcar
    spec = arload.load_1dspec(spec_file, extract='box')
    assert isinstance(spec, XSpectrum1D)
    # By objname
    spec2 = arload.load_1dspec(spec_file, objname='O473-S5473-D01-I0008')
    assert isinstance(spec2, XSpectrum1D)


