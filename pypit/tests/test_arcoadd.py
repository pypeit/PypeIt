# Module to run tests on arcoadd

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


def test_find_standard():
    from pypit import arutils as arut
    from pypit import arflux as arflx
    # Dummy self
    slf = arut.dummy_self()
    # G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:52:01.0'
    # Grab
    std_dict = arflx.find_standard_file(slf._argflag, (std_ra, std_dec))
    # Test
    assert std_dict['name'] == 'G191B2B'
    assert std_dict['file'] == '/data/standards/calspec/g191b2b_mod_005.fits'
    assert std_dict['fmt'] == 1
    # Fail to find
    # near G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:22:01.0'
    std_dict = arflx.find_standard_file(slf._argflag, (std_ra,std_dec))
    assert std_dict is None


def test_load_extinction():
    from pypit import arflux as arflx
    from pypit import arutils as arut
    # Dummy self
    slf = arut.dummy_self()
    slf._spect['mosaic']['latitude'] = 37.3413889
    slf._spect['mosaic']['longitude'] = 121.6428
    # Load
    extinct = arflx.load_extinction_data(slf)
    np.testing.assert_allclose(extinct['wave'][0], 3200.)
    assert extinct['wave'].unit == u.AA
    np.testing.assert_allclose(extinct['mag_ext'][0], 1.084)
    # Fail
    slf._spect['mosaic']['latitude'] = 37.3413889
    slf._spect['mosaic']['longitude'] = 0.
    #
    extinct = arflx.load_extinction_data(slf)
    assert extinct is None


def test_extinction_correction():
    from pypit import arflux as arflx
    from pypit import arutils as arut
    # Dummy self
    slf = arut.dummy_self()
    slf._spect['mosaic']['latitude'] = 37.3413889
    slf._spect['mosaic']['longitude'] = 121.6428
    # Load
    extinct = arflx.load_extinction_data(slf)
    # Correction
    wave = np.arange(3000.,10000.)*u.AA
    AM=1.5
    flux_corr = arflx.extinction_correction(wave,AM,extinct)
    # Test
    np.testing.assert_allclose(flux_corr[0], 4.47095192)
