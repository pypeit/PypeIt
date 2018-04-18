# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import os
import sys

import pdb
import numpy as np
import pytest

try:
    tsterror = FileExistsError
except NameError:
    FileExistsError = OSError

from astropy import units
import linetools.utils

from pypit import arparse as settings
from pypit import arflux
from pypit import arload
from pypit import arutils
from pypit import arsciexp
from pypit import armasters

#from xastropy.xutils import afits as xafits
#from xastropy.xutils import xdebug as xdb

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_bspline_fit():
    # Testing the bspline works ok (really testing bkspace)
    fit_dict = linetools.utils.loadjson(data_path('flux_data.json'))
    wave = np.array(fit_dict['wave'])
    magfunc = np.array(fit_dict['magf'])
    logivar = np.array(fit_dict['logiv'])
    kwargs = dict(bkspace=fit_dict['bkspec'])
    mask, tck = arutils.robust_polyfit(wave, magfunc, 3, function='bspline', weights=np.sqrt(logivar), **kwargs)


def test_gen_sensfunc():
    # Load a random spectrum for the sensitivity function
    sfile = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    specobjs = arload.load_specobj(sfile)
    # Settings, etc.
    settings.dummy_settings()
    settings.argflag['run']['spectrograph'] = 'shane_kast_blue'
    settings.argflag['reduce']['masters']['setup'] = 'C_01_aa'
    settings.spect['arc'] = {}
    settings.spect['arc']['index'] = [[0]]
    fitsdict = arutils.dummy_fitsdict()
    slf = arsciexp.dummy_self()
    slf._msstd[0]['RA'] = '05:06:36.6'
    slf._msstd[0]['DEC'] = '52:52:01.0'
    # Generate
    slf._sensfunc = arflux.generate_sensfunc(slf, 4, [specobjs], fitsdict)
    # Save
    try:
        os.mkdir('MF_shane_kast_blue')
    except FileExistsError:
        pass
    armasters.save_sensfunc(slf, 'C_01_aa')
    # Test
    assert isinstance(slf._sensfunc, dict)
    assert isinstance(slf._sensfunc['wave_min'], units.Quantity)


def test_find_standard():
    # G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:52:01.0'
    # Grab
    std_dict = arflux.find_standard_file((std_ra, std_dec))
    # Test
    assert std_dict['name'] == 'G191B2B'
    assert std_dict['file'] == '/data/standards/calspec/g191b2b_mod_005.fits'
    assert std_dict['fmt'] == 1
    # Fail to find
    # near G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:22:01.0'
    std_dict = arflux.find_standard_file((std_ra,std_dec))
    assert std_dict is None


def test_load_extinction():
    # Dummy self
    settings.spect['mosaic']['latitude'] = 37.3413889
    settings.spect['mosaic']['longitude'] = 121.6428
    # Load
    extinct = arflux.load_extinction_data()
    np.testing.assert_allclose(extinct['wave'][0], 3200.)
    assert extinct['wave'].unit == units.AA
    np.testing.assert_allclose(extinct['mag_ext'][0], 1.084)
    # Fail
    settings.spect['mosaic']['latitude'] = 37.3413889
    settings.spect['mosaic']['longitude'] = 0.
    #
    extinct = arflux.load_extinction_data()
    assert extinct is None


def test_extinction_correction():
    # Dummy self
    settings.spect['mosaic']['latitude'] = 37.3413889
    settings.spect['mosaic']['longitude'] = 121.6428
    # Load
    extinct = arflux.load_extinction_data()
    # Correction
    wave = np.arange(3000.,10000.)*units.AA
    AM=1.5
    flux_corr = arflux.extinction_correction(wave,AM,extinct)
    # Test
    np.testing.assert_allclose(flux_corr[0], 4.47095192)

