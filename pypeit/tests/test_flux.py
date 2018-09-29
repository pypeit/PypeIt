# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import os
import sys

import numpy as np
import pytest

try:
    tsterror = FileExistsError
except NameError:
    FileExistsError = OSError

from astropy import units
import linetools.utils

from pypeit.core import flux
from pypeit.core import load
from pypeit import utils
from pypeit import metadata
from pypeit import telescopes
from pypeit.spectrographs.util import load_spectrograph

#from xastropy.xutils import afits as xafits
#from xastropy.xutils import xdebug as xdb

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

# JFH This test is defunct
#def test_bspline_fit():
#    # Testing the bspline works ok (really testing bkspace)
#    fit_dict = linetools.utils.loadjson(data_path('flux_data.json'))
#    wave = np.array(fit_dict['wave'])
#    magfunc = np.array(fit_dict['magf'])
#    logivar = np.array(fit_dict['logiv'])
#    bspline_par = dict(bkspace=fit_dict['bkspec'])
#    mask, tck = utils.robust_polyfit(wave, magfunc, 3, function='bspline',
#                                       weights=np.sqrt(logivar), bspline_par=bspline_par)


def test_gen_sensfunc():

    kastb = load_spectrograph('shane_kast_blue')

    # Load a random spectrum for the sensitivity function
    sfile = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    specobjs = load.load_specobj(sfile)
    telescope = telescopes.ShaneTelescopePar()
    fitstbl = metadata.dummy_fitstbl()
    RA = '05:06:36.6'
    DEC = '52:52:01.0'

    # Get the sensitivity function
    sens_dict = flux.generate_sensfunc(specobjs[0][0].boxcar['WAVE'],
                                      specobjs[0][0].boxcar['COUNTS'],
                                      specobjs[0][0].boxcar['COUNTS_IVAR'],
                                      fitstbl['airmass'][4], fitstbl['exptime'][4], kastb,
                                      ra=RA, dec=DEC)

    # Test
    assert isinstance(sens_dict, dict)
    assert isinstance(sens_dict['wave_min'], units.Quantity)


def test_find_standard():
    # G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:52:01.0'
    # Grab
    std_dict = flux.find_standard_file(std_ra, std_dec)
    # Test
    assert std_dict['name'] == 'G191B2B'
    assert std_dict['calibfile'] == '/data/standards/calspec/g191b2b_mod_005.fits'
    assert std_dict['fmt'] == 1
    # Fail to find
    # near G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:22:01.0'
    std_dict = flux.find_standard_file(std_ra,std_dec)
    assert std_dict is None


def test_load_extinction():
    # Load
    extinct = flux.load_extinction_data(121.6428, 37.3413889)
    np.testing.assert_allclose(extinct['wave'][0], 3200.)
    assert extinct['wave'].unit == units.AA
    np.testing.assert_allclose(extinct['mag_ext'][0], 1.084)
    # Fail
    extinct = flux.load_extinction_data(0., 37.3413889)
    assert extinct is None


def test_extinction_correction():
    # Load
    extinct = flux.load_extinction_data(121.6428, 37.3413889)
    # Correction
    wave = np.arange(3000.,10000.)*units.AA
    AM=1.5
    flux_corr = flux.extinction_correction(wave, AM, extinct)
    # Test
    np.testing.assert_allclose(flux_corr[0], 4.47095192)

