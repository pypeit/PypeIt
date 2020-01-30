"""
Module to run tests on simple fitting routines for arrays
"""
import os
import sys

import numpy as np
import pytest

from astropy import units

from pypeit.core import flux_calib
from pypeit.core import load
from pypeit.spectrographs.util import load_spectrograph
from pypeit import specobjs
from pypeit.tests.tstutils import dummy_fitstbl
from pypeit.pypmsgs import PypeItError


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


# JFH This test is defunct
# TODO: Can it be repurposed to test the relevant functionality? If
# not, delete it.
#def test_bspline_fit():
#    # Testing the bspline works ok (really testing bkspace)
#    fit_dict = linetools.utils.loadjson(data_path('flux_data.json'))
#    wave = np.array(fit_dict['wave'])
#    magfunc = np.array(fit_dict['magf'])
#    logivar = np.array(fit_dict['logiv'])
#    bspline_par = dict(bkspace=fit_dict['bkspec'])
#    mask, tck = utils.robust_polyfit(wave, magfunc, 3, function='bspline',
#                                       weights=np.sqrt(logivar), bspline_par=bspline_par)

# TODO: This needs to be replaced with new tests of SensFunc!!
#def test_gen_sensfunc():
#
#    kastr = load_spectrograph('shane_kast_red')
#
#    # Load a random spectrum for the sensitivity function
#    sfile = data_path('spec1d_r153-J0025-0312_KASTr_2015Jan23T025323.850.fits')
#    sobjs = specobjs.SpecObjs.from_fitsfile(sfile)
##    telescope = telescopes.ShaneTelescopePar()
#    fitstbl = dummy_fitstbl()
#    RA = '05:06:36.6'
#    DEC = '52:52:01.0'
#
#    # Get the sensitivity function
#    sens_dict = flux_calib.generate_sensfunc(sobjs[0].BOX_WAVE,
#                                             sobjs[0].BOX_COUNTS,
#                                             sobjs[0].BOX_COUNTS_IVAR,
#                                             fitstbl['airmass'][4], fitstbl['exptime'][4],
#                                             kastr.telescope['longitude'],
#                                             kastr.telescope['latitude'],
#                                             ra=RA, dec=DEC)
#
#    # Test
#    assert isinstance(sens_dict, dict)
#    assert isinstance(sens_dict['wave_min'], units.Quantity)


def test_find_standard():
    # G191b2b
    std_ra = '05:06:30.6'
    std_dec = '52:49:51.0'
    # Grab
    std_dict = flux_calib.find_standard_file(std_ra, std_dec) 
    # Test
    assert std_dict['name'] == 'G191B2B'
    assert os.path.split(std_dict['cal_file'])[1] == 'g191b2b_stisnic_002.fits'
    assert std_dict['std_source'] == 'calspec'
    # Fail to find
    # near G191b2b
    std_ra = '05:06:36.6'
    std_dec = '52:22:01.0'
    with pytest.raises(PypeItError):
        std_dict = flux_calib.find_standard_file(std_ra, std_dec)


def test_load_extinction():
    # Load
    extinct = flux_calib.load_extinction_data(121.6428, 37.3413889)
    np.testing.assert_allclose(extinct['wave'][0], 3200.)
    assert extinct['wave'].unit == units.AA
    np.testing.assert_allclose(extinct['mag_ext'][0], 1.084)
    # Fail
    extinct = flux_calib.load_extinction_data(0., 37.3413889)
    assert extinct is None


def test_extinction_correction():
    # Load
    extinct = flux_calib.load_extinction_data(121.6428, 37.3413889)
    # Correction
    wave = np.arange(3000.,10000.)*units.AA
    AM=1.5
    flux_corr = flux_calib.extinction_correction(wave, AM, extinct)
    # Test
    np.testing.assert_allclose(flux_corr[0], 4.47095192)


