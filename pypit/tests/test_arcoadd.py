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

@pytest.fixture
def dummy_spectra(s2n=10., seed=1234):
    """ Generate a set of normalized spectra with varying wavelength
    and noise

    Parameters
    ----------
    s2n : float, optional

    Returns
    -------
    dspec : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    from linetools.spectra.utils import collate
    wvmnx = [[5000., 6000.],
            [4000.5, 5800.5],
            [4500.8, 6300.8],
            ]
    npix = [1000, 1001, 1100]
    slist = []
    for ii, inpix in enumerate(npix):
        wave = np.linspace(wvmnx[ii][0], wvmnx[ii][1], inpix)
        flux = np.ones_like(wave)
        sig = np.ones_like(wave) / s2n
        spec = XSpectrum1D.from_tuple((wave,flux,sig))
        # Noise and append
        slist.append(spec.add_noise(seed=seed))
    # Collate
    dspec = collate(slist)
    #
    return dspec

def test_load():
    from pypit import arcoadd as arco
    files = [data_path('spec1d_J2202p1708_KASTb_2015Nov06T024436.08.fits'),
             data_path('spec1d_J2202p1708_KASTb_2015Nov06T031500.09.fits'),
             data_path('spec1d_J2202p1708_KASTb_2015Nov06T034520.54.fits')]
    spectra = arco.load_spec(files)
    assert spectra.nspec == 3
    # Boxcar too
    spectra = arco.load_spec(files, extract='box')


def test_new_wave_grid():
    from pypit import arcoadd as arco
    # Dummy spectrum
    dspec = dummy_spectra()
    # iref [default]
    iref_wave = arco.new_wave_grid(dspec.data['wave'])
    np.testing.assert_allclose(iref_wave[0], 5000.)
    np.testing.assert_allclose(iref_wave[-1], 6000.)
    # Concatenate
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    np.testing.assert_allclose(cat_wave[0], 4000.5)
    np.testing.assert_allclose(cat_wave[-1], 6300.8)
    # Velocity
    vel_wave = arco.new_wave_grid(dspec.data['wave'], method='velocity')
    np.testing.assert_allclose(vel_wave[0], 4000.5)
    np.testing.assert_allclose(vel_wave[-1], 6302.7837748108632)
    # Pixel
    pix_wave = arco.new_wave_grid(dspec.data['wave'], method='pixel', pix_size=2.5)
    np.testing.assert_allclose(pix_wave[0], 4000.5)
    np.testing.assert_allclose(pix_wave[-1], 6303.15)


def test_sn_weight():
    """ Test sn_weight method """
    from pypit import arcoadd as arco
    #  Low S/N first
    dspec = dummy_spectra(s2n=3.)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True)
    sn2, weights = arco.sn_weight(cat_wave, rspec.data['flux'], rspec.data['sig']**2)
    np.testing.assert_allclose(sn2[0], 8.85, atol=0.1)  # Noise is random
    #  High S/N now
    dspec2 = dummy_spectra(s2n=10.)
    cat_wave = arco.new_wave_grid(dspec2.data['wave'], method='concatenate')
    rspec2 = dspec2.rebin(cat_wave*u.AA, all=True, do_sig=True)
    sn2, weights = arco.sn_weight(cat_wave, rspec2.data['flux'], rspec2.data['sig']**2)
    np.testing.assert_allclose(sn2[0], 98.3, atol=0.1)  # Noise is random


def test_grow_mask():
    """ Test grow_mask method"""
    from pypit import arcoadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    mask = np.ma.getmaskarray(dspec.data['wave'])
    mask[:] = False
    # Set some
    mask[0, 100] = True
    mask[1, 0] = True
    mask[2, -1] = True
    # Grow
    new_mask = arco.grow_mask(mask, n_grow=1)
    # Test
    badp = np.where(new_mask[0,:])[0]
    assert np.all(badp == np.array([99,100,101]))
    badpb = np.where(new_mask[1,:])[0]
    assert np.all(badpb == np.array([0,1]))
    badpc = np.where(new_mask[2,:])[0]
    assert np.all(badpc == np.array([1098,1099]))
    # Grow 2
    new_mask2 = arco.grow_mask(mask, n_grow=2)
    badp2 = np.where(new_mask2[0,:])[0]
    assert np.all(badp2 == np.array([98,99,100,101,102]))


def test_sigma_clip():
    """ Test sigma_clip method """
    from pypit import arcoadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    # Rebin
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True)
    sn2, weights = arco.sn_weight(cat_wave, rspec.data['flux'], rspec.data['sig']**2)
    # Here we go
    rspec.data['flux'][0, 700] = 999.
    final_mask = arco.sigma_clip(rspec.data['flux'], rspec.data['sig']**2, sn2=sn2)


