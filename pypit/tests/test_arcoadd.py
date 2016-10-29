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
def dummy_spectrum(s2n=10., seed=1234, wave=None):
    """
    Parameters
    ----------
    s2n
    seed
    wave

    Returns
    -------
    spec : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    if wave is None:
        wave = np.linspace(4000., 5000., 2000)
    # Create
    flux = np.ones_like(wave)
    sig = np.ones_like(wave) / s2n
    ispec = XSpectrum1D.from_tuple((wave,flux,sig))
    # Noise and append
    spec = ispec.add_noise(seed=seed)
    return spec

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
    from linetools.spectra.utils import collate
    wvmnx = [[5000., 6000.],
            [4000.5, 5800.5],
            [4500.8, 6300.8],
            ]
    npix = [1000, 1001, 1100]
    slist = []
    for ii, inpix in enumerate(npix):
        wave = np.linspace(wvmnx[ii][0], wvmnx[ii][1], inpix)
        #flux = np.ones_like(wave)
        #sig = np.ones_like(wave) / s2n
        #spec = XSpectrum1D.from_tuple((wave,flux,sig))
        ## Noise and append
        #slist.append(spec.add_noise(seed=seed))
        slist.append(dummy_spectrum(wave=wave, s2n=s2n, seed=seed))
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


def test_median_flux():
    """ Test median flux algorithm """
    from pypit import arcoadd as arco
    spec = dummy_spectrum(s2n=10)
    med_flux, std_flux = arco.median_flux(spec)
    # Put in a bad pixel
    spec.data['flux'][0,500] = 0.
    med_flux, std_flux = arco.median_flux(spec)
    np.testing.assert_allclose(med_flux, 1.0, atol=0.05)  # Noise is random
    np.testing.assert_allclose(std_flux, 0.095, atol=0.004)  # Noise is random


def test_sn_weight():
    """ Test sn_weight method """
    from pypit import arcoadd as arco
    #  Very low S/N first
    dspec = dummy_spectra(s2n=0.3, seed=1234)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sn2, weights = arco.sn_weight(rspec)
    np.testing.assert_allclose(sn2[0], 0.095, atol=0.1)  # Noise is random
    #  Low S/N first
    dspec = dummy_spectra(s2n=3., seed=1234)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sn2, weights = arco.sn_weight(rspec)
    np.testing.assert_allclose(sn2[0], 8.6, atol=0.1)  # Noise is random
    #  High S/N now
    dspec2 = dummy_spectra(s2n=10., seed=1234)
    cat_wave = arco.new_wave_grid(dspec2.data['wave'], method='concatenate')
    rspec2 = dspec2.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sn2, weights = arco.sn_weight(rspec2)
    np.testing.assert_allclose(sn2[0], 97.9, atol=0.1)  # Noise is random


def test_scale():
    """ Test scale algorithms """
    from pypit import arcoadd as arco
    # Hand
    dspec = dummy_spectra(s2n=10.)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sv_high = rspec.copy()
    sn2, weights = arco.sn_weight(rspec)
    _, _ = arco.scale_spectra(rspec, sn2, hand_scale=[3., 5., 10.], method='hand')
    np.testing.assert_allclose(np.median(rspec.flux.value), 3., atol=0.01)  # Noise is random
    # Median
    rspec = sv_high.copy()
    sn2, weights = arco.sn_weight(rspec)
    _, mthd = arco.scale_spectra(rspec, sn2, method='median')
    assert mthd == 'median'
    np.testing.assert_allclose(np.median(rspec.flux.value), 1., atol=0.01)  # Noise is random
    #  Auto-none
    dspec = dummy_spectra(s2n=0.1)
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sn2, weights = arco.sn_weight(rspec)
    an_scls, an_mthd = arco.scale_spectra(rspec, sn2)
    assert an_mthd == 'none_SN'
    #  Auto-median
    dspec = dummy_spectra(s2n=1.5)
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    rspec.data['flux'][1,:] *= 10.
    rspec.data['sig'][1,:] *= 10.
    sn2, weights = arco.sn_weight(rspec)
    am_scls, am_mthd = arco.scale_spectra(rspec, sn2)
    assert am_mthd == 'median'
    np.testing.assert_allclose(am_scls[1], 0.1, atol=0.01)


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


def test_1dcoadd():
    """ Test 1dcoadd method"""
    from pypit import arcoadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    sn2, weights = arco.sn_weight(rspec)
    # Coadd
    spec1d = arco.one_d_coadd(rspec, weights)
    assert spec1d.npix == 1740

def test_cleancr():
    """ Test clean CR method"""
    from pypit import arcoadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    cat_wave = arco.new_wave_grid(dspec.data['wave'], method='concatenate')
    rspec = dspec.rebin(cat_wave*u.AA, all=True, do_sig=True, masking='none')
    arco.clean_cr(rspec)

def test_coadd():
    """ Test full coadd method"""
    from pypit import arcoadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    arco.coadd_spectra(dspec, wave_grid_method='concatenate')


'''
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


'''
