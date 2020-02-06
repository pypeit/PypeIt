"""
Module to run tests on arcoadd
"""
import os

import pytest
import numpy as np

from astropy import units
from linetools.spectra.utils import collate
from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit.core import coadd
from pypeit.spectrographs.util import load_spectrograph
from pypeit import msgs
from pypeit import utils
from IPython import embed

kast_blue = load_spectrograph('shane_kast_blue')

import warnings
warnings.simplefilter("ignore", UserWarning)

# TODO: Need to rewrite the test for coadd1d. FW commented out most tests at this moment.

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

#@pytest.fixture -- Crashes in pytest 4.0.0 as we call this from inside the test, not the call to the test
def dummy_spectrum(s2n=10., rstate=None, seed=1234, wave=None):
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
    if rstate is None:
        rstate=np.random.RandomState(seed)
    if wave is None:
        wave = np.linspace(4000., 5000., 2000)
    # Create
    flux = np.ones_like(wave)
    sig = np.ones_like(wave) / s2n
    ispec = XSpectrum1D.from_tuple((wave,flux,sig))
    # Noise and append
    spec = ispec.add_noise(rstate=rstate)
    flux, sig, mask = spec.data['flux'], spec.data['sig'], spec.data['flux'].mask
    ivar = utils.inverse(sig**2)
    return flux, ivar, mask

#@pytest.fixture
def dummy_spectra(s2n=10., seed=1234, wvmnx=None, npix=None):
    """ Generate a set of normalized spectra with varying wavelength
    and noise

    Parameters
    ----------
    s2n : float, optional

    Returns
    -------
    dspec : XSpectrum1D

    """
    rstate=np.random.RandomState(seed)
    if wvmnx is None:
        wvmnx = [[5000., 6000.],
                [4000.5, 5800.5],
                [4500.8, 6300.8],
                ]
    if npix is None:
        npix = [1000, 1001, 1100]
    slist = []
    for ii, inpix in enumerate(npix):
        wave = np.linspace(wvmnx[ii][0], wvmnx[ii][1], inpix)
        #flux = np.ones_like(wave)
        #sig = np.ones_like(wave) / s2n
        #spec = XSpectrum1D.from_tuple((wave,flux,sig))
        ## Noise and append
        #slist.append(spec.add_noise(seed=seed))
        slist.append(dummy_spectrum(wave=wave, s2n=s2n, rstate=rstate))
    # Collate
    dspec = collate(slist, masking='edges')
    #
    return dspec

def test_qa():
    """ Test QA """
    if os.getenv('RUN_ALL_PYPIT_TESTS') is None:
        assert True
        return
    # Setup
    #wvmnx = [[5000., 6000.],
    #         [5000.5, 6000.5],
    #         [5000.8, 6000.8],
    #         ]
    #npix = [1000, 1000, 1000]
    dspec = dummy_spectra(s2n=10.)#, wvmnx=wvmnx, npix=npix)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    #TODO rewrite this test
    #coadd.coadd_spectra(dspec, wave_method='concatenate', qafile='tst.pdf')


'''  THIS ARE OLD MODELS
def test_load():
    files = [data_path('spec1d_J2202p1708_KASTb_2015Nov06T024436.08.fits'),
             data_path('spec1d_J2202p1708_KASTb_2015Nov06T031500.09.fits'),
             data_path('spec1d_J2202p1708_KASTb_2015Nov06T034520.54.fits')]
    spectra = coadd.load_spec(files, flux=False)  # Need to add in fluxed test files
    assert spectra.nspec == 3
    # Boxcar too
    spectra = coadd.load_spec(files, extract='box')
'''

'''
Need to rewrite the test for coadd1d.
def test_new_wave_grid():
    # Dummy spectrum
    dspec = dummy_spectra()
    # iref [default]
    iref_wave = coadd1d.new_wave_grid(dspec.data['wave'])
    np.testing.assert_allclose(iref_wave[0], 5000.)
    np.testing.assert_allclose(iref_wave[-1], 6000.)
    # Concatenate
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    np.testing.assert_allclose(cat_wave[0], 4000.5)
    np.testing.assert_allclose(cat_wave[-1], 6300.8)
    # Velocity
    vel_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='velocity')
    np.testing.assert_allclose(vel_wave[0], 4000.5)
    np.testing.assert_allclose(vel_wave[-1], 6300.25691664)
    vel_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='velocity', v_pix=100.)
    np.testing.assert_allclose(vel_wave[0], 4000.5)
    np.testing.assert_allclose(vel_wave[-1], 6300.6820934900243)
    # Pixel
    pix_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='pixel', A_pix=2.5)
    np.testing.assert_allclose(pix_wave[0], 4000.5)
    np.testing.assert_allclose(pix_wave[-1], 6303.0)



def test_median_flux():
    """ Test median flux algorithm """
    from pypeit import coadd as arco
    spec = dummy_spectrum(s2n=10)
    # Put in a bad pixel
    spec.data['flux'][0,500] = 0.
    med_flux, std_flux = coadd.median_flux(spec)
    np.testing.assert_allclose(med_flux, 1.0, atol=0.05)  # Noise is random
    np.testing.assert_allclose(std_flux, 0.095, atol=0.004)  # Noise is random
'''

'''
Need to rewrite the test for coadd1d.

def test_median_ratio_flux():
    """ Test median ratio flux algorithm """
    # Setup
    spec1 = dummy_spectrum(s2n=10)
    spec2 = dummy_spectrum(s2n=10)
    spec2.data['flux'][0,:] *= 2.
    spec = collate([spec1,spec2])
    smask = spec.data['sig'].filled(0.) > 0.
    # Put in a bad pixel
    med_flux = coadd1d.median_ratio_flux(spec, smask, 1, 0)
    np.testing.assert_allclose(med_flux, 0.5, atol=0.05)


def test_sn_weight():
    """ Test sn_weight method """
    #  Very low S/N first
    dspec = dummy_spectra(s2n=0.3, seed=1234)
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    smask = rspec.data['sig'].filled(0.) > 0.
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    ivars = utils.inverse(sigs)
    rms_sn, weights = coadd1d.sn_weights(wave, fluxes, ivars, mask=smask)
    np.testing.assert_allclose(rms_sn[0], 0.318, atol=0.1)  # Noise is random
    #  Low S/N first
    dspec = dummy_spectra(s2n=3., seed=1234)
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    smask = rspec.data['sig'].filled(0.) > 0.
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    ivars = utils.inverse(sigs)
    rms_sn, weights = coadd1d.sn_weights(wave, fluxes, ivars, mask=smask)
    np.testing.assert_allclose(rms_sn[0], 2.934, atol=0.1)  # Noise is random
    #  High S/N now
    dspec2 = dummy_spectra(s2n=10., seed=1234)
    cat_wave = coadd1d.new_wave_grid(dspec2.data['wave'], wave_method='concatenate')
    rspec2 = dspec2.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    smask = rspec2.data['sig'].filled(0.) > 0.
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec2)
    ivars = utils.inverse(sigs)
    rms_sn, weights = coadd1d.sn_weights(wave, fluxes, ivars, mask=smask)
    np.testing.assert_allclose(rms_sn[0], 9.904, atol=0.1)  # Noise is random


def test_scale():
    """ Test scale algorithms """
    # Hand
    dspec = dummy_spectra(s2n=10.)
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    smask = rspec.data['sig'].filled(0.) > 0.
    sv_high = rspec.copy()
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    rms_sn, weights = coadd1d.sn_weights(fluxes, sigs, smask, wave)
    _, _ = coadd1d.scale_spectra(rspec, smask, rms_sn, hand_scale=[3., 5., 10.], scale_method='hand')
    np.testing.assert_allclose(np.median(rspec.flux.value[rspec.sig>0.]), 3., atol=0.01)  # Noise is random
    # Median
    rspec = sv_high.copy()
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    rms_sn, weights = coadd1d.sn_weights(fluxes, sigs, smask, wave)
    _, mthd = coadd1d.scale_spectra(rspec, smask, rms_sn, scale_method='median')
    assert mthd == 'median_flux'
    np.testing.assert_allclose(np.median(rspec.flux.value[rspec.sig>0.]), 1., atol=0.01)  # Noise is random
    #  Auto-none
    dspec = dummy_spectra(s2n=0.1)
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    rms_sn, weights = coadd1d.sn_weights(fluxes, sigs, smask, wave)
    an_scls, an_mthd = coadd1d.scale_spectra(rspec, smask, rms_sn)
    assert an_mthd == 'none_SN'
    #  Auto-median
    dspec = dummy_spectra(s2n=1.5)
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    rspec.data['flux'][1,:] *= 10.
    rspec.data['sig'][1,:] *= 10.
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    rms_sn, weights = coadd1d.sn_weights(fluxes, sigs, smask, wave)
    am_scls, am_mthd = coadd1d.scale_spectra(rspec, smask, rms_sn, scale_method='median')
    assert am_mthd == 'median_flux'
    np.testing.assert_allclose(am_scls[1], 0.1, atol=0.01)


def test_grow_mask():
    """ Test grow_mask method.  Now works on 1d spectra"""
    # Setup
    dspec = dummy_spectrum(s2n=10.)
    mask = np.ma.getmaskarray(dspec.data['wave'][0,:])
    mask[:] = True
    # Set some
    mask[0] = False
    mask[50] = False
    mask[500] = False
    # Grow
    new_mask = coadd1d.grow_mask(mask, n_grow=1)
    # Test
    badp = np.where(new_mask == False)[0]
    assert np.all(badp == np.array([0,1,49,50,51,499,500,501]))
    #badpb = np.where(new_mask[1,:])[0]
    #assert np.all(badpb == np.array([0,1]))
    #badpc = np.where(new_mask[2,:])[0]
    #assert np.all(badpc == np.array([1098,1099]))
    # Grow 2
    new_mask2 = coadd1d.grow_mask(mask, n_grow=2)
    badp2 = np.where(new_mask2 == False)[0]
    assert len(badp2) == 13


def test_1dcoadd():
    """ Test 1dcoadd method"""
    # Setup
    dspec = dummy_spectra(s2n=10.)
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    smask = rspec.data['sig'].filled(0.) > 0.
    fluxes, sigs, wave = coadd1d.unpack_spec(rspec)
    rms_sn, weights = coadd1d.sn_weights(fluxes, sigs, smask, wave)
    # Coadd
    spec1d = coadd1d.one_d_coadd(rspec, smask, weights)
    assert spec1d.npix == 1740

def test_cleancr():
    """ Test clean CR method"""
    # Setup
    dspec = dummy_spectra(s2n=10.)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True, masking='none')
    #
    smask = rspec.data['sig'].filled(0.) <= 0.
    coadd1d.clean_cr(rspec, smask)


def test_coadd():
    """ Test full coadd method"""
    # Setup
    dspec = dummy_spectra(s2n=10.)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    spec1d = coadd1d.coadd_spectra(kast_blue, None, dspec, wave_grid_method='concatenate')
    assert np.isclose(np.median(spec1d.flux.value), 1., atol=0.003)


# TODO: This needs to be fixed.
def test_coadd_with_fluxing():
    """ Test full coadd method with flux scaling"""
    scale_dict = dict(filter='DECAM-R', mag=19.0, mag_type='AB')
    # Setup
    dspec = dummy_spectra(s2n=10.)
    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
    dspec.data['sig'][0, 700] *= 500.
    spec1d = coadd1d.coadd_spectra(kast_blue, None, dspec, wave_grid_method='concatenate', flux_scale=scale_dict)
    # Test
    assert np.median(spec1d.flux.value) > 6.61


def test_coadd_qa():
    if os.getenv('PYPIT') is None:
        assert True
        return
'''

'''
def test_sigma_clip():
    """ Test sigma_clip method """
    from pypeit import coadd as arco
    # Setup
    dspec = dummy_spectra(s2n=10.)
    cat_wave = coadd1d.new_wave_grid(dspec.data['wave'], wave_method='concatenate')
    # Rebin
    rspec = dspec.rebin(cat_wave*units.AA, all=True, do_sig=True)
    sn2, weights = coadd1d.sn_weight(cat_wave, rspec.data['flux'], rspec.data['sig']**2)
    # Here we go
    rspec.data['flux'][0, 700] = 999.
    final_mask = coadd1d.sigma_clip(rspec.data['flux'], rspec.data['sig']**2, sn2=sn2)


'''
