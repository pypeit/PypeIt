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
from pypeit.core.datacube import coadd_cube
from pypeit import msgs
from pypeit import utils
from IPython import embed
from pypeit.tests.tstutils import cooked_required, data_path

import warnings
warnings.simplefilter("ignore", UserWarning)

# TODO: Need to rewrite the test for coadd1d. FW commented out most tests at this moment.


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

# TODO: This test needs to be re-written.  It currently does nothing, so I
# removed it.
#def test_qa():
#    """ Test QA """
#    if os.getenv('RUN_ALL_PYPIT_TESTS') is None:
#        assert True
#        return
#    # Setup
#    #wvmnx = [[5000., 6000.],
#    #         [5000.5, 6000.5],
#    #         [5000.8, 6000.8],
#    #         ]
#    #npix = [1000, 1000, 1000]
#    dspec = dummy_spectra(s2n=10.)#, wvmnx=wvmnx, npix=npix)
#    dspec.data['flux'][0, 700] *= 1000.  # One bad pixel
#    dspec.data['sig'][0, 700] *= 500.
#    #TODO rewrite this test
#    #coadd.coadd_spectra(dspec, wave_method='concatenate', qafile='tst.pdf')
#

@cooked_required
def test_coadd_datacube():
    """ Test the coaddition of spec2D files into datacubes """
    droot = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')
    files = [os.path.join(droot,
                          'spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits'),
             os.path.join(droot,
                          'spec2d_KB.20191219.57662-BB1245p4238_KCWI_20191219T160102.755.fits')]
    coadd_cube(files, overwrite=True)
    os.remove('datacube.fits')


