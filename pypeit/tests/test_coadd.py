"""
Module to run tests on arcoadd
"""
import os

import pytest
import numpy as np

from astropy.table import Table

from linetools.spectra.utils import collate
from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit.core import coadd
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.datacube import coadd_cube
from pypeit import utils
from pypeit.tests.tstutils import cooked_required, data_path
from pypeit import inputfiles

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

def test_input_coadd1d_file():
    """ Test I/O of Coadd1D input files """

    # Generate an input file
    coadd1d_input_file = data_path('test.coadd1d')
    if os.path.isfile(coadd1d_input_file):
        os.remove(coadd1d_input_file)

    cfg_lines = ['[coadd1d]']
    cfg_lines += ['  coaddfile = YOUR_OUTPUT_FILE_NAME # Please set your output file name']
    cfg_lines += ['  sensfuncfile = YOUR_SENSFUNC_FILE # Please set your SENSFUNC file name. Only required for Echelle']

    # These files need to be in tests/files/
    data = Table()
    data['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data['obj_id'] = 'SPAT0059-SLIT0175-DET01' # Make believe
    # 
    paths = [data_path('')]

    coadd1dFile = inputfiles.Coadd1DFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data)
    # Write
    coadd1dFile.write(coadd1d_input_file)

    # Read
    coadd1dFile2 = inputfiles.Coadd1DFile.from_file(coadd1d_input_file)
    assert np.all(coadd1dFile2.data['filename'] == data['filename'])

    # Test path
    assert coadd1dFile2.file_paths[0] == paths[0]
    assert coadd1dFile2.filenames[0] == os.path.join(paths[0], data['filename'][0])

    # Other obj_id approach
    data3 = Table()
    data3['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data3['obj_id'] = ['SPAT0059-SLIT0175-DET01', '']

    coadd1dFile3 = inputfiles.Coadd1DFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data3)
    assert coadd1dFile3.objids[1] == data3['obj_id'][0]

def test_input_coadd2d_file():
    """ Test I/O of Coadd2D input files """

    # Generate an input file
    coadd2d_input_file = data_path('test.coadd2d')
    if os.path.isfile(coadd2d_input_file):
        os.remove(coadd2d_input_file)

    cfg_lines = ['[coadd2d]']
    cfg_lines += ['  use_slits4wvgrid = True']

    # These files need to be in tests/files/
    #  The ones below are bogus (i.e. not spec2d files)
    data = Table()
    data['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    # 
    paths = [data_path('')]

    coadd2dFile = inputfiles.Coadd2DFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data)
    # Write
    coadd2dFile.write(coadd2d_input_file)

    # Read
    coadd2dFile2 = inputfiles.Coadd2DFile.from_file(coadd2d_input_file)
    assert np.all(coadd2dFile2.data['filename'] == data['filename'])

    # Test path
    assert coadd2dFile2.file_paths[0] == paths[0]
    assert coadd2dFile2.filenames[0] == os.path.join(paths[0], data['filename'][0])

@cooked_required
def test_coadd_datacube():
    """ Test the coaddition of spec2D files into datacubes """
    droot = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')
    files = [os.path.join(droot,
                          'spec2d_KB.20191219.56886-BB1245p4238_KCWI_20191219T154806.538.fits'),
             os.path.join(droot,
                          'spec2d_KB.20191219.57662-BB1245p4238_KCWI_20191219T160102.755.fits')]
    output_filename = "BB1245p4238_KCWI_20191219.fits"
    # Grab the spectrograph and parset
    spec = load_spectrograph("keck_kcwi")
    parset = spec.default_pypeit_par()
    parset['reduce']['cube']['output_filename'] = output_filename
    parset['reduce']['cube']['combine'] = True
    parset['reduce']['cube']['astrometric'] = False
    parset['reduce']['cube']['grating_corr'] = False
    coadd_cube(files, parset=parset, overwrite=True)
    # Now test the fluxing
    flux_files = [files[0]]
    output_fileflux = "BB1245p4238_KCWI_20191219_fluxing.fits"
    parset['reduce']['cube']['output_filename'] = output_fileflux
    parset['reduce']['cube']['combine'] = False
    parset['reduce']['cube']['standard_cube'] = output_filename
    parset['reduce']['cube']['astrometric'] = False
    parset['reduce']['cube']['grating_corr'] = False
    coadd_cube(flux_files, parset=parset, overwrite=True)
    # Check the files exist
    assert(os.path.exists(output_filename))
    assert(os.path.exists(output_fileflux))
    # Remove the created files
    os.remove(output_filename)
    os.remove(output_fileflux)
