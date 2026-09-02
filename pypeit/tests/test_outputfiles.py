"""
Module to run tests on pypeit.outputfiles functions.
"""
from pathlib import Path

from pypeit import outputfiles

_MJD = 58000.5
_TARGET = 'M31'
_CAMERA = 'CAM'
_ROOT = f'{_TARGET}_{_CAMERA}_20170904T120000.000'

def test_strip_raw_extension_fits():
    assert outputfiles.strip_raw_extension('b27.fits', ['.fits']) == 'b27'

def test_strip_raw_extension_fits_gz():
    assert outputfiles.strip_raw_extension('b27.fits.gz', ['.fits', '.fits.gz']) == 'b27'

def test_strip_raw_extension_fits_bz2():
    # Gemini GMOS-style compressed raw file extension
    assert outputfiles.strip_raw_extension('S0183.fits.bz2', ['.fits', '.fits.bz2']) == 'S0183'

def test_strip_raw_extension_non_fits_extension():
    # SOAR Goodman-style raw file extension that does not contain '.fits'
    assert outputfiles.strip_raw_extension('b27.fz', ['.fz']) == 'b27'

def test_strip_raw_extension_path_input():
    path = Path('/some/dir/b27.fits.gz')
    assert outputfiles.strip_raw_extension(path, ['.fits', '.fits.gz']) == 'b27'

def test_strip_raw_extension_unrecognized_extension():
    # No extension in allowed_extensions matches; the full name is kept and a
    # warning is issued, but the function should not raise an exception.
    assert outputfiles.strip_raw_extension('b27.dat', ['.fits']) == 'b27.dat'

def test_construct_basename_fits():
    basename = outputfiles.construct_basename('b27.fits', _TARGET, _CAMERA, _MJD, ['.fits'])
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_fits_gz():
    basename = outputfiles.construct_basename('b27.fits.gz', _TARGET, _CAMERA, _MJD,
                                               ['.fits', '.fits.gz'])
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_non_fits_extension():
    # SOAR Goodman-style raw file extension that does not contain '.fits'
    basename = outputfiles.construct_basename('b27.fz', _TARGET, _CAMERA, _MJD, ['.fz'])
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_path_input():
    basename = outputfiles.construct_basename(Path('/some/dir/b27.fits.gz'), _TARGET, _CAMERA,
                                               _MJD, ['.fits', '.fits.gz'])
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_unrecognized_extension():
    # No extension in allowed_extensions matches; the full name is kept and a
    # warning is issued, but the function should not raise an exception.
    basename = outputfiles.construct_basename('b27.dat', _TARGET, _CAMERA, _MJD, ['.fits'])
    assert basename == f'b27.dat-{_ROOT}'
