"""
Module to run tests on pypeit.outputfiles functions.
"""
from pathlib import Path

from pypeit import outputfiles

_MJD = 58000.5
_TARGET = 'M31'
_CAMERA = 'CAM'
_ROOT = f'{_TARGET}_{_CAMERA}_20170904T120000.000'

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
