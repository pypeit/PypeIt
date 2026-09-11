"""
Module to run tests on pypeit.outputfiles functions.
"""
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

from pypeit import inputfiles
from pypeit import outputfiles
from pypeit.spectrographs.util import load_spectrograph

_MJD = 58000.5
_TARGET = 'M31'
_CAMERA = 'CAM'
_ROOT = f'{_TARGET}_{_CAMERA}_20170904T120000.000'


def _write_spec2d(path, target='J0750+6927'):
    hdr = fits.Header()
    hdr['TARGET'] = target
    hdr['PYP_SPEC'] = 'keck_kcrm'
    hdr['PYPELINE'] = 'SlicerIFU'
    fits.PrimaryHDU(header=hdr).writeto(path)

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
    basename = outputfiles.construct_basename('b27.fits', _CAMERA, ['.fits'], target=_TARGET,
                                               mjd=_MJD)
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_fits_gz():
    basename = outputfiles.construct_basename('b27.fits.gz', _CAMERA, ['.fits', '.fits.gz'],
                                               target=_TARGET, mjd=_MJD)
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_non_fits_extension():
    # SOAR Goodman-style raw file extension that does not contain '.fits'
    basename = outputfiles.construct_basename('b27.fz', _CAMERA, ['.fz'], target=_TARGET,
                                               mjd=_MJD)
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_path_input():
    basename = outputfiles.construct_basename(Path('/some/dir/b27.fits.gz'), _CAMERA,
                                               ['.fits', '.fits.gz'], target=_TARGET, mjd=_MJD)
    assert basename == f'b27-{_ROOT}'

def test_construct_basename_unrecognized_extension():
    # No extension in allowed_extensions matches; the full name is kept and a
    # warning is issued, but the function should not raise an exception.
    basename = outputfiles.construct_basename('b27.dat', _CAMERA, ['.fits'], target=_TARGET,
                                               mjd=_MJD)
    assert basename == f'b27.dat-{_ROOT}'


def test_spec2d_target_prefers_target_over_targname_and_object(tmp_path):
    path = tmp_path / 'spec2d_test.fits'
    hdr = fits.Header()
    hdr['TARGET'] = 'M31'
    hdr['TARGNAME'] = 'not-this-one'
    hdr['OBJECT'] = 'not-this-one-either'
    fits.PrimaryHDU(header=hdr).writeto(path)
    assert outputfiles.spec2d_target(path) == 'M31', \
        'TARGET should take precedence over TARGNAME and OBJECT'


def test_spec2d_target_falls_back_to_targname_then_object(tmp_path):
    path = tmp_path / 'spec2d_test.fits'
    hdr = fits.Header()
    hdr['OBJECT'] = 'M31'
    fits.PrimaryHDU(header=hdr).writeto(path)
    assert outputfiles.spec2d_target(path) == 'M31', \
        'should fall back to OBJECT when TARGET and TARGNAME are absent'


def test_spec2d_target_none_when_no_keyword_present(tmp_path):
    path = tmp_path / 'spec2d_test.fits'
    fits.PrimaryHDU().writeto(path)
    assert outputfiles.spec2d_target(path) is None, \
        'should return None when no TARGET/TARGNAME/OBJECT keyword is present'


def test_find_reduced_spec2d_exact_match_primary_path(tmp_path):
    # The primary path should find the file from the row's own metadata alone, with
    # no header read -- so a deliberately WRONG TARGET header must not prevent a match.
    spec = load_spectrograph('keck_kcrm')
    row = Table(
        {'filename': ['kr260610_00054.fits.gz'], 'target': ['J0750+6927'], 'mjd': [59742.123456]}
    )[0]
    expected_basename = outputfiles.construct_basename(
        row['filename'], spec.camera, spec.allowed_extensions, target=row['target'],
        mjd=row['mjd']
    )
    spec2d_path = tmp_path / f'spec2d_{expected_basename}.fits'
    _write_spec2d(spec2d_path, target='not-the-requested-target')

    assert outputfiles.find_reduced_spec2d(tmp_path, row, spec) == spec2d_path, \
        'exact-match path should find the file by reconstructed name, ignoring the header TARGET'


def test_find_reduced_spec2d_fallback_path(tmp_path):
    # Name the spec2d file with a compound raw extension (.fits.gz) and a timestamp
    # token that does NOT match construct_basename's reconstruction from `row['mjd']`
    # (simulating an mjd-precision miss), so only the raw-stem-prefix glob, verified
    # against the spec2d header's TARGET, can find it.
    spec = load_spectrograph('keck_kcrm')
    row = Table(
        {'filename': ['kr260610_00058.fits.gz'], 'target': ['J0750+6927'], 'mjd': [59742.654321]}
    )[0]
    raw_stem = outputfiles.strip_raw_extension(row['filename'], spec.allowed_extensions)
    spec2d_path = tmp_path / f'spec2d_{raw_stem}-J0750+6927_KCRM_notarealtimestamp.fits'
    _write_spec2d(spec2d_path, target='J0750+6927')

    assert outputfiles.find_reduced_spec2d(tmp_path, row, spec) == spec2d_path, \
        'exact match should miss (timestamp mismatch); should be found via the glob+header fallback'

    # A row with no mjd column at all still falls back correctly (this is the
    # original bug: a naive `Path(filename).stem` would have left a stray `.fits` on
    # the glob prefix for a `.fits.gz` raw file).
    no_mjd_row = Table({'filename': ['kr260610_00058.fits.gz'], 'target': ['J0750+6927']})[0]
    assert outputfiles.find_reduced_spec2d(tmp_path, no_mjd_row, spec) == spec2d_path, \
        'a row with no mjd column should go straight to the fallback and still find the file'


def test_find_reduced_spec2d_missing_returns_none(tmp_path):
    spec = load_spectrograph('keck_kcrm')
    row = Table(
        {'filename': ['kr260610_99999.fits.gz'], 'target': ['J0750+6927'], 'mjd': [59742.0]}
    )[0]
    assert outputfiles.find_reduced_spec2d(tmp_path, row, spec) is None, \
        'should return None when no candidate spec2d file exists on disk'


def test_existing_spec2d_files_found_and_missing(tmp_path):
    pypeit_file = inputfiles.PypeItFile(
        config={'rdx': {'spectrograph': 'keck_kcrm'}},
        file_paths=['/tmp/raw'],
        data_table=Table({
            'filename': ['kr260610_00054.fits.gz', 'kr260610_00058.fits.gz'],
            'frametype': ['tilt, science', 'tilt, science'],
            'target': ['J0750+6927', 'J0750+6927'],
            'mjd': [59742.123456, 59742.654321],
            'comb_id': [1, 2],
        }),
        setup={'Setup A': {'binning': '2,2'}},
    )
    spec = load_spectrograph('keck_kcrm')

    # Only the first exposure has been reduced so far.
    expected_basename = outputfiles.construct_basename(
        'kr260610_00054.fits.gz', spec.camera, spec.allowed_extensions, target='J0750+6927',
        mjd=59742.123456
    )
    spec2d_path = tmp_path / f'spec2d_{expected_basename}.fits'
    _write_spec2d(spec2d_path, target='J0750+6927')

    files, missing, target_name = outputfiles.existing_spec2d_files(
        pypeit_file, 'J0750+6927', tmp_path, spec
    )
    assert files == [spec2d_path], 'the one reduced spec2d file should be found'
    assert missing == ['kr260610_00058'], \
        'the second comb_id group has no reduced spec2d file yet and should be reported missing'
    assert target_name == 'J0750+6927', 'should return the literal target string from the data table'
