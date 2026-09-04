from pathlib import Path

from astropy.io import fits
from astropy.table import Table
import pytest

from pypeit import PypeItError
from pypeit import inputfiles
from pypeit import outputfiles
from pypeit.scripts import setup_datacube
from pypeit.scripts.setup_datacube import SetupDataCube
from pypeit.spectrographs.util import load_spectrograph


def _write_spec2d(path, target='J0750+6927'):
    hdr = fits.Header()
    hdr['TARGET'] = target
    hdr['PYP_SPEC'] = 'keck_kcrm'
    hdr['PYPELINE'] = 'SlicerIFU'
    fits.PrimaryHDU(header=hdr).writeto(path)


def _write_pypeit_file(path):
    path.write_text(
        '\n'.join([
            '[rdx]',
            '    spectrograph = keck_kcrm',
            '',
            'setup read',
            'Setup A:',
            '  binning: 2,2',
            'setup end',
            '',
            'data read',
            ' path /tmp/raw',
            ' filename | frametype | target | comb_id',
            ' kr260610_00054.fits | tilt, science | J0750+6927 | 1',
            ' kr260610_00055.fits | tilt, science | J0750+6927 | 1',
            ' kr260610_00058.fits | tilt, science | J0750+6927 | 2',
            ' kr260610_00059.fits | tilt, science | J0750+6927 | 2',
            ' kr260610_00062.fits | tilt, science | J0913+6007 | 3',
            'data end',
            ''
        ])
    )


def test_setup_datacube_write_and_append(tmp_path):
    science_dir = tmp_path / 'Science'
    science_dir.mkdir()
    sensfile = tmp_path / 'sens_gd71_000.fits'
    fits.PrimaryHDU().writeto(sensfile)
    pypeit_file = tmp_path / 'kcrm_jun10_hizqso.pypeit'
    _write_pypeit_file(pypeit_file)

    first_spec2d = science_dir / 'spec2d_kr260610_00054-J0750+6927_KCRM_test.fits'
    _write_spec2d(first_spec2d)

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927'])
    SetupDataCube.main(args)

    source_dir = tmp_path / 'sources' / 'J0750+6927'
    coadd3d_file = source_dir / 'J0750+6927.coadd3d'
    extract_file = source_dir / 'J0750+6927.extract'
    coadd3d_text = coadd3d_file.read_text()
    extract_text = extract_file.read_text()

    assert 'output_filename = J0750+6927' in coadd3d_text
    assert 'whitelight_range = None,None' in coadd3d_text
    assert 'whitelight_range = None,None' in extract_text
    assert 'sensfile =' not in coadd3d_text
    assert '# weights_init_obj_pos = x:y' in coadd3d_text
    assert first_spec2d.name in coadd3d_text
    assert 'kr260610_00058' not in coadd3d_text
    assert 'opt_prof_method = fit_gauss' in extract_text
    assert 'manual =' not in extract_text

    # Simulate a user edit and a later-reduced second comb_id product.
    edited_coadd3d = coadd3d_text.replace('weight_method = uniform', 'weight_method = auto')
    coadd3d_file.write_text(edited_coadd3d)
    edited_extract = extract_text + '# user edit\n'
    extract_file.write_text(edited_extract)
    second_spec2d = science_dir / 'spec2d_kr260610_00058-J0750+6927_KCRM_test.fits'
    _write_spec2d(second_spec2d)

    args = SetupDataCube.parse_args(['--append', str(pypeit_file), 'J0750+6927'])
    SetupDataCube.main(args)

    appended_text = coadd3d_file.read_text()
    assert 'weight_method = auto' in appended_text
    assert appended_text.count(first_spec2d.name) == 1
    assert appended_text.count(second_spec2d.name) == 1
    assert extract_file.read_text() == edited_extract

    args = SetupDataCube.parse_args([
        str(pypeit_file), 'J0750+6927', '--wl_range', '9400,10000'
    ])
    assert args.whitelight_range == '9400,10000'

    args = SetupDataCube.parse_args([
        str(pypeit_file), 'J0750+6927', '--sensfile', str(sensfile), '-o'
    ])
    SetupDataCube.main(args)
    assert 'sensfile = ' + str(sensfile.absolute()) in coadd3d_file.read_text()
    assert '# user edit' not in extract_file.read_text()

    alias_rows = setup_datacube.matching_science_rows(
        inputfiles.PypeItFile.from_file(str(pypeit_file)), 'J0750p6927'
    )
    assert len(alias_rows) == 4


def test_find_reduced_spec2d_exact_match_primary_path(tmp_path):
    # The primary path should find the file from the row's own metadata alone, with
    # no header read -- so a deliberately WRONG TARGET header must not prevent a match.
    spec = load_spectrograph('keck_kcrm')
    row = Table(
        {'filename': ['kr260610_00054.fits.gz'], 'target': ['J0750+6927'], 'mjd': [59742.123456]}
    )[0]
    expected_basename = outputfiles.construct_basename(
        row['filename'], row['target'], spec.camera, row['mjd'], spec.allowed_extensions
    )
    spec2d_path = tmp_path / f'spec2d_{expected_basename}.fits'
    _write_spec2d(spec2d_path, target='not-the-requested-target')

    assert setup_datacube.find_reduced_spec2d(tmp_path, row, spec) == spec2d_path


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

    assert setup_datacube.find_reduced_spec2d(tmp_path, row, spec) == spec2d_path

    # A row with no target match in the header, and no mjd column at all, still falls
    # back correctly (this is the original bug: a naive `Path(filename).stem` would
    # have left a stray `.fits` on the glob prefix for a `.fits.gz` raw file).
    no_mjd_row = Table({'filename': ['kr260610_00058.fits.gz'], 'target': ['J0750+6927']})[0]
    assert setup_datacube.find_reduced_spec2d(tmp_path, no_mjd_row, spec) == spec2d_path


def test_setup_datacube_exact_match_end_to_end(tmp_path):
    # Full CLI-level regression test for the bug that motivated this integration: a
    # compound (.fits.gz) raw extension, matched via the primary exact-match path. The
    # spec2d header's TARGET is deliberately wrong, so a successful run here can only be
    # explained by the exact-match path (the fallback would reject this file).
    science_dir = tmp_path / 'Science'
    science_dir.mkdir()
    mjd = 59742.123456
    pypeit_file = tmp_path / 'kcrm_jun10_hizqso.pypeit'
    pypeit_file.write_text(
        '\n'.join([
            '[rdx]',
            '    spectrograph = keck_kcrm',
            '',
            'setup read',
            'Setup A:',
            '  binning: 2,2',
            'setup end',
            '',
            'data read',
            ' path /tmp/raw',
            ' filename | frametype | target | mjd | comb_id',
            f' kr260610_00054.fits.gz | tilt, science | J0750+6927 | {mjd} | 1',
            'data end',
            ''
        ])
    )

    spec = load_spectrograph('keck_kcrm')
    expected_basename = outputfiles.construct_basename(
        'kr260610_00054.fits.gz', 'J0750+6927', spec.camera, mjd, spec.allowed_extensions
    )
    spec2d_path = science_dir / f'spec2d_{expected_basename}.fits'
    _write_spec2d(spec2d_path, target='not-the-requested-target')

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927'])
    SetupDataCube.main(args)

    coadd3d_text = (tmp_path / 'sources' / 'J0750+6927' / 'J0750+6927.coadd3d').read_text()
    assert spec2d_path.name in coadd3d_text


def test_setup_datacube_manual_validation(tmp_path):
    extract_file = tmp_path / 'J0750+6927.extract'

    setup_datacube.write_extract_file(
        extract_file, 'J0750+6927', 'None,None', manual='9.8:13.6'
    )
    extract_text = extract_file.read_text()
    assert 'manual = 9.8:13.6' in extract_text
    assert 'opt_prof_method = user_gauss' in extract_text

    with pytest.raises(PypeItError, match='colon-separated x:y'):
        setup_datacube.write_extract_file(
            extract_file, 'J0750+6927', 'None,None', manual='9.8,13.6'
        )
