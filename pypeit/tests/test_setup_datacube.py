from astropy.io import fits
import pytest

from pypeit import PypeItError
from pypeit.scripts.setup_datacube import SetupDataCube, _parse_whitelight_range


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


def test_parse_whitelight_range():
    assert _parse_whitelight_range('9400,10000') == [9400.0, 10000.0], \
        'should parse a numeric min,max pair into a float list'
    assert _parse_whitelight_range('None,None') == [None, None], \
        'a "none" entry (case-insensitive) should parse to None'
    assert _parse_whitelight_range('9400,None') == [9400.0, None], \
        'the two entries should be parsed independently'
    with pytest.raises(ValueError):
        _parse_whitelight_range('9400')  # only one entry


def test_setup_datacube_append_missing_coadd3d_file_raises(tmp_path):
    science_dir = tmp_path / 'Science'
    science_dir.mkdir()
    pypeit_file = tmp_path / 'kcrm_jun10_hizqso.pypeit'
    _write_pypeit_file(pypeit_file)
    first_spec2d = science_dir / 'spec2d_kr260610_00054-J0750+6927_KCRM_test.fits'
    _write_spec2d(first_spec2d)

    # --append before any .coadd3d file has been written should raise a clear,
    # append-specific error, not a generic FileNotFoundError from deep inside
    # inputfiles.InputFile.from_file.
    args = SetupDataCube.parse_args(['--append', str(pypeit_file), 'J0750+6927'])
    with pytest.raises(PypeItError, match='Cannot append to missing .coadd3d file'):
        SetupDataCube.main(args)


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

    assert 'output_filename = J0750+6927' in coadd3d_text, \
        'the target should be used as the coadd3d output_filename'
    assert 'whitelight_range = None, None' in coadd3d_text, \
        'the default whitelight_range should be written as an unset [None, None] pair'
    assert 'whitelight_range = None, None' in extract_text, \
        'the extract file should carry the same default whitelight_range'
    assert 'sensfile =' not in coadd3d_text, \
        'no sensfile line should be written when --sensfile is not given'
    assert '# weights_init_obj_pos = x:y' in coadd3d_text, \
        'the commented-out weights_init_obj_pos example should still be included'
    assert first_spec2d.name in coadd3d_text, \
        'the one reduced spec2d file should be listed in the coadd3d data block'
    assert 'kr260610_00058' not in coadd3d_text, \
        'a comb_id group with no reduced spec2d file yet should not appear'
    assert 'opt_prof_method = fit_gauss' in extract_text, \
        'opt_prof_method should default to fit_gauss when --manual is not given'
    assert 'manual =' not in extract_text, \
        'no manual line should be written when --manual is not given'

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
    assert 'weight_method = auto' in appended_text, \
        'the user edit to weight_method should survive --append'
    assert appended_text.count(first_spec2d.name) == 1, \
        'the original filename should not be duplicated by --append'
    assert appended_text.count(second_spec2d.name) == 1, \
        'the newly reduced filename should be appended'
    assert extract_file.read_text() == edited_extract, \
        '--append should leave the .extract file completely untouched'

    args = SetupDataCube.parse_args([
        str(pypeit_file), 'J0750+6927', '--wl_range', '9400,10000'
    ])
    assert args.whitelight_range == [9400.0, 10000.0], \
        '--wl_range should parse into a two-element float list'

    args = SetupDataCube.parse_args([
        str(pypeit_file), 'J0750+6927', '--sensfile', str(sensfile), '-o'
    ])
    SetupDataCube.main(args)
    assert 'sensfile = ' + str(sensfile.absolute()) in coadd3d_file.read_text(), \
        '--sensfile should be written into the refreshed .coadd3d file'
    assert '# user edit' not in extract_file.read_text(), \
        '-o/--overwrite should replace the .extract file'


def test_setup_datacube_manual_validation(tmp_path):
    science_dir = tmp_path / 'Science'
    science_dir.mkdir()
    pypeit_file = tmp_path / 'kcrm_jun10_hizqso.pypeit'
    _write_pypeit_file(pypeit_file)
    first_spec2d = science_dir / 'spec2d_kr260610_00054-J0750+6927_KCRM_test.fits'
    _write_spec2d(first_spec2d)

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927', '--manual', '9.8:13.6'])
    SetupDataCube.main(args)
    extract_text = (tmp_path / 'sources' / 'J0750+6927' / 'J0750+6927.extract').read_text()
    assert 'manual = 9.8:13.6' in extract_text, \
        '--manual value should be written into the .extract file'
    assert 'opt_prof_method = user_gauss' in extract_text, \
        'providing --manual should switch opt_prof_method to user_gauss'

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927', '--manual', '9.8,13.6'])
    with pytest.raises(PypeItError, match='colon-separated x:y'):
        SetupDataCube.main(args)

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927', '--manual', '9.8:13.6:1.0'])
    with pytest.raises(PypeItError, match='x and y'):
        SetupDataCube.main(args)

    args = SetupDataCube.parse_args([str(pypeit_file), 'J0750+6927', '--manual', 'abc:def'])
    with pytest.raises(PypeItError, match='numeric'):
        SetupDataCube.main(args)
