from pathlib import Path

from astropy.io import fits

from pypeit.scripts.setup_datacube import SetupDataCube


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

    args = SetupDataCube.parse_args([
        str(pypeit_file), 'J0750+6927', '9400,10000', str(sensfile)
    ])
    SetupDataCube.main(args)

    source_dir = tmp_path / 'sources' / 'J0750p6927'
    coadd3d_file = source_dir / 'J0750p6927.coadd3d'
    extract_file = source_dir / 'J0750p6927.extract'
    coadd3d_text = coadd3d_file.read_text()
    extract_text = extract_file.read_text()

    assert 'output_filename = J0750p6927' in coadd3d_text
    assert 'whitelight_range = 9400,10000' in coadd3d_text
    assert 'sensfile = ' + str(sensfile.absolute()) in coadd3d_text
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

    args = SetupDataCube.parse_args([
        '--append', str(pypeit_file), 'J0750+6927', '9400,10000', str(sensfile)
    ])
    SetupDataCube.main(args)

    appended_text = coadd3d_file.read_text()
    assert 'weight_method = auto' in appended_text
    assert appended_text.count(first_spec2d.name) == 1
    assert appended_text.count(second_spec2d.name) == 1
    assert extract_file.read_text() == edited_extract
