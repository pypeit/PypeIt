"""
Tests on io module
"""
from pathlib import Path
import tempfile
import textwrap

from astropy.table import Table
from IPython import embed
import pytest
from pypeit import PypeItError
import numpy as np
from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.par.pypeitpar import PypeItPar

from pypeit import dataPaths
from pypeit import inputfiles
from pypeit.tests import tstutils


def _pypeitfile_components():
    """
    Bits needed to generate a PypeIt file
    """
    confdict = {'rdx': {'spectrograph': 'keck_hires'}}

    data = Table()
    data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
    data['frametype'] = ['arc', 'science']
    data['exptime'] = [1., 10.]

    file_paths = [tstutils.data_output_path('')]

    setup_dict = {'Setup A': ' '}

    # Return
    return confdict, data, file_paths, setup_dict


def test_grab_rawfiles():

    tst_file = Path(tstutils.data_output_path('test.rawfiles')).absolute()
    if tst_file.exists():
        tst_file.unlink()

    # Download and move all the b*fits.gz files into the local package
    # installation
    tstutils.install_shane_kast_blue_raw_data()

    root = Path(tstutils.data_output_path('')).absolute()
    raw_files = [root / 'b11.fits.gz', root / 'b12.fits.gz']
    assert all([f.exists() for f in raw_files]), 'Files missing'

    tbl = Table()
    tbl['filename'] = [r.name for r in raw_files]
    inputfiles.RawFiles(file_paths=[str(root)], data_table=tbl).write(tst_file)

    _raw_files = inputfiles.grab_rawfiles(file_of_files=str(tst_file))
    assert all(f == Path(_f) for f, _f in zip(raw_files, _raw_files)), \
        'Bad file_of_files read'

    _raw_files = inputfiles.grab_rawfiles(list_of_files=tbl['filename'], raw_paths=[str(root)])
    assert all(f == Path(_f) for f, _f in zip(raw_files, _raw_files)), \
        'Bad list_of_files read'

    _raw_files = inputfiles.grab_rawfiles(raw_paths=[str(root)], extension='.fits.gz')
    assert len(_raw_files) == 9, 'Found the wrong number of files'
    assert all([str(root / f) in _raw_files for f in tbl['filename']]), \
        'Missing expected files'

    tst_file.unlink()


def test_instantiate_pypeitfile():
    # Test of instantiation
    confdict, data, file_paths, setup_dict = _pypeitfile_components()
    pypeItFile = inputfiles.PypeItFile(
        config=confdict, file_paths=file_paths, data_table=data, setup=setup_dict
    )
    # Data files                                    
    filenames = pypeItFile.filenames
    assert 'b1.fits.gz' in filenames[0], "First filename should contain 'b1.fits.gz'"

    # Frame types
    frame_type_dict = pypeItFile.frametypes
    assert frame_type_dict['b1.fits.gz'] == 'arc', "Frametype for 'b1.fits.gz' should be 'arc'"

    # More tests
    assert pypeItFile.setup_name == 'A', "Setup name should be 'A'"


def test_read_pypeitfile():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file.pypeit')
    pypeItFile = inputfiles.PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict), 'Configuration should be parsed into a dict'


def test_read_backwards_pypeitfile():
    # Read the PypeIt file (backwards compatibility)
    ifile = dataPaths.tests.get_file_path('example_pypeit_file_backwards.pypeit')
    pypeItFile = inputfiles.PypeItFile.from_file(ifile)
    assert isinstance(pypeItFile.config, dict), 'Backwards-format PypeIt file should parse to a dict'


def test_write_pypeitfile():
    # Test writing a PypeIt file
    outfile = Path(tstutils.data_output_path('tmp_file.pypeit')).absolute()
    if outfile.is_file():
        outfile.unlink()

    # Instantiate
    confdict, data, file_paths, setup_dict = _pypeitfile_components()
    pypeItFile = inputfiles.PypeItFile(
        config=confdict, file_paths=file_paths, data_table=data, setup=setup_dict,
    )
    # Write
    pypeItFile.write(outfile)
    assert outfile.is_file(), 'PypeIt file not written'

    # Let's read it too
    pypeItFile2 = inputfiles.PypeItFile.from_file(outfile)
    assert dict(pypeItFile2.config) == confdict, 'Config dict not the same after read/write'
    assert all(pypeItFile2.data == data), 'Data table not the same after read/write'
    assert (str(pypeItFile2.file_paths[0]) == file_paths[0]), \
        'File paths not the same after read/write'
    assert pypeItFile2.setup == setup_dict, \
        'Setup dict not the same after read/write'

    # Clean up
    outfile.unlink()


def test_path_and_files_commented_and_skip_blank():
    # Prepare a small test directory and files
    root = Path(tstutils.data_output_path('')).absolute()

    f1 = root / 'good.fits'
    f2 = root / 'masked.fits'
    # Ensure clean state
    if f1.exists():
        f1.unlink()
    if f2.exists():
        f2.unlink()
    f1.touch()
    f2.touch()

    tbl = Table()
    tbl['filename'] = ['good.fits', '#masked.fits', '', 'none', 'None']
    tbl['frametype'] = ['science'] * len(tbl['filename'])

    # Instantiate without vetting to avoid requiring full valid config
    pfile = inputfiles.PypeItFile(
        config={'rdx': {}},
        file_paths=[str(root)],
        data_table=tbl,
        setup={'Setup A': ' '},
        vet=False,
    )

    # skip_blank should remove '', 'none', 'None'; include_commented_out False
    # should omit commented entry
    res = pfile.path_and_files(key='filename', skip_blank=True)
    assert len(res) == 1 and Path(res[0]) == f1, \
        'Expected only the un-commented, existing file to be returned'

    # include_commented_out True should include the masked file (without leading '#')
    res2 = pfile.path_and_files(key='filename', skip_blank=True, include_commented_out=True)
    assert any(Path(r) == f2 for r in res2), \
        'include_commented_out should return the masked file without # prefix'

    # Clean up
    f1.unlink()
    f2.unlink()


def test_path_and_files_check_exists_flag():
    root = Path(tstutils.data_output_path('')).absolute()

    f_exists = root / 'exists1.fits'
    if f_exists.exists():
        f_exists.unlink()
    f_exists.touch()

    tbl = Table()
    tbl['filename'] = ['exists1.fits', 'no_such.fits']
    tbl['frametype'] = ['science', 'science']

    pfile = inputfiles.PypeItFile(
        config={'rdx': {}}, file_paths=[str(root)], data_table=tbl, setup={'Setup A': ' '},
        vet=False,
    )

    # With check_exists=False, both entries are returned regardless of on-disk presence
    res = pfile.path_and_files(key='filename', check_exists=False)
    assert any('exists1.fits' in r for r in res), \
        'Existing file should appear in results when check_exists=False'
    assert any('no_such.fits' in r for r in res), \
        'Non-existent filename should be returned when check_exists=False'

    # With check_exists=True, missing file should raise
    with pytest.raises(FileNotFoundError):
        pfile.path_and_files(key='filename')

    # Clean up
    f_exists.unlink()


def test_vet_raises_when_missing_datablock():
    # datablock_required for PypeItFile is True; vet should raise when data is None
    root = Path(tstutils.data_output_path('')).absolute()
    pfile = inputfiles.PypeItFile(
        config={'rdx': {}}, file_paths=[str(root)], data_table=None, setup={'Setup A': ' '},
        vet=False,
    )
    with pytest.raises(PypeItError):
        pfile.vet()


def test_vet_raises_when_required_columns_missing():
    # Missing required column 'frametype' should trigger vet() to raise
    root = Path(tstutils.data_output_path('')).absolute()
    tbl = Table()
    tbl['filename'] = ['a.fits']
    pfile = inputfiles.PypeItFile(
        config={'rdx': {}}, file_paths=[str(root)], data_table=tbl, setup={'Setup A': ' '},
        vet=False,
    )
    with pytest.raises(PypeItError):
        pfile.vet()


def test_remove_comments_and_blanks():
    lines = np.array([
        '',
        '   ',
        '# full comment',
        '  data line 1  # inline comment',
        '\tdata line 2',
        'data#another'
    ])
    out = inputfiles.InputFile.remove_comments_and_blanks(lines)
    assert isinstance(out, np.ndarray), 'Output should be a numpy array'
    assert list(out) == ['data line 1', 'data line 2', 'data'], \
        'Comments and blanks should be removed and inline comments stripped'


def test_get_spectrograph_success_and_missing():
    # Success case
    confdict, data, file_paths, setup_dict = _pypeitfile_components()
    pfile = inputfiles.PypeItFile(
        config=confdict, file_paths=file_paths, data_table=data, setup=setup_dict,
    )
    spec = pfile.get_spectrograph()
    assert isinstance(spec, Spectrograph), 'get_spectrograph should return a Spectrograph instance'

    # Missing spectrograph should raise
    pfile2 = inputfiles.PypeItFile(
        config={'rdx': {}}, file_paths=file_paths, data_table=data, setup=setup_dict, vet=False,
    )
    with pytest.raises(PypeItError):
        pfile2.get_spectrograph()


def test_find_block():
    lines = np.array([
        '# comment',
        'setup read',
        'Setup A: something',
        'setup end',
        'data read',
        'col1 | col2',
        'row1 | row2',
        'data end'
    ])
    s, e = inputfiles.InputFile.find_block(lines, 'setup')
    assert s == 2 and e == 3, \
        'find_block should identify the start and end indices of the setup block'
    s2, e2 = inputfiles.InputFile.find_block(lines, 'data')
    assert s2 == 5 and e2 == 7, \
        'find_block should identify the start and end indices of the data block'


# NOTE: This is a reasonable, AI-generated test of the get_pypeitpar method, but
# it faults because of the ham-fisted way we define the config_specific_file in
# that function.  We should fix that, and then we can uncomment this test.  We
# could also move this test to the dev-suite unit tests so that it has a real
# raw file to go from.
#def test_get_pypeitpar_selects_science_file():
#    # Create example files and a PypeIt table with an arc and a science frame
#    root = Path(tstutils.data_output_path('')).absolute()
#    sci = root / 'sci1.fits'
#    arc = root / 'arc1.fits'
#    if sci.exists():
#        sci.unlink()
#    if arc.exists():
#        arc.unlink()
#    sci.touch()
#    arc.touch()
#
#    confdict, _, _, setup_dict = _pypeitfile_components()
#    tbl = Table()
#    tbl['filename'] = ['arc1.fits', 'sci1.fits']
#    tbl['frametype'] = ['arc', 'science']
#    tbl['exptime'] = [1.0, 10.0]
#
#    pfile = inputfiles.PypeItFile(
#        config=confdict, file_paths=[str(root)], data_table=tbl, setup=setup_dict, vet=False,
#    )
#
#    spec, par, csf = pfile.get_pypeitpar()
#    assert isinstance(spec, Spectrograph), \
#        'get_pypeitpar should return a Spectrograph instance as first element'
#    assert isinstance(par, PypeItPar), \
#        'get_pypeitpar should return a PypeItPar as second element'
#    assert csf.endswith('sci1.fits'), \
#        'get_pypeitpar should select the first science file as config_specific_file'
#
#    # Clean up
#    sci.unlink()
#    arc.unlink()

def test_preserve_comments_in_filenames_property():
    root = Path(tstutils.data_output_path('')).absolute()
    f1 = root / 'good.fits'
    f2 = root / 'masked.fits'
    if f1.exists():
        f1.unlink()
    if f2.exists():
        f2.unlink()
    f1.touch()
    f2.touch()

    tbl = Table()
    tbl['filename'] = ['good.fits', '#masked.fits']
    tbl['frametype'] = ['science', 'science']

    pfile = inputfiles.PypeItFile(
        config={'rdx': {}}, file_paths=[str(root)], data_table=tbl, setup={'Setup A': ' '},
        vet=False,
    )

    # When preserve_comments is True, filenames should include the masked file
    pfile.preserve_comments = True
    fnames = pfile.filenames
    assert any('masked.fits' in f for f in fnames), \
        'When preserve_comments=True, commented-out filenames should be included'

    # When preserve_comments is False, masked file should be omitted
    pfile.preserve_comments = False
    fnames2 = pfile.filenames
    assert all('masked.fits' not in f for f in fnames2), \
        'When preserve_comments=False, commented-out filenames should be omitted'

    # Clean up
    f1.unlink()
    f2.unlink()


def test_readlines_tab_and_missing():
    # Create a temporary file with tabs and trailing spaces
    tf = tempfile.NamedTemporaryFile(delete=False, mode='w')

    try:
        tf.write('Line\twith\ttab  \nSecond\tline  \n')
        tf.close()
        out = inputfiles.InputFile.readlines(tf.name)
        assert isinstance(out, np.ndarray), 'readlines should return a numpy array'
        assert out[0] == 'Line with tab', \
            'Tabs should be replaced by spaces and trailing whitespace stripped'
        # Non-existent file raises
        with pytest.raises(FileNotFoundError):
            inputfiles.InputFile.readlines(tf.name + '.nope')
    finally:
        Path.unlink(tf.name)


def test_from_file_parses_paths_and_table():
    # Build a small PypeIt-format file with path entries and a table
    with tempfile.TemporaryDirectory() as td:
        p = Path(td)
        # Create a dummy data file referenced by the table
        (p / 'a.fits').touch()

        content = textwrap.dedent(
            f"""
            # Example header
            setup read
            Setup A:
                foo: bar
            setup end

            data read
            path {p}
            filename | frametype
            a.fits   | science  
            data end
            """
        )

        fp = p / 'tmp.pypeit'
        fp.write_text(content)

        pfile = inputfiles.PypeItFile.from_file(str(fp), vet=False)
        assert len(pfile.file_paths) == 1 and str(p) == str(pfile.file_paths[0]), \
            'from_file should capture path entries from the data block'
        assert pfile.data is not None and 'filename' in pfile.data.colnames, \
            'Data table should be parsed and contain filename column'


def test_parse_setup_lines_colon_and_multi_error():
    # Missing ':' should be handled
    setups, sdict = inputfiles.InputFile._parse_setup_lines(np.array(['Setup A']))
    assert setups == ['A'] and 'Setup A' in sdict, \
        '_parse_setup_lines should add missing : and parse single Setup'

    # Multiple Setups should raise
    with pytest.raises(PypeItError):
        inputfiles.InputFile._parse_setup_lines(np.array(['Setup A:', 'Setup B:']))


def test_read_data_file_table_paths_and_preserve_comments():
    # Simulate a data block with two path entries and a commented-out row
    with tempfile.TemporaryDirectory() as td1, tempfile.TemporaryDirectory() as td2:
        lines = [
            f'path {td1}',
            f'path {td2}',
            'filename | frametype',
            'good.fits | science',
            '# masked.fits | arc'
        ]

        # When preserve_comments is False the commented row should be ignored
        paths, tbl = inputfiles.InputFile._read_data_file_table(
            np.array(lines), preserve_comments=False
        )
        assert str(td1) in paths and str(td2) in paths, \
            '_read_data_file_table should collect path entries'
        assert len(tbl) == 1, \
            'Commented table rows should be ignored when preserve_comments=False'

        # When preserve_comments is True the commented row should be preserved
        paths2, tbl2 = inputfiles.InputFile._read_data_file_table(
            np.array(lines), preserve_comments=True
        )
        assert len(tbl2) == 2, \
            'Commented table rows should be preserved when preserve_comments=True'


def test_path_and_files_no_file_paths():
    # If file_paths is None, absolute filenames should be accepted
    with tempfile.TemporaryDirectory() as td:
        p = Path(td)

        f = p / 'abs.fits'
        f.touch()

        tbl = Table()
        tbl['filename'] = [str(f)]
        tbl['frametype'] = ['science']

        pfile = inputfiles.PypeItFile(
            config={'rdx': {}}, file_paths=None, data_table=tbl, setup={'Setup A': ' '}, vet=False
        )

        res = pfile.filenames
        assert res is not None and Path(res[0]).name == 'abs.fits', \
            'When file_paths is None absolute filenames should be returned unchanged'


def test_write_preserve_comments_true_and_cfg_lines():
    # Create a PypeIt text with a config comment and preserve_comments=True
    with tempfile.TemporaryDirectory() as td:
        p = Path(td)

        content = textwrap.dedent(
            """
            # KEEPME This comment should be preserved
            setup read
            Setup A:
            setup end

            data read
            | filename | frametype |
            | foo.fits | science   |
            data end
            """
        )
        inp = p / 'keep.pypeit'
        inp.write_text(content)

        pfile = inputfiles.PypeItFile.from_file(str(inp), vet=False, preserve_comments=True)

        out = p / 'out.pypeit'
        pfile.write(out)

        txt = out.read_text()
        assert 'KEEPME' in txt, \
            'Writing with preserve_comments=True should keep original config comments'

        # cfg_lines when config is None
        pfile2 = inputfiles.PypeItFile(
            config=None, file_paths=[str(p)], data_table=Table(), setup={'Setup B': ' '},
            vet=False
        )
        assert pfile2.cfg_lines is None, 'cfg_lines should be None when no config provided'
        assert pfile2.setup_name == 'B', 'setup_name should return the setup letter'


def test_sensfile_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    sfile = inputfiles.SensFile(config={'rdx': {}}, file_paths=[str(root)], data_table=None, setup=None, vet=False)
    sfile.vet()


def test_extractfile_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    efile = inputfiles.ExtractFile(config={'rdx': {}}, file_paths=[str(root)], data_table=None, setup=None, vet=False)
    efile.vet()


def test_fluxfile_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tbl = Table()
    tbl['filename'] = ['f1.fits']
    tbl['frametype'] = ['science']
    ffile = inputfiles.FluxFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl, setup={'Setup A': ' '}, vet=False)
    ffile.vet()
    assert 'sensfile' in ffile.data.colnames and all(v == '' for v in ffile.data['sensfile']), \
        'FluxFile.vet should add an empty sensfile column when not provided'


def test_input_flux_file():
    """Tests for generating and reading fluxing input files
    """
    # Generate an input file
    flux_input_file = Path(tstutils.data_output_path('test.flux')).absolute()
    if flux_input_file.is_file():
        flux_input_file.unlink()

    cfg_lines = ['[fluxcalib]']
    cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
    cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']

    # These files need to be in tests/files/
    data = Table()
    data['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data['sensfile'] = 'sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits'
    # 
    paths = [Path(tstutils.data_output_path('')).absolute()]
    # If pulling from the cache, make sure there are symlinks at the expected path
    for f in data['filename']:
        dataPaths.tests.get_file_path(f, to_pkg='symlink')
    dataPaths.tests.get_file_path(data['sensfile'][0], to_pkg='symlink')

    fluxFile = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data)
    # Write
    fluxFile.write(flux_input_file)

    # Read
    fluxFile2 = inputfiles.FluxFile.from_file(flux_input_file)
    assert np.all(fluxFile2.data['filename'] == data['filename'])

    # Test path
    assert fluxFile2.file_paths[0] == paths[0]
    assert fluxFile2.filenames[0] == str(paths[0] / data['filename'][0])

    # #################
    # Tickle the other ways to do sensfiles
    data3 = Table()
    data3['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data3['sensfile'] = ['sens_cN20170331S0206-HIP62745_GNIRS_20170331T083351.681.fits',
                         '']

    fluxFile3 = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data3)
    assert fluxFile3.sensfiles[1] == str(paths[0] / data['sensfile'][0])
    
    data4 = Table()
    data4['filename'] = ['spec1d_cN20170331S0216-pisco_GNIRS_20170331T085412.181.fits',
                        'spec1d_cN20170331S0217-pisco_GNIRS_20170331T085933.097.fits']
    data4['sensfile'] = ''

    fluxFile4 = inputfiles.FluxFile(config=cfg_lines, 
                        file_paths=paths,
                        data_table=data4)
    assert len(fluxFile4.sensfiles) == 0

    # Clean up
    flux_input_file.unlink()


def test_coadd1dfile_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tbl2 = Table()
    tbl2['filename'] = ['a.fits']
    with pytest.raises(PypeItError):
        inputfiles.Coadd1DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl2, setup=None, vet=True)
    tbl2['obj_id'] = [1]
    c1 = inputfiles.Coadd1DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl2, setup=None, vet=False)
    c1.vet()


def test_coadd2d_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tbl3 = Table()
    tbl3['filename'] = ['b.fits']
    with pytest.raises(PypeItError):
        inputfiles.Coadd2DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=True)
    c2 = inputfiles.Coadd2DFile(config={'rdx': {'spectrograph': 'keck_hires'}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=False)
    c2.vet()


def test_coadd3d_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tbl3 = Table()
    tbl3['filename'] = ['b.fits']
    with pytest.raises(PypeItError):
        inputfiles.Coadd3DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=True)
    c3 = inputfiles.Coadd3DFile(config={'rdx': {'spectrograph': 'keck_hires'}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=False)
    c3.vet()


def test_telluric_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tfile = inputfiles.TelluricFile(config={'rdx': {}}, file_paths=[str(root)], data_table=None, setup=None, vet=False)
    tfile.vet()


def test_flexure_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    tbl3 = Table()
    tbl3['filename'] = ['b.fits']
    with pytest.raises(PypeItError):
        inputfiles.FlexureFile(config={'rdx': {}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=True)
    flex = inputfiles.FlexureFile(config={'rdx': {'spectrograph': 'keck_hires'}}, file_paths=[str(root)], data_table=tbl3, setup=None, vet=False)
    flex.vet()


def test_collate1d_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    with pytest.raises(PypeItError):
        inputfiles.Collate1DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=Table(), setup=None, vet=True)
    rtbl = Table()
    rtbl['filename'] = ['r1.fits']
    coll = inputfiles.Collate1DFile(config={'rdx': {}}, file_paths=[str(root)], data_table=rtbl, setup=None, vet=False)
    coll.vet()


def test_rawfiles_basic():
    root = Path(tstutils.data_output_path('')).absolute()
    with pytest.raises(PypeItError):
        inputfiles.RawFiles(config={'rdx': {}}, file_paths=[str(root)], data_table=Table(), setup=None, vet=True)
    rtbl = Table()
    rtbl['filename'] = ['r1.fits']
    raw = inputfiles.RawFiles(config={'rdx': {}}, file_paths=[str(root)], data_table=rtbl, setup=None, vet=False)
    raw.vet()

