"""
Tests on io module
"""
from pathlib import Path

from IPython import embed

from astropy.table import Table

from pypeit import io
from pypeit import inputfiles
from pypeit.tests.tstutils import data_path


def test_remove_suffix():
    assert io.remove_suffix('unzipped_file.txt') == 'unzipped_file', 'bad unzipped removal'
    assert io.remove_suffix('/path/to/unzipped_file.fits') == 'unzipped_file', \
            'bad path and/or suffix removal for unzipped file'
    assert io.remove_suffix('dot.separated.file.name.txt') == 'dot.separated.file.name', \
            'bad many suffix removal'
    assert io.remove_suffix('gzipped_file.fits.gz') == 'gzipped_file', 'bad gzipped removal'


def test_grab_rawfiles():

    tst_file = Path(data_path('test.rawfiles')).resolve()
    if tst_file.exists():
        tst_file.unlink()

    root = Path(data_path('')).resolve()
    raw_files = [root / 'b11.fits.gz', root / 'b12.fits.gz']
    assert all([f.exists() for f in raw_files]), 'Files missing'

    tbl = Table()
    tbl['filename'] = [r.name for r in raw_files]
    inputfiles.RawFiles(file_paths=[str(root)], data_table=tbl).write(tst_file)

    _raw_files = inputfiles.grab_rawfiles(file_of_files=str(tst_file))
    assert [str(f) for f in raw_files] == _raw_files, 'Bad file_of_files read'

    _raw_files = inputfiles.grab_rawfiles(list_of_files=tbl['filename'], raw_paths=[str(root)])
    assert [str(f) for f in raw_files] == _raw_files, 'Bad list_of_files read'

    _raw_files = inputfiles.grab_rawfiles(raw_paths=[str(root)], extension='.fits.gz')
    assert len(_raw_files) == 8, 'Found the wrong number of files'
    assert all([str(root / f) in _raw_files for f in tbl['filename']]), 'Missing expected files'

