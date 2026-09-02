"""
Tests on io module
"""
from IPython import embed

from pypeit import io


def test_remove_suffix():
    assert io.remove_suffix('unzipped_file.txt') == 'unzipped_file', 'bad unzipped removal'
    assert io.remove_suffix('/path/to/unzipped_file.fits') == 'unzipped_file', \
            'bad path and/or suffix removal for unzipped file'
    assert io.remove_suffix('dot.separated.file.name.txt') == 'dot.separated.file.name', \
            'bad many suffix removal'
    assert io.remove_suffix('gzipped_file.fits.gz') == 'gzipped_file', 'bad gzipped removal'
    assert io.remove_suffix('bz2_file.fits.bz2') == 'bz2_file', 'bad bz2 removal'
    assert io.remove_suffix('fpacked_file.fits.fz') == 'fpacked_file', 'bad fz removal'


