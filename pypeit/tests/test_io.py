"""
Tests on io module
"""
from IPython import embed

from pypeit import io
from pypeit.scripts.versions import Versions


def test_remove_suffix():
    assert io.remove_suffix('unzipped_file.txt') == 'unzipped_file', 'bad unzipped removal'
    assert io.remove_suffix('/path/to/unzipped_file.fits') == 'unzipped_file', \
            'bad path and/or suffix removal for unzipped file'
    assert io.remove_suffix('dot.separated.file.name.txt') == 'dot.separated.file.name', \
            'bad many suffix removal'
    assert io.remove_suffix('gzipped_file.fits.gz') == 'gzipped_file', 'bad gzipped removal'


def test_runtime_versions_in_header():
    hdr = io.initialize_header()
    for _, keyword, package_version, comment in io.runtime_versions():
        assert hdr[keyword] == package_version
        assert hdr.comments[keyword] == comment


def test_versions_script(capsys):
    Versions.main(None)
    output = capsys.readouterr().out.splitlines()
    assert output == [f'{package}: {package_version}'
                      for package, _, package_version, _ in io.runtime_versions()]

