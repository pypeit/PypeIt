"""
Module to run tests on sensfilearchive code.
"""

import pytest
from functools import partial
from pathlib import Path

from pypeit.sensfilearchive import SensFileArchive
from astropy.io import fits
from pypeit.pypmsgs import PypeItError

def test_getinstance():
    # Test success
    sfa = SensFileArchive.get_instance("keck_deimos")

    assert sfa is not None

    # Test failure
    with pytest.raises(ValueError, match = "No SensFileArchive found for invalid_spec"):
        sfa = SensFileArchive.get_instance("invalid_spec")

def test_supported_spectrgraphs():
    assert SensFileArchive.supported_spectrographs() == ['keck_deimos']

def test_get_archived_sensfile(monkeypatch):
    sfa = SensFileArchive.get_instance("keck_deimos")
    
    def mock_get_header(file, mock_return):
        return {"DISPNAME": mock_return}

    # Test success by making sure the file with the right name is returned, and that the file exists
    with monkeypatch.context() as m:
        monkeypatch.setattr(fits, "getheader", partial(mock_get_header, mock_return='600ZD'))
        # Note: `symlink_in_pkgdir=True` is needed to appease the `assert` below
        sensfile = Path(sfa.get_archived_sensfile("test_file.fits", symlink_in_pkgdir=True))
        assert sensfile.exists() is True and sensfile.name == "keck_deimos_600ZD_sensfunc.fits"

    # Test failure
    with monkeypatch.context() as m:
        monkeypatch.setattr(fits, "getheader", partial(mock_get_header, mock_return='should_fail'))
        with pytest.raises(PypeItError):
            sensfile = Path(sfa.get_archived_sensfile("test_file.fits"))
