"""Tests for Binospec static bad pixel masks."""
import numpy as np
import pytest
from astropy.io import fits

from pypeit import dataPaths
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECSpectrograph


@pytest.mark.parametrize('det', [1, 2])
def test_bpm_file_exists(det):
    """Both static BPM files exist in static_calibs."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    assert bpm_path.exists(), f'BPM file not found: {bpm_path}'


@pytest.mark.parametrize('det', [1, 2])
def test_bpm_shape(det):
    """Each BPM file has the correct NumPy shape and dtype."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    bpm = fits.getdata(bpm_path)
    assert bpm.shape == (4096, 4112)
    assert bpm.dtype == np.int8


@pytest.mark.parametrize('det,min_bad,max_bad', [
    (1, 30000, 80000),
    (2, 15000, 50000),
])
def test_bpm_nonzero_count(det, min_bad, max_bad):
    """Each BPM has a reasonable number of masked pixels."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    bpm = fits.getdata(bpm_path)
    n_bad = np.sum(bpm > 0)
    assert min_bad < n_bad < max_bad, \
        f'det={det}: {n_bad} bad pixels outside expected range [{min_bad}, {max_bad}]'


@pytest.mark.parametrize('det', [1, 2])
def test_bpm_applied(det):
    """bpm() returns a mask that includes the static BPM pixels."""
    spec = MMTBINOSPECSpectrograph()
    # Use shape parameter to avoid needing a raw data file
    bpm_img = spec.bpm(filename=None, det=det, shape=(4096, 4112))

    # Load the static file directly for comparison
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    static_bpm = fits.getdata(bpm_path)

    # Every pixel marked in the static file should be marked in bpm_img
    assert np.all(bpm_img[static_bpm > 0] > 0)
