"""
Test the Keck/MOSFIRE read noise from the readout mode.
"""
import numpy as np
import pytest

from astropy.io import fits

from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs.keck_mosfire import mosfire_read_noise, MOSFIRE_READ_NOISE, KeckMOSFIRESpectrograph


def _hdu(**cards):
    hdr = fits.Header()
    for key, value in cards.items():
        hdr[key] = value
    return fits.HDUList([fits.PrimaryHDU(header=hdr)])


def test_mosfire_read_noise_table():
    # CDS is one read pair, whatever NUMREADS says
    assert mosfire_read_noise(2, 1) == 21.0
    assert mosfire_read_noise(2, 16) == 21.0
    # Tabulated MCDS values are returned exactly
    for n, rn in MOSFIRE_READ_NOISE.items():
        if n > 1:
            assert mosfire_read_noise(3, n) == pytest.approx(rn)
    # Linear in log2(N) between entries: MCDS-2 is halfway (in log2) between CDS and MCDS-4,
    # MCDS-12 is at log2(12) between MCDS-8 and MCDS-16
    assert mosfire_read_noise(3, 2) == pytest.approx(0.5 * (21.0 + 10.8))
    f = np.log2(12) - 3
    assert mosfire_read_noise(3, 12) == pytest.approx(7.7 + f * (5.8 - 7.7))
    # Held at the end values outside the table
    assert mosfire_read_noise(3, 256) == pytest.approx(3.0)
    # Single, UTR, non-positive or missing reads are not tabulated
    assert mosfire_read_noise(1, 1) is None
    assert mosfire_read_noise(4, 16) is None
    assert mosfire_read_noise(3, 0) is None
    assert mosfire_read_noise(None, 16) is None
    assert mosfire_read_noise(3, None) is None


def test_mosfire_detector_ronoise():
    spec = load_spectrograph('keck_mosfire')
    # No header: the MCDS-16 default
    assert spec.get_detector_par(1).ronoise[0] == pytest.approx(5.8)
    # CDS (e.g. the 2022-04-09 dome flats) and MCDS-16 (science and standard frames)
    assert spec.get_detector_par(1, hdu=_hdu(SAMPMODE=2, NUMREADS=1)).ronoise[0] == pytest.approx(21.0)
    assert spec.get_detector_par(1, hdu=_hdu(SAMPMODE=3, NUMREADS=16)).ronoise[0] == pytest.approx(5.8)
    assert spec.get_detector_par(1, hdu=_hdu(SAMPMODE=3, NUMREADS=4)).ronoise[0] == pytest.approx(10.8)
    # Untabulated modes and a header without the cards fall back to the default
    assert spec.get_detector_par(1, hdu=_hdu(SAMPMODE=4, NUMREADS=16)).ronoise[0] == pytest.approx(5.8)
    assert spec.get_detector_par(1, hdu=_hdu()).ronoise[0] == pytest.approx(5.8)


LONG2POS_BARS = [(26, 0.7, -45.556), (25, 4.0, -44.998), (24, 0.7, -44.44), (23, 4.0, 0.0),
                 (22, 0.7, 44.438), (21, 4.0, 44.996), (20, 0.7, 45.554)]


def test_long2pos_match_bars():
    match = KeckMOSFIRESpectrograph.long2pos_match_bars
    # 2014-06-01 J: every bar traced (CRPIX2 1035, PSCALE 0.1799)
    assert match(LONG2POS_BARS, [902, 945, 989, 1033, 1078, 1122, 1165], 1035, 0.1799) == [0, 1, 2, 3, 4, 5, 6]
    # 2014-06-01 H: bar 26 not traced
    assert match(LONG2POS_BARS, [943, 987, 1031, 1076, 1120, 1163], 1035, 0.1799) == [None, 0, 1, 2, 3, 4, 5]
    # 2015-09-17 H: a spurious 847-pixel slit above the mask
    assert match(LONG2POS_BARS, [903, 943, 987, 1032, 1076, 1120, 1164, 1610], 1035, 0.1799) \
        == [0, 1, 2, 3, 4, 5, 6]
    # a slit serves one bar only, and nothing matches far from the mask
    assert match(LONG2POS_BARS, [1033], 1035, 0.1799) == [None, None, None, 0, None, None, None]
    assert match(LONG2POS_BARS, [], 1035, 0.1799) == [None] * 7
